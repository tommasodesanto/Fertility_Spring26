"""Bounded nine-coordinate refinement with beta_annual free on [.94, .99].

The controller only orchestrates frozen scored-candidate evaluations.  It never
changes targets, weights, numerical gates, or the original economic source.
"""
import argparse
import concurrent.futures as cf
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time

for _n in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_n] = '1'
import numpy as np
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_profile as support

APPROVED_OBJECTIVE = support.APPROVED_OBJECTIVE
ALL_NAMES = support.ALL_NAMES
BETA_LOWER, BETA_UPPER = .94, .99
WRAPPER_SECONDS = 2100


def read(p): return support.read(p)
def write(p, v): return support.write(p, v)
def sha(p): return support.sha(p)
def numeric_signature(score): return support.numeric_signature(score)
def parameters(score): return support.parameters(score)
def residual(score): return support.residual(score)
def transform(name, value): return support.transform(name, value)
def inverse(name, value): return support.inverse(name, value)


def metadata():
    return dict(beta_annual_status='estimated_capped', beta_lower=BETA_LOWER, beta_upper=BETA_UPPER,
                free_parameter_count=9, free_parameters=list(ALL_NAMES), scored_moment_count=12,
                separate_normalization_target=2.1, objective_canonical_sha256=APPROVED_OBJECTIVE,
                production_eligible=False, identification_established=False)


def validate_plan(p):
    required = ('source_root', 'smoke_summary_sha256', 'smoke_contract_sha256', 'resume_score_path',
                'resume_score_sha256', 'resume_summary_sha256', 'resume_proposal', 'controller_sha256',
                'file_sha256', 'output_dir', 'seed_manifest_path', 'seed_manifest_sha256', 'beta_upper', 'beta_lower', 'workers', 'rounds',
                'maximum_search_cases', 'search_seconds')
    for key in required:
        if key not in p: raise ValueError('Missing plan field: ' + key)
    if p['beta_lower'] != BETA_LOWER or p['beta_upper'] != BETA_UPPER: raise ValueError('Unexpected beta cap')
    if p['workers'] != 18 or p['rounds'] != 3 or p['maximum_search_cases'] != 90: raise ValueError('Plan exceeds bounded controller')
    if p['search_seconds'] != 7800: raise ValueError('Search deadline must reserve final verification')
    if set(p['resume_proposal'].get('parameters', {})) != set(ALL_NAMES): raise ValueError('Incomplete resume proposal')


def validate_seed_manifest(plan):
    path = Path(plan['seed_manifest_path'])
    if sha(path) != plan['seed_manifest_sha256']: raise ValueError('Seed manifest changed')
    manifest = read(path); pins = manifest.get('file_sha256', {})
    if len(pins) != 5: raise ValueError('Seed manifest must pin five evidence files')
    for filename, pin in pins.items():
        if sha(filename) != pin: raise ValueError('Seed evidence changed: ' + filename)
    by_name = {Path(k).name: Path(k) for k in pins}
    # The score is distinguished by its scored-repetition parent; the other four
    # names are unique in the five-file packet.
    score_path = next(Path(k) for k in pins if Path(k).parent.name == 'scored_repetition_01')
    preflight = read(by_name['preflight.json']); initial = read(by_name['initial_contract.json']); run = read(by_name['run_contract.json'])
    summary = read(by_name['summary.json']); score = read(score_path)
    if (preflight.get('status') != 'verified' or preflight.get('objective_canonical_sha256') != APPROVED_OBJECTIVE
            or summary.get('status') != 'verified_scored_candidate' or summary.get('repetitions', 1) != 1
            or summary.get('objective_canonical_sha256') != APPROVED_OBJECTIVE or summary.get('loss') != score.get('loss')):
        raise ValueError('Seed receipt is not the required one-repeat verified score')
    if (preflight.get('run_contract_sha256') != sha(by_name['run_contract.json'])
            or preflight.get('initial_solve_contract_sha256') != sha(by_name['initial_contract.json'])
            or run.get('initial_solve_contract', {}).get('sha256') != sha(by_name['initial_contract.json'])):
        raise ValueError('Seed contract hash mismatch')
    proposal = plan['resume_proposal']
    if initial.get('structural_candidate') != proposal['parameters'] or initial.get('initial_psi') != proposal['initial_psi']:
        raise ValueError('Seed initial input differs from proposal')
    return score


def bounds(names, restrictions, x, radius):
    lo, hi = [], []
    for i, name in enumerate(names):
        r = restrictions[name]
        # beta's log(-log(beta)) coordinate is decreasing in beta.
        lower, upper = ((transform(name, r['upper']), transform(name, r['lower'])) if name == 'beta_annual'
                        else (-math.inf if r['lower'] == 0 else transform(name, r['lower']), transform(name, r['upper'])))
        lo.append(max(-radius, lower - x[i])); hi.append(min(radius, upper - x[i]))
    return np.array(lo), np.array(hi)


def feasible_steps(lower, upper):
    return [(sign, float(step)) for sign, step in ((-1, lower), (1, upper)) if abs(step) >= 1e-10]


def clamp(name, value, restrictions):
    r = restrictions[name]; lo, hi = float(r['lower']), float(r['upper'])
    if name == 'beta_annual': lo, hi = max(lo, BETA_LOWER), min(hi, BETA_UPPER)
    eps = 1e-12 * max(1., abs(lo), abs(hi))
    if value < lo - eps or value > hi + eps: raise ValueError('Proposal outside bound: ' + name)
    return min(hi, max(lo, value))


def validate_proposal(item, restrictions):
    p = item['parameters']
    if set(p) != set(ALL_NAMES): raise ValueError('Incomplete structural proposal')
    for name in ALL_NAMES:
        if not math.isfinite(float(p[name])) or clamp(name, float(p[name]), restrictions) != float(p[name]):
            raise ValueError('Proposal violates bounds: ' + name)


def validate_score(score, item, restrictions):
    validate_proposal(item, restrictions)
    got = parameters(score)
    if any(got[n] != item['parameters'][n] for n in ALL_NAMES): raise ValueError('Scored parameters differ from proposal')
    if score.get('contract_sha256') != APPROVED_OBJECTIVE or score.get('free_parameter_count') != 9: raise ValueError('Mixed score contract')
    if len(score['target_fit']) != 13 or residual(score).shape != (12,): raise ValueError('Incomplete target system')
    norms = [r for r in score['target_fit'] if not r['scored']]
    if len(norms) != 1 or float(norms[0]['target']) != 2.1 or abs(float(norms[0]['gap'])) > 5e-4: raise ValueError('Normalization gate failed')
    if not math.isclose(float(residual(score) @ residual(score)), float(score['loss']), rel_tol=1e-12, abs_tol=1e-12): raise ValueError('Loss mismatch')


def derivative_column(samples, center): return support.derivative_column(samples, center)


def derivative_proposals(center, restrictions, round_index):
    v = parameters(center['score']); x = np.array([transform(n, v[n]) for n in ALL_NAMES]); lo, hi = bounds(ALL_NAMES, restrictions, x, .02)
    ans = []
    for j, name in enumerate(ALL_NAMES):
        for sign, step in feasible_steps(lo[j], hi[j]):
            q = {n: v[n] for n in ALL_NAMES}; q[name] = clamp(name, inverse(name, x[j] + step), restrictions)
            ans.append(dict(case_id=f'r{round_index}_d{j}_{sign:+d}', parameters=q, initial_psi=v['psi_child'], column=j, step=transform(name, q[name])-x[j]))
    return ans


def jacobian(results, center):
    J = np.empty((12, 9)); c = residual(center['score'])
    for j in range(9):
        samples = [(r['proposal']['step'], residual(r['score'])) for r in results if r['status'] == 'verified' and r['proposal']['column'] == j]
        J[:, j] = derivative_column(samples, c)
    return J


def joint_proposals(center, J, restrictions, round_index):
    if J.shape != (12, 9): raise ValueError('Expected complete 12 by 9 Jacobian')
    v = parameters(center['score']); x = np.array([transform(n, v[n]) for n in ALL_NAMES]); r = residual(center['score'])
    penalty = np.diag(np.maximum(np.linalg.norm(J, axis=0), 1e-8)); ans = []
    for radius in (.1, .2, .35, .5):
        lower, upper = bounds(ALL_NAMES, restrictions, x, radius)
        for ridge in (.01, .1, 1.):
            fit = lsq_linear(np.vstack((J, math.sqrt(ridge)*penalty)), np.r_[-r, np.zeros(9)], bounds=(lower, upper), tol=1e-12, max_iter=300)
            if not fit.success: raise RuntimeError('Linear proposal failed')
            q = {n: clamp(n, inverse(n, x[i]+fit.x[i]), restrictions) for i, n in enumerate(ALL_NAMES)}
            ans.append(dict(case_id=f'r{round_index}_joint_{len(ans):02d}', parameters=q, initial_psi=v['psi_child'], radius=radius, ridge=ridge))
    return ans


def _batch(items, worker, completed, workers):
    if len(items) > workers: raise ValueError('A capped stage must fit one 18-worker wave')
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        results = []
        for future in cf.as_completed([pool.submit(worker, item) for item in items]):
            q = future.result(); completed(q); results.append(q)
    return results


def search(plan, seed_score, restrictions, run_case, out):
    """Complete synthetic-testable loop; run_case is the frozen source wrapper in production."""
    out = Path(out); started = time.monotonic(); deadline = started + plan['search_seconds']; records=[]; best=None
    state=dict(phase='seed_exact_repetitions', completed=0, failed=0, status='running'); stop=threading.Event(); attempted=0
    def persist():
        write(out/'best_so_far.json', None if best is None else {k:v for k,v in best.items() if k != 'score'})
        write(out/'latest_completed.json', dict(state, latest=records[-1] if records else None)); write(out/'cases.json', records)
    def beat():
        while not stop.wait(30): write(out/'heartbeat.json', dict(state, elapsed_seconds=time.monotonic()-started))
    def completed(q):
        nonlocal best
        state['completed'] += 1; state['failed'] += q['status'] != 'verified'
        if q['status'] == 'verified':
            validate_score(q['score'], q['proposal'], restrictions)
            if q['loss'] != q['score']['loss']: raise ValueError('Worker loss mismatch')
            if best is None or q['loss'] < best['loss']: best=q
        records.append({k:v for k,v in q.items() if k != 'score'}); persist()
        print(json.dumps(dict(completed=state['completed'], case=q['case_id'], status=q['status'], best_loss=None if best is None else best['loss'])), flush=True)
    def stage(items, label):
        nonlocal attempted
        state['phase']=label; attempted += len(items)
        for item in items: validate_proposal(item, restrictions)
        result = _batch(items, run_case, completed, plan['workers'])
        if any(q['status'] == 'failed' for q in result): raise RuntimeError('Unknown failure: stopped for review')
        if sum(q['status'] == 'rejected_equilibrium' for q in result) > len(result)/2: raise RuntimeError('Housing rejection majority: stopped for review')
        if sum(q['status'] == 'rejected_mass_gate' for q in records) > 1: raise RuntimeError('Repeated mass rejection: stopped for review')
        return result
    persist(); write(out/'heartbeat.json', state); threading.Thread(target=beat, daemon=True).start(); status='completed'; reason='round_cap'; final_ok=False
    pinned = dict(case_id='pinned_seed', parameters=dict(plan['resume_proposal']['parameters']), initial_psi=plan['resume_proposal']['initial_psi'])
    validate_score(seed_score, pinned, restrictions)
    try:
        # The inherited trial is evidence only. Two fresh single repeats must agree with it before it enters search.
        repeats=[]
        for i in (1, 2):
            q=run_case(dict(pinned, case_id=f'seed_exact_repeat_{i:02d}')); completed(q)
            if q['status'] != 'verified' or q.get('receipt_status') != 'verified_scored_candidate' or numeric_signature(q['score']) != numeric_signature(seed_score): raise RuntimeError('Fresh exact seed repetition failed')
            repeats.append(q)
        best=repeats[0]; persist()
        for round_index in range(plan['rounds']):
            if attempted + 30 > plan['maximum_search_cases'] or time.monotonic()+2*WRAPPER_SECONDS > deadline: reason='stage_deadline_or_case_budget'; break
            center=copy.deepcopy(best); probes=derivative_proposals(center, restrictions, round_index)
            if len(probes)>18: raise RuntimeError('Derivative cap exceeded')
            results=stage(probes, f'round_{round_index}_derivatives'); J=jacobian(results, center)
            cond=np.linalg.cond(J); write(out/f'jacobian_round_{round_index}.json', dict(names=list(ALL_NAMES), weighted_jacobian=J.tolist(), rank=int(np.linalg.matrix_rank(J)), condition=float(cond) if np.isfinite(cond) else None, derivative_cases=[dict(case_id=q['case_id'], status=q['status'], column=q['proposal']['column'], step=q['proposal']['step']) for q in results]))
            joint=joint_proposals(center,J,restrictions,round_index)
            if attempted + len(joint) > plan['maximum_search_cases'] or time.monotonic()+WRAPPER_SECONDS > deadline: reason='stage_deadline_or_case_budget'; break
            write(out/f'proposals_round_{round_index}.json',joint); stage(joint,f'round_{round_index}_joint')
        selection=copy.deepcopy(best); q=run_case(dict(case_id='selected_exact_repetitions', parameters={n:parameters(selection['score'])[n] for n in ALL_NAMES}, initial_psi=selection['proposal']['initial_psi'], repetitions=2)); completed(q)
        if q['status'] != 'verified' or q.get('receipt_status') != 'verified_scored_candidate' or numeric_signature(q['score']) != numeric_signature(selection['score']) or not q.get('second_signature_equal', False): raise RuntimeError('Final repetitions failed')
        best=q; final_ok=True
    except Exception as e:
        status='stopped_for_review'; reason='failure'; write(out/'failure.json',dict(error=str(e),type=type(e).__name__))
    finally:
        stop.set(); state.update(status=status,phase='finished'); persist(); write(out/'heartbeat.json',dict(state,elapsed_seconds=time.monotonic()-started))
        summary=dict(metadata(),status=status,stop_reason=reason,elapsed_seconds=time.monotonic()-started,attempted_cases=len(records),search_attempted_cases=attempted,new_case_evaluations=len(records),maximum_stationary_solves=8*sum(r['proposal'].get('repetitions',1) for r in records),best_case=None if best is None else best['case_id'],best_loss=None if best is None else best['loss'],selected_exact_repetitions_verified=final_ok)
        write(out/'summary.json',summary)
    return summary,best


def save_tables(out,best,restrictions):
    if best is None:return
    rows=[]
    for row in best['score']['parameters']:
        q=dict(row); n=q['parameter']
        if n in ALL_NAMES:
            lo,hi=restrictions[n]['lower'],restrictions[n]['upper']
            q.update(status='estimated_capped' if n=='beta_annual' else q.get('status','estimated'), lower=lo, upper=hi,
                     near_bound=min(q['estimate']-lo,hi-q['estimate']) <= .01*(hi-lo))
            if n=='beta_annual':q.update(capped_beta_near_bound_fraction=.01, capped_beta_lower=BETA_LOWER, capped_beta_upper=BETA_UPPER)
        rows.append(q)
    support.write_table(Path(out)/'selected_parameters.csv',rows); support.write_table(Path(out)/'selected_target_fit.csv',best['score']['target_fit'])
    write(Path(out)/'selected_capped_score.json',dict(metadata(),loss=best['score']['loss'],parameters=rows,target_fit=best['score']['target_fit']))
    raw=sorted((Path(best['output'])/'raw').glob('repetition_*'))
    if raw: shutil.copytree(raw[-1]/'standard_diagnostics',Path(out)/'selected_standard_diagnostics')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--plan',type=Path,required=True); args=ap.parse_args(); p=read(args.plan); validate_plan(p); base=args.plan.parent
    for n,pin in p['file_sha256'].items():
        if sha(base/n)!=pin:raise ValueError('Launch file changed: '+n)
    if sha(__file__)!=p['controller_sha256']:raise ValueError('Controller changed')
    manifest_seed=validate_seed_manifest(p)  # validate all five remote receipts before source imports
    source=Path(p['source_root']); old=source/'batches/md_exact_loop'; sys.path.insert(0,str(old/'inputs'))
    import run_scored_candidate as wrapper
    import score_initial as scorer
    smoke=source/'output/md_exact_loop_17370427/summary.json'
    if sha(smoke)!=p['smoke_summary_sha256'] or not read(smoke)['exact_loss_equality']:raise ValueError('Smoke changed')
    wrapper.preflight(old/'contract.json',p['smoke_contract_sha256']); template=read(old/'contract.json'); initial=read(old/'initial_contract.json'); restrictions={r['parameter']:dict(r) for r in read(old/'inputs/working_contract.json')['parameter_restrictions']}; restrictions['beta_annual'].update(lower=BETA_LOWER,upper=BETA_UPPER)
    seed_path=Path(p['resume_score_path']); seed=read(seed_path)
    if sha(seed_path)!=p['resume_score_sha256'] or sha(seed_path.parent.parent/'summary.json')!=p['resume_summary_sha256']:raise ValueError('Resume evidence changed')
    if numeric_signature(seed) != numeric_signature(manifest_seed):raise ValueError('Plan score and seed-manifest score differ')
    receipt=read(seed_path.parent.parent/'summary.json')
    if receipt.get('status')!='verified_scored_candidate' or receipt.get('repetitions',1)!=1:raise ValueError('Resume must be one verified scored trial')
    out=Path(p['output_dir']); out.mkdir(parents=True,exist_ok=False); write(out/'capped_beta_contract.json',metadata())
    def run_case(item):
        validate_proposal(item,restrictions); folder=out/'cases'/item['case_id']; folder.mkdir(parents=True,exist_ok=False); ic=copy.deepcopy(initial); ic.update(case_id=item['case_id'],structural_candidate=item['parameters'],initial_psi=item['initial_psi'],repetitions=item.get('repetitions',1),maximum_GE_solves=8*item.get('repetitions',1),round_id='capped_beta_refinement',scope='Nine free structural coordinates; beta capped at .94,.99') ; ic['run_input_fingerprint']=scorer.fingerprint(item); write(folder/'initial_contract.json',ic)
        c=copy.deepcopy(template);c['case_id']=item['case_id'];
        for k in ('working_objective','scorer','validator'):c[k]['path']=str(old/template[k]['path'])
        for k,e in c['objective_source_files'].items():e['path']=str(old/template['objective_source_files'][k]['path'])
        c['initial_solve_contract']=dict(path=str(folder/'initial_contract.json'),sha256=sha(folder/'initial_contract.json'));write(folder/'run_contract.json',c);dest=folder/'evaluation'
        with (folder/'wrapper.log').open('w') as log: code=subprocess.call([sys.executable,str(old/'inputs/run_scored_candidate.py'),'run','--contract',str(folder/'run_contract.json'),'--contract-sha256',sha(folder/'run_contract.json'),'--output',str(dest)],cwd=source,env=dict(os.environ,PYTHONOPTIMIZE='0',NUMBA_DISABLE_JIT='0',MPLCONFIGDIR=str(folder/'mpl')),stdout=log,stderr=subprocess.STDOUT)
        if code or not (dest/'summary.json').exists():
            detail=read(dest/'raw/failure.json') if (dest/'raw/failure.json').exists() else None
            return dict(case_id=item['case_id'],status=support.failure_status(dest),output=str(dest),proposal=item,failure_detail=detail)
        rec=read(dest/'summary.json'); score=read(dest/'scored_repetition_01/score.json'); validate_score(score,item,restrictions)
        if rec.get('status') != 'verified_scored_candidate' or rec.get('objective_canonical_sha256') != APPROVED_OBJECTIVE or rec.get('loss') != score['loss'] or rec.get('repetitions') != item.get('repetitions',1):raise ValueError('Scored receipt mismatch')
        q=dict(case_id=item['case_id'],status='verified',loss=score['loss'],output=str(dest),proposal=item,score=score,receipt_status=rec.get('status'))
        if item.get('repetitions')==2:
            second=read(dest/'scored_repetition_02/score.json');validate_score(second,item,restrictions);q['second_signature_equal']=bool(rec.get('exact_loss_equality')) and numeric_signature(score)==numeric_signature(second)
        return q
    summary,best=search(p,seed,restrictions,run_case,out);save_tables(out,best,restrictions);print(json.dumps(summary),flush=True)
    if summary['status']!='completed':raise SystemExit(2)

if __name__=='__main__':main()
