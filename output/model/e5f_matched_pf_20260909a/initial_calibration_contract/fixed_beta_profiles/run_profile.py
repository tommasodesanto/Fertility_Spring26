"""Bounded fixed-beta profiles around the unchanged frozen scored-candidate wrapper.

--plan accepts the existing extended-refinement pins plus fixed_beta, workers<=8,
rounds<=2 and maximum_search_cases<=56 (excludes the new fixed seed and two final
repetitions). --seed-only exercises the exact first solve/normalization/score gate.
--smoke exercises one COMPLETE adaptive round and the final repetition pair.
No raw objective, scientific source or unrestricted scorer metadata is modified.
"""
import argparse
import concurrent.futures as cf
import copy
import csv
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

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
APPROVED_OBJECTIVE = 'c0e266d3a0d430343c469d780d1aedb45fa87f8763c9c938889e0c37daa31de2'
ALL_NAMES = ('beta_annual', 'kappa_fert', 'kappa_fert_continuation', 'chi', 'H0',
             'theta0', 'theta1', 'first_birth_fixed_cost', 'h_P')
FREE_NAMES = tuple(n for n in ALL_NAMES if n != 'beta_annual')
WRAPPER_SECONDS = 2100


def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p, value):
    p = Path(p); q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'); q.replace(p)
def write_table(path, rows):
    if not rows: return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys); writer.writeheader(); writer.writerows(rows)
def parameters(score): return {r['parameter']: float(r['estimate']) for r in score['parameters']}
def residual(score):
    return np.array([float(r['gap']) * math.sqrt(float(r['actual_weight']))
                     for r in score['target_fit'] if r['scored']])
def transform(name, value): return math.log(-math.log(value)) if name == 'beta_annual' else math.log(value)
def inverse(name, value): return math.exp(-math.exp(value)) if name == 'beta_annual' else math.exp(value)


def validate_plan(p):
    required = ('source_root', 'fixed_beta', 'workers', 'rounds', 'search_seconds',
                'maximum_search_cases', 'controller_sha256', 'file_sha256',
                'smoke_summary_sha256', 'smoke_contract_sha256', 'resume_score_path',
                'resume_score_sha256', 'resume_summary_sha256', 'resume_proposal')
    for key in required:
        if key not in p: raise ValueError('Missing plan field: ' + key)
    if p['fixed_beta'] not in (.98, .99): raise ValueError('Only author-approved beta profiles allowed')
    if p['workers'] != 8: raise ValueError('These complete derivative stages require exactly eight workers')
    if type(p['rounds']) is not int or not 0 <= p['rounds'] <= 2: raise ValueError('At most two rounds')
    if type(p['maximum_search_cases']) is not int or not 0 <= p['maximum_search_cases'] <= 56:
        raise ValueError('At most 56 new search cases')
    if not math.isfinite(p['search_seconds']) or p['search_seconds'] <= 0: raise ValueError('Invalid deadline')
    if set(p['resume_proposal']['parameters']) != set(ALL_NAMES): raise ValueError('Incomplete seed proposal')


def profile_metadata(beta):
    return dict(fixed_beta_annual=beta, free_parameter_count=8, free_parameters=list(FREE_NAMES),
                fixed_parameter_restrictions={'beta_annual': dict(value=beta, status='author_fixed_profile_restriction')},
                scored_moment_count=12, separate_normalization_target=2.1,
                objective_canonical_sha256=APPROVED_OBJECTIVE,
                raw_scorer_free_parameter_count=9,
                raw_scorer_metadata_note='Frozen scorer describes the unrestricted nine-coordinate contract; this profile imposes the additional author-fixed beta restriction externally.',
                production_eligible=False, identification_established=False)


def full_parameters(values, beta):
    return {n: beta if n == 'beta_annual' else float(values[n]) for n in ALL_NAMES}


def validate_score(score, proposal, beta):
    if set(proposal['parameters']) != set(ALL_NAMES) or proposal['parameters']['beta_annual'] != beta:
        raise ValueError('Candidate violates complete fixed-beta input contract')
    v = parameters(score)
    if any(v[n] != proposal['parameters'][n] for n in ALL_NAMES):
        raise ValueError('Scored structural parameters differ from exact proposal')
    if score.get('contract_sha256') != APPROVED_OBJECTIVE: raise ValueError('Mixed score objective')
    if score.get('free_parameter_count') != 9: raise ValueError('Frozen raw scorer metadata changed')
    if len(score['target_fit']) != 13 or residual(score).shape != (12,): raise ValueError('Incomplete target system')
    norms = [r for r in score['target_fit'] if not r['scored']]
    if len(norms) != 1 or float(norms[0]['target']) != 2.1 or abs(float(norms[0]['gap'])) > 5e-4:
        raise ValueError('Separate normalization gate failed')
    if not math.isclose(float(residual(score) @ residual(score)), score['loss'], rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError('Output loss mismatch')


def numeric_signature(score):
    # Check every numerical target cell, all parameter metadata and normalization;
    # checkpoint provenance legitimately differs between exact repetitions.
    fits = [{k: v for k, v in row.items() if isinstance(v, (int, float, bool)) or k == 'restriction_id'}
            for row in score['target_fit']]
    normalization = {k: score['normalization'][k] for k in
                     ('target', 'psi_child', 'completed_fertility', 'absolute_gap', 'status')
                     if k in score['normalization']}
    return dict(loss=score['loss'], target_fit=fits, parameters=score['parameters'], normalization=normalization)


def bounds(names, restrictions, x, radius):
    lo = []; hi = []
    for i, n in enumerate(names):
        r = restrictions[n]
        lower = -math.inf if r['lower'] == 0 else transform(n, r['lower'])
        upper = transform(n, r['upper'])
        lo.append(max(-radius, lower - x[i])); hi.append(min(radius, upper - x[i]))
    return np.array(lo), np.array(hi)


def feasible_steps(lower, upper):
    return [(sign, float(step)) for sign, step in [(-1, lower), (1, upper)] if abs(step) >= 1e-10]


def derivative_column(samples, center):
    if len(samples) == 2:
        samples = sorted(samples, key=lambda a: a[0])
        return (samples[1][1] - samples[0][1]) / (samples[1][0] - samples[0][0])
    if len(samples) == 1: return (samples[0][1] - center) / samples[0][0]
    raise RuntimeError('No valid derivative observation; no invented derivative')


def proposals(center, restrictions, beta, round_index):
    v = parameters(center['score']); x = np.array([transform(n, v[n]) for n in FREE_NAMES])
    lo, hi = bounds(FREE_NAMES, restrictions, x, .02); result = []
    for j, n in enumerate(FREE_NAMES):
        for sign, step in feasible_steps(lo[j], hi[j]):
            candidate = full_parameters(v, beta)
            candidate[n] = min(restrictions[n]['upper'], max(restrictions[n]['lower'], inverse(n, x[j] + step)))
            result.append(dict(case_id=f'r{round_index}_d{j}_{sign:+d}', parameters=candidate,
                               initial_psi=v['psi_child'], column=j, step=transform(n, candidate[n]) - x[j]))
    return result


def jacobian(results, center):
    J = np.empty((12, 8))
    for j in range(8):
        samples = [(q['proposal']['step'], residual(q['score'])) for q in results
                   if q['proposal']['column'] == j and q['status'] == 'verified']
        J[:, j] = derivative_column(samples, residual(center['score']))
    return J


def joint_proposals(center, J, restrictions, beta, round_index):
    if J.shape != (12, 8): raise ValueError('Fixed-beta Jacobian must be 12 by 8')
    v = parameters(center['score']); x = np.array([transform(n, v[n]) for n in FREE_NAMES])
    r = residual(center['score']); penalty = np.diag(np.maximum(np.linalg.norm(J, axis=0), 1e-8)); joint = []
    for radius in (.1, .2, .35, .5):
        lower, upper = bounds(FREE_NAMES, restrictions, x, radius)
        for ridge in (.01, .1, 1.):
            fit = lsq_linear(np.vstack([J, math.sqrt(ridge) * penalty]), np.r_[-r, np.zeros(8)],
                             bounds=(lower, upper), tol=1e-12, max_iter=300)
            if not fit.success: raise RuntimeError('Linear proposal failed')
            candidate = full_parameters(v, beta)
            for j, n in enumerate(FREE_NAMES):
                candidate[n] = min(restrictions[n]['upper'], max(restrictions[n]['lower'], inverse(n, x[j] + fit.x[j])))
            joint.append(dict(case_id=f'r{round_index}_joint_{len(joint):02d}', parameters=candidate,
                              initial_psi=v['psi_child'], radius=radius, ridge=ridge))
    return joint


def batch(items, worker, completed, workers):
    # Explicit waves provide a strict <=2-wave stage bound, even after a rejection.
    if len(items) > 2 * workers: raise ValueError('Stage exceeds two worker waves')
    results = []
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        for offset in range(0, len(items), workers):
            futures = [pool.submit(worker, item) for item in items[offset:offset + workers]]
            for future in cf.as_completed(futures):
                result = future.result(); completed(result); results.append(result)
    return results


def failure_status(dest):
    preflight = Path(dest) / 'preflight.json'; raw = Path(dest) / 'raw/failure.json'
    if preflight.exists() and raw.exists():
        f = read(raw)
        if (f.get('error_type') == 'RuntimeError' and f.get('phase') == 'stationary_equilibrium'
                and f.get('error') == 'Initial housing equilibrium failed its unchanged strict gate'):
            return 'rejected_equilibrium'
    return 'failed'


def search(p, seed_score, restrictions, run_case, out, *, seed_only=False, smoke=False):
    """Production loop, injected worker also used by complete deterministic smoke tests."""
    beta = p['fixed_beta']; started = time.monotonic(); deadline = started + p['search_seconds']
    best = None; records = []; state = dict(phase='fixed_beta_seed', completed=0, failed=0, status='running')
    stop = threading.Event(); search_attempts = 0; fixed_seed_loss = None; verified_repeat = False
    def persist():
        write(out / 'best_so_far.json', None if best is None else {k: v for k, v in best.items() if k != 'score'})
        write(out / 'latest_completed.json', dict(state, latest=records[-1] if records else None))
        write(out / 'cases.json', records)
    def heartbeat():
        while not stop.wait(30): write(out / 'heartbeat.json', dict(state, elapsed_seconds=time.monotonic() - started))
    def completed(result):
        nonlocal best
        state['completed'] += 1
        if result['status'] != 'verified': state['failed'] += 1
        else:
            validate_score(result['score'], result['proposal'], beta)
            if result['loss'] != result['score']['loss']: raise ValueError('Worker loss mismatch')
            if best is None or result['loss'] < best['loss']: best = result
        records.append({k: v for k, v in result.items() if k != 'score'}); persist()
        print(json.dumps(dict(completed=state['completed'], case=result['case_id'], status=result['status'],
                              best_loss=None if best is None else best['loss'])), flush=True)
    def run_stage(items, label):
        nonlocal search_attempts
        state['phase'] = label; search_attempts += len(items)
        results = batch(items, run_case, completed, p['workers'])
        if any(r['status'] == 'failed' for r in results): raise RuntimeError('Unexpected/code/source/observer failure: stopped for review')
        if sum(r['status'] == 'rejected_equilibrium' for r in results) > len(results) / 2:
            raise RuntimeError('Majority of proposals fail equilibrium: stopped for review')
        return results
    def can_search(count):
        waves = math.ceil(count / p['workers'])
        return waves <= 2 and time.monotonic() + waves * WRAPPER_SECONDS <= deadline and search_attempts + count <= p['maximum_search_cases']
    persist(); write(out / 'heartbeat.json', state); threading.Thread(target=heartbeat, daemon=True).start()
    final_status = 'completed'; stop_reason = 'round_cap'
    try:
        # The unrestricted selected point is provenance/initialization only. It can
        # never enter best, even if its loss is lower than every fixed-beta point.
        seed = dict(case_id='fixed_beta_seed', parameters=full_parameters(parameters(seed_score), beta),
                    initial_psi=p.get('seed_initial_psi', p['resume_proposal']['initial_psi']))
        fixed = run_case(seed); completed(fixed)
        if fixed['status'] == 'rejected_equilibrium':
            # Predeclared recovery is limited to two H0 alternatives. These are
            # newly normalized fixed-beta solves, never the unrestricted score.
            alternatives = []
            for factor in (.98, 1.02):
                candidate = copy.deepcopy(seed); candidate['case_id'] = f'fixed_beta_seed_H0_{factor:.2f}'
                candidate['parameters']['H0'] *= factor
                if restrictions['H0']['lower'] <= candidate['parameters']['H0'] <= restrictions['H0']['upper']:
                    alternatives.append(candidate)
            state['phase'] = 'fixed_beta_seed_recovery'
            recovered = batch(alternatives, run_case, completed, p['workers'])
            if any(r['status'] == 'failed' for r in recovered): raise RuntimeError('Unexpected fixed-seed recovery failure')
        elif fixed['status'] != 'verified':
            raise RuntimeError('Unexpected fixed-seed failure')
        if best is None: raise RuntimeError('No valid fixed-beta center after bounded recovery; unrestricted seed cannot be selected')
        fixed_seed_loss = best['loss']
        if seed_only:
            final_status = 'seed_smoke_passed'; stop_reason = 'seed_only'
        else:
            for round_index in range(min(p['rounds'], 1) if smoke else p['rounds']):
                center = copy.deepcopy(best); probes = proposals(center, restrictions, beta, round_index)
                if not can_search(len(probes)):
                    stop_reason = 'stage_deadline_or_case_budget'; break
                results = run_stage(probes, f'round_{round_index}_feasible_derivatives')
                J = jacobian(results, center)
                condition = np.linalg.cond(J)
                write(out / f'jacobian_round_{round_index}.json', dict(center=center['case_id'], names=list(FREE_NAMES),
                      fixed_beta=beta, weighted_jacobian=J.tolist(), rank=int(np.linalg.matrix_rank(J)),
                      condition=float(condition) if np.isfinite(condition) else None,
                      derivative_cases=[dict(case_id=q['case_id'], status=q['status'], column=q['proposal']['column'], step=q['proposal']['step']) for q in results]))
                joint = joint_proposals(center, J, restrictions, beta, round_index)
                if not can_search(len(joint)):
                    stop_reason = 'stage_deadline_or_case_budget'; break
                write(out / f'proposals_round_{round_index}.json', joint)
                run_stage(joint, f'round_{round_index}_12_joint_proposals')
            state['phase'] = 'selected_exact_repetitions'; selection = copy.deepcopy(best)
            # Reuse the selected INPUT psi, not its final normalized psi: otherwise
            # the approximate normalization stopping point can change on replay.
            item = dict(case_id='selected_exact_repetitions', parameters=full_parameters(parameters(selection['score']), beta),
                        initial_psi=selection['proposal']['initial_psi'], repetitions=2)
            result = run_case(item); completed(result)
            if result['status'] != 'verified' or numeric_signature(result['score']) != numeric_signature(selection['score']):
                raise RuntimeError('Selected fixed-beta candidate does not exactly reproduce')
            verified_repeat = True; best = result
    except Exception as exc:
        final_status = 'stopped_for_review'; stop_reason = 'failure'
        write(out / 'failure.json', dict(error=str(exc), type=type(exc).__name__))
    finally:
        stop.set(); state.update(status=final_status, phase='finished'); persist()
        write(out / 'heartbeat.json', dict(state, elapsed_seconds=time.monotonic() - started))
        summary = dict(profile_metadata(beta), status=final_status, stop_reason=stop_reason,
                       elapsed_seconds=time.monotonic() - started, attempted_cases=len(records),
                       search_attempted_cases=search_attempts,
                       maximum_stationary_solves=8 * sum(r['proposal'].get('repetitions', 1) for r in records),
                       best_case=None if best is None else best['case_id'], best_loss=None if best is None else best['loss'],
                       unrestricted_seed_loss_for_reference_only=seed_score['loss'], fixed_beta_seed_loss=fixed_seed_loss,
                       selected_exact_repetitions_verified=verified_repeat)
        write(out / 'summary.json', summary)
    return summary, best


def save_profile_tables(out, best, beta):
    if best is None: return
    selected = Path(best['output']); score = best['score']; rows = []
    for original in score['parameters']:
        row = copy.deepcopy(original); name = row['parameter']
        row.update(original_unrestricted_status=row.get('status'), profile_free_parameter=name in FREE_NAMES,
                   profile_fixed_beta_annual=beta, profile_free_parameter_count=8,
                   profile_restriction_status='reoptimized' if name in FREE_NAMES else 'inherited_fixed_or_derived')
        if name == 'beta_annual':
            row.update(status='author_fixed_profile_restriction', profile_restriction_status='author_fixed',
                       profile_fixed_value=beta, profile_lower=beta, profile_upper=beta)
        rows.append(row)
    write_table(out / 'selected_parameters.csv', rows)
    write_table(out / 'selected_target_fit.csv', [dict(r, profile_fixed_beta_annual=beta, profile_free_parameter_count=8)
                                                for r in score['target_fit']])
    write(out / 'selected_profile_score.json', dict(profile_metadata(beta), loss=score['loss'], target_fit=score['target_fit'],
          parameters=rows, raw_unchanged_score_path=str(selected / 'scored_repetition_01/score.json')))
    raw_reps = sorted((selected / 'raw').glob('repetition_*'))
    if raw_reps:
        shutil.copytree(raw_reps[-1] / 'standard_diagnostics', out / 'selected_standard_diagnostics')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=HERE / 'plan.json')
    group = parser.add_mutually_exclusive_group(); group.add_argument('--smoke', action='store_true'); group.add_argument('--seed-only', action='store_true')
    args = parser.parse_args(); p = read(args.plan); validate_plan(p); base = args.plan.resolve().parent
    for name, pin in p['file_sha256'].items():
        if sha(base / name) != pin: raise ValueError('Launch file changed: ' + name)
    if sha(__file__) != p['controller_sha256']: raise ValueError('Controller changed')
    source = Path(p['source_root']); old = source / 'batches/md_exact_loop'
    sys.path.insert(0, str(old / 'inputs'))
    import run_scored_candidate as wrapper
    import score_initial as scorer
    if tuple(scorer.PARAMETERS) != ALL_NAMES or wrapper.APPROVED_OBJECTIVE != APPROVED_OBJECTIVE:
        raise ValueError('Original scorer/wrapper interface changed')
    smoke_path = source / 'output/md_exact_loop_17370427/summary.json'; smoke_receipt = read(smoke_path)
    if sha(smoke_path) != p['smoke_summary_sha256'] or not smoke_receipt['exact_loss_equality']:
        raise ValueError('Original exact-loop smoke changed/incomplete')
    if smoke_receipt['objective_canonical_sha256'] != APPROVED_OBJECTIVE: raise ValueError('Smoke objective mismatch')
    wrapper.preflight(old / 'contract.json', p['smoke_contract_sha256'])
    template = read(old / 'contract.json'); initial = read(old / 'initial_contract.json')
    objective = read(old / 'inputs/working_contract.json')
    restrictions = {r['parameter']: r for r in objective['parameter_restrictions']}
    if not restrictions['beta_annual']['lower'] <= p['fixed_beta'] <= restrictions['beta_annual']['upper']:
        raise ValueError('Fixed beta violates original bounds')
    seed_path = Path(p['resume_score_path'])
    if sha(seed_path) != p['resume_score_sha256']: raise ValueError('Resume score changed')
    seed_score = read(seed_path); prior_path = seed_path.parent.parent / 'summary.json'
    if sha(prior_path) != p['resume_summary_sha256']: raise ValueError('Selected repetition receipt changed')
    prior = read(prior_path)
    if not (prior['exact_loss_equality'] and prior['repetitions'] == 2 and prior['loss'] == seed_score['loss']
            and prior['objective_canonical_sha256'] == APPROVED_OBJECTIVE): raise ValueError('Unverified source seed')
    if full_parameters(parameters(seed_score), parameters(seed_score)['beta_annual']) != p['resume_proposal']['parameters']:
        raise ValueError('Pinned resume proposal differs from seed parameters')
    out = Path(p.get('output_dir', base / ('seed_smoke_results' if args.seed_only else 'smoke_results' if args.smoke else 'results')))
    out.mkdir(parents=True, exist_ok=False); write(out / 'profile_contract.json', profile_metadata(p['fixed_beta']))
    def run_case(item):
        if item['parameters'] != full_parameters(item['parameters'], p['fixed_beta']): raise ValueError('Nonfixed proposal')
        case = item['case_id']; folder = out / 'cases' / case; folder.mkdir(parents=True, exist_ok=False)
        ic = copy.deepcopy(initial)
        ic.update(case_id=case, structural_candidate=item['parameters'], initial_psi=item['initial_psi'],
                  repetitions=item.get('repetitions', 1), maximum_GE_solves=8 * item.get('repetitions', 1),
                  round_id='fixed_beta_profile', scope='Author-fixed annual beta; complete unchanged 12-row objective and separate normalization')
        ic['run_input_fingerprint'] = scorer.fingerprint(item); write(folder / 'initial_contract.json', ic)
        c = copy.deepcopy(template); c['case_id'] = case
        for key in ('working_objective', 'scorer', 'validator'): c[key]['path'] = str(old / template[key]['path'])
        for key, entry in c['objective_source_files'].items(): entry['path'] = str(old / template['objective_source_files'][key]['path'])
        c['initial_solve_contract'] = dict(path=str(folder / 'initial_contract.json'), sha256=sha(folder / 'initial_contract.json'))
        write(folder / 'run_contract.json', c)
        env = dict(os.environ, PYTHONOPTIMIZE='0', NUMBA_DISABLE_JIT='0', MPLCONFIGDIR=str(folder / 'mpl')); dest = folder / 'evaluation'
        with (folder / 'wrapper.log').open('w') as log:
            code = subprocess.call([sys.executable, str(old / 'inputs/run_scored_candidate.py'), 'run', '--contract',
                   str(folder / 'run_contract.json'), '--contract-sha256', sha(folder / 'run_contract.json'), '--output', str(dest)],
                   cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
        if code or not (dest / 'summary.json').exists():
            return dict(case_id=case, status=failure_status(dest), returncode=code, output=str(dest), proposal=item)
        receipt = read(dest / 'summary.json')
        if receipt['objective_canonical_sha256'] != APPROVED_OBJECTIVE: raise ValueError('Mixed output objective')
        score = read(dest / 'scored_repetition_01/score.json'); validate_score(score, item, p['fixed_beta'])
        if score['loss'] != receipt['loss']: raise ValueError('Receipt loss mismatch')
        if item.get('repetitions', 1) == 2:
            second = read(dest / 'scored_repetition_02/score.json'); validate_score(second, item, p['fixed_beta'])
            if not (receipt['repetitions'] == 2 and receipt['exact_loss_equality'] and numeric_signature(score) == numeric_signature(second)):
                raise ValueError('Incomplete/nonexact final repetition pair')
        return dict(case_id=case, status='verified', loss=score['loss'], output=str(dest), score=score, proposal=item)
    summary, best = search(p, seed_score, restrictions, run_case, out, seed_only=args.seed_only, smoke=args.smoke)
    save_profile_tables(out, best, p['fixed_beta']); print(json.dumps(summary), flush=True)
    if summary['status'] == 'stopped_for_review': sys.exit(2)


if __name__ == '__main__': main()
