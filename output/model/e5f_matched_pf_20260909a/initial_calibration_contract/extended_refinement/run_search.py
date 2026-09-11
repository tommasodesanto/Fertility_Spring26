"""Bounded orchestration around the unchanged, exactly reproduced scored solver."""
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

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.optimize import lsq_linear

HERE=Path(__file__).resolve().parent
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    p=Path(p);q=p.with_suffix(p.suffix+'.tmp')
    q.write_text(json.dumps(v,indent=2,sort_keys=True,allow_nan=False)+'\n');q.replace(p)
def transform(name,value): return math.log(-math.log(value)) if name=='beta_annual' else math.log(value)
def inverse(name,value): return math.exp(-math.exp(value)) if name=='beta_annual' else math.exp(value)
def bounds(names,restrictions,x,radius):
    lo=[];hi=[]
    for i,n in enumerate(names):
        r=restrictions[n]
        lower,upper=(transform(n,r['upper']),transform(n,r['lower'])) if n=='beta_annual' else (-math.inf if r['lower']==0 else transform(n,r['lower']),transform(n,r['upper']))
        lo.append(max(-radius,lower-x[i]));hi.append(min(radius,upper-x[i]))
    return np.array(lo),np.array(hi)
def feasible_steps(lower,upper):
    """At a bound, omit the zero outward step and retain the feasible inward probe."""
    return [(sign,float(step)) for sign,step in [(-1,lower),(1,upper)] if abs(step)>=1e-10]

def derivative_column(samples,center):
    if len(samples)==2:
        samples=sorted(samples,key=lambda a:a[0])
        return (samples[1][1]-samples[0][1])/(samples[1][0]-samples[0][0])
    if len(samples)==1:return (samples[0][1]-center)/samples[0][0]
    raise RuntimeError('No valid derivative observation; no invented derivative')

def residual(score):
    return np.array([float(r['gap'])*math.sqrt(float(r['actual_weight'])) for r in score['target_fit'] if r['scored']])
def parameters(score): return {r['parameter']:float(r['estimate']) for r in score['parameters']}
def batch(items,worker,completed,workers):
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(worker,item):item for item in items}
        results=[]
        for f in cf.as_completed(futures):
            result=f.result();completed(result);results.append(result)
    return results

def failure_status(dest):
    preflight=Path(dest)/'preflight.json';raw=Path(dest)/'raw/failure.json'
    if preflight.exists() and raw.exists():
        f=read(raw)
        if (f.get('error_type')=='RuntimeError' and f.get('phase')=='stationary_equilibrium'
            and f.get('error')=='Initial housing equilibrium failed its unchanged strict gate'):
            return 'rejected_equilibrium'
    return 'failed'


def main():
    p=read(HERE/'plan.json')
    for name,pin in p['file_sha256'].items():
        if sha(HERE/name)!=pin:raise ValueError('Launch file changed: '+name)
    if sha(__file__)!=p['controller_sha256']:raise ValueError('Controller changed')
    source=Path(p['source_root']);old=source/'batches/md_exact_loop'
    sys.path.insert(0,str(old/'inputs'))
    import run_scored_candidate as wrapper
    import score_initial as scorer
    smoke=read(source/'output/md_exact_loop_17370427/summary.json')
    if sha(source/'output/md_exact_loop_17370427/summary.json')!=p['smoke_summary_sha256'] or not smoke['exact_loss_equality']:
        raise ValueError('Exact-loop smoke changed or incomplete')
    if smoke['objective_canonical_sha256']!=wrapper.APPROVED_OBJECTIVE:raise ValueError('Smoke objective mismatch')
    template=read(old/'contract.json');initial=read(old/'initial_contract.json')
    # Full source/observer/target/checkpoint verification before any new solve.
    wrapper.preflight(old/'contract.json',p['smoke_contract_sha256'])
    objective=read(old/'inputs/working_contract.json')
    names=list(scorer.PARAMETERS);restrictions={r['parameter']:r for r in objective['parameter_restrictions']}
    out=HERE/'results';out.mkdir(exist_ok=False)
    started=time.monotonic();search_deadline=started+p['search_seconds']
    state={'phase':'starting','completed':0,'failed':0,'status':'running'}
    stop=threading.Event()
    def heartbeat():
        while not stop.wait(30):write(out/'heartbeat.json',dict(state,elapsed_seconds=time.monotonic()-started))
    threading.Thread(target=heartbeat,daemon=True).start()
    seed_path=Path(p['resume_score_path'])
    if sha(seed_path)!=p['resume_score_sha256']:raise ValueError('Resume score changed')
    seed_score=read(seed_path)
    prior=seed_path.parent.parent/'summary.json'
    if sha(prior)!=p['resume_summary_sha256']:raise ValueError('Selected repetition receipt changed')
    prior_receipt=read(prior)
    if not (prior_receipt['exact_loss_equality'] and prior_receipt['repetitions']==2 and prior_receipt['loss']==seed_score['loss'] and prior_receipt['objective_canonical_sha256']==wrapper.APPROVED_OBJECTIVE):raise ValueError('Prior selected point was not exactly reproduced under this objective')
    best={'case_id':'verified_seed','loss':seed_score['loss'],'score':seed_score,
          'output':str(seed_path.parent.parent),'inherited':True,'proposal':p['resume_proposal']}
    records=[];attempts=0
    def persist():
        write(out/'best_so_far.json',{k:v for k,v in best.items() if k!='score'})
        write(out/'latest_completed.json',dict(state,latest=records[-1] if records else None))
        write(out/'cases.json',records)
    persist();write(out/'heartbeat.json',state)
    def run_case(item):
        case=item['case_id'];folder=out/'cases'/case;folder.mkdir(parents=True,exist_ok=False)
        ic=copy.deepcopy(initial)
        ic.update(case_id=case,structural_candidate=item['parameters'],initial_psi=item['initial_psi'],
                  repetitions=item.get('repetitions',1),maximum_GE_solves=8*item.get('repetitions',1),
                  round_id='three_hour_refinement',scope='Fixed complete early working objective; no production promotion')
        ic['run_input_fingerprint']=scorer.fingerprint(item)
        write(folder/'initial_contract.json',ic)
        c=copy.deepcopy(template);c['case_id']=case
        for key in ('working_objective','scorer','validator'):
            c[key]['path']=str(old/template[key]['path'])
        for key,entry in c['objective_source_files'].items():entry['path']=str(old/template['objective_source_files'][key]['path'])
        c['initial_solve_contract']={'path':str(folder/'initial_contract.json'),'sha256':sha(folder/'initial_contract.json')}
        write(folder/'run_contract.json',c)
        env=dict(os.environ,PYTHONOPTIMIZE='0',NUMBA_DISABLE_JIT='0',MPLCONFIGDIR=str(folder/'mpl'))
        dest=folder/'evaluation'
        with (folder/'wrapper.log').open('w') as log:
            # The unchanged wrapper enforces its own 2100-second cap and reaps model children.
            code=subprocess.call([sys.executable,str(old/'inputs/run_scored_candidate.py'),'run',
                '--contract',str(folder/'run_contract.json'),'--contract-sha256',sha(folder/'run_contract.json'),
                '--output',str(dest)],cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT)
        if code or not (dest/'summary.json').exists():
            return dict(case_id=case,status=failure_status(dest),returncode=code,output=str(dest),proposal=item)
        s=read(dest/'summary.json')
        if s['objective_canonical_sha256']!=wrapper.APPROVED_OBJECTIVE:raise ValueError('Mixed output objective')
        score=read(dest/'scored_repetition_01/score.json')
        if score['loss']!=s['loss'] or not math.isclose(float(residual(score)@residual(score)),s['loss'],rel_tol=1e-12):
            raise ValueError('Output score mismatch')
        return dict(case_id=case,status='verified',loss=s['loss'],output=str(dest),score=score,proposal=item)
    def completed(result):
        nonlocal best,attempts
        attempts+=1;state['completed']+=1
        if result['status']!='verified':state['failed']+=1
        elif result['loss']<best['loss']:best=result
        records.append({k:v for k,v in result.items() if k!='score'});persist()
        print(json.dumps({'completed':attempts,'case':result['case_id'],'status':result['status'],'best_loss':best['loss']}),flush=True)
    def can_search(additional=18):return time.monotonic()+2100<=search_deadline and attempts+additional<=p['maximum_search_cases']
    def run_stage(items,label):
        state['phase']=label
        results=batch(items,run_case,completed,p['workers'])
        # Only the exact, source-verified housing-gate failure is an inadmissible proposal.
        # Unexpected errors still stop; no rejected proposal receives a fabricated score.
        if any(r['status']=='failed' for r in results):raise RuntimeError('Unexpected/code/source/observer failure: stopped for review')
        if sum(r['status']=='rejected_equilibrium' for r in results)>len(results)/2:raise RuntimeError('Majority of proposals fail equilibrium: stopped for review')
        return results
    final_status='completed';verified_repeat=False
    try:
        # The selected seed already has pinned, exact repetitions; do not repeat those again.
        for round_index in range(p['rounds']):
            if not can_search():break
            center=copy.deepcopy(best);v=parameters(center['score']);x=np.array([transform(n,v[n]) for n in names])
            lo,hi=bounds(names,restrictions,x,.02)
            probes=[]
            for j,n in enumerate(names):
                for sign,d in feasible_steps(lo[j],hi[j]):
                    params={k:v[k] for k in names};params[n]=inverse(n,x[j]+d)
                    probes.append(dict(case_id=f'r{round_index}_d{j}_{sign:+d}',parameters=params,
                        initial_psi=v['psi_child'],column=j,step=float(d)))
            results=run_stage(probes,f'round_{round_index}_feasible_derivatives')
            J=np.empty((12,9))
            for j in range(9):
                pair=sorted([r for r in results if r['proposal']['column']==j and r['status']=='verified'],key=lambda r:r['proposal']['step'])
                J[:,j]=derivative_column([(q['proposal']['step'],residual(q['score'])) for q in pair],residual(center['score']))
            write(out/f'jacobian_round_{round_index}.json',dict(center=center['case_id'],names=names,
                weighted_jacobian=J.tolist(),derivative_cases=[dict(case_id=q['case_id'],status=q['status'],column=q['proposal']['column'],step=q['proposal']['step']) for q in results],rank=int(np.linalg.matrix_rank(J)),condition=float(np.linalg.cond(J)) if np.isfinite(np.linalg.cond(J)) else None))
            if not can_search(12):break
            r=residual(center['score']);penalty=np.diag(np.maximum(np.linalg.norm(J,axis=0),1e-8));joint=[]
            for radius in (.1,.2,.35,.5):
                lower,upper=bounds(names,restrictions,x,radius)
                for ridge in (.01,.1,1.):
                    fit=lsq_linear(np.vstack([J,math.sqrt(ridge)*penalty]),np.r_[-r,np.zeros(9)],bounds=(lower,upper),tol=1e-12,max_iter=300)
                    if not fit.success:raise RuntimeError('Linear proposal failed')
                    joint.append(dict(case_id=f'r{round_index}_joint_{len(joint):02d}',parameters={n:inverse(n,x[j]+fit.x[j]) for j,n in enumerate(names)},
                        initial_psi=v['psi_child'],radius=radius,ridge=ridge))
            write(out/f'proposals_round_{round_index}.json',joint)
            run_stage(joint,f'round_{round_index}_12_joint_proposals')
        state['phase']='selected_exact_repetitions';v=parameters(best['score']);selection=copy.deepcopy(best)
        result=run_case(dict(case_id='selected_exact_repetitions',parameters={n:v[n] for n in names},
                            initial_psi=selection['proposal']['initial_psi'] if 'proposal' in selection else initial['initial_psi'],repetitions=2))
        completed(result)
        if result['status']!='verified' or result['loss']!=selection['loss']:raise RuntimeError('Selected candidate does not exactly reproduce')
        verified_repeat=True;best=result
    except Exception as exc:
        final_status='stopped_for_review';write(out/'failure.json',dict(error=str(exc),type=type(exc).__name__))
    finally:
        stop.set();state.update(status=final_status,phase='finished');persist()
        selected=Path(best['output']);score_folder=selected/'scored_repetition_01'
        shutil.copy2(score_folder/'target_fit.csv',out/'selected_target_fit.csv')
        shutil.copy2(score_folder/'parameters.csv',out/'selected_parameters.csv')
        raw_reps=sorted((selected/'raw').glob('repetition_*'))
        shutil.copytree(raw_reps[-1]/'standard_diagnostics',out/'selected_standard_diagnostics')
        summary=dict(status=final_status,elapsed_seconds=time.monotonic()-started,attempted_cases=attempts,
            maximum_stationary_solves=8*sum(r['proposal'].get('repetitions',1) for r in records),best_case=best['case_id'],best_loss=best['loss'],
            initial_loss=seed_score['loss'],selected_exact_repetitions_verified=verified_repeat,
            objective_canonical_sha256=wrapper.APPROVED_OBJECTIVE,production_eligible=False)
        write(out/'summary.json',summary);print(json.dumps(summary),flush=True)
    if final_status!='completed':sys.exit(2)

if __name__=='__main__':main()
