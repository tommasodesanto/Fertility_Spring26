"""Canonical cached-policy budget/purchase and aggregate audit; no Bellman call."""
from pathlib import Path
from types import SimpleNamespace
import argparse,json,sys,time,hashlib,traceback,os,signal,resource
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
from model_policy_tools import aggregate_solution

def encode(v):
    if isinstance(v,np.ndarray): return v.tolist()
    if isinstance(v,np.generic): return v.item()
    return str(v)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('case',type=Path);ap.add_argument('out',type=Path);args=ap.parse_args()
    args.case=args.case.resolve();out=args.out.resolve();out.mkdir(exist_ok=False)
    if not out.is_relative_to(HERE):raise RuntimeError('Audit output must remain in owned audit directory')
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    prep=ROOT/'output/model/credit_mechanism_20261004/diagnostics'/('phi095_preparation.json' if args.case.name.startswith('phi') else 'price110_preparation.json')
    contract=json.loads(prep.read_text());snapshot=ROOT/'output/model/credit_mechanism_20261004/production_source_snapshot.json'
    source_mismatches=[k for k,v in contract['source_sha256'].items() if sha(ROOT/k)!=v]
    if source_mismatches:raise RuntimeError('Reached-source identity differs: '+str(source_mismatches))
    status=json.loads((args.case/'status.json').read_text())
    assert sha(args.case/'solution_arrays.npz')==status['arrays_sha256']
    assert sha(args.case/'executed_P.json')==status['parameters_sha256']
    def alarm(*_):raise TimeoutError('300-second cached native audit budget')
    signal.signal(signal.SIGALRM,alarm);signal.setitimer(signal.ITIMER_REAL,300)
    started=time.monotonic()
    with np.load(args.case/'solution_arrays.npz',allow_pickle=False) as z:
        sol=SimpleNamespace(**{k:z[k].copy() for k in z.files if not k.startswith('shared.')})
        sd=SimpleNamespace(**{k.removeprefix('shared.'):z[k].copy() for k in z.files if k.startswith('shared.')})
    p=json.loads((args.case/'executed_P.json').read_text());P=SimpleNamespace(**{k:np.asarray(v) if isinstance(v,list) else v for k,v in p.items()})
    result={'case':str(args.case),'audit_out':str(out),'source_contract':str(prep),'source_contract_sha256':sha(prep),'production_source_snapshot':str(snapshot),'production_source_snapshot_sha256':sha(snapshot),'source_pin_count':len(contract['source_sha256']),'arrays_sha256':status['arrays_sha256'],'parameters_sha256':status['parameters_sha256'],'threads':{k:os.environ.get(k) for k in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}}
    audit_stage='initializing'
    try:
        from production.reporting import build_context
        ctx=build_context(P,sol.b_grid,out,price_start=float(sol.p_eq[0]),deadline=time.time()+90,max_lifecycle=0,closure='fixed_h0')
        audit_stage='computed_checks'
        rt=ctx['prepared'].rt;cal=rt['primitive'].pf.calendar
        policy=cal.policy_from_solution(sol,sol.p_eq,P,sol.b_grid,sd)
        pre,recon=cal.reconstruct_stationary_pre_fertility(sol,policy,P,sol.b_grid,sd)
        supply=cal.HousingSupplyRule('static-elastic',float(sol.p_eq[0]),float(P.H0[0]*(P.user_cost_rate*sol.p_eq[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
        counter=cal.SolveCounter()
        ev=cal.evaluate_period(sol.p_eq,pre,P,sol.b_grid,sd,counter,supply_rule=supply,supplied_policy=policy)
        result['counter']=vars(counter)
        if counter.total!=0:raise RuntimeError('Cached audit unexpectedly performed lifecycle solve')
        result['reconstruction']=recon
        result['canonical_dated_budget']=rt['primitive'].dated_budget(ev,P,sd,sol.b_grid,float(P.user_cost_rate*sol.p_eq[0]))
        result['canonical_purchase']=rt['accounting'].audit_purchase_accounting(ev,P,sd,sol.b_grid,rt['model'])
        # Exact unchanged canonical gate thresholds from fixed_price_v1.gates.
        gate=sys.modules['e5f_evening_calibration_runtime'].require_abs_gate
        gate(result['canonical_dated_budget']['budget_excess_mass'],2e-10,'Household budget')
        gate(result['canonical_purchase']['maximum_occupied_transaction_wealth_error'],1e-9,'Transaction wealth')
        for k,v in result['canonical_purchase'].items():
            if k.endswith('violation_mass') or k in ('transaction_outside_grid_mass','negative_estate_exposure_mass','saving_outside_grid_mass'):gate(v,2e-10,k)
        result['canonical_budget_purchase_status']='passed'
    except Exception as exc:
        result['canonical_budget_purchase_status']='unexecuted_initialization_blocked' if audit_stage=='initializing' else 'computed_gate_failed';result['canonical_budget_purchase_error']=str(exc)
        result['canonical_budget_purchase_missing_path']=getattr(exc,'filename',None)
        result['exception_type']=type(exc).__name__
        result['traceback']=traceback.format_exc()
    result['elapsed_seconds']=time.monotonic()-started
    result['peak_rss']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['source_mismatches_after']=[k for k,v in contract['source_sha256'].items() if sha(ROOT/k)!=v]
    signal.setitimer(signal.ITIMER_REAL,0)
    (out/'receipt.json').write_text(json.dumps(result,indent=2,default=encode)+'\n')
    print(json.dumps({k:v for k,v in result.items() if 'status' in k or 'error' in k},indent=2))

if __name__=='__main__':main()
