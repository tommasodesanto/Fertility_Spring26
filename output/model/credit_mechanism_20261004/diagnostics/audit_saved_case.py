"""Canonical cached-policy budget/purchase and aggregate audit; no Bellman call."""
from pathlib import Path
from types import SimpleNamespace
import argparse,json,sys,time
import numpy as np
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
from model_policy_tools import aggregate_solution

def encode(v):
    if isinstance(v,np.ndarray): return v.tolist()
    if isinstance(v,np.generic): return v.item()
    return str(v)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('case',type=Path);ap.add_argument('--aggregate-only',action='store_true');args=ap.parse_args();out=args.case/'cached_audits';out.mkdir(exist_ok=True)
    with np.load(args.case/'solution_arrays.npz',allow_pickle=False) as z:
        sol=SimpleNamespace(**{k:z[k].copy() for k in z.files if not k.startswith('shared.')})
        sd=SimpleNamespace(**{k.removeprefix('shared.'):z[k].copy() for k in z.files if k.startswith('shared.')})
    p=json.loads((args.case/'executed_P.json').read_text());P=SimpleNamespace(**{k:np.asarray(v) if isinstance(v,list) else v for k,v in p.items()})
    result={}
    for name in ('g','g_beginning_distribution','g_stay_distribution'):
        a=getattr(sol,name);result[name]=dict(mass=float(a.sum()),minimum=float(a.min()),finite=bool(np.isfinite(a).all()))
    result['probabilities']={}
    for name in ('fert_probs','fert2_probs','tenure_probs','loc_probs'):
        a=getattr(sol,name);result['probabilities'][name]=dict(minimum=float(a.min()),maximum=float(a.max()),finite=bool(np.isfinite(a).all()))
    try:
        aggregate=aggregate_solution(sol,houses=P.H_own,age_start=int(P.age_start),period_years=int(P.period_years))
        result['aggregate_policy_status']='passed';result['aggregate_policy_overall']=aggregate['overall']
    except Exception as exc:result['aggregate_policy_status']='failed';result['aggregate_policy_error']=repr(exc)
    if args.aggregate_only:
        result['canonical_budget_purchase_status']='unexecuted'
        result['canonical_budget_purchase_reason']='Reporting initialization blocked in the .95 audit; no broader repair or repeat requested.'
        (out/'receipt.json').write_text(json.dumps(result,indent=2,default=encode)+'\n')
        print(json.dumps({k:v for k,v in result.items() if 'status' in k or 'error' in k},indent=2));return
    audit_stage='initializing'
    try:
        from production.reporting import build_context
        ctx=build_context(P,sol.b_grid,out,price_start=float(sol.p_eq[0]),deadline=time.time()+90,max_lifecycle=1,closure='fixed_h0')
        audit_stage='computed_checks'
        rt=ctx['prepared'].rt;cal=rt['primitive'].pf.calendar
        policy=cal.policy_from_solution(sol,sol.p_eq,P,sol.b_grid,sd)
        pre,recon=cal.reconstruct_stationary_pre_fertility(sol,policy,P,sol.b_grid,sd)
        supply=cal.HousingSupplyRule('static-elastic',float(sol.p_eq[0]),float(P.H0[0]*(P.user_cost_rate*sol.p_eq[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
        counter=cal.SolveCounter()
        ev=cal.evaluate_period(sol.p_eq,pre,P,sol.b_grid,sd,counter,supply_rule=supply,supplied_policy=policy)
        result['counter']=vars(counter)
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
    (out/'receipt.json').write_text(json.dumps(result,indent=2,default=encode)+'\n')
    print(json.dumps({k:v for k,v in result.items() if 'status' in k or 'error' in k},indent=2))

if __name__=='__main__':main()
