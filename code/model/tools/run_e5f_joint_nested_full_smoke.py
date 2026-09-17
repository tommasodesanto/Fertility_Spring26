"""Bounded lifecycle/operator smoke at inherited prices; no calibration claim."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ.setdefault(name,'1')
import argparse,copy,json,time,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import numpy as np
import run_e5f_independent_numerical_audit as audit
from intergen_eqscale_seq_optimized import solver as model

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--checkpoint',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--scale',type=float,default=2.);ap.add_argument('--lambda',dest='lam',type=float,default=.8);ap.add_argument('--reference',action='store_true')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    packet=audit.load_checkpoint(args.checkpoint);P=copy.deepcopy(packet['parameters']);bg=packet['b_grid'];price=packet['evaluation'].policy.price
    start=time.time()
    if args.reference:
        P.joint_nested_choice=False
    else:
        P.joint_nested_choice=True;P.tenure_choice_kappa=args.scale;P.joint_nest_lambda=args.lam
    print('starting_full_lifecycle',flush=True)
    sd=model.precompute_shared(P,bg)
    if args.reference:
        reference=packet['evaluation'].policy
        expected_second=np.array(P._fert2_probs,copy=True)
        result=model.solve_bellman_full_markov_income(P.user_cost_rate*price,price,P,bg,sd,continuation_V=reference.V)
        fields=('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value')
        checked=[]
        for field,value in zip(fields,result[:9]):
            if not np.array_equal(value,getattr(reference,field)):
                raise RuntimeError(f'Default-off baseline does not reproduce {field}')
            checked.append(field)
        if not np.array_equal(P._fert2_probs,expected_second):
            raise RuntimeError('Default-off continuation births do not reproduce')
        receipt=dict(status='exact',arrays=checked+['fert2_probs'],elapsed_seconds=time.time()-start)
        (args.output/'baseline_reference.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps(receipt),flush=True)
        return
    sol=model.solve_markov_income_at_prices(price,P,bg,SD=sd)
    result=dict(elapsed_seconds=time.time()-start,reference=args.reference,scale=P.tenure_choice_kappa,
        lambda_=getattr(P,'joint_nest_lambda',None),mass=sol.total_mass,births=sol.total_births_kfe,
        ownership=sol.own_rate,mean_parity=sol.mean_parity,market_residual=sol.best_max_abs_rel_excess,
        timings=sol.timings)
    if not args.reference:
        audit.transition.configure_sequential_model()
        pol=audit.calendar.policy_from_solution(sol,price,P,bg,sd)
        pre,recon=audit.calendar.reconstruct_stationary_pre_fertility(sol,pol,P,bg,sd)
        ev=audit.calendar.evaluate_period(price,pre,P,bg,sd,audit.calendar.SolveCounter(),supplied_policy=pol)
        result['reconstruction']=recon
        result['stationary_current_l1']=float(np.abs(ev.g_current-sol.g).sum())
        result['stationary_birth_gap']=float(ev.births-sol.total_births_kfe)
        if result['stationary_current_l1']>2e-8 or abs(result['stationary_birth_gap'])>2e-10:raise RuntimeError(result)
        modified=dict(packet,parameters=P,shared=sd,evaluation=ev)
        audit.standard_diagnostics(modified,args.output,validate_production_young=False)
        audit.budget_audit(modified,args.output)
        audit.policy_array_audit(modified,args.output)
    (args.output/'smoke.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
