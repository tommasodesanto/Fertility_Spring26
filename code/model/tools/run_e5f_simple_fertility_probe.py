"""Bounded, unchanged-grid fixed-price lifecycle probe; never a calibration.

No figures: author explicitly requested no unsolicited illustrations.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ.setdefault(key,'1')
import argparse,copy,json,time,sys,threading,hashlib,faulthandler
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/model'),str(ROOT/'code/model/tools')]
import numpy as np
import run_e5f_independent_numerical_audit as audit
from intergen_eqscale_seq_optimized import solver as model

def main():
    audit.transition.configure_sequential_model()
    audit.calendar.apply_fertility = audit.transition.apply_sequential_fertility
    audit.calendar.advance_calendar_distribution = audit.transition.advance_sequential_calendar_distribution
    audit.calendar.distribution_rows = audit.transition.independent_child_distribution_rows
    ap=argparse.ArgumentParser();ap.add_argument('--checkpoint',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--seconds',type=int,default=600)
    ap.add_argument('--skip-reference',action='store_true',help='Local exact-reference receipt must already exist')
    ap.add_argument('--cases',nargs='+',default=['old_fixed_price','sequential_exhaustive_fixed_price','fertility_nest_fixed_price'])
    ap.add_argument('--local-reference',type=Path,help='Pristine parent arrays generated in this same runtime')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();state={'phase':'loading','completed':[]}
    def save(name,obj):audit.save_json(args.output/name,obj)
    def heartbeat():
        while True:
            elapsed=time.monotonic()-start
            save('heartbeat.json',dict(state,elapsed_seconds=elapsed))
            if elapsed>args.seconds:
                save('timeout.json',dict(state,elapsed_seconds=elapsed));os._exit(124)
            time.sleep(30)
    threading.Thread(target=heartbeat,daemon=True).start()
    packet=audit.load_checkpoint(args.checkpoint);P0=packet['parameters'];bg=packet['b_grid'];price=packet['evaluation'].policy.price
    save('contract.json',dict(status='diagnostic_only_no_calibration_no_market_clearing',checkpoint_sha256=audit.digest(args.checkpoint),
        local_reference_sha256=audit.digest(args.local_reference) if args.local_reference else None,
        price=price,grid=list(packet['evaluation'].policy.V.shape),time_budget_seconds=args.seconds,
        cases=['default_off_reference','old_fixed_price','sequential_exhaustive_fixed_price','fertility_nest_fixed_price'],
        shocks='simple fertility GEV over complete contingent housing plans',
        housing_scale=P0.tenure_choice_kappa,first_fertility_scale=P0.kappa_fert,later_fertility_scale=P0.kappa_fert_continuation,
        plan_menu='wait products and Cartesian success/failure products; one inner scale; no committed tenure',
        figures='suppressed at explicit author request',source_hashes={str(p.relative_to(ROOT)):audit.digest(p) for p in (ROOT/'code/model/intergen_eqscale_seq_optimized').glob('*.py')}))
    if not args.skip_reference:
        # Baseline equivalence check uses the checkpoint's exact continuation.
        P=copy.deepcopy(P0);P.joint_nested_choice=False;P.two_shock_choice=False;P.fertility_nest_choice=False;P.exhaustive_saving_control=False
        sd=model.precompute_shared(P,bg);ref=packet['evaluation'].policy
        state['phase']='default_off_reference';print(state['phase'],flush=True)
        expected_second=P._fert2_probs.copy()
        result=model.solve_bellman_full_markov_income(P.user_cost_rate*price,price,P,bg,sd,continuation_V=ref.V)
        fields=('V','c_pol','hR_pol','bp_pol','tenure_choice','tenure_probs','loc_probs','fert_probs','fert_value')
        local=np.load(args.local_reference) if args.local_reference else None
        differences={}
        for field,value in zip(fields,result[:9]):
            differences[field]=float(np.max(abs(value-getattr(ref,field))))
            expected=local[field] if local is not None else getattr(ref,field)
            if not np.array_equal(value,expected):raise RuntimeError('Baseline changed: '+field)
        expected=local['fert2_probs'] if local is not None else expected_second
        if not np.array_equal(P._fert2_probs,expected):raise RuntimeError('Baseline continuation fertility changed')
        save('baseline_reference.json',dict(status='exact_same_runtime_parent' if local is not None else 'exact_checkpoint',checkpoint_max_abs_differences=differences,arrays=list(fields)+['fert2_probs'],elapsed_seconds=time.monotonic()-start))
        state['completed'].append('default_off_reference')
    cases=[]
    faulthandler.dump_traceback_later(180,repeat=True)
    modes={'old_fixed_price':(False,False),'sequential_exhaustive_fixed_price':(False,True),'fertility_nest_fixed_price':(True,False)}
    for name in args.cases:
        active,matched=modes[name]
        state['phase']=name;print(name,flush=True);begin=time.monotonic()
        P=copy.deepcopy(P0);P.joint_nested_choice=active;P.two_shock_choice=False;P.fertility_nest_choice=active;P.exhaustive_saving_control=matched
        sd=model.precompute_shared(P,bg)
        sol=model.solve_markov_income_at_prices(price,P,bg,SD=sd)
        row=dict(case=name,elapsed_seconds=time.monotonic()-begin,mass=sol.total_mass,births=sol.total_births_kfe,
            ownership=sol.own_rate,mean_parity=sol.mean_parity,market_residual=sol.best_max_abs_rel_excess,timings=sol.timings)
        if active or matched:
            audit.transition.configure_sequential_model()
            state['phase']=name+'_policy_bundle';print(state['phase'],flush=True)
            pol=audit.calendar.policy_from_solution(sol,price,P,bg,sd)
            state['phase']=name+'_reconstruction';print(state['phase'],flush=True)
            pre,recon=audit.calendar.reconstruct_stationary_pre_fertility(sol,pol,P,bg,sd)
            state['phase']=name+'_evaluate';print(state['phase'],flush=True)
            ev=audit.calendar.evaluate_period(price,pre,P,bg,sd,audit.calendar.SolveCounter(),supplied_policy=pol)
            row['stationary_current_l1']=float(np.abs(ev.g_current-sol.g).sum())
            row['stationary_birth_gap']=float(ev.births-sol.total_births_kfe)
            row['reconstruction']=recon
            if row['stationary_current_l1']>2e-8 or abs(row['stationary_birth_gap'])>2e-10:raise RuntimeError(row)
            modified=dict(packet,parameters=P,shared=sd,evaluation=ev)
            folder=args.output/name;folder.mkdir(exist_ok=True)
            state['phase']=name+'_audits';print(state['phase'],flush=True)
            budget=audit.budget_audit(modified,folder);policy=audit.policy_array_audit(modified,folder)
            if budget['budget_excess_mass']>2e-10 or policy['occupied_negative_steps']:
                raise RuntimeError(dict(budget=budget,policy=policy))
        cases.append(row);state['completed'].append(name)
        save('latest_completed.json',row);save('summary.json',dict(status='in_progress',cases=cases))
        np.savez_compressed(args.output/(name+'_policies.npz'),V=sol.V,bp_pol=sol.bp_pol,g=sol.g)
        print(json.dumps(row),flush=True)
    state['phase']='complete';save('summary.json',dict(status='complete',cases=cases,elapsed_seconds=time.monotonic()-start))

if __name__=='__main__':
    try:
        main()
    except Exception as error:
        if '--output' in sys.argv:
            out=Path(sys.argv[sys.argv.index('--output')+1]);out.mkdir(parents=True,exist_ok=True)
            audit.save_json(out/'failure.json',dict(error=repr(error)))
        raise
