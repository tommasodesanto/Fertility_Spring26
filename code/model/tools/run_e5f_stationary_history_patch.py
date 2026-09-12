"""Fit diagnostic stationary fertility points; never reset the carried history."""
from pathlib import Path
import argparse,copy,csv,gzip,hashlib,json,os,pickle,sys,threading,time
import numpy as np

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix('.tmp')
    q.write_text(json.dumps(d,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x),indent=2)+'\n');q.replace(p)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True)
    ap.add_argument('--arm',type=int,default=0);ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    p=read(args.plan);root=Path(p['source_root']);batch=Path(p['batch'])
    out=batch/('smoke' if args.smoke else f'fit_{args.arm}');out.mkdir(exist_ok=False)
    start=time.monotonic();deadline=start+(900 if args.smoke else 1800)
    state=dict(phase='preflight');stop=threading.Event()
    def heartbeat():
        while not stop.wait(60):
            save(out/'heartbeat.json',dict(state,elapsed=time.monotonic()-start))
            if time.monotonic()>deadline:save(out/'failure.json',dict(error='Explicit stationary-patch budget exhausted',**state));os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    for path,pin in p['file_sha256'].items():
        if sha(path)!=pin:raise ValueError('Changed pinned source/input: '+path)
    sys.path[:0]=[str(root/'code/model/tools'),str(root/'code/model')]
    import run_e5f_candidate_terminal as td
    import run_e5f_matched_pf_smoke as primitive
    import run_e5f_transition_calibration as fertility
    from e5f_stationary_paygo import solve_balanced_initial_equilibrium,certify_initial_pension
    from e5f_parenthood_utility import validate_parenthood_utility
    import run_e5f_independent_numerical_audit as diagnostic
    c=p['terminal_template'];td.validate_contract(c);td.verify_sources(c)
    for key in ('initial_checkpoint','initial_summary','initial_contract'):td.verify(c[key]['path'],c[key]['sha256'])
    score=read(p['initial_score_path'])
    if score['contract_sha256']!=p['target_fingerprint']:raise ValueError('Changed initial target/weight contract')
    with gzip.open(c['initial_checkpoint']['path'],'rb') as f:seed=pickle.load(f)
    initial=seed['parameters'];validate_parenthood_utility(initial)
    if initial.joint_nested_choice or not initial.exhaustive_saving_control:raise ValueError('Wrong choice/optimizer specification')
    grid=seed['b_grid'];supply=seed['supply_rule'];psi0=float(initial.psi_child)
    chain,model=primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility=primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution=primitive.pf.transition.advance_sequential_calendar_distribution
    empirical=list(csv.DictReader(Path(p['empirical_blocks']).open()))
    target=float(empirical[args.arm]['period_tfr_arithmetic_mean']);records=[]
    def evaluate(psi,label):
        if time.monotonic()+60>deadline:raise TimeoutError('Insufficient remaining case budget')
        folder=out/label;folder.mkdir();state.update(phase='stationary_solve',psi=psi,case=label)
        Q=copy.deepcopy(initial);Q.psi_child=float(psi);tick=time.monotonic()
        sol,Q,price,fiscal=solve_balanced_initial_equilibrium(model=model,parameters=Q,b_grid=grid,
            initial_prices=seed['solution'].p_eq,payroll_tax=.179,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        shared=model.precompute_shared(Q,grid);Q._fert2_probs=np.asarray(sol.fert2_probs).copy()
        policy=primitive.pf.calendar.policy_from_solution(sol,price,Q,grid,shared)
        pre,reconstruction=primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol,policy,Q,grid,shared)
        gates=primitive.pf.transition.operator_gates(sol,policy,pre,Q,grid,shared);gates.update(reconstruction)
        for key in ('stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1','mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):
            if not np.isfinite(gates[key]) or abs(gates[key])>5e-9:raise RuntimeError('Stationary operator gate: '+key)
        if abs(gates['zero_entry_mass_accounting_residual'])>2e-8 or gates['stationary_feasibility_projection_mass']>1e-6:raise RuntimeError('Stationary population/feasibility gate')
        e=primitive.pf.calendar.evaluate_period(price,pre,Q,grid,shared,primitive.pf.calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
        if not np.isfinite(e.relative_market_residual) or abs(e.relative_market_residual)>2e-4:raise RuntimeError('Inherited supply does not clear')
        fiscal=certify_initial_pension(e.g_current,Q,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        budget=primitive.dated_budget(e,Q,shared,grid,float(Q.user_cost_rate*price[0]))
        measurements=fertility.period_fertility_diagnostics(e,Q)
        snapshot=dict(parameters=Q,b_grid=grid,evaluation=e,shared=shared,supply_rule=supply,solution=sol,stationary_g_pre=pre)
        checkpoint=folder/'state.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as stream:pickle.dump(snapshot,stream,protocol=5)
        r=dict(psi=psi,model=measurements['period_tfr_topcode_adjusted'],target=target,price=float(price[0]),
            seconds=time.monotonic()-tick,checkpoint=str(checkpoint),checkpoint_sha256=sha(checkpoint),
            market=e.relative_market_residual,fiscal=fiscal,operator=gates,budget=budget,
            all_stationary_moments=chain.extract_moments(sol,Q),period_fertility=measurements,
            stationary_demographic_composition_fixed=True,historical_state_carried=False,production_eligible=False)
        r['gap']=r['model']-target;records.append(r);save(folder/'receipt.json',r)
        save(out/'latest_completed.json',r);save(out/'best_so_far.json',min(records,key=lambda q:abs(q['gap'])))
        return r,snapshot
    try:
        if args.smoke:
            first,_=evaluate(psi0,'seed_replay')
            old=fertility.period_fertility_diagnostics(seed['evaluation'],initial)['period_tfr_topcode_adjusted']
            if abs(first['model']-old)>1e-6 or abs(first['price']-float(seed['solution'].p_eq[0]))>1e-6:raise RuntimeError('Pinned seed reproduction failed')
            second,snapshot=evaluate(psi0-.02,'lower_preference')
            diagnostic.standard_diagnostics(snapshot,out/'graphs',validate_production_young=False)
            if len(list((out/'graphs/standard_diagnostics').glob('*.png')))!=17:raise RuntimeError('Missing diagnostic graphs')
            save(out/'summary.json',dict(status='passed_exact_two_point_loop',points=[first,second],seed_period_fertility=old,seconds=time.monotonic()-start))
        else:
            smoke=read(batch/'smoke/summary.json')
            if smoke['status']!='passed_exact_two_point_loop':raise ValueError('Exact loop smoke must pass first')
            from run_e5f_successive_surprises_overnight import next_psi
            trials=[(r['psi'],r['model']-target) for r in smoke['points']]
            for attempt in range(8):
                psi=next_psi(trials,psi0,-.02,(psi0-.15,psi0+.02))
                result,snapshot=evaluate(psi,f'case_{attempt:02d}');trials.append((psi,result['gap']))
                if abs(result['gap'])<=.005:
                    diagnostic.standard_diagnostics(snapshot,out/'graphs',validate_production_young=False)
                    save(out/'summary.json',dict(status='stationary_fertility_target_matched',decision_year=int(empirical[args.arm]['decision_year']),selected=result,
                        carry_forward_still_required=True,post2023_pf_still_required=True));break
            else:raise RuntimeError('Stationary fit exhausted eight candidates without tolerance')
    except Exception as exc:
        save(out/'failure.json',dict(error=str(exc),error_type=type(exc).__name__,**state));raise
    finally:stop.set()

if __name__=='__main__':main()
