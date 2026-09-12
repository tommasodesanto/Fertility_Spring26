"""Conditional post-history policy roots with separate PAYGO and rebate budgets."""
from pathlib import Path
import argparse,copy,gzip,json,os,pickle,subprocess,sys,time
from types import SimpleNamespace as NS
import numpy as np


def policy_residual(demand,supply,accounts,revenue,transfer_outlays,rebated):
    if not np.isfinite(supply) or supply<=0:raise ValueError('Positive supply required')
    r=[(demand-supply)/supply,200*(accounts['payroll_tax_revenue']-accounts['pension_outlays'])/max(abs(accounts['payroll_tax_revenue']),abs(accounts['pension_outlays']),1e-12)]
    if rebated:r.append(200*(revenue-transfer_outlays)/max(abs(revenue),abs(transfer_outlays),1e-12))
    return r


def run_policies(*,inherited,psi,old,demographics,plan,out,deadline):
    from concurrent.futures import ThreadPoolExecutor
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    inputs=out/'inputs.pkl.gz'
    with gzip.open(inputs,'wb',compresslevel=1) as f:pickle.dump(dict(inherited=inherited,psi=psi,old=old,demographics=demographics,plan=plan),f,protocol=5)
    def run(case):
        with (out/(case+'.log')).open('w') as log:
            try:
                result=subprocess.run([sys.executable,'-B',__file__,'--inputs',str(inputs),'--case',case,'--output',str(out/case),'--seconds',str(max(60,int(deadline-time.monotonic()-60)))],stdout=log,stderr=subprocess.STDOUT,timeout=max(60,deadline-time.monotonic()))
                return dict(case=case,exit=result.returncode)
            except subprocess.TimeoutExpired:return dict(case=case,exit=124)
    cases=['baseline','equal-rebate-1pct','equal-rebate-2pct']
    with ThreadPoolExecutor(max_workers=3) as pool:statuses=list(pool.map(run,cases))
    (out/'cases.json').write_text(json.dumps(statuses,indent=2)+'\n')


def execute(packet,case,out,seconds):
    from run_e5f_successive_surprises_overnight import save
    plan=packet['plan'];root=Path(plan['source_root']);sys.path[:0]=[str(root/'code/model/tools'),str(root/'code/model')]
    import e5f_balanced_history as baseline
    import e5f_balanced_terminal as balanced
    from e5f_social_security import bind_social_security_income,fiscal_accounts
    from e5f_matched_pf_endpoint import evaluate_endpoint,EndpointControls
    from e5f_matched_pf_path_root import solve_price_path
    import run_e5f_transition_calibration as fertility
    joined,primitive,checks,rent_domain=baseline._runtime();model,calendar=primitive.model,primitive.calendar
    old=packet['old'];P=copy.deepcopy(old.parameters);P.psi_child=packet['psi'];tax=.08 if case=='equal-rebate-2pct' else .04;rebated=case!='baseline'
    P.tau_H=tax;P.user_cost_rate=P.q+P.delta+tax
    grid=old.b_grid;supply=old.supply_rule;demographics=packet['demographics'];start=time.monotonic();deadline=start+seconds
    controls=EndpointControls(**plan['terminal_template']['endpoint_controls']);audit=balanced.TerminalAuditControls(**plan['terminal_template']['audit_controls'])
    out.mkdir(parents=True,exist_ok=False);last=None;coordinates=None
    def root_solve(initial,evaluate,project,folder,initial_jacobian=None):
        def progress(record):
            save(folder/'latest_completed.json',record)
            if record.get('new_best'):save(folder/'best_so_far.json',record)
        result=solve_price_path(initial_prices=np.asarray(initial),evaluate=evaluate,project=project,slope=1.63,market_tolerance=2e-4,
            max_log_step=.2,damping=1.,max_evaluations=8,deadline_monotonic=min(deadline-60,time.monotonic()+3600),
            max_condition_number=1e10,worsening_factor=1.5,final_reproduction_tolerance=2e-10,callback=progress,initial_jacobian=initial_jacobian,
            default_jacobian=np.diag(np.r_[np.full(len(initial)//(3 if rebated else 2),-1.63),np.full(len(initial)-len(initial)//(3 if rebated else 2),-200.)]))
        save(folder/'root_receipt.json',result);return result
    def residual(demand,supply,accounts,revenue,transfer_outlays):
        return policy_residual(demand,supply,accounts,revenue,transfer_outlays,rebated)
    def terminal_evaluate(x):
        nonlocal last,coordinates
        Q=copy.deepcopy(P);Q.property_tax_lump_sum_transfer=float(x[2]) if rebated else 0.;bind_social_security_income(Q,pension_period=float(x[1]),payroll_tax=.179)
        shared=model.precompute_shared(Q,grid);price=np.array([x[0]])
        solution=model.solve_markov_income_at_prices(price,Q,grid,SD=shared);policy=calendar.policy_from_solution(solution,price,Q,grid,shared)
        seed,reconstruction=calendar.reconstruct_stationary_pre_fertility(solution,policy,Q,grid,shared);del solution
        endpoint=evaluate_endpoint(parameters=Q,b_grid=grid,policy=policy,asset_price=float(x[0]),transfer=Q.property_tax_lump_sum_transfer,psi_child=Q.psi_child,
            demographic_primitives=demographics,initial_g_pre=seed,supply_rule=supply,fiscal_regime='fixed_transfer',controls=controls)
        e=calendar.evaluate_period(price,endpoint.fixed_point.g_pre,Q,grid,shared,calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
        d,gates=balanced._household_checks(e,Q,shared,grid,float(Q.user_cost_rate*x[0]),primitive,audit)
        accounts=fiscal_accounts(e.g_current,Q);revenue=endpoint.residuals['tax_revenue'];outlays=Q.property_tax_lump_sum_transfer*float(e.g_current.sum())
        gates['seed_reconstruction']=all(abs(reconstruction[k])<=audit.reconstruction_tolerance for k in ('stationary_post_fertility_nesting_l1','stationary_post_fertility_nesting_max_abs'))
        gates['seed_projection']=reconstruction['stationary_feasibility_projection_mass']<=audit.feasibility_projection_tolerance
        last=NS(P=Q,policy=policy,endpoint=endpoint,evaluation=e,shared=shared);coordinates=np.asarray(x).copy()
        return dict(residual=residual(endpoint.residuals['housing_demand'],endpoint.residuals['housing_supply'],accounts,revenue,outlays),mapping_valid=bool(endpoint.mapping_valid and all(gates.values())),payload=dict(accounts=accounts,household_gates=gates,rebate_revenue=revenue,rebate_outlays=outlays))
    guess=[float(old.policy.price[0]),float(P.pension)]+([.25] if rebated else [])
    r=None
    for attempt in range(3):
        r=root_solve(guess,terminal_evaluate,lambda x:np.clip(x,.000001,10),out/f'terminal_{attempt}',None if r is None else r['final_jacobian'])
        if r['converged'] and r['final'] is not None and np.array_equal(r['final']['prices'],coordinates):break
        f=r.get('final') or r.get('best')
        if f is None:raise RuntimeError('No valid policy terminal')
        guess=f['prices']
    if not r['converged'] or not np.array_equal(r['final']['prices'],coordinates):raise RuntimeError('Policy terminal did not clear all budgets with its retained evaluation')
    terminal=last;terminal_x=coordinates.copy();nfirst=6
    reference=NS(parameters=terminal.P,asset_price=float(terminal_x[0]),renter_price=float(terminal.P.user_cost_rate*terminal_x[0]),
        equal_transfer=float(terminal_x[2]) if rebated else 0.,psi_child=P.psi_child,state=joined.person_pf.PersonPFState(terminal.endpoint.fixed_point.g_pre,terminal.endpoint.fixed_point.persons))
    history=[];last=None;previous=None;snapshot={}
    for count in (6,28,56):
        if time.monotonic()+180>deadline:break
        guess=np.r_[np.full(count,terminal_x[0]),np.full(count,terminal_x[1]),np.full(count,terminal_x[2]) if rebated else []]
        if previous is not None:
            for j in range(3 if rebated else 2):guess[j*count:j*count+min(count,previous_count)]=previous[j*previous_count:(j+1)*previous_count][:count]
        def project(x):
            x=np.clip(x,1e-6,10);x[:count]=rent_domain.project_price_path_to_positive_rents(x[:count],terminal=reference,minimum_rent_share=1e-6)[0];return x
        def path_evaluate(x):
            nonlocal last,coordinates
            prices=x[:count];benefits=x[count:2*count];transfers=x[2*count:] if rebated else np.zeros(count);accounts=[];all_gates=[];diagnostics=[]
            rents=joined.pf.rents_from_asset_prices(prices,reference.asset_price,P)
            def observe(i,e,Q,g,shared):
                d,gates=balanced._household_checks(e,Q,shared,g,float(rents[i]),primitive,audit)
                if not all(gates.values()):raise RuntimeError('Policy household gate failed')
                accounts.append(fiscal_accounts(e.g_current,Q));all_gates.append(gates);diagnostics.append(fertility.period_fertility_diagnostics(e,Q))
                if i==0:snapshot.clear();snapshot.update(parameters=Q,b_grid=g,evaluation=e,shared=shared,supply_rule=supply,actual_renter_price=float(rents[i]))
            path=joined.person_pf.evaluate_path_at_prices_person_demography(prices=prices,psi_path=np.full(count,P.psi_child),transfer_path=transfers,
                terminal_price=reference.asset_price,terminal_V=terminal.policy.V,base_parameters=P,b_grid=grid,initial_state=packet['inherited'].households,
                demographic_primitives=demographics,supply_rule=supply,pension_path=benefits,payroll_tax_path=np.full(count,.179),observer=observe)
            limits={'maximum_person_identity_error':2e-9,'maximum_head_identity_error':2e-9,'maximum_household_person_head_gap':2e-9,'maximum_age_head_gap':2e-9,'maximum_policy_reproduction_error':2e-10,'maximum_feasibility_projection_mass':1e-6}
            if any(not np.isfinite(getattr(path,k)) or getattr(path,k)>v for k,v in limits.items()):raise RuntimeError('Policy population/feasibility gate failed')
            blocks=[residual(row['housing_demand'],row['housing_supply'],a,row['property_tax_revenue'],row['equal_transfer_outlays']) for row,a in zip(path.rows,accounts)]
            distance=checks.terminal_convergence_diagnostics(path,terminal=reference,psi_path=np.full(count,P.psi_child));distance['pension_gap']=abs(float(benefits[-1])-float(terminal.P.pension))/float(terminal.P.pension)
            last=path;coordinates=x.copy()
            return dict(residual=np.asarray(blocks).T.reshape(-1),mapping_valid=True,payload=dict(terminal_distance=distance,fertility=diagnostics,accounts=accounts))
        r=None
        for attempt in range(3):
            r=root_solve(guess,path_evaluate,project,out/f'path_{count}_{attempt}',None if r is None else r['final_jacobian'])
            f=r.get('final') or r.get('best')
            if f is None:break
            guess=f['prices'];previous=guess.copy();previous_count=count
            if r['converged'] and np.array_equal(f['prices'],coordinates):break
        if not r['converged'] or not np.array_equal(r['final']['prices'],coordinates):
            if count==6:raise RuntimeError('Native policy smoke failed; no long policy solve')
            continue
        joined.pf.write_csv(out/f'policy_path_{count}.csv',last.rows);save(out/f'fertility_{count}.json',r['final']['payload']['fertility'])
        import run_e5f_independent_numerical_audit as audit_writer
        from unittest.mock import patch
        original=audit_writer.write_diagnostics
        def writer(stats,Q,destination):stats.owner_user_cost=np.array([snapshot['actual_renter_price']]);return original(stats,Q,destination)
        with patch.object(audit_writer,'write_diagnostics',writer):audit_writer.standard_diagnostics(snapshot,out/f'graphs_{count}',validate_production_young=False)
        if len(list((out/f'graphs_{count}'/'standard_diagnostics').glob('*.png')))!=17:raise RuntimeError('Incomplete stable diagnostic packet')
        distance=r['final']['payload']['terminal_distance'];history.append(dict(count=count,finite_converged=True,terminal_distance=distance))
        save(out/'summary.json',dict(case=case,annual_property_tax=tax/4,equal_rebate=rebated,paths=history,production_eligible=False))
        if distance['all_checks_pass'] and distance['pension_gap']<=.01:break
    if not history:raise RuntimeError('No converged policy path within budget')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--case',required=True,choices=['baseline','equal-rebate-1pct','equal-rebate-2pct']);p.add_argument('--output',type=Path,required=True);p.add_argument('--seconds',type=int,required=True);a=p.parse_args()
    with gzip.open(a.inputs,'rb') as f:packet=pickle.load(f)
    try:execute(packet,a.case,a.output,a.seconds)
    except Exception as e:
        a.output.mkdir(parents=True,exist_ok=True);(a.output/'failure.json').write_text(json.dumps(dict(error=str(e),type=type(e).__name__))+'\n');raise
