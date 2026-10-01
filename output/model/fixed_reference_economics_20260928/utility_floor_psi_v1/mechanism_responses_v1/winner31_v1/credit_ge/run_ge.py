#!/usr/bin/env python3
"""Bounded support-limited native-credit GE; no refit or full-support certificate."""
from __future__ import annotations
import argparse,copy,gzip,hashlib,importlib.util,json,math,os,pickle,signal,subprocess,sys,time,traceback
from pathlib import Path
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
HERE=Path(__file__).resolve().parent  # Isolated winner-specific GE controller
ROOT=HERE.parents[6]
MECHANISM=HERE.parent
OLD_GE=ROOT/'output/model/fixed_reference_economics_20260928/credit_ge_v1/run_ge.py'
RENEWAL_TOL=1e-6
PAYGO_TOL=1e-6

def require(ok,message):
    if not ok:raise RuntimeError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n');temp.replace(path)
def progress(out,status,**extra):write(Path(out)/'latest.json',dict(status=status,checked_epoch=time.time(),pid=os.getpid(),**extra))
def old_math():
    spec=importlib.util.spec_from_file_location('frozen_monday_ge_closure',OLD_GE)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def verify_plan(path,*,zeroLC=False):
    p=read(path)
    if not zeroLC:require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm required')
    require(p['schema']=='utility_floor_native_credit_ge_v1' and p['candidate_case']=='chain_7/0173_nm','Candidate/schema drift')
    require(p['maximum_lifecycle_solves']==11 and p['prior_attempts_reserved']==1 and p['reserve_exact_repeat_solves']==1 and p['case_seconds']==300 and p['total_seconds']==2400 and p['threads']==1 and p['memory_gib']==24,'Budget drift')
    require(p['lower_price_factors']==[.95,.90,.80] and p['upper_price_factors']==[1.05,1.10,1.20,1.35]
            and p['lower_price_search_authorized'] and p['upper_price_search_authorized']
            and p['lower_price_domain']==[.8,1.] and p['upper_price_domain']==[1.,1.35],
            'Reviewed directional price domains drift')
    require(p['credit_mode']=='existing_native_solvency_support_limited' and not p['natural_support_certified'] and not p['refit'] and p['hbar_first_child_jump']==2.3 and p['hbar_child_rooms']==0. and p['psi_child']==0.17156192800028292 and p['candidate_q0']==0.719168368828958,'Economic contract drift')
    for relative,digest in p['files'].items():require(sha(ROOT/relative)==digest,'Source/input drift: '+relative)
    require(p['files'][str(Path(__file__).resolve().relative_to(ROOT))]==sha(__file__),'Driver missing pin')
    return p

def candidate_runtime(out):
    sys.path.insert(0,str(MECHANISM))
    import fixed_price_responses as d
    auth=d.authenticate_candidate(out)
    P=auth['natural']
    require(P.hbar_first_child_jump==2.3 and P.hbar_child_rooms==0. and P.child_room_floor,'First-child floor drift')
    require(P.psi_child==0.17156192800028292 and P.native_solvency_credit and not P.native_due_stayer_credit and P.unsecured_credit_limit is None,'Credit/preferences drift')
    require(len(auth['actual_parameters'])==31 and len(auth['grid'])==120 and P.Nz==9,'Parameters/grid drift')
    require(auth['entry']['conditional_sha256']=='bf976b033678629c4f9e91df4037092d6e35ab113eebaf3fcf63c59510d2ca4f','Entry drift')
    return d,auth

def solve_case(args,p):
    import numpy as np
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    d,auth=candidate_runtime(out/'runtime')
    from refactor_lab.inputs import serialized
    P=copy.deepcopy(auth['natural']);grid=auth['grid'];q=p['candidate_q0']*args.factor
    P.native_inherited_distribution_evidence_dir=str(out/'inherited_state_diagnostics')
    public_before=serialized({k:v for k,v in vars(P).items() if not k.startswith('_')})
    sd=auth['solver'].precompute_shared(P,grid)
    require(time.time()<args.deadline_epoch,'Case deadline before lifecycle')
    progress(out,'lifecycle_claimed',price=q,factor=args.factor,lifecycle_solves=1,deadline_epoch=args.deadline_epoch)
    started=time.monotonic()
    import refactor_lab.engine.household as household
    from credit_mode import trace_native_support,lower_grid_diagnostic
    with trace_native_support(household) as calls:
        sol=auth['solver'].solve_markov_income_at_prices(np.asarray([q]),P,grid,SD=sd,verbose=False,fast_stats=False)
    elapsed=time.monotonic()-started
    require(time.time()<args.deadline_epoch,'Case deadline during lifecycle')
    require(float(getattr(P,'_entry_censored_mass',0.))<=auth['credit'].DEAD_MASS_TOL,'Entry censoring')
    prepared=auth['context']['prepared'];rt=prepared.rt;cal=rt['primitive'].pf.calendar
    P._fert2_probs=sol.fert2_probs.copy()
    policy=cal.policy_from_solution(sol,np.asarray([q]),P,grid,sd)
    pre,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,sd)
    runtime=auth['context']['runtime']
    runtime.require_abs_gate(reconstruction['stationary_post_fertility_nesting_l1'],5e-9,'Cohort reconstruction')
    runtime.require_abs_gate(reconstruction['stationary_feasibility_projection_mass'],0.,'Cohort projection')
    supply=cal.HousingSupplyRule('static-elastic',q,float(P.H0[0]*(P.user_cost_rate*q/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
    ev=cal.evaluate_period(np.asarray([q]),pre,P,grid,sd,cal.SolveCounter(),supply_rule=supply,supplied_policy=policy)
    packet=dict(parameters=P,b_grid=grid,shared=sd,solution=sol,evaluation=ev,stationary_g_pre=pre,supply_rule=supply,demographic_seed=None)
    support=lower_grid_diagnostic(sol,grid,calls,P,realized_distribution=ev.g_current)
    write(out/'support_diagnostic.json',support)
    require(support['status']=='occupied_support_pass_unoccupied_alternatives_unverified','Occupied support failed; no GE result')
    gates=d.case_gates(auth,packet,out,'lifetime_repayment_only',stationary=True)
    fiscal=float(gates['fiscal']['scaled_pension_budget_residual'])
    require(abs(fiscal)<=PAYGO_TOL and gates['fiscal_certificate']['fiscal_gate'] and gates['fiscal_certificate']['marginal_gate'],'Actual PAYGO/age-income margin gate')
    adjusted=float(rt['primitive'].pf.transition.calendar_topcode_birth_accounting(ev.g_pre,ev.g_post_fertility,float(ev.births),P)['topcode_adjusted_birth_children'])
    old=old_math()
    closure=old.closed_accounting(float(sol.entry_rate),adjusted,float(np.asarray(ev.demand_by_loc).sum()),float(np.asarray(ev.supply_by_loc).sum()),q)
    require(float(P.xi_supply[0])==.63 and abs(closure['absolute_housing_residual'])<=1e-12,'Absolute supply closure')
    selected=abs(closure['renewal_residual'])<=RENEWAL_TOL
    queue=old.scaled_native_step(packet,prepared,closure) if selected else None
    fertility={kind:rt['observe_initial_fertility'](ev,P,age_projection=kind) for kind in ('uniform_birth_time','constant_post_cell')}
    housing=rt['observe_initial_housing_wealth'](ev,P,grid,sd,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
    recent=rt['observe_recent_parent_flow'](ev,P,diagnostic_enabled=True,snapshot=rt['SNAPSHOT'],age_projection=rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,input_provenance=dict(case_id=args.case,candidate_case=p['candidate_case']))
    completed=float(rt['chain'].extract_moments(sol,P)['tfr'])
    fits=runtime.score_targets(auth['context']['objective'],fertility,housing,recent['model_value'],completed)
    d.table(out/'target_fit.csv',fits)
    fits=auth['native'].readtable(out/'target_fit.csv')
    require(len(fits)==14 and auth['native'].target_identity(fits)==auth['native'].PLAN['target_contract'],'All14 original targets/weights')
    auth['native'].residual(fits)
    actual=auth['context']['fp'].actual_parameters(prepared,P,grid)
    auth['ge'].validate_parameter_estimates(auth['context'],auth['params_rows'],actual)
    require(actual==auth['actual_parameters'],'Effective31 parameter drift')
    params=copy.deepcopy(auth['params_rows'])
    for row in params:row['status']='Fixed candidate value; support-limited credit GE, no refit'
    d.table(out/'target_fit.csv',fits);d.table(out/'parameters.csv',params)
    write(out/'observers.json',cal.jsonable(dict(fertility=fertility,housing_wealth=housing,recent_parent=recent)))
    plots=[]
    if selected or args.role=='seed':
        report_ev=copy.copy(ev)
        if selected:
            report_ev.supply_by_loc=np.asarray(ev.supply_by_loc)/closure['population_scale']
            report_ev.relative_market_residual=float(np.max(np.abs((ev.demand_by_loc-report_ev.supply_by_loc)/report_ev.supply_by_loc)))
        rt['audit'].standard_diagnostics(dict(packet,evaluation=report_ev),out,validate_production_young=False)
        plots=sorted(x.name for x in (out/'standard_diagnostics').glob('*.png'))
        require(plots==sorted(auth['context']['manifest']['standard_diagnostic_names']),'Standard17 plot names')
        write(out/'reporting_units.json',dict(economic_absolute_supply=closure['absolute_housing_supply'],reporting_supply_per_household=float(report_ev.supply_by_loc.sum()),population_scale=closure['population_scale'],economic_H0_unchanged=True,supply_divided_by_population_for_plot=selected,q0_not_cleared=not selected))
    plot_hashes={name:sha(out/'standard_diagnostics'/name) for name in plots}
    repeat=old.compare_repeat(packet,fits,params,out,auth['context']['fp']) if args.role=='repeat' else None
    if args.role=='repeat':
        chosen=read(out.parent/'selected.json');first=out.parent/chosen['case']
        first_receipt=read(first/'receipt.json')
        first_names=sorted(x.name for x in (first/'standard_diagnostics').glob('*.png'))
        first_hashes={name:sha(first/'standard_diagnostics'/name) for name in first_names}
        require(len(first_hashes)==len(plot_hashes)==17 and first_hashes==first_receipt['standard_plot_sha256']==plot_hashes,'Selected/repeat17 PNG hashes differ')
        repeat['standard_plot_hashes_exact']=17
    with gzip.open(out/'conditional_cohort_state.pkl.gz','wb',compresslevel=1) as stream:pickle.dump(packet,stream,protocol=5)
    checkpoint=dict(path=str(out/'conditional_cohort_state.pkl.gz'),sha256=sha(out/'conditional_cohort_state.pkl.gz'))
    require(serialized({k:v for k,v in vars(P).items() if not k.startswith('_')})==public_before,'Public parameter mutation')
    verify_plan(args.plan)
    require(time.time()<args.deadline_epoch,'Case deadline during reporting')
    receipt=dict(status='passed_support_limited_case',case=args.case,role=args.role,price=q,factor=args.factor,lifecycle_solves=1,lifecycle_seconds=elapsed,renewal_residual=closure['renewal_residual'],paygo_residual=fiscal,closure=closure,native_scaled_step=queue,repeat_check=repeat,support_diagnostic=support,gates=gates,reconstruction=reconstruction,checkpoint=checkpoint,standard_plot_count=len(plots),standard_plot_sha256=plot_hashes,target_fit_sha256=sha(out/'target_fit.csv'),parameters_sha256=sha(out/'parameters.csv'),candidate_case=p['candidate_case'],plan_sha256=sha(args.plan),natural_support_certified=False,production_adoption=False,estate_closure='Candidate provisional net-estate funding/residual sink retained; counterparties unresolved')
    write(out/'receipt.json',cal.jsonable(receipt));progress(out,'case_completed',lifecycle_solves=1,receipt_sha256=sha(out/'receipt.json'))

def controller(args,p,*,evaluate=None,now=time.time,mock=False):
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    started=now();deadline=min(float(args.deadline_epoch),started+p['total_seconds'])
    records=[];attempts=0;old=old_math()
    # The six-price experiment has already solved the exact q0 expanded-credit
    # cell. Reuse its signed renewal residual and certificate as a search seed.
    # This is not counted as a new lifecycle solve or a GE root certificate.
    q0_case=Path(args.q0_credit_case)
    q0_receipt=read(q0_case/'receipt.json');q0_closure=read(q0_case/'closure.json')
    require(q0_receipt['status']=='completed_support_limited_diagnostic'
            and q0_receipt['regime']=='lifetime_repayment_only'
            and float(q0_receipt['price_factor'])==1.0
            and q0_receipt['source_binding_sha256']==sha(MECHANISM/'source_binding.json'),
            'Exact q0 expanded-credit fixed-price receipt is required')
    require(float(q0_closure['candidate_base_loss'])==float(p['candidate_loss'])
            and float(q0_closure['candidate_price_q0'])==float(p['candidate_q0'])
            and float(q0_closure['price'])==float(p['candidate_q0'])
            and float(q0_closure['candidate_psi_child'])==float(p['psi_child'])
            and int(q0_closure['grid_nodes'])==120
            and int(q0_closure['standard_plot_count'])==17,
            'q0 fixed-price seed candidate or price differs')
    require(q0_closure['support_diagnostic']['status']=='occupied_support_pass_unoccupied_alternatives_unverified'
            and q0_closure['natural_support_certified'] is False,
            'q0 fixed-price occupied-support certificate is insufficient')
    require(q0_receipt['target_fit_sha256']==sha(q0_case/'target_fit.csv')
            and q0_receipt['parameters_sha256']==sha(q0_case/'parameters.csv'),
            'q0 credit tables differ from their receipt')
    sys.path.insert(0,str(ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_round2_v1'))
    import runner as floor_runner
    q0_targets=floor_runner.readtable(q0_case/'target_fit.csv')
    require(len(q0_targets)==14 and floor_runner.target_identity(q0_targets)==floor_runner.PLAN['target_contract'],
            'q0 credit target identity differs from original 14-row contract')
    q0_parameters=floor_runner.readtable(q0_case/'parameters.csv')
    winner_parameters=floor_runner.readtable(ROOT/'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/ROOT/parameters.csv')
    q0_effective={r['parameter']:float(r['estimate']) for r in q0_parameters}
    winner_effective={r['parameter']:float(r['estimate']) for r in winner_parameters}
    require(len(q0_parameters)==len(winner_parameters)==31 and q0_effective==winner_effective,
            'q0 credit effective parameters differ from verified winner')
    q0_seed=dict(case='fixed_price_q0_seed',role='fixed_price_seed',factor=1.,
        price=float(q0_closure['price']),
        residual=float(q0_closure['renewal_residual_reported_not_imposed']),
        population=None,receipt_sha256=sha(q0_case/'receipt.json'),checkpoint=None)
    require(math.isfinite(q0_seed['residual']), 'q0 renewal residual is not finite')
    direction='root_at_q0' if abs(q0_seed['residual'])<=RENEWAL_TOL else ('lower' if q0_seed['residual']<0 else 'upper')
    factors=[] if direction=='root_at_q0' else p[direction+'_price_factors']
    records.append(q0_seed)
    write(out/'launch.json',dict(started_epoch=started,deadline_epoch=deadline,plan_sha256=sha(args.plan),plan=p,slurm_job=os.environ.get('SLURM_JOB_ID'),mock=mock,lifecycle_cap=11,prior_attempts_reserved=1,natural_support_certified=False))
    write(out/'latest_completed.json',dict(completed=[],lifecycle_attempts=0));write(out/'best_so_far.json',dict(status='no_case_completed'))
    def run_case(name,factor,role):
        nonlocal attempts
        if not mock:require(verify_plan(args.plan)==p,'Source/plan drift during controller')
        require(attempts<11-(1 if role!='repeat' else 0),'11-new-attempt cap/repeat reserve plus one prior reservation')
        require(now()<deadline and (role=='repeat' or deadline-now()>300),'Global deadline or repeat time reserve')
        case_end=min(deadline if role=='repeat' else deadline-300.,now()+300);attempts+=1
        progress(out,'case_dispatch',case=name,factor=factor,role=role,lifecycle_attempts=attempts,deadline_epoch=case_end)
        if evaluate is not None:receipt=evaluate(name,factor,role,case_end,out)
        else:
            command=[sys.executable,str(Path(__file__).resolve()),'--plan',str(args.plan),'--out',str(out/name),'--case',name,'--factor',repr(factor),'--role',role,'--deadline-epoch',repr(case_end)]
            with (out/(name+'.log')).open('w') as log:
                process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1'),start_new_session=True)
                try:
                    while process.poll() is None:
                        progress(out,'case_running',case=name,child_pid=process.pid,lifecycle_attempts=attempts,deadline_epoch=case_end)
                        if now()>=case_end:raise TimeoutError('300-second case or global deadline: '+name)
                        time.sleep(2)
                    require(process.returncode==0,'Case failed without retry: '+name)
                finally:
                    if process.poll() is None:
                        os.killpg(process.pid,signal.SIGTERM)
                        try:process.wait(timeout=min(3.,max(.001,deadline-now())))
                        except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
            receipt=read(out/name/'receipt.json')
        require(receipt['status']=='passed_support_limited_case' and receipt['lifecycle_solves']==1 and receipt['plan_sha256']==sha(args.plan),'Case receipt identity')
        row=dict(case=name,role=role,factor=factor,price=receipt['price'],residual=receipt['renewal_residual'],population=receipt['closure']['population_scale'],receipt_sha256=sha(out/name/'receipt.json'),checkpoint=receipt['checkpoint'])
        records.append(row);write(out/'latest_completed.json',dict(completed=records,lifecycle_attempts=attempts))
        best=min(records,key=lambda r:abs(r['residual']))
        write(out/'best_so_far.json',dict(status='fixed_price_seed_only' if best['role']=='fixed_price_seed' else 'completed_support_limited_case',**best))
        return row,receipt
    try:
        selected=None
        if abs(q0_seed['residual'])<=RENEWAL_TOL:
            # A fixed-price seed at the root still needs one native GE check
            # before its selected-point repeat can be certified.
            q0,qreceipt=run_case('q0_selected',1.,'seed')
            require(qreceipt['standard_plot_count']==17
                    and abs(q0['residual'])<=RENEWAL_TOL
                    and abs(q0['residual']-q0_seed['residual'])<=1e-9,
                    'q0 native GE check disagrees with fixed-price seed')
            selected=q0
        bracket=None
        if selected is None:
            for i,factor in enumerate(factors):
                row,receipt=run_case(direction+'_'+str(int(round(100*factor))),factor,'search')
                bracket=old.choose_bracket(records)
                if abs(row['residual'])<=RENEWAL_TOL:selected=row;break
                if bracket is not None:break
        require(selected is not None or bracket is not None,
                'No signed renewal bracket inside reviewed '+direction+' domain; no result')
        if selected is None:
            left,right=bracket
            for i in range(10-attempts):
                row,receipt=run_case('root_'+str(i+1).zfill(2),old.next_log_factor(left,right),'search')
                if abs(row['residual'])<=RENEWAL_TOL:selected=row;break
                if row['residual']*left['residual']>0:left=row
                else:right=row
                write(out/'bracket.json',dict(left=left,right=right))
        require(selected is not None,'Solve cap without1e-6 renewal root; no result')
        chosen=read(out/selected['case']/'receipt.json')
        require(chosen['native_scaled_step']['status']=='passed' and chosen['standard_plot_count']==len(chosen['standard_plot_sha256'])==17,'Selected root queue/plots incomplete')
        write(out/'selected.json',selected)
        repeat,repeated=run_case('selected_repeat',selected['factor'],'repeat')
        require(repeated['repeat_check']['status']=='passed' and repeat['residual']==selected['residual'] and abs(repeat['residual'])<=RENEWAL_TOL and repeated['native_scaled_step']['status']=='passed' and repeated['standard_plot_count']==17,'Fresh selected repeat gates')
        require(repeated['repeat_check']['standard_plot_hashes_exact']==17 and repeated['standard_plot_sha256']==chosen['standard_plot_sha256'],'Fresh repeat PNG hashes differ')
        require(now()<deadline,'Global deadline during final reporting')
        result=dict(status='completed_support_limited_stationary_ge_diagnostic',direction=direction,selected=selected,repeat=repeat,total_lifecycle_attempts=attempts,natural_support_certified=False,production_adoption=False,elapsed_seconds=now()-started,deadline_epoch=deadline,interpretation='Native-credit GE closure passes; full natural support and grid convergence unverified; no transition or refit')
        write(out/'completed.json',result);progress(out,'complete',lifecycle_attempts=attempts)
        return result
    except BaseException as exc:
        failure=dict(status='failed_no_ge_result',error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),lifecycle_attempts=attempts,completed=records,no_retry=True,natural_support_certified=False)
        write(out/'failure.json',failure);write(out/'completed.json',failure);raise

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--deadline-epoch',type=float,required=True);ap.add_argument('--q0-credit-case',type=Path);ap.add_argument('--case');ap.add_argument('--factor',type=float);ap.add_argument('--role',choices=('seed','search','repeat'))
    args=ap.parse_args();args.out=args.out.resolve();args.plan=args.plan.resolve();p=verify_plan(args.plan)
    require(time.time()<args.deadline_epoch,'Absolute deadline already reached')
    if args.case:
        require(args.factor is not None and .8<=args.factor<=1.35 and args.role and args.deadline_epoch<=time.time()+300.,'Invalid bounded child')
        try:solve_case(args,p)
        except BaseException as exc:
            if args.out.exists():write(args.out/'failure.json',dict(status='fatal_case_failure',error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc(),no_retry=True))
            raise
    else:
        require(args.factor is None and args.role is None,'Controller owns price and role')
        require(args.q0_credit_case is not None,'Exact completed q0 expanded-credit seed required')
        controller(args,p)
if __name__=='__main__':main()
