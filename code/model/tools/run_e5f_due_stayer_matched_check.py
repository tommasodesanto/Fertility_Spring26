#!/usr/bin/env python3
"""Pinned same-price DUE diagnostic; preparation never authorizes a solve.

This separate diagnostic reports market residuals and never calls them an
 equilibrium. It preserves fixed de_0093 parameters and strict household gates.
"""
from __future__ import annotations
import argparse,copy,csv,gzip,hashlib,json,os,pickle,signal,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
MAIN=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
PORTABLE=MAIN/'tmp/e5f_overnight_local_20260927/portable'
REFERENCE=PORTABLE/'night_launch_v4/primary_continuation/search/de_0093/case'
CONTRACT=PORTABLE/'night_launch_v4/primary_continuation/production_contract.json'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def table(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def verify(plan):
    assert plan['runner_sha256']==sha(__file__)
    assert plan['contract_sha256']==sha(plan['contract'])
    assert Path(plan['source_root']).resolve()==ROOT
    for relative,h in plan['source_files'].items():assert sha(ROOT/relative)==h,relative
    for name,h in plan['reference_files'].items():assert sha(Path(plan['reference'])/name)==h,name
    assert plan['audit_contract']['origin_specific'] is True
    for name,h in plan['audit_contract']['files'].items():assert sha(ROOT/name)==h,name
    assert plan['seconds_per_arm']<=300 and plan['total_seconds']<=1200
    return plan

def prepare(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'code/model').rglob('*.py')) if '__pycache__' not in p.parts}
    audit_files=['code/model/tools/e5f_overnight_estate_audit.py','code/model/tools/audit_e5f_estate_resource_account.py','code/model/tools/run_e5f_matched_pf_smoke.py','code/model/tools/run_dynamic_population_transition.py']
    plan=dict(status='prepared_not_authorized',execution_authorized=False,runner_sha256=sha(__file__),source_root=str(ROOT),source_files=files,contract=str(CONTRACT),contract_sha256=sha(CONTRACT),reference=str(REFERENCE),reference_files={n:sha(REFERENCE/n) for n in ['receipt.json','initial_state.pkl.gz','target_fit.csv','parameters.csv']},audit_contract=dict(origin_specific=True,files={n:files[n] for n in audit_files},purchase_audit_path=None,purchase_audit_sha256=None,purchase_audit_capability='SUPPORTS_NATIVE_DUE_STAYER_CREDIT'),seconds_per_arm=300,total_seconds=1200,absolute_end_epoch=None,arms=[dict(id='baseline',due=False,price_multiplier=1.,inherited_distribution=False),dict(id='due',due=True,price_multiplier=1.,inherited_distribution=False),dict(id='baseline_pricefall',due=False,price_multiplier=.9,inherited_distribution=True),dict(id='due_pricefall',due=True,price_multiplier=.9,inherited_distribution=True)],pricefall_status='Experimental 10 percent diagnostic; lead chooses and records concrete future-price path within existing authority',economic_changes=['DUE arm alone changes existing-owner borrowing rule; purchases retain baseline rule'],remaining_blockers=['Lead numerical diff review and test receipt','Explicit authenticated origin-aware purchase audit','Lead specification of diagnostic future-price path and dated observer; no additional user permission required','Immutable source snapshot including native source and audit closure','Explicit execution authorization and absolute deadline'])
    write(output/'plan_draft.json',plan)
    for rel in files:
        q=output/'source_snapshot'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes((ROOT/rel).read_bytes());assert sha(q)==files[rel]
    return plan

def run(plan_path,arm_id,output):
    import numpy as np
    import e5f_current_transition_runtime as native
    plan=verify(read(plan_path));assert plan['execution_authorized'] is True
    assert plan['absolute_end_epoch'] and time.time()<plan['absolute_end_epoch']
    arm=next(r for r in plan['arms'] if r['id']==arm_id)
    assert 0<arm['price_multiplier']<=1
    if arm['inherited_distribution']:
        raise RuntimeError('Pricefall preparation only: dated expectations and normalization-row observer need lead-approved specification; no solve')
    audit=plan['audit_contract'];audit_path=audit['purchase_audit_path']
    # Old origin-blind purchase accounting must never be admitted by accident.
    if not audit_path or not audit['purchase_audit_sha256']:raise RuntimeError('Authenticated origin-aware purchase audit missing; zero solves')
    assert sha(audit_path)==audit['purchase_audit_sha256']
    purchase=native.load('due_explicit_purchase_audit',audit_path)
    assert getattr(purchase,audit['purchase_audit_capability'],False) is True
    estate_path=ROOT/'code/model/tools/e5f_overnight_estate_audit.py'
    estate=native.load('due_explicit_estate_audit',estate_path)
    assert getattr(estate,'SUPPORTS_NATIVE_DUE_STAYER_CREDIT',False) is True
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    deadline=min(time.time()+plan['seconds_per_arm'],plan['absolute_end_epoch'])
    def timeout(*_):raise TimeoutError('Matched-check arm deadline')
    previous=signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,deadline-time.time())
    try:
        prepared=native.setup(out/'preparation',contract=Path(plan['contract']),reference=Path(plan['reference']),fixed_reference_price=True)
        for source_key,h in prepared['preparation']['current_source_files'].items():
            source_path=Path(source_key)
            rel=str(source_path.resolve().relative_to(ROOT)) if source_path.is_absolute() else str(source_path)
            assert plan['source_files'].get(rel)==h,source_key
        rt=prepared['runtime'];model=rt['model'];P=copy.deepcopy(prepared['parameters']);selected=prepared['selected'];grid=selected['b_grid']
        P.native_due_stayer_credit=bool(arm['due'])
        assert not getattr(P,'native_solvency_credit',False)
        rt['accounting']=purchase
        prepared['objective'].estate=native.EstateAuditContract(estate)
        original=prepared['parameters'];price=np.asarray(selected['solution'].p_eq)*arm['price_multiplier']
        start=time.monotonic();sol=model.solve_markov_income_at_prices(price,P,grid,verbose=False,fast_stats=False);elapsed=time.monotonic()-start
        for key,value in vars(original).items():
            if key.startswith('_') or key in ('native_due_stayer_credit','eq_iter'):continue
            other=getattr(P,key,None)
            try:equal=np.array_equal(value,other,equal_nan=True)
            except (TypeError,ValueError):equal=repr(value)==repr(other)
            assert equal,'Parameter changed: '+key
        shared=model.precompute_shared(P,grid);P._fert2_probs=sol.fert2_probs.copy()
        cal=rt['primitive'].pf.calendar;policy=cal.policy_from_solution(sol,price,P,grid,shared)
        if arm['inherited_distribution']:
            pre=np.asarray(selected['stationary_g_pre']).copy();original_pre=pre.copy()
            # Strict gate before calendar's projection machinery: never move a household.
            for age in range(P.J):
                model._gate_dead_mass_at_age(pre[:,:,:,age,:,:,:],policy.V[:,:,:,age,:,:,:],age,'matched_pricefall_inherited',P.user_cost_rate*price,price,P,grid,shared,markov_income=True)
            np.testing.assert_array_equal(pre,original_pre)
            operator={'status':'dated inherited distribution; stationarity not imposed'}
        else:
            pre,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
            operator=rt['primitive'].pf.transition.operator_gates(sol,policy,pre,P,grid,shared);operator.update(reconstruction)
            native.write(out/'operator_diagnostics.json',prepared['tax'].finite_json(cal.jsonable(operator)))
            for name in ('stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1','mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):assert abs(operator[name])<=5e-9,name
            assert abs(operator['zero_entry_mass_accounting_residual'])<=2e-8
            assert operator['stationary_feasibility_projection_mass']==0.
        supply=cal.HousingSupplyRule('static-elastic',float(price[0]),float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))
        evaluation=cal.evaluate_period(price,pre,P,grid,shared,cal.SolveCounter(),supply_rule=supply,supplied_policy=policy)
        if arm['inherited_distribution']:
            assert evaluation.feasibility_projection_mass==0.,'Inherited distribution projection prohibited'
            np.testing.assert_array_equal(pre,original_pre)
        fiscal=rt['certify_initial_pension'](evaluation.g_current,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        budget=rt['primitive'].dated_budget(evaluation,P,shared,grid,float(P.user_cost_rate*price[0]))
        purchasing=purchase.audit_purchase_accounting(evaluation,P,shared,grid,model)
        estates=prepared['objective'].estate.audit(evaluation,P,grid)
        packet=dict(parameters=P,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,solution=sol,stationary_g_pre=pre,demographic_seed=selected.get('demographic_seed'))
        checkpoint=out/'initial_state.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as f:pickle.dump(packet,f,protocol=5)
        arrays=rt['audit'].policy_array_audit(packet,out)
        assert arrays['occupied_negative_steps']==0
        assert all(not r['nonfinite'] and r['minimum']>=0 and r['maximum']<=1 for r in arrays['probabilities'].values())
        fertility={p:rt['observe_initial_fertility'](evaluation,P,age_projection=p) for p in ('uniform_birth_time','constant_post_cell')}
        housing=rt['observe_initial_housing_wealth'](evaluation,P,grid,shared,diagnostic_enabled=True,age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=True)
        recent=rt['observe_recent_parent_flow'](evaluation,P,diagnostic_enabled=True,snapshot=rt['SNAPSHOT'],age_projection=rt['AGE_PROJECTION'],diagnostic_allow_residence_proxy=True,input_provenance=dict(case_id=arm_id,checkpoint_sha256=sha(checkpoint)))
        reporter=sys.modules[prepared['objective'].__class__.__module__]
        fit=reporter.score_targets(prepared['objective_definition'],fertility,housing,recent['model_value'],float(rt['chain'].extract_moments(sol,P)['tfr']))
        table(out/'target_fit.csv',fit)
        params=rows(Path(plan['reference'])/'parameters.csv')
        for row in params:row['status']='Fixed de_0093 value; no recalibration or fertility normalization'
        assert len(fit)==14 and len(params)==31;table(out/'parameters.csv',params)
        comparison=native.compare_arrays(selected,packet);write(out/'reference_array_comparison.json',comparison)
        if arm_id=='baseline':
            # Preserve the authored initial array exactly. Removing the old tiny
            # projection changes only three evaluation distributions; core
            # policies, shared arrays and stationary distribution stay exact.
            np.testing.assert_array_equal(evaluation.g_pre,selected['stationary_g_pre'])
            evaluation_mass={'evaluation.g_pre','evaluation.g_post_fertility','evaluation.g_current'}
            for name,row in comparison['arrays'].items():
                if name in evaluation_mass:
                    assert row.get('finite') and row['l1']<=2*model.DEAD_MASS_TOL,(name,row)
                else:
                    assert row.get('exact',False),(name,row)
            old={r['moment']:r for r in rows(Path(plan['reference'])/'target_fit.csv')}
            for row in fit:
                for key in ('target','model','gap','weight','loss_contribution'):
                    if key=='loss_contribution' and row[key]!='':
                        assert float(row[key])==float(row['weight'])*float(row['gap'])**2
                        continue
                    assert str(row[key])==old[row['moment']][key] or (row[key]!='' and abs(float(row[key])-float(old[row['moment']][key])) <= (1e-12 if key in ('model','gap','loss_contribution') else 0.)),(row['moment'],key)
        rt['audit'].standard_diagnostics(packet,out,validate_production_young=False)
        assert len(list((out/'standard_diagnostics').glob('*.png')))==17
        verify(plan)
        result=dict(status='matched_partial_equilibrium_diagnostic',arm=arm,market_residual=float(evaluation.relative_market_residual),equilibrium_certified=False,solve_seconds=elapsed,parameters_fixed=True,operator=operator,fiscal=fiscal,budget=budget,purchase=purchasing,estate_funding=estates,policy_array_audit=arrays,plan_sha256=sha(plan_path),checkpoint_sha256=sha(checkpoint),source_files=plan['source_files'])
        native.write(out/'complete.json',prepared['tax'].finite_json(cal.jsonable(result)));return result
    except BaseException as exc:
        write(out/'failure.json',dict(status='failed_no_gate_relaxation',error_type=type(exc).__name__,error=str(exc),arm=arm))
        raise
    finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous)

def main():
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--plan',type=Path);p.add_argument('--arm');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.prepare:prepare(a.output)
    else:
        assert a.plan and a.arm
        for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
        run(a.plan,a.arm,a.output)
if __name__=='__main__':main()
