#!/usr/bin/env python3
"""Saved-packet verification and summaries; no model solves or changed policies."""
import argparse,ast,csv,gzip,hashlib,inspect,json,pickle
from pathlib import Path
import numpy as np
import run_e5f_housing_fertility_cost_diagnostic as driver
import run_e5f_soft_housing_probe as common


def consolidate_saved(out):
    """Recover a complete partial-run manifest after the numerical time cutoff."""
    if (out/'complete.json').exists():return
    history=list(csv.DictReader((out/'solve_history.csv').open()))
    design=json.loads((out/'design.json').read_text());records=[]
    for path in out.glob('*/receipt.json'):
        record=json.loads(path.read_text())
        if record.get('phase') in ('fixed','matched'):records.append(record)
    order={arm:i for i,(arm,_) in enumerate(driver.ARMS)}
    records.sort(key=lambda r:(r['phase'],order[r['arm']],r['price_factor']))
    matches={}
    for arm,_ in driver.ARMS:
        base=next((r for r in records if r['phase']=='matched' and r['arm']==arm and r['price_factor']==1),None)
        shock=next((r for r in records if r['phase']=='matched' and r['arm']==arm and r['price_factor']==1.1),None)
        searches=[r for r in history if r['arm']==arm and r['case'].startswith('matchsearch')]
        if base and shock:
            matches[arm]=dict(status='matched',b=base['b'],fertility=base['fertility'],search_evaluations=len(searches))
        else:
            last=searches[-1] if searches else None
            matches[arm]=dict(status='pending_time_budget',last_fertility=float(last['fertility']) if last else 0.,last_b=float(last['b']) if last else None,search_evaluations=len(searches))
    common.write(out/'complete.json',dict(status='partial_numerical_budget_stop',records=records,matching_status=matches,evaluations=len(history)+1,elapsed_seconds=float(history[-1]['elapsed_seconds']),job=design['job'],stop_reason='Recovered from saved case receipts after numerical time cutoff; no new solves.'))


def all_income_grid(packet, reference_g):
    P=packet['parameters'];policy=packet['evaluation'].policy;grid=packet['b_grid'];rows=[]
    for j in range(P.J):
        age=P.age_start+j*P.da
        if age not in (22,26,30):continue
        for z in range(reference_g.shape[4]):
            for n,m in ((0,0),(1,1)):
                pair=policy.fert_probs[:,0,0,j,z,:2] if n==0 else policy.fert2_probs[:,0,0,j,z,:,n-1,m]
                for ib,b in enumerate(grid):
                    p0,p1=map(float,pair[ib]);finite=p0>0 and p1>0
                    rows.append(dict(age=age,income_index=z,income_value=float(P.z_grid[z]),birth_number=n+1,wealth=float(b),common_mass=float(reference_g[ib,0,0,j,z,n,m]),attempt_probability=p1,gap_over_shock_scale=float(np.log(p1)-np.log(p0)) if finite else '',choice_status='finite' if finite else ('unavailable' if p0==p1==0 else 'endpoint'),conditional_renter_housing=float(policy.hR_pol[ib,0,0,j,z,n,m]),conditional_renter_consumption=float(policy.c_pol[ib,0,0,j,z,n,m])))
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();out=a.output
    consolidate_saved(out)
    rt,tax,obj,selected=driver.setup(out/'audit_context')
    binary_source=inspect.getsource(rt['model'].solve_bellman_full_markov_income)
    assert 'lf = Vfa / P.kappa_fert' in binary_source
    assert 'ls, pr = logsumexp(lf, axis=3)' in binary_source
    assert 'logsumexp(V2 / kf_cont, axis=3)' in binary_source
    from e5f_social_security import fiscal_accounts
    completed=json.loads((out/'complete.json').read_text());records=completed['records']
    from e5f_housing_fertility_cost_resume import fingerprint
    parent=Path(completed['original_source']) if 'original_source' in completed else None
    parent_before=fingerprint(parent) if parent else None
    refP=selected['packet']['parameters'];refg=selected['packet']['stationary_g_pre']
    j26=int((26-refP.age_start)/refP.da);income_mass=refg[:,0,0,j26,:,0,0].sum(axis=0)
    income_values=np.asarray(refP.z_grid);order=np.argsort(income_values);cdf=np.cumsum(income_mass[order]);cdf/=cdf[-1]
    indices=[int(order[min(int(np.searchsorted(cdf,q)),len(order)-1)]) for q in (.1,.5,.9)]
    assert len(set(indices))==3,indices
    common.write(out/'plot_income_states.json',dict(indices=indices,income_values=income_values[indices].tolist(),quantiles=[.1,.5,.9],definition='Income-state weighted deciles among control childless current renters age26; fixed across all panels.'))
    checks=[];sensitivity=[];states=[];parameter_changes=[];grid_summaries=[];common_gates=[]
    for rec in records:
        folder=out/rec['case'];fits=list(csv.DictReader((folder/'target_fit.csv').open()));params=list(csv.DictReader((folder/'parameters.csv').open()))
        assert len(fits)==13 and len(params)==32
        for row in fits:
            np.testing.assert_allclose(float(row['gap']),float(row['model'])-float(row['target']),rtol=0,atol=1e-12)
            if row['weight']!='':np.testing.assert_allclose(float(row['loss_contribution']),float(row['weight'])*float(row['gap'])**2,rtol=1e-12,atol=1e-10)
        with gzip.open(folder/'state.pkl.gz','rb') as f:packet=pickle.load(f)
        P=packet['parameters'];fiscal=fiscal_accounts(packet['evaluation'].g_current,P)
        for row in params:
            if row['status']=='experimental free coordinate':row['status']='Inherited search coordinate; fixed in this diagnostic'
            if row['parameter']=='pension_period':row['status']='Held at authenticated baseline pension; no rebalancing'
            if row['parameter']=='payroll_tax':row['status']='Held at authenticated baseline tax rate'
            if row['parameter']=='psi_child':row['estimate']=float(P.psi_child)
            if row['parameter']=='psi_child' and rec['phase']=='matched':
                lo,hi=driver.B_BRACKET;row.update(lower=lo,upper=hi,near_bound=min(P.psi_child-lo,hi-P.psi_child)<=.02*(hi-lo))
        common.table(folder/'parameters.csv',params)
        values={r['parameter']:float(r['estimate']) for r in params}
        expected=tax.actual_parameters(P)
        expected.update(theta1=P.theta1,psi_child=P.psi_child,payroll_tax=P.tau_pay,
            pension_period=P.pension,housing_supply_elasticity=P.xi_supply[0],
            tenure_choice_kappa=P.tenure_choice_kappa,alpha_cons=P.alpha_cons,sigma=P.sigma,
            selling_cost=P.psi,financed_share=P.phi[0],annual_depreciation=1-(1-P.delta)**(1/P.da),
            period_depreciation=P.delta,annual_property_tax=P.tau_H/P.da,period_property_tax=P.tau_H,
            income_process=len(P.z_grid),entrant_conversion_factor=P.entrant_conversion_factor,
            adult_entry_birth_to_household_conversion=1/2.1,
            child_benefit_exponent=P.utility_child_benefit_exponent,utility_reference_rent=P.utility_reference_rent,
            child_benefit_curvature=1-P.utility_child_benefit_exponent,
            child_benefit_CRRA_coefficient=P.psi_child*P.utility_child_benefit_exponent,
            h_P=P.hbar_first_child_jump)
        # The pension target is an external restriction, not the unweighted income-array mean.
        reference_values={r['parameter']:float(r['estimate']) for r in csv.DictReader((selected['case']/'parameters.csv').open())}
        expected['pension_to_gross_worker_earnings']=reference_values['pension_to_gross_worker_earnings']
        assert set(values)<=set(expected),(set(values)-set(expected))
        for name,value in values.items():np.testing.assert_allclose(value,expected[name],rtol=1e-12,atol=1e-14,err_msg=name)
        assert values['psi_child']==P.psi_child
        np.testing.assert_allclose(values['child_benefit_CRRA_coefficient'],.86*P.psi_child,rtol=0,atol=1e-15)
        from intergen_eqscale_seq_optimized.parameters import readiness_gate_active
        assert P.sequential_births and not getattr(P,'joint_nested_choice',False)
        assert not readiness_gate_active(P)
        assert P.kappa_fert>0 and P.kappa_fert_continuation>0
        assert P.utility_reference_rent==selected['packet']['parameters'].utility_reference_rent
        assert P.pension==selected['packet']['parameters'].pension
        assert P.tau_pay==selected['packet']['parameters'].tau_pay
        np.testing.assert_array_equal(P.income,selected['packet']['parameters'].income)
        common.write(folder/'fiscal_accounts.json',fiscal)
        cal=rt['primitive'].pf.calendar;policy=packet['evaluation'].policy
        reference_g=selected['packet']['stationary_g_pre'];grid=packet['b_grid'];sd=packet['shared']
        common.table(folder/'common_grid.csv',all_income_grid(packet,reference_g))
        for tenure,label in ((0,'current_renters'),(1,'current_owners')):
            for low in (False,True):
                weights=np.zeros_like(reference_g)
                for j in range(P.J):
                    age=P.age_start+j*P.da
                    if age>34 or not P.A_f_start<=j+1<=P.A_f_end:continue
                    if tenure==0:weights[:,0,:,j,:,0,0]=reference_g[:,0,:,j,:,0,0]
                    else:weights[:,1:,:,j,:,0,0]=reference_g[:,1:,:,j,:,0,0]
                if low:weights[grid>1.]=0
                ev=cal.evaluate_period(policy.price,weights,P,grid,sd,cal.SolveCounter(),supplied_policy=policy)
                unavailable=(policy.fert_probs[...,0]==0)&(policy.fert_probs[...,1]==0)
                unavailable_mass=float(weights[...,0,0][unavailable].sum())
                assert unavailable_mass<1e-12 and ev.feasibility_projection_mass<1e-6
                budget=rt['primitive'].dated_budget(ev,P,sd,grid,float(P.user_cost_rate*policy.price[0]))
                assert budget['budget_excess_mass']<=2e-10
                common_gates.append(dict(case=rec['case'],group=label,wealth_le_one=low,unavailable_mass=unavailable_mass,projection_mass=float(ev.feasibility_projection_mass),budget=budget))
        checks.append(dict(case=rec['case'],checkpoint_sha256=hashlib.sha256((folder/'state.pkl.gz').read_bytes()).hexdigest(),targets=len(fits),parameters=len(params),all_parameter_rows_verified=True,external_restrictions=['pension_to_gross_worker_earnings','adult_entry_birth_to_household_conversion'],standard_figures=len(list((folder/'standard_diagnostics').glob('*.png'))),pension_budget_scaled_residual=fiscal['scaled_pension_budget_residual'],b=float(P.psi_child),reference_rent=float(P.utility_reference_rent)))
        if rec['price_factor']==1:
            baseP=P
            pair=next((r for r in records if r['phase']==rec['phase'] and r['arm']==rec['arm'] and r['price_factor']==1.1),None)
            if pair:
                with gzip.open(out/pair['case']/'state.pkl.gz','rb') as f:shock=pickle.load(f)
                common.table(out/pair['case']/'common_grid.csv',all_income_grid(shock,reference_g))
                for key in vars(baseP):
                    if key.startswith('_'):continue
                    x,y=getattr(baseP,key),getattr(shock['parameters'],key)
                    try:same=np.array_equal(x,y,equal_nan=True)
                    except (TypeError,ValueError):same=str(x)==str(y)
                    if not same:parameter_changes.append(dict(case=rec['case'],field=key))
                bm={v['moment']:v['model'] for v in rec['target_fit']};sm={v['moment']:v['model'] for v in pair['target_fit']}
                sensitivity.append(dict(phase=rec['phase'],arm=rec['arm'],b=rec['b'],base_fertility=rec['fertility'],shock_fertility=pair['fertility'],fertility_change=pair['fertility']-rec['fertility'],fertility_percent_change=100*(pair['fertility']/rec['fertility']-1),fertility_log_elasticity=np.log(pair['fertility']/rec['fertility'])/np.log(1.1),base_birth_age=bm['nchs_mean_age'],shock_birth_age=sm['nchs_mean_age'],birth_age_change=sm['nchs_mean_age']-bm['nchs_mean_age'],base_childlessness=bm['cps_childlessness'],shock_childlessness=sm['cps_childlessness'],base_birth_rooms=bm['first_birth_rooms'],shock_birth_rooms=sm['first_birth_rooms']))
                base={r['group']:r for r in csv.DictReader((folder/'common_states.csv').open())};sh={r['group']:r for r in csv.DictReader((out/pair['case']/'common_states.csv').open())}
                for group,row in base.items():
                    d=dict(phase=rec['phase'],arm=rec['arm'],group=group)
                    for key in ('birth_probability','realized_ownership','realized_housing','realized_nonhousing'):
                        d[key+'_base']=float(row[key]);d[key+'_shock']=float(sh[group][key]);d[key+'_change']=float(sh[group][key])-float(row[key])
                    states.append(d)
                bgrid=list(csv.DictReader((folder/'common_grid.csv').open()));sgrid=list(csv.DictReader((out/pair['case']/'common_grid.csv').open()))
                assert len(bgrid)==len(sgrid)
                for birth in (1,2):
                    for low in (False,True):
                        vals=[(b,s) for b,s in zip(bgrid,sgrid) if int(b['birth_number'])==birth and (not low or float(b['wealth'])<=1)]
                        finite=[(b,s) for b,s in vals if b['gap_over_shock_scale']!='' and s['gap_over_shock_scale']!='']
                        total=sum(float(b['common_mass']) for b,s in vals);mass=sum(float(b['common_mass']) for b,s in finite)
                        grid_summaries.append(dict(phase=rec['phase'],arm=rec['arm'],birth_number=birth,wealth_le_one=low,common_mass=total,finite_gap_mass=mass,weighted_gap_change=sum(float(b['common_mass'])*(float(s['gap_over_shock_scale'])-float(b['gap_over_shock_scale'])) for b,s in finite)/mass if mass else None))
    assert not parameter_changes,parameter_changes
    if parent:assert parent_before==fingerprint(parent),'Original run modified by saved-result audit'
    common.table(out/'price_sensitivity.csv',sensitivity);common.table(out/'common_state_changes.csv',states);common.table(out/'common_grid_gap_changes.csv',grid_summaries)
    history=list(csv.DictReader((out/'solve_history.csv').open()));search=[]
    for arm,lam in driver.ARMS:
        rows=sorted([r for r in history if r['arm']==arm and r['case'].startswith('matchsearch')],key=lambda r:float(r['b']))
        search.append(dict(arm=arm,evaluations=len(rows),fertility_nondecreasing_in_b=all(float(y['fertility'])>=float(x['fertility'])-1e-10 for x,y in zip(rows,rows[1:])),full_audit_scope='Final reported points only; intermediate bracket/search points have native solve and reconstruction checks.'))
    common.write(out/'common_state_gate_checks.json',common_gates)
    numerical_names=('setup','loop_smoke','solve','validate','common_states','main')
    original=Path(driver.__file__).with_name('e5f_housing_fertility_cost_executed.py')
    def function_asts(path):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef)}
    old_functions,new_functions=function_asts(original),function_asts(Path(driver.__file__))
    assert all(old_functions[n]==new_functions[n] for n in numerical_names)
    receipt=dict(status='passed',cases=checks,parameter_changes_under_shock=parameter_changes,search_checks=search,binary_logit_source=inspect.getsourcefile(rt['model'].solve_bellman_full_markov_income),binary_logit_function_sha256=hashlib.sha256(binary_source.encode()).hexdigest(),driver_final_sha256=hashlib.sha256(Path(driver.__file__).read_bytes()).hexdigest(),audit_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),evaluations=completed['evaluations'],closure='Fixed-price PE; no market clearing imposed. Actual pension residual measured; estate and income rules fixed.')
    common.write(out/'final_verification.json',receipt)
    common.write(out/'finalization_checks.json',dict(original_run_unchanged=parent is not None,numerical_function_asts_unchanged=list(numerical_names),all_parameter_rows_verified=True,parameter_external_restrictions=['pension_to_gross_worker_earnings','adult_entry_birth_to_household_conversion']))
    print(json.dumps(dict(status='passed',cases=len(checks),price_sensitivity=sensitivity)),flush=True)

if __name__=='__main__':main()
