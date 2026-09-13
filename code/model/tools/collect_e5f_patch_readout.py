"""Replay one saved patch forecast at its solved coordinates; extract dated profiles."""
from pathlib import Path
import copy,csv,gzip,hashlib,json,pickle,sys,time
import numpy as np
import argparse

ROOT=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
BATCH=ROOT/'batches/stationary_2019_pf_test_20260912'
OUT=ROOT/'batches/patch_readout_20260912'
sys.path[:0]=[str(BATCH),str(ROOT/'batches/matched_continuation_20260912'),str(ROOT/'code/model/tools'),str(ROOT/'code/model')]
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    with gzip.open(p,'rb') as f:return pickle.load(f)
def save(p,v):Path(p).write_text(json.dumps(v,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x),indent=2)+'\n')

def profile(e,P,grid):
    g=e.g_current;rows=[];large=[]
    for j in range(P.J):
        x=g[:,:,:,j,:,:,:];mass=float(x.sum());td=x.sum(axis=(0,2,3,4))
        rooms=float(np.sum(x[:,0]*e.policy.hR_pol[:,0,:,j,:,:,:]));cap=float(np.sum(x[:,0]*np.minimum(e.policy.hR_pol[:,0,:,j,:,:,:],9)))
        for t,h in enumerate(P.H_own,1):rooms+=float(x[:,t].sum())*h;cap+=float(x[:,t].sum())*min(h,9)
        slots=[t for t,h in enumerate(P.H_own,1) if h>=6]
        large.append(dict(age=float(P.age_start+j*P.da),age_width=float(P.da),without_children=float(td[slots,0].sum()),with_children=float(td[slots,1:].sum())))
        rows.append(dict(age=float(P.age_start+j*P.da),age_width=float(P.da),households=mass,owners=float(x[:,1:].sum()),rooms=rooms,capped_rooms=cap,with_children=float(x[:,:,:,:,:,1:].sum()),
            consumption=float(np.sum(x*e.policy.c_pol[:,:,:,j,:,:,:])),wealth=float(np.sum(x*grid[:,None,None,None,None,None]))))
    counts=g.sum(axis=(0,1,2,3,4,6));totals={k:sum(r[k] for r in rows) for k in ('households','owners','rooms','capped_rooms','with_children','consumption','wealth')}
    assert abs(sum(counts)-totals['households'])<1e-10
    return dict(rows=rows,large_owner_age_cells=large,number_children_mass=counts,totals=totals,child_observer='Model dependent count; ACS resident own minor children are an approximation')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--batch',type=Path,default=BATCH,help='Pinned source batch (default: prior patch batch).')
    ap.add_argument('--out',type=Path,default=OUT,help='Collector output directory.')
    ap.add_argument('--forecast-stage',type=Path,default=None,help='Selected forecast_6_* directory; defaults to the prior forecast.')
    ap.add_argument('--plan',type=Path,default=None,help='Plan carrying the approved initial and stationary source receipts.')
    ap.add_argument('--inherited-checkpoint',type=Path,help='Opt-in actual carried2019state, instead of the stationary patch restart.')
    ap.add_argument('--inherited-sha256',help='Required pin for the carried state.')
    ap.add_argument('--all-dates',action='store_true',help='Also observe fertility stocks and capped housing at every saved forecast date.')
    args=ap.parse_args()
    if args.all_dates and not args.inherited_checkpoint:
        raise ValueError('All-date observations require the actual carried state')
    batch=args.batch; out=args.out
    out.mkdir(parents=True,exist_ok=True)
    forecast_stage=args.forecast_stage or (batch/'results/arm_0/trial_00_2019/forecast_6_2')
    plan_path=args.plan or (batch/'plan.json')
    import run_e5f_candidate_terminal as td
    import run_e5f_transition_calibration as fert
    import e5f_successive_surprises as surprise
    import e5f_balanced_terminal as term
    from e5f_approved_initial_state import build_approved_initial_state
    from run_e5f_matched_pf_history import pf
    import run_e5f_matched_pf_smoke as primitive
    primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility=primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution=primitive.pf.transition.advance_sequential_calendar_distribution
    plan=read(plan_path);c=plan['terminal_template'];td.validate_contract(c);td.verify_sources(c)
    for k in ('initial_checkpoint','initial_summary','initial_contract'):td.verify(c[k]['path'],c[k]['sha256'])
    seed=load(c['initial_checkpoint']['path']);initial_summary=read(c['initial_summary']['path'])
    old=build_approved_initial_state(packet=seed,normalization=initial_summary['normalization'],outside_origin_entry_share=plan['outside_origin_entry_share'],preference_change_2023=0.,fertility_tolerance=5e-4)
    demographics=seed['demographic_seed']
    save(out/'initial.json',dict(profile=profile(seed['evaluation'],seed['parameters'],seed['b_grid']),fertility=fert.period_fertility_diagnostics(seed['evaluation'],seed['parameters']),checkpoint_sha256=c['initial_checkpoint']['sha256']))
    if args.inherited_checkpoint:
        if not args.inherited_sha256:raise ValueError('Carried state requires an explicit hash')
        td.verify(args.inherited_checkpoint,args.inherited_sha256);inherited=load(args.inherited_checkpoint)
        if inherited.year!=2019:raise ValueError('Collector requires the actual2019inherited state')
        if plan['resume_fitted_prefix']['checkpoint']!=dict(path=str(args.inherited_checkpoint),sha256=args.inherited_sha256):
            raise ValueError('Carried state differs from solved forecast plan')
    else:
        hist=[]
        static=ROOT/'batches/stationary_history_patch_20260912'
        for i in range(3):
            s=read(static/f'fit_{i}/summary.json');a=s['selected'];td.verify(a['checkpoint'],a['checkpoint_sha256']);q=load(a['checkpoint'])
            hist.append(dict(window_end=s['decision_year']+4,model=a['model'],data=a['target'],psi=a['psi'],price=a['price'],profile=profile(q['evaluation'],q['parameters'],q['b_grid']),checkpoint_sha256=a['checkpoint_sha256']))
        save(out/'stationary_history.json',hist)
        s=read(plan['stationary_restart_2019']['fit_summary']);a=s['selected'];q=load(a['checkpoint']);Q=q['parameters'];sol=q['solution']
        birth=fert.closure.topcode_consistent_renewal_accounting(sol,Q)
        inherited=surprise.InheritedState(2019,pf.PFInitialState(q['stationary_g_pre'].copy(),[float(birth['topcode_adjusted_birth_children'])/2.1]*4,[float(sol.total_births_kfe)/2.1]*4))
    stage=forecast_stage; trial=stage.parent; receipt=read(stage/'vintage/2019/root_receipt.json');assert receipt['finite_horizon_market_fiscal_converged']
    t=load(trial/'terminal/terminal_state.pkl.gz');tr=read(trial/'terminal/root_receipt.json');payload=tr['final']['payload']
    terminal=term.BalancedTerminalEndpoint(t['parameters'],t['b_grid'],t['policy'],t['endpoint'],t['social_security'],payload['diagnostics'],payload['household_gates'])
    observed={};branch=None;measurement_errors={};observed_dates=[]
    def observer(i,e,P,grid,shared):
        nonlocal branch
        save(out/'latest_date.json',dict(year=2019+4*i))
        if args.all_dates:
            from e5f_initial_fertility_observer import observe_initial_fertility
            dated=dict(calendar_year=2019+4*i,profile=profile(e,P,grid),
                fertility=fert.period_fertility_diagnostics(e,P),
                fertility_stock_timing=observe_initial_fertility(e,P,age_projection='uniform_birth_time'))
            observed_dates.append(dated)
        if args.inherited_checkpoint and i==0:
            try:branch=fert.begin_dated_first_birth_housing_branch(e,P,grid,shared,origin_period=0)
            except (ValueError,RuntimeError) as exc:measurement_errors['first_birth_origin']=str(exc)
        if i==1:
            observed.update(profile=profile(e,P,grid),fertility=fert.period_fertility_diagnostics(e,P),calendar_year=2023)
            if args.inherited_checkpoint:
                from e5f_initial_fertility_observer import observe_initial_fertility
                from e5f_initial_housing_observer import observe_initial_housing_wealth
                from e5f_recent_parent_flow_observer import observe_recent_parent_flow,SNAPSHOT,AGE_PROJECTION
                # Reuse the same cross-sectional arithmetic without asserting
                # stationarity; the same-policy birth-response option stays off.
                for name,call in (
                    ('fertility_stock_timing',lambda:observe_initial_fertility(e,P,age_projection='uniform_birth_time')),
                    ('housing_wealth',lambda:observe_initial_housing_wealth(e,P,grid,shared,diagnostic_enabled=True,
                        age_projection='uniform_within_age_cell',diagnostic_allow_family_proxies=True,include_wealth=True,include_birth_response=False)),
                    ('recent_parent',lambda:observe_recent_parent_flow(e,P,diagnostic_enabled=True,snapshot=SNAPSHOT,
                        age_projection=AGE_PROJECTION,diagnostic_allow_residence_proxy=True)),
                    ('dated_first_birth_rooms',lambda:fert.finish_dated_first_birth_housing_branch(branch,e,P,grid,shared,destination_period=1) if branch is not None else None)):
                    try:observed[name]=call()
                    except (ValueError,RuntimeError) as exc:measurement_errors[name]=str(exc)
                save(out/'moments_2023_partial.json',dict(observed,measurement_errors=measurement_errors,
                    observer_scope='Dated cross-sectional reuse; no stationary certification; first-birth rooms uses2019and2023policies'))
                save(out/'dependent_age_2023.json',dict(calendar_year=2023,
                    g_age_m_post_fertility=e.g_post_fertility.sum(axis=(0,1,2,4,5)),
                    g_age_m_current=e.g_current.sum(axis=(0,1,2,4,5)),
                    survival_probs=P.survival_probs,child_state_mode=P.child_state_mode,
                    Pi_child=getattr(P,'Pi_child',None),child_exit_prob=getattr(P,'child_exit_prob',None),
                    period_years=P.period_years,age_start=P.age_start,da=P.da,births=e.births,
                    source_receipt_sha256=sha(stage/'vintage/2019/root_receipt.json')))
                with gzip.open(out/'state_2023.pkl.gz','wb',compresslevel=1) as stream:
                    pickle.dump(dict(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=old.supply_rule),stream,protocol=5)
    final=receipt['final'];start=time.monotonic()
    path=surprise.evaluate_forecast(inherited=inherited,old_state=old,demographics=demographics,prices=final['prices'],pensions=final['fiscal_values'],psi=receipt['psi'],terminal=terminal,observer=observer)
    prior=list(csv.DictReader((stage/'vintage/2019/expected_transition.csv').open()))
    keys=('asset_price','renter_price','housing_demand','housing_supply','owner_rate','birth_children_topcode_adjusted','pension_period_units','payroll_tax_revenue','pension_outlays')
    maximum=max(abs(float(a[k])-float(b[k])) for a,b in zip(prior,path.rows) for k in keys)
    if len(path.rows)!=len(prior) or maximum>2e-10:raise RuntimeError(f'Patch replay failed:{maximum}')
    assert observed['calendar_year']==2023
    agg=observed['profile']['totals'];r=path.rows[1]
    assert abs(agg['rooms']-r['housing_demand'])<2e-10 and abs(agg['owners']/agg['households']-r['owner_rate'])<2e-10
    if args.all_dates:
        assert len(observed_dates)==len(path.rows)
        for dated,row in zip(observed_dates,path.rows):
            assert dated['calendar_year']==int(row['calendar_year'])
            totals=dated['profile']['totals']
            assert abs(totals['rooms']-row['housing_demand'])<2e-10
            assert 0<=totals['capped_rooms']<=totals['rooms']+2e-10
        dated2023=next(d for d in observed_dates if d['calendar_year']==2023)
        normalize=lambda obj:json.dumps(obj,sort_keys=True,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x))
        assert normalize(dated2023['profile'])==normalize(observed['profile'])
        assert dated2023['fertility_stock_timing']==observed['fertility_stock_timing']
        save(out/'observed_dates.json',dict(dates=observed_dates,
            forecast_receipt_sha256=sha(stage/'vintage/2019/root_receipt.json'),
            replay_maximum_abs=maximum,finite_converged=True,horizon_verified=False))
    # Keep the selected forecast's immutable aggregate rows and receipt beside
    # the derived observer files, making the source packet reproducible offline.
    import shutil
    shutil.copy2(stage/'vintage/2019/expected_transition.csv',out/'expected_transition.csv')
    shutil.copy2(stage/'vintage/2019/root_receipt.json',out/'root_receipt.json')
    save(out/'model_2023.json',dict(observed,forecast_receipt_sha256=sha(stage/'vintage/2019/root_receipt.json'),finite_converged=True,horizon_verified=False))
    if args.inherited_checkpoint:save(out/'measurement_verification.json',dict(errors=measurement_errors,inherited_checkpoint=str(args.inherited_checkpoint),inherited_sha256=args.inherited_sha256,stationary_restart=False))
    save(out/'verification.json',dict(status='PASS',replay_maximum_abs=maximum,seconds=time.monotonic()-start,year=2023,initial_checkpoint_sha256=c['initial_checkpoint']['sha256'],root_receipt_sha256=sha(stage/'vintage/2019/root_receipt.json'),horizon_verified=False,historical_fit_complete=False,forecast_stage=str(stage)))
    print(json.dumps(read(out/'verification.json')))
if __name__=='__main__':main()
