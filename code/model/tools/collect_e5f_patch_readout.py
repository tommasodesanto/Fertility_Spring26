"""Replay one saved patch forecast at its solved coordinates; extract dated profiles."""
from pathlib import Path
import copy,csv,gzip,hashlib,json,pickle,sys,time
import numpy as np

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
    plan=read(BATCH/'plan.json');c=plan['terminal_template'];td.validate_contract(c);td.verify_sources(c)
    for k in ('initial_checkpoint','initial_summary','initial_contract'):td.verify(c[k]['path'],c[k]['sha256'])
    seed=load(c['initial_checkpoint']['path']);initial_summary=read(c['initial_summary']['path'])
    old=build_approved_initial_state(packet=seed,normalization=initial_summary['normalization'],outside_origin_entry_share=plan['outside_origin_entry_share'],preference_change_2023=0.,fertility_tolerance=5e-4)
    demographics=seed['demographic_seed']
    save(OUT/'initial.json',dict(profile=profile(seed['evaluation'],seed['parameters'],seed['b_grid']),fertility=fert.period_fertility_diagnostics(seed['evaluation'],seed['parameters']),checkpoint_sha256=c['initial_checkpoint']['sha256']))
    hist=[]
    static=ROOT/'batches/stationary_history_patch_20260912'
    for i in range(3):
        s=read(static/f'fit_{i}/summary.json');a=s['selected'];td.verify(a['checkpoint'],a['checkpoint_sha256']);q=load(a['checkpoint'])
        hist.append(dict(window_end=s['decision_year']+4,model=a['model'],data=a['target'],psi=a['psi'],price=a['price'],profile=profile(q['evaluation'],q['parameters'],q['b_grid']),checkpoint_sha256=a['checkpoint_sha256']))
    save(OUT/'stationary_history.json',hist)
    s=read(plan['stationary_restart_2019']['fit_summary']);a=s['selected'];q=load(a['checkpoint']);Q=q['parameters'];sol=q['solution']
    birth=fert.closure.topcode_consistent_renewal_accounting(sol,Q)
    inherited=surprise.InheritedState(2019,pf.PFInitialState(q['stationary_g_pre'].copy(),[float(birth['topcode_adjusted_birth_children'])/2.1]*4,[float(sol.total_births_kfe)/2.1]*4))
    trial=BATCH/'results/arm_0/trial_00_2019';stage=trial/'forecast_6_2';receipt=read(stage/'vintage/2019/root_receipt.json');assert receipt['finite_horizon_market_fiscal_converged']
    t=load(trial/'terminal/terminal_state.pkl.gz');tr=read(trial/'terminal/root_receipt.json');payload=tr['final']['payload']
    terminal=term.BalancedTerminalEndpoint(t['parameters'],t['b_grid'],t['policy'],t['endpoint'],t['social_security'],payload['diagnostics'],payload['household_gates'])
    observed={}
    def observer(i,e,P,grid,shared):
        save(OUT/'latest_date.json',dict(year=2019+4*i))
        if i==1:observed.update(profile=profile(e,P,grid),fertility=fert.period_fertility_diagnostics(e,P),calendar_year=2023)
    final=receipt['final'];start=time.monotonic()
    path=surprise.evaluate_forecast(inherited=inherited,old_state=old,demographics=demographics,prices=final['prices'],pensions=final['fiscal_values'],psi=receipt['psi'],terminal=terminal,observer=observer)
    prior=list(csv.DictReader((stage/'vintage/2019/expected_transition.csv').open()))
    keys=('asset_price','renter_price','housing_demand','housing_supply','owner_rate','birth_children_topcode_adjusted','pension_period_units','payroll_tax_revenue','pension_outlays')
    maximum=max(abs(float(a[k])-float(b[k])) for a,b in zip(prior,path.rows) for k in keys)
    if len(path.rows)!=len(prior) or maximum>2e-10:raise RuntimeError(f'Patch replay failed:{maximum}')
    assert observed['calendar_year']==2023
    agg=observed['profile']['totals'];r=path.rows[1]
    assert abs(agg['rooms']-r['housing_demand'])<2e-10 and abs(agg['owners']/agg['households']-r['owner_rate'])<2e-10
    save(OUT/'model_2023.json',dict(observed,forecast_receipt_sha256=sha(stage/'vintage/2019/root_receipt.json'),finite_converged=True,horizon_verified=False))
    save(OUT/'verification.json',dict(status='PASS',replay_maximum_abs=maximum,seconds=time.monotonic()-start,year=2023,initial_checkpoint_sha256=c['initial_checkpoint']['sha256'],root_receipt_sha256=sha(stage/'vintage/2019/root_receipt.json'),horizon_verified=False,historical_fit_complete=False))
    print(json.dumps(read(OUT/'verification.json')))
if __name__=='__main__':main()
