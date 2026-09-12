"""Inspect carried households under fitted stationary policies; never certify GE."""
from pathlib import Path
from dataclasses import replace
import argparse,gzip,json,pickle,sys,time
import numpy as np
from run_e5f_stationary_history_patch import read,save,sha

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    p=read(a.plan);root=Path(p['source_root']);batch=Path(p['batch']);out=batch/'carry';out.mkdir(exist_ok=False)
    for path,pin in p['file_sha256'].items():
        if sha(path)!=pin:raise ValueError('Pinned input changed: '+path)
    sys.path[:0]=[str(root/'code/model/tools'),str(root/'code/model')]
    import e5f_balanced_history as balanced
    joined,primitive,checks,rents=balanced._runtime();pf=joined.pf
    import run_e5f_transition_calibration as fertility
    from e5f_approved_initial_state import build_approved_initial_state
    from e5f_social_security import fiscal_accounts
    from e5f_stationary_paygo import certify_initial_pension
    c=p['terminal_template']
    with gzip.open(c['initial_checkpoint']['path'],'rb') as f:seed=pickle.load(f)
    old=build_approved_initial_state(packet=seed,normalization=read(c['initial_summary']['path'])['normalization'],
        outside_origin_entry_share=p['outside_origin_entry_share'],preference_change_2023=0.,fertility_tolerance=5e-4)
    inherited=old.initial_state;rows=[];start=time.monotonic()
    try:
        for i,year in enumerate((2007,2011,2015,2019)):
            s=read(batch/f'fit_{i}/summary.json');assert s['status']=='stationary_fertility_target_matched'
            selected=s['selected'];path=Path(selected['checkpoint']);assert sha(path)==selected['checkpoint_sha256']
            with gzip.open(path,'rb') as f:stationary=pickle.load(f)
            Q=stationary['parameters'];price=float(stationary['solution'].p_eq[0]);snapshot={}
            def observer(index,e,P,grid,shared):
                snapshot.update(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=old.supply_rule)
            conditioning=replace(old.historical_conditioning,start_year=year,next_age_targets={1:year+4},observer=observer)
            save(out/'heartbeat.json',dict(phase='one_period_stationary_policy_carry',year=year,elapsed=time.monotonic()-start))
            result=pf.evaluate_path_at_prices(prices=[price],psi_path=[Q.psi_child],terminal_price=price,
                terminal_V=stationary['evaluation'].policy.V,base_parameters=Q,b_grid=old.b_grid,initial_state=inherited,
                supply_rule=old.supply_rule,birth_to_entry_conversion=1/2.1,historical_conditioning=conditioning,
                transfer_path=[0.],pension_path=[Q.pension],payroll_tax_path=[.179])
            e=snapshot['evaluation'];P=snapshot['parameters'];f=fiscal_accounts(e.g_current,P)
            measured=fertility.period_fertility_diagnostics(e,P)
            row=dict(year=year,psi=Q.psi_child,target=selected['target'],stationary_model=selected['model'],
                carried_model=measured['period_tfr_topcode_adjusted'],path=result.rows[0],fiscal=f,
                household_budget=primitive.dated_budget(e,P,snapshot['shared'],old.b_grid,float(P.user_cost_rate*price)),
                mass_error=result.maximum_mass_accounting_error,replay_error=result.maximum_policy_reproduction_error,
                approximation='Stationary future values/prices/pension; optimal policies applied to carried households. No current price/pension reclearing.',
                equilibrium_certified=False)
            rows.append(row);save(out/'latest_completed.json',row);save(out/'carried_history.json',rows)
            inherited=result.terminal_state
            with gzip.open(out/f'state_{year+4}.pkl.gz','wb',compresslevel=1) as stream:pickle.dump(inherited,stream,protocol=5)
        save(out/'summary.json',dict(status='diagnostic_carry_completed',rows=rows,market_fiscal_reclearing_required=True,
            fitted_historical_equilibrium=False,post2023_pf_started=False,seconds=time.monotonic()-start))
    except Exception as exc:save(out/'failure.json',dict(error=str(exc),type=type(exc).__name__,completed_years=[r['year'] for r in rows]));raise

if __name__=='__main__':main()
