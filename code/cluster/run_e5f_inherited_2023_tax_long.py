"""Matched inherited-2023 tax continuations with verified original-queue endpoints.

Standard-library preparation submits two independent compute-node pipelines.
No native scientific kernel is modified. The policy extension changes only the
annual property tax and its user-cost identity, retaining the fixed asset-price
housing supply curve and every native household, fiscal and population gate.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import shlex
import subprocess
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch

BASE=Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a')
PYTHON='/share/apps/anaconda3/2025.06/bin/python'

def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)

def prepare(batch):
    if (batch/'manifest.json').exists(): raise ValueError('Refuse duplicate experiment')
    recovery=BASE/'batches/permanent_2023_profile_20260913_retry2'
    rm=read(recovery/'manifest.json');ro=Path(rm['output'])
    refs=[Path(__file__).resolve(),Path(rm['spec']),ro/'native_2023_snapshot.pkl.gz',ro/'verification.json',ro/'rows.json',Path(rm['endpoint']),Path(rm['endpoint_receipt']),Path(rm['snapshot'])/'latest_completed.json',ro/'continuation_2027.pkl.gz']
    proof=read(ro/'verification.json')
    if proof['status']!='PASS' or proof['observed_year']!=2023 or proof['root_solves']!=0: raise ValueError('Saved 2023 recovery is not verified')
    m=dict(spec=rm['spec'], recovery=str(ro), endpoint=rm['endpoint'], endpoint_receipt=rm['endpoint_receipt'], seed=str(Path(rm['snapshot'])/'latest_completed.json'), output=str(batch), file_sha256={str(p):sha(p) for p in refs},
        annual_taxes=[.01,.02],start_year=2023,periods=100,psi=.09221854783921073,
        max_path_evaluations=16,terminal_max_evaluations=24,terminal_seconds=3600,smoke_seconds=1500,path_seconds=8*3600,total_seconds=10*3600,
        expected_mapping_seconds=1900,expected_hours_per_arm=[4,9],
        closure=dict(preference='Held permanently at fitted final level; no new preference shock',information='Unexpected tax reform in 2023, then known permanently',initial_distribution='Saved one-permanent-preference-shock iteration3 pre-choice2023, unconverged history',population='Original four-vintage birth queue, births/2.1; no immigration or rescaling',fiscal='Equal property-tax rebate to household heads; balanced PAYGO at fixed payroll .179',housing='Fixed calibrated supply curve in asset prices, elasticity .63',production_eligible=False),
        objects=dict(structural_parameters='estimated, retained 2007 calibration',preference_levels='estimated in earlier historical exercise, fixed here',housing_supply_scale='empirically normalized at calibration, unchanged',housing_supply_elasticity='externally fixed .63',birth_to_entry_conversion='externally fixed 1/2.1',entry_distribution='retained calibrated entry distribution',migration='externally fixed zero',annual_property_tax='externally fixed policy values .01/.02',payroll_tax='externally fixed .179',pension='endogenous balanced amount',rebate='endogenous equal-per-head balanced amount',inherited_history='outstanding: not converged',terminal_approach_and_horizon='outstanding until assessed'))
    batch.mkdir(parents=True,exist_ok=True);save(batch/'manifest.json',m);jobs=[]
    for tax in m['annual_taxes']:
        label=f'tax{round(tax*100)}';folder=batch/label;folder.mkdir()
        argv=[PYTHON,str(Path(__file__).resolve()),'--run',str(batch/'manifest.json'),'--tax',str(tax)]
        script=folder/'run.sbatch';script.write_text('\n'.join(['#!/bin/bash',f'#SBATCH --job-name=e5f_2023_{label}_100','#SBATCH --account=torch_pr_570_general','#SBATCH --cpus-per-task=1','#SBATCH --mem=32G','#SBATCH --time=610',f'#SBATCH --output={folder}/slurm_%j.log','set -euo pipefail','export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1',shlex.join(argv),'']))
        r=subprocess.run(['sbatch','--parsable',str(script)],capture_output=True,text=True,check=True)
        receipt=dict(job_id=int(r.stdout.strip().split(';')[0]),annual_tax=tax,script_sha256=sha(script),manifest_sha256=sha(batch/'manifest.json'))
        jobs.append(receipt);save(folder/'submission.json',receipt)
    save(batch/'submission.json',dict(jobs=jobs));print(json.dumps(dict(jobs=jobs)),flush=True)

def policy_old(c,tax):
    import numpy as np
    old=copy.copy(c.old);old.parameters=copy.deepcopy(c.old.parameters);P=old.parameters
    if tax not in (.01,.02): raise ValueError('Only approved 1%/2% comparison')
    if not np.isclose(P.tau_H,.04,atol=1e-15,rtol=0): raise ValueError('Initial tax must be 1% annually')
    if tax==.01:
        return old  # Preserve the calibrated floating-point user cost exactly.
    P.tau_H=float(P.period_years)*tax
    P.user_cost_rate=float(P.R_gross)+float(P.delta)+float(P.tau_H)-1.
    return old

def run(manifest,tax):
    import numpy as np
    m=read(manifest)
    for path,digest in m['file_sha256'].items():
        if sha(path)!=digest: raise ValueError('Pinned input changed: '+path)
    out=Path(m['output'])/f'tax{round(tax*100)}'
    if (out/'startup.json').exists(): raise ValueError('Refuse duplicate run')
    deadline=time.monotonic()+m['total_seconds'];stop=threading.Event()
    save(out/'startup.json',dict(tax=tax,start_year=2023,count=100,manifest_sha256=sha(manifest),phase='loading_verified_state'))
    def heartbeat():
        while not stop.wait(60):
            save(out/'heartbeat.json',dict(remaining_seconds=deadline-time.monotonic()))
            if time.monotonic()>=deadline:
                save(out/'failure.json',dict(error='Ten-hour hard deadline'));os._exit(124)
    threading.Thread(target=heartbeat,daemon=True).start()
    spec=read(m['spec']);sys.path.insert(0,str(Path(spec['batch'])/'source'))
    import run_e5f_original_queue_experiments as runner
    c=runner.load_context(m['spec']);old=policy_old(c,tax)
    import e5f_original_queue_terminal as terminal
    import run_e5f_transition_calibration as fertility
    from e5f_balanced_terminal import _household_checks
    pf=c.joined.pf;ro=Path(m['recovery']);rows=read(ro/'rows.json')
    with gzip.open(ro/'native_2023_snapshot.pkl.gz','rb') as f: snapshot=pickle.load(f)
    saved=snapshot['evaluation'];P23=snapshot['parameters']
    assert snapshot['calendar_year']==2023 and [x['calendar_year'] for x in rows]==[2007,2011,2015,2019,2023]
    assert np.isclose(P23.psi_child,m['psi'],rtol=0,atol=1e-14)
    np.testing.assert_array_equal(snapshot['b_grid'],c.old.b_grid)
    queue=[r['birth_children_topcode_adjusted']/2.1 for r in rows[:4]]
    raw=[r['birth_children']/2.1 for r in rows[:4]]
    inherited=c.rebated.InheritedState(2023,pf.PFInitialState(g_pre=np.asarray(saved.g_pre).copy(),scheduled_entries=queue,scheduled_raw_entries=raw))
    assert abs(queue[0]-rows[4]['effective_mature_entrant_flow_B'])<2e-10
    assert abs(raw[0]-rows[4]['raw_state_scheduled_mature_entrant_flow_B'])<2e-10
    # Reevaluate only aggregation with the saved household policy: no Bellman.
    proof=pf.calendar.evaluate_period(saved.policy.price,inherited.households.g_pre,P23,c.old.b_grid,snapshot['shared'],pf.calendar.SolveCounter(),supply_rule=c.old.supply_rule,supplied_policy=saved.policy)
    np.testing.assert_array_equal(proof.g_current,saved.g_current)
    assert abs(float(proof.demand_by_loc[0])-rows[4]['housing_demand'])<2e-10
    assert abs(float(proof.births)-rows[4]['birth_children'])<2e-10
    # The dated row's rent includes expected capital gains and is not steady-state user cost.
    diag,gates=_household_checks(saved,P23,snapshot['shared'],c.old.b_grid,rows[4]['renter_price'],c.primitive,c.audit)
    if not all(gates.values()): raise ValueError('Inherited dated household checks failed')
    c.driver.save(out/'initial_state_verification.json',dict(passed=True,calendar_year=2023,scheduled_entries=queue,scheduled_raw_entries=raw,native_gates=gates,source_history_converged=False,initial_mass=float(saved.g_pre.sum()),source_sha256=sha(ro/'native_2023_snapshot.pkl.gz')))
    with gzip.open(m['endpoint'],'rb') as f: reference_endpoint=pickle.load(f)
    # Validate baseline supply with the existing validator. For the tax policy,
    # check the intended tax/user-cost change and preserve the identical explicit
    # supply curve instead of silently recalibrating H0 or r_bar under new tax.
    native_validate=terminal._validate_supply
    native_validate(c.old)
    def validate_policy(candidate):
        if candidate.supply_rule is not c.old.supply_rule: raise ValueError('Supply object changed')
        Q=candidate.parameters
        if not np.isclose(Q.tau_H,4*tax,rtol=0,atol=1e-15): raise ValueError('Tax conversion differs')
        if not np.isclose(Q.user_cost_rate,Q.R_gross+Q.delta+Q.tau_H-1,rtol=0,atol=1e-15): raise ValueError('User-cost identity fails')
        for name in ['H0','r_bar','xi_supply']:
            np.testing.assert_array_equal(getattr(Q,name),getattr(c.old.parameters,name))
        # Every other baseline validation is retained on a copy with the two
        # explicitly changed policy coordinates restored to their baseline values.
        baseline=copy.copy(candidate);baseline.parameters=copy.deepcopy(Q)
        baseline.parameters.tau_H=c.old.parameters.tau_H
        baseline.parameters.user_cost_rate=c.old.parameters.user_cost_rate
        native_validate(baseline)
    seed=np.asarray(read(m['seed'])['prices'],float).reshape(3,int(m.get('seed_horizon',100)))
    start=np.asarray(reference_endpoint.coordinates,float).copy()
    if tax==.02:
        start[0]*=c.old.parameters.user_cost_rate/old.parameters.user_cost_rate
        start[2]*=2
    def stage(name): save(out/'phase.json',dict(phase=name,remaining_seconds=deadline-time.monotonic()))
    try:
        with c.queue.original_queue_adapter(),c.cache.policy_cache(pf,max_bytes=12*1024**3),patch.object(terminal,'_validate_supply',validate_policy):
            stage('terminal_equilibrium')
            if tax==.01 or m.get('shared_tax_endpoint'):
                if tax==.01:
                    endpoint=reference_endpoint
                else:
                    with gzip.open(m['shared_tax_endpoint'],'rb') as f:endpoint=pickle.load(f)
                    if not endpoint.verified or not read(m['shared_tax_endpoint_receipt']).get('verified'):raise ValueError('Shared policy endpoint is not verified')
                    if not np.isclose(endpoint.parameters.tau_H,4*tax,rtol=0,atol=1e-15):raise ValueError('Shared endpoint tax differs')
                fresh=terminal._evaluate_trial(old=old,psi=m['psi'],coordinates=np.asarray(endpoint.coordinates),audit=c.audit,deadline=min(deadline,time.monotonic()+3600),trial=1)
                audit=terminal._one_step_audit(fresh,old,m['psi'])
                if not fresh.mapping_valid or audit.get('status')!='passed' or not all(audit['checks'].values()): raise ValueError('Baseline endpoint failed fresh audit')
                np.testing.assert_allclose(fresh.policy.V,endpoint.policy.V,rtol=0,atol=2e-10)
                c.driver.save(out/'endpoint_verification.json',dict(passed=True,fresh_audit=audit,source=m['endpoint'] if tax==.01 else m['shared_tax_endpoint'],tax=tax))
            else:
                controls=dict(c.controls,max_evaluations=m['terminal_max_evaluations'])
                endpoint=terminal.solve_terminal(old=old,psi=m['psi'],audit=c.audit,controls=controls,start=start,deadline=min(deadline,time.monotonic()+m['terminal_seconds']),folder=out/'endpoint')
                if not endpoint.verified: raise RuntimeError('No verified 2% terminal equilibrium')
            boundary=NS(parameters=endpoint.parameters,policy=endpoint.policy,asset_price=endpoint.asset_price)
            # Two dates run the exact same joint root/observer/save pipeline as
            # the full path; a native household-only continuation at fixed initial
            # prices avoids imposing the long-run capital loss in this test.
            stage('exact_loop_smoke')
            Q=copy.deepcopy(old.parameters);Q.psi_child=m['psi'];Q.pension=P23.pension
            Q.property_tax_lump_sum_transfer=P23.property_tax_lump_sum_transfer*(tax/.01)
            from e5f_social_security import bind_social_security_income
            bind_social_security_income(Q,pension_period=P23.pension,payroll_tax=.179)
            shared=pf.calendar.model.precompute_shared(Q,c.old.b_grid);q=rows[4]['asset_price']
            objects=pf.calendar.model.solve_bellman_full_markov_income(np.array([Q.user_cost_rate*q]),np.array([q]),Q,c.old.b_grid,shared)
            smoke_policy=pf.policy_from_objects(objects,q,Q,c.old.b_grid,shared)
            smoke_boundary=NS(parameters=Q,policy=smoke_policy,asset_price=q)
            smoke_guess=np.repeat([q,Q.pension,Q.property_tax_lump_sum_transfer],2).reshape(3,2)
            solve_path(c,old,inherited,smoke_boundary,endpoint,m,out/'smoke',smoke_guess,min(deadline,time.monotonic()+m['smoke_seconds']),2)
            smoke_files=list((out/'smoke/mappings').glob('*/rows.json'))
            if not smoke_files: raise RuntimeError('No valid exact-loop smoke mapping')
            save(out/'smoke_passed.json',dict(passed=True,mappings=len(smoke_files),terminal_household_only=True,full_path_has_verified_ge_endpoint=True,manifest_sha256=sha(manifest)))
            stage('full_100_period_transition')
            # Inherit the saved 2023--2403 numerical guess and extend four dates.
            tail=seed[:,4:]
            if tail.shape[1] not in (96,100):raise ValueError('Expected 100/104-date historical seed')
            padding=100-tail.shape[1]
            guess=np.concatenate([tail,np.repeat(np.asarray(endpoint.coordinates)[:,None],padding,axis=1)],axis=1)
            if tax==.02:
                weight=np.linspace(0,1,100)
                guess[0]*=np.exp(weight*np.log(endpoint.asset_price/reference_endpoint.asset_price))
                guess[1]*=np.exp(weight*np.log(endpoint.coordinates[1]/reference_endpoint.coordinates[1]))
                guess[2]*=(tax/.01)*np.exp(weight*np.log(endpoint.coordinates[2]/(reference_endpoint.coordinates[2]*(tax/.01))))
                if padding:guess[:,-padding:]=np.asarray(endpoint.coordinates)[:,None]
            path_deadline=min(deadline,time.monotonic()+m['path_seconds'])
            resume=None
            for round_number in (1,2):
                result=solve_path(c,old,inherited,boundary,endpoint,m,out/'transition'/f'round_{round_number:02d}',guess,path_deadline,8,resume=resume)
                root=result.root_receipt
                c.driver.save(out/'transition/latest_round.json',dict(round=round_number,root_receipt=root))
                if root.get('finite_horizon_market_fiscal_converged') or time.monotonic()+2100>=path_deadline:break
                best=root.get('best') or root.get('final')
                if not best or not best.get('mapping_valid'):break
                guess=np.asarray(best['prices'],float).reshape(3,100)
                resume=root
            c.driver.save(out/'transition/root_receipt.json',root)
            stage('complete');save(out/'complete.json',dict(completed=True,production_eligible=False))
    except BaseException as exc:
        save(out/'failure.json',dict(error_type=type(exc).__name__,error=str(exc)));raise
    finally: stop.set()


def solve_path(c,old,inherited,boundary,endpoint,m,folder,guess,deadline,max_evaluations,resume=None):
    import numpy as np
    import run_e5f_transition_calibration as fertility
    count=guess.shape[1];folder.mkdir(parents=True,exist_ok=True)
    controls=dict(c.controls)
    for k in ['automatic_fiscal_polish','fiscal_tolerance','fiscal_slope','initial_jacobian']:controls.pop(k,None)
    controls['slope']=controls.pop('market_slope',1.63);controls['max_evaluations']=max_evaluations
    if resume is not None:
        if resume.get('final_jacobian') is not None:controls['initial_jacobian']=resume['final_jacobian']
        if resume.get('final_damping') is not None:controls['damping']=resume['final_damping']
    observations=[];snap={};counter=[0]
    def observer(i,e,P,grid,shared):
        observations.append(dict(period=i,calendar_year=2023+4*i,**fertility.period_fertility_diagnostics(e,P)))
        if i==0:snap.update(parameters=P,b_grid=grid,evaluation=e,shared=shared,supply_rule=old.supply_rule)
        if i%10==0:save(folder/'date_progress.json',dict(period=i,total=count))
    native=c.rebated.evaluate_forecast
    def evaluate(**kwargs):
        observations.clear();snap.clear();path=native(**kwargs);counter[0]+=1
        here=folder/'mappings'/f'mapping_{counter[0]:02d}'
        c.driver.save(here/'rows.json',path.rows);c.driver.save(here/'fertility.json',observations)
        actual=path.person_tail.terminal_state;target=endpoint.state
        gaps=dict(population_relative_gap=float(actual.g_pre.sum()/target.g_pre.sum()-1),distribution_relative_l1=float(np.abs(actual.g_pre-target.g_pre).sum()/target.g_pre.sum()),queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_entries)/target.scheduled_entries-1))),raw_queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_raw_entries)/target.scheduled_raw_entries-1))),terminal_verified=True,horizon_comparison_passed=False)
        c.driver.save(here/'terminal_distance.json',gaps)
        c.driver.save(here/'native_gates.json',dict(mass=path.maximum_mass_accounting_error,policy=path.maximum_policy_reproduction_error,feasibility=path.maximum_feasibility_projection_mass))
        save(folder/'latest_completed_mapping.json',dict(mapping=counter[0],path=str(here)))
        plot(here)
        if count==100 and counter[0]==1:
            from run_e5f_successive_surprises_overnight import standard_graphs
            standard_graphs(snap,NS(path=path),here/'graphs')
        return path
    def progress(row):
        c.driver.save(folder/'latest_completed.json',row)
        if row.get('new_best'):c.driver.save(folder/'best_so_far.json',row)
    with patch.object(c.rebated,'evaluate_forecast',evaluate):
        result=c.rebated.solve_rebated_forecast(inherited=inherited,psi=m['psi'],old_state=old,terminal=boundary,demographic_primitives=None,count=count,initial_prices=guess[0],initial_pensions=guess[1],initial_transfers=guess[2],audit_controls=c.audit,root_controls=controls,deadline_monotonic=deadline,callback=progress,observer=observer)
    c.driver.save(folder/'root_receipt.json',result.root_receipt)
    return result


def plot(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=read(folder/'rows.json');fert=read(folder/'fertility.json');fig,ax=plt.subplots(2,3,figsize=(12,6),constrained_layout=True)
    years=[r['calendar_year'] for r in rows]
    specs=[('asset_price','House price'),('renter_price','Rent'),('adult_population','Household population'),('pension_period_units','Pension'),('equal_transfer_period_units','Equal rebate')]
    for a,(key,title) in zip(ax.flat,specs):a.plot(years,[r[key] for r in rows]);a.set_title(title)
    ax.flat[5].plot(years,[r['housing_demand'] for r in rows],label='Demand');ax.flat[5].plot(years,[r['housing_supply'] for r in rows],ls='--',label='Supply');ax.flat[5].set_title('Housing');ax.flat[5].legend()
    fig.suptitle('2023 tax continuation — price-path evaluation; check equilibrium receipt')
    fig.savefig(folder/'path.png',dpi=150);fig.savefig(folder/'path.pdf');plt.close(fig)


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',type=Path);g.add_argument('--run',type=Path);p.add_argument('--tax',type=float);a=p.parse_args()
    if a.prepare:prepare(a.prepare)
    else:
        try:run(a.run,a.tax)
        except BaseException as exc:
            save(Path(read(a.run)['output'])/f'tax{round(a.tax*100)}'/'failure.json',dict(error_type=type(exc).__name__,error=str(exc)))
            raise
if __name__=='__main__':main()
