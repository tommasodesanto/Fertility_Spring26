#!/usr/bin/env python3
"""Pinned permanent-price-fall household/operator diagnostic, not equilibrium.

No stationary distribution, calibration moments, or normalization is computed.
A parent controller must authenticate its supervision receipt before execution.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import json
import os
import pickle
import signal
import time
from pathlib import Path
import run_e5f_due_stayer_matched_check as matched

ROOT = Path(__file__).resolve().parents[3]
THREAD_VARS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
               'NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')


def array_sha(array):
    return hashlib.sha256(array.tobytes(order='C')).hexdigest()


def validate_spec(plan):
    if plan['arms'] != ['baseline', 'due']:
        raise ValueError('Exactly two independent baseline/DUE arms required')
    if plan['price_multiplier'] != .9 or plan['expectations'] != 'permanent_constant_price':
        raise ValueError('This diagnostic is a permanent ten-percent price fall')
    if plan['fiscal_closure'] != 'fixed_tax_fixed_pension_partial_equilibrium':
        raise ValueError('Fixed tax and pension must be explicitly labeled')
    if not (0 < plan['seconds_per_arm'] <= 300 and 0 < plan['total_seconds'] <= 900):
        raise ValueError('Case/total budgets exceed authorized bounds')
    if plan['initial_state'] != 'exact_original_stationary_g_pre':
        raise ValueError('Initial distribution cannot be replaced or projected')
    if plan['target_reporting'] != 'none_dated_distribution':
        raise ValueError('Stationary target reporting prohibited')


def verify(plan):
    validate_spec(plan)
    if plan['runner_sha256'] != matched.sha(__file__):
        raise ValueError('Runner source changed')
    if Path(plan['source_root']).resolve() != ROOT:
        raise ValueError('Wrong source root')
    for relative, digest in plan['source_files'].items():
        if matched.sha(ROOT/relative) != digest:
            raise ValueError('Source changed: '+relative)
        if matched.sha(Path(plan['snapshot_root'])/relative) != digest:
            raise ValueError('Snapshot changed: '+relative)
    if matched.sha(plan['contract']) != plan['contract_sha256']:
        raise ValueError('Contract changed')
    for name, digest in plan['reference_files'].items():
        if matched.sha(Path(plan['reference'])/name) != digest:
            raise ValueError('Reference changed: '+name)


def verify_supervision(plan, now=None):
    now = time.time() if now is None else now
    if plan.get('execution_authorized') is not True:
        raise ValueError('Prepared plan is not execution authorization')
    receipt_path = plan.get('parent_supervisor_receipt')
    if not receipt_path or matched.sha(receipt_path) != plan.get('parent_supervisor_sha256'):
        raise ValueError('Authenticated parent supervision is required')
    receipt = matched.read(receipt_path)
    if (receipt.get('status') != 'supervising' or receipt.get('maximum_active_children') != 1
            or receipt.get('seconds_per_arm') != plan['seconds_per_arm']
            or receipt.get('absolute_end_epoch') != plan['absolute_end_epoch']):
        raise ValueError('Supervision receipt does not match bounded plan')
    start, end = float(receipt['start_epoch']), float(plan['absolute_end_epoch'])
    if not start <= now < end or end-start > plan['total_seconds']:
        raise ValueError('Expired or extended global deadline')
    if int(receipt['pid']) <= 1:
        raise ValueError('Invalid supervisor PID')
    os.kill(int(receipt['pid']), 0)
    return receipt


def prepare(output):
    output=Path(output)
    # Reuse only the zero-solve source snapshot builder; replace its draft spec.
    old=matched.prepare(output)
    plan=dict(schema='due_permanent_price_fall_v1', source_root=str(ROOT),
        runner_sha256=matched.sha(__file__),source_files=old['source_files'],
        snapshot_root=str(output.resolve()/'source_snapshot'),
        contract=old['contract'],contract_sha256=old['contract_sha256'],
        reference=old['reference'],reference_files=old['reference_files'],
        arms=['baseline','due'],price_multiplier=.9,expectations='permanent_constant_price',
        fiscal_closure='fixed_tax_fixed_pension_partial_equilibrium',
        initial_state='exact_original_stationary_g_pre',target_reporting='none_dated_distribution',
        seconds_per_arm=300,total_seconds=900,execution_authorized=False,
        absolute_end_epoch=None,parent_supervisor_receipt=None,parent_supervisor_sha256=None,
        expected_household_solves=2,expected_stationary_distribution_solves=0,
        economic_changes=['Both arms: experimental permanent 10% house-price fall; rent equals unchanged user-cost rate times new price',
          'DUE arm only: grandfather existing-owner debt with separate current-price net-estate solvency at possible death',
          'All 31 original de_0093 parameters, grid, earnings, entry distribution, tax and pension fixed'],
        missing_production_closures=['Housing market clearing','PAYGO pension adjustment','Equilibrium transition and horizon checks'],
        plot_status='No stationary diagnostic packet: dated policy payload retained; supporting plots require a dated semantic review')
    validate_spec(plan)
    matched.write(output/'plan_draft.json',plan)
    return plan


def audit_gates(budget, purchase, estate, mass_residual):
    try:
        matched.require_no_negative_estates(estate)
        estate_ok=True
    except RuntimeError:
        estate_ok=False
    gates={'budget':float(budget['budget_excess_mass'])<=2e-10,
           'mass':abs(float(mass_residual))<=2e-8,
           'estate_funded':estate['status']=='funded',
           'dated_entry':estate['audit_id']=='estate_funded_dated_entry_provisional_net_v1',
           'negative_estates':estate_ok,
           'transaction_wealth':float(purchase['maximum_occupied_transaction_wealth_error'])<=1e-9}
    for key,value in purchase.items():
        if key.endswith('violation_mass') or key in ('transaction_outside_grid_mass','negative_estate_exposure_mass','saving_outside_grid_mass'):
            gates[key]=float(value)<=2e-10
    return gates


def run(plan_path, arm, output):
    for name in THREAD_VARS: os.environ[name]='1'
    plan=matched.read(plan_path);verify(plan);verify_supervision(plan)
    if arm not in plan['arms']:raise ValueError('Unknown arm')
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    deadline=min(time.time()+plan['seconds_per_arm'],float(plan['absolute_end_epoch']))
    def timed_out(*_): raise TimeoutError('Permanent-price-fall arm deadline')
    previous=signal.signal(signal.SIGALRM,timed_out)
    signal.setitimer(signal.ITIMER_REAL,max(.001,deadline-time.time()))
    stage='setup';g0_hash=None
    def progress(name):
        nonlocal stage
        stage=name;matched.write(out/'heartbeat.json',dict(stage=stage,epoch=time.time(),arm=arm,deadline=deadline))
    try:
        import numpy as np
        import e5f_current_transition_runtime as native
        progress('setup')
        prepared=native.setup(out/'preparation',contract=Path(plan['contract']),reference=Path(plan['reference']),fixed_reference_price=True)
        rt=prepared['runtime'];model=rt['model'];pf=rt['primitive'].pf;cal=pf.calendar
        selected=prepared['selected'];P=copy.deepcopy(prepared['parameters']);original=copy.deepcopy(P)
        P.native_due_stayer_credit=(arm=='due')
        P.native_inherited_distribution_evidence_dir=str(out/'support')
        if getattr(P,'native_solvency_credit',False):raise ValueError('Natural credit must remain off')
        grid=selected['b_grid'];shared=model.precompute_shared(P,grid)
        g0=np.asarray(selected['stationary_g_pre']).copy();g0_hash=array_sha(g0)
        price=float(np.asarray(selected['solution'].p_eq)[0])*plan['price_multiplier']
        rent=float(P.user_cost_rate)*price
        # No dated continuation override: finite-age recursion uses the same
        # permanent lower price/rent at every future age, all other inputs fixed.
        progress('permanent_price_household_policy')
        objects=model.solve_bellman_full_markov_income(np.array([rent]),np.array([price]),P,grid,shared)
        policy=pf.policy_from_objects(objects,price,P,grid,shared)
        for key,value in vars(original).items():
            if key.startswith('_') or key in ('native_due_stayer_credit','native_inherited_distribution_evidence_dir','eq_iter'):continue
            try:equal=np.array_equal(value,getattr(P,key),equal_nan=True)
            except (TypeError,ValueError):equal=repr(value)==repr(getattr(P,key))
            if not equal:raise ValueError('Parameter changed: '+key)
        params=matched.rows(Path(plan['reference'])/'parameters.csv')
        if len(params)!=31:raise ValueError('Expected all31 original parameters')
        for row in params:row['status']='Fixed original de_0093 parameter; dated experiment, no recalibration'
        matched.table(out/'parameters.csv',params)
        checkpoint=out/'policy_and_original_state.pkl.gz'
        with gzip.open(checkpoint,'wb',compresslevel=1) as stream:
            pickle.dump(dict(parameters=P,policy=policy,b_grid=grid,shared=shared,g_pre=g0,price=price,rent=rent),stream,protocol=5)
        matched.write(out/'initial_state_identity.json',dict(reference_checkpoint_sha256=plan['reference_files']['initial_state.pkl.gz'],g_pre_sha256=g0_hash,shape=list(g0.shape),mass=float(g0.sum()),unchanged=True))
        progress('strict_original_state_support')
        # This writes its own census on rejection or accepted numerical tails.
        cal._require_exact_inherited_distribution(g0,policy,P,grid)
        if array_sha(g0)!=g0_hash:raise ValueError('Initial state modified')
        progress('dated_current_operator')
        evaluation=cal.evaluate_period(np.array([price]),g0,P,grid,shared,cal.SolveCounter(),supply_rule=selected['supply_rule'],supplied_policy=policy)
        if evaluation.feasibility_projection_mass!=0 or array_sha(evaluation.g_pre)!=g0_hash:
            raise ValueError('Original distribution was projected or changed')
        initial=pf.stationary_initial_state(g0,float(g0[:,:,:,0].sum()),float(selected['evaluation'].births),P,1/2.1)
        accounting=pf.transition.calendar_topcode_birth_accounting(evaluation.g_pre,evaluation.g_post_fertility,float(evaluation.births),P)
        due_entry,next_queue=pf.transition.advance_adult_entry_clock(initial.scheduled_entries,float(accounting['topcode_adjusted_birth_children']),1/2.1,pf.entry_clock_timing(P))
        if abs(float(np.asarray(P.entry_shares).sum())-1.)>1e-12:raise ValueError('Entry shares do not sum to one')
        next_cohort=cal.entrant_cohort(float(due_entry)*np.asarray(P.entry_shares),P,grid)
        nxt,_,deaths,_=pf.transition.advance_sequential_calendar_distribution(evaluation,np.zeros(int(P.I)),P,grid,shared)
        nxt[:,:,:,0]=next_cohort
        mass_residual=float(nxt.sum()-evaluation.g_post_fertility.sum()+deaths-next_cohort.sum())
        if not np.isfinite(nxt).all() or np.any(nxt<0):raise ValueError('Invalid next distribution')
        progress('dated_accounts')
        purchase=native.load('due_pricefall_purchase',ROOT/'code/model/tools/e5f_due_purchase_audit.py')
        estate=native.load('due_pricefall_estate',ROOT/'code/model/tools/e5f_overnight_estate_audit.py')
        if not purchase.SUPPORTS_NATIVE_DUE_STAYER_CREDIT or not estate.SUPPORTS_NATIVE_DUE_STAYER_CREDIT:raise ValueError('Origin-aware audits required')
        budget=rt['primitive'].dated_budget(evaluation,P,shared,grid,rent)
        purchasing=purchase.audit_purchase_accounting(evaluation,P,shared,grid,model)
        try:funding=estate.audit(evaluation,P,grid,next_entrant_cohort=next_cohort)
        except estate.EstateFundingShortfall as exc:funding=exc.audit
        fiscal=pf.social_security.fiscal_accounts(evaluation.g_current,P)
        gates=audit_gates(budget,purchasing,funding,mass_residual)
        result=dict(status='dated_partial_equilibrium_diagnostic' if all(gates.values()) else 'dated_audit_failed',arm=arm,price=price,rent=rent,
            scope='One dated operator from exact original state under permanent lower-price expectations; fixed tax and pension',
            equilibrium_certified=False,initial_g_pre_sha256=g0_hash,parameters_fixed=True,
            housing_market_residual=float(evaluation.relative_market_residual),fiscal_accounts=fiscal,
            budget=budget,purchase=purchasing,estate=funding,gates=gates,births=float(evaluation.births),
            deaths=float(deaths),actual_next_entrants=float(next_cohort.sum()),mass_residual=mass_residual,
            queue_before=pf.birth_queue_values(initial.scheduled_entries).tolist(),queue_after=pf.birth_queue_values(next_queue).tolist(),
            plan_sha256=matched.sha(plan_path),checkpoint_sha256=matched.sha(checkpoint),target_reporting='not_applicable_to_dated_distribution',
            plots='not_generated_stationary_templates_not_semantically_valid_for_this_payload')
        native.write(out/'dated_accounts.json',prepared['tax'].finite_json(cal.jsonable(result)))
        with gzip.open(out/'dated_operator.pkl.gz','wb',compresslevel=1) as stream:
            pickle.dump(dict(evaluation=evaluation,next_g_pre=nxt,next_entrant_cohort=next_cohort,parameters=P,shared=shared,b_grid=grid),stream,protocol=5)
        verify(plan)
        if not all(gates.values()):raise RuntimeError('Dated household/accounting gates failed; see dated_accounts.json')
        matched.write(out/'complete.json',dict(status=result['status'],arm=arm,plan_sha256=matched.sha(plan_path),completed_epoch=time.time(),equilibrium_certified=False))
        return result
    except BaseException as exc:
        evidence=getattr(exc,'audit',None)
        matched.write(out/'failure.json',dict(status='failed_no_gate_relaxation',arm=arm,stage=stage,error_type=type(exc).__name__,error=str(exc),initial_g_pre_sha256=g0_hash,support_evidence=evidence,plan_sha256=matched.sha(plan_path)))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');parser.add_argument('--plan',type=Path);parser.add_argument('--arm',choices=['baseline','due']);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.prepare:prepare(args.output)
    else:
        if not args.plan or not args.arm:parser.error('--plan and --arm are required')
        run(args.plan,args.arm,args.output)
if __name__=='__main__':main()
