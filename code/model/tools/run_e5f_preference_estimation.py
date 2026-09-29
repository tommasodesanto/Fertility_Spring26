#!/usr/bin/env python3
"""Estimate four successive surprises, or one permanent shock fitted to 2023.

Default is plan inspection. Numerical execution is Torch-only and separately
enabled in a pinned plan. Preparation never estimates historical shocks.
"""
from __future__ import annotations
import argparse
import copy
import csv
import gzip
import json
import math
import os
from pathlib import Path
import pickle
import shutil
import signal
import sys
import threading
import time
from types import SimpleNamespace

import run_e5f_preference_transition as inner

SOURCES=inner.SOURCE_NAMES+('e5f_preference_shock_fit.py','run_e5f_preference_estimation.py')
EMPIRICAL=inner.ROOT/'output/model/e5f_matched_pf_20260909a/current_candidate_transition/inputs/empirical_blocks.csv'
ANNUAL=inner.ROOT/'output/model/e5f_matched_pf_20260909a/path_pilot_20260910/fertility_data/annual_fertility_2007_2023.csv'


class CandidateRejected(RuntimeError):
    """A completed but uncertified numerical proposal, never a fertility residual."""


def accepted(condition,message):
    if not condition:raise CandidateRejected(message)


def table(path,rows):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def target_contract(blocks=EMPIRICAL,annual=ANNUAL):
    """Recompute the retained four targets from all sixteen published annual rates."""
    with Path(blocks).open() as stream:rows=list(csv.DictReader(stream))
    with Path(annual).open() as stream:rates=list(csv.DictReader(stream))
    result=[]
    for i,year in enumerate((2007,2011,2015,2019)):
        block=[r for r in rows if int(r['decision_year'])==year]
        sample=[r for r in rates if year<int(r['year'])<=year+4]
        inner.require(len(block)==1 and len(sample)==4 and len({r['year'] for r in sample})==4,
                      'Exactly four unique annual observations per fertility window required')
        b=block[0];value=sum(float(r['period_tfr_births_per_woman']) for r in sample)/4
        inner.require(int(b['birth_year_start'])==year+1 and int(b['birth_year_end'])==year+4 and
            abs(value-float(b['period_tfr_arithmetic_mean']))<=1e-12,'Historical target construction changed')
        inner.require(all(r['status']=='verified_published_final' for r in sample),'Unverified annual data')
        result.append(dict(moment='period_tfr_'+str(year+1)+'_'+str(year+4),decision_year=year,
            period=i,birth_year_start=year+1,birth_year_end=year+4,target=value,
            estimator='equal-weight arithmetic mean of four published NCHS annual TFRs',
            sample='US annual birth-registration rates; published female exposure',
            fixed_effects='not applicable',clustering='not applicable',uncertainty=None,
            sources=sorted({r['source_url'] for r in sample})))
    return dict(schema='retained_nchs_four_windows_v1',rows=result,
        blocks=dict(path=str(blocks),sha256=inner.sha(blocks)),annual=dict(path=str(annual),sha256=inner.sha(annual)),
        model_measurement='period_tfr_topcode_adjusted; sum of four-year birth flows divided by own-age household mass',
        caveat='Retained household-rate analogue, not literal female-exposure TFR; no new age-support target adopted')


def draft_plan(kind='four_successive'):
    inner.require(kind in ('one_permanent','four_successive'),'Unknown shock experiment')
    return dict(schema='block0506_surprise_estimation_v1',execution_enabled=False,kind=kind,
        reference_manifest_sha256=inner.MANIFEST_SHA,housing='fixed_stock',
        credit='saved_reference_unchanged',expectations='current_shock_permanent_until_next_surprise',
        fiscal='fixed_payroll_tax_endogenous_pension',outside_entry=0,retention=1,property_rebate=0,
        initial_level='saved_reference_level',search_bound_ratios=[.01,2.],
        target_contract=None,source_pins={},readiness_receipt=None,
        acceleration=dict(seed_horizon=10,perturbed_date=5,log_step=1e-5,seed_receipt=None),
        horizons=None,horizon_comparison_periods=1 if kind=='four_successive' else 4,horizon_relative_tolerance=1e-3,
        budget=dict(total_seconds=None,candidate_seconds=None,endpoint_seconds=None,
            mapping_seconds=None,path_seconds=None,jacobian_seconds=None,maximum_policy_calls=None),
        fit=dict(max_evaluations=12,log_difference_step=.01,
            fertility_tolerance=.005,max_log_step=.15,damping=.7,max_condition_number=1e8,
            worsening_factor=1.5,reproduction_tolerance=1e-8),
        endpoint=dict(max_evaluations=24,price_bound_ratios=[.05,20.],max_log_step=.15,
            damping=.7,slope=1.,renewal_tolerance=1e-6,reproduction_tolerance=1e-10),
        path=dict(max_evaluations=16,price_bound_ratios=[.05,20.],pension_bound_ratios=[.05,20.],
            cache_max_bytes=2*1024**3,market_tolerance=2e-4,fiscal_tolerance=1e-6,
            market_slope=1.,fiscal_slope=1.,max_log_step=.15,damping=.7,
            final_reproduction_tolerance=1e-10,
            terminal_tolerances={k:1e-3 for k in inner.TERMINAL_KEYS},raw_queue_relative_tolerance=1e-3),
        labels=dict(bounds='numerical search domain, report every boundary hit',
            one_shock_objective='final 2020–2023 window fitted; earlier three validation',
            four_shock_objective='fit 2007/2011/2015/2019 sequentially to the next four-year window; no future shocks in any forecast',
            endpoint='price enforces closed demographic renewal at each candidate psi; psi is never renormalized',
            estate='provisional inherited estate settlement remains outstanding'))


def validate_plan(plan,launching=False):
    inner.require(plan['schema']=='block0506_surprise_estimation_v1' and
                  plan['reference_manifest_sha256']==inner.MANIFEST_SHA,'Wrong estimator reference')
    inner.require(plan['kind'] in ('one_permanent','four_successive') and
        plan['credit']=='saved_reference_unchanged' and plan['housing'] in ('fixed_stock','elastic_reference') and
        plan['expectations']=='current_shock_permanent_until_next_surprise' and
        plan['fiscal']=='fixed_payroll_tax_endogenous_pension' and plan['outside_entry']==0 and
        plan['retention']==1 and plan['property_rebate']==0 and
        plan['initial_level']=='saved_reference_level','Historical/economic contract changed')
    missing=[key for key in ('horizons','target_contract','readiness_receipt','source_pins') if not plan[key]]
    missing += ['budget.'+key for key,value in plan['budget'].items() if value is None]
    if launching:
        inner.require(plan['execution_enabled'] is True,'Shock estimation is disabled')
        inner.require(not missing,'Unresolved numerical launch settings: '+', '.join(missing))
        H=plan['horizons'];inner.require(len(H)>=2 and all(type(h) is int and h>=6 for h in H) and
            all(a<b for a,b in zip(H,H[1:])),'At least two increasing horizons required')
        inner.require(plan['horizon_comparison_periods']==(1 if plan['kind']=='four_successive' else 4) and
                      0<plan['horizon_relative_tolerance']<=1e-3,'Every implemented period requires horizon verification')
        for value in plan['budget'].values():inner.require(math.isfinite(value) and value>0,'Finite positive budgets required')
        lower,upper=plan['search_bound_ratios'];inner.require(0<lower<1<upper and math.isfinite(upper),'Bounds must bracket reference start')
        f,p,e=plan['fit'],plan['path'],plan['endpoint']
        a=plan['acceleration']
        inner.require(type(a['seed_horizon']) is int and a['seed_horizon']>=6 and
            type(a['perturbed_date']) is int and 0<a['perturbed_date']<a['seed_horizon']-1 and
            math.isfinite(a['log_step']) and 0<a['log_step']<.01 and
            (a['seed_receipt'] is None or (set(a['seed_receipt'])=={'path','sha256'} and
             isinstance(a['seed_receipt']['path'],str) and isinstance(a['seed_receipt']['sha256'],str) and
             len(a['seed_receipt']['sha256'])==64)),
            'Explicit measured cross-date initialization or pinned seed receipt required')
        inner.require(type(f['max_evaluations']) is int and f['max_evaluations']>=5 and
                      0<f['log_difference_step']<.25,'Finite derivative and replay budget required')
        inner.require(0<f['fertility_tolerance']<=.005 and 0<=f['reproduction_tolerance']<=1e-8,
                      'Fertility fit/replay gates cannot be loosened')
        inner.require(0<p['market_tolerance']<=2e-4 and 0<p['fiscal_tolerance']<=1e-6 and
            0<=p['final_reproduction_tolerance']<=1e-10 and 0<e['renewal_tolerance']<=1e-6 and
            0<=e['reproduction_tolerance']<=1e-10,'Equilibrium gates cannot be loosened')
        inner.require(0<=p['cache_max_bytes']<=2*1024**3,'Cache cap exceeds verified allocation')
        for section in (p,e):
            inner.require(type(section['max_evaluations']) is int and section['max_evaluations']>=2,
                          'Root budget must reserve fresh replay')
        inner.require(set(p['terminal_tolerances'])==inner.TERMINAL_KEYS and
            all(math.isfinite(v) and 0<v<=1e-3 for v in p['terminal_tolerances'].values()) and
            math.isfinite(p['raw_queue_relative_tolerance']) and 0<p['raw_queue_relative_tolerance']<=1e-3,
            'Complete finite terminal tolerances required')
        stages=4 if plan['kind']=='four_successive' else 1
        # Each of the five residual mappings solves both the backward Bellman
        # problem and the forward distribution, for the full seed horizon.
        # Keep this ceiling even when a pinned receipt avoids the work.
        expected=10*a['seed_horizon']+stages*f['max_evaluations']*(e['max_evaluations']+2+sum(2*h*p['max_evaluations'] for h in H))+8
        inner.require(plan['budget']['maximum_policy_calls']==expected,'Conservative solve count must include all roots and final diagnostics')
    return missing


class NativeEstimator:
    """One authenticated model, bounded endpoint reuse, and fresh dated roots."""
    def __init__(self,plan,output,manifest,packet,evaluator):
        self.plan=plan;self.out=Path(output);self.m=manifest;self.reference=packet;self.rt=evaluator
        self.deadline=time.monotonic()+plan['budget']['total_seconds'];self.candidate_deadline=self.deadline
        self.endpoint_index={};self.endpoint_attempts=0;self.warm={};self.count=0;self.latest=None
        self.inherited=None;self.stage=0;self.year=2007;self.realized=[];self.realized_fertility=[];self.parameters=[]
        self.seed_receipt=None;self.seed_matrices={}

    def guarded(self,seconds,call):
        remaining=min(seconds,self.deadline-time.monotonic(),self.candidate_deadline-time.monotonic())
        if remaining<=0:raise TimeoutError('Estimation/candidate time budget exhausted')
        def timeout(*_):raise TimeoutError('Bounded numerical call exhausted its time budget')
        previous=signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,remaining)
        try:return call()
        finally:signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous)

    def stationary(self,psi,price,output):
        """Fixed-psi solve; only price adjusts renewal, with pension verified on actual mass."""
        import numpy as np
        from e5f_stationary_paygo import bind_initial_balanced_pension
        import e5f_overnight_estate_audit as estate
        rt=self.rt.rt;pf=rt['primitive'].pf;cal=pf.calendar;grid=self.reference['b_grid']
        P=copy.deepcopy(self.reference['parameters']);P.psi_child=float(psi)
        P.native_inherited_distribution_evidence_dir=str(output/'inherited_state_evidence')
        P,pension_rule=bind_initial_balanced_pension(P,payroll_tax=P.tau_pay)
        output.mkdir(parents=True,exist_ok=False)
        sol=self.guarded(self.plan['budget']['mapping_seconds'],lambda:rt['model'].solve_markov_income_at_prices(
            np.array([price]),P,grid,verbose=False,fast_stats=False))
        inner.check_endpoint_primitives(self.reference['parameters'],P)
        inner.require(np.array_equal(sol.b_grid,grid),'Endpoint changed the reference grid')
        fiscal=rt['certify_initial_pension'](sol.g,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        shared=rt['model'].precompute_shared(P,grid);P._fert2_probs=sol.fert2_probs.copy()
        policy=cal.policy_from_solution(sol,np.array([price]),P,grid,shared)
        g,reconstruction=cal.reconstruct_stationary_pre_fertility(sol,policy,P,grid,shared)
        op=pf.transition.operator_gates(sol,policy,g,P,grid,shared);op.update(reconstruction)
        demand=float(np.asarray(sol.housing_demand).sum())
        unit_supply=cal.HousingSupplyRule('fixed-stock',price,demand,0.)
        ev=cal.evaluate_period(np.array([price]),g,P,grid,shared,cal.SolveCounter(),
                               supply_rule=unit_supply,supplied_policy=policy)
        E=float(g[:,:,:,0].sum());account=pf.transition.calendar_topcode_birth_accounting(g,ev.g_post_fertility,float(ev.births),P)
        renewal=float(account['topcode_adjusted_birth_children'])/(2.1*E)-1
        supply=inner.supply_rule(self.reference,pf,self.plan['housing'])
        scale=float(supply.quantity([price])[0])/float(ev.demand_by_loc[0])
        budget=rt['primitive'].dated_budget(ev,P,shared,grid,P.user_cost_rate*price)
        purchase=rt['accounting'].audit_purchase_accounting(ev,P,shared,grid,rt['model'])
        funding=estate.audit(ev,P,grid,next_entrant_cohort=cal.entrant_cohort(np.array([E]),P,grid))
        packet=dict(parameters=P,solution=sol,evaluation=ev,b_grid=grid,shared=shared,
                    stationary_g_pre=g,supply_rule=unit_supply)
        arrays=rt['audit'].policy_array_audit(packet,output)
        gates=dict(unit_mass=abs(float(g.sum())-1)<=1e-9,housing=abs(ev.relative_market_residual)<=2e-4,
            population=math.isfinite(scale) and scale>0,budget=budget['budget_excess_mass']<=2e-10,
            purchase=abs(purchase['maximum_occupied_transaction_wealth_error'])<=1e-9,
            funded=funding['status']=='funded',negative_estates=funding['estate']['totals']['net_negative']<=1e-10,
            occupied_values=arrays['occupied_negative_steps']==0,
            probabilities=all(not v['nonfinite'] and 0<=v['minimum']<=v['maximum']<=1 for v in arrays['probabilities'].values()),
            projection=op['stationary_feasibility_projection_mass']==0,
            mass=abs(op['zero_entry_mass_accounting_residual'])<=2e-8)
        for key in ('stationary_post_fertility_nesting_l1','one_step_constant_path_nesting_l1',
                    'mature_flow_abs_error','birth_flow_abs_error','topcode_adjusted_birth_flow_abs_error'):
            gates[key]=abs(op[key])<=5e-9
        for key,value in purchase.items():
            if key.endswith('violation_mass') or key in ('transaction_outside_grid_mass','saving_outside_grid_mass'):
                gates[key]=abs(float(value))<=2e-10
        record=dict(price=price,psi_child=psi,pension=P.pension,population_scale=scale,
            renewal_residual=renewal,gates=gates,pension_accounts=fiscal,estate=funding,operator=op,
            absolute_housing_demand=scale*float(ev.demand_by_loc[0]),absolute_housing_supply=float(supply.quantity([price])[0]))
        inner.write(output/'point.json',record)
        return packet,record

    def prepare_jacobian(self):
        """Measure once at the unchanged reference, before estimating any shock."""
        import numpy as np
        a=self.plan['acceleration'];P=self.reference['parameters'];H=a['seed_horizon']
        q=float(self.reference['solution'].p_eq[0]);folder=self.out/'jacobian_seed';count=0
        if a['seed_receipt'] is not None:
            pinned=inner.pinned(a['seed_receipt']);receipt=inner.read(pinned)
            self._validate_seed_receipt(receipt)
            stored=inner.pinned(receipt['matrix']);matrix=np.load(stored,allow_pickle=False)
            reconstructed=self._reconstruct_seed(receipt,H)
            inner.require(matrix.shape==(2*H,2*H) and np.isfinite(matrix).all() and
                np.array_equal(matrix,reconstructed),'Pinned seed matrix does not reconstruct its finite measured Jacobian')
            self.seed_receipt=receipt
            inner.write(folder/'reused_receipt.json',dict(status='reused_pinned_measured_seed',
                seed_receipt=dict(path=str(pinned),sha256=a['seed_receipt']['sha256']),
                matrix=receipt['matrix'],reference_manifest_sha256=inner.MANIFEST_SHA,
                psi_child=P.psi_child,scientific_validation=receipt['scientific_validation']))
            return
        previous_deadline=self.candidate_deadline
        self.candidate_deadline=min(self.deadline,time.monotonic()+self.plan['budget']['jacobian_seconds'])
        def evaluate(prices,pensions):
            nonlocal count
            count+=1
            result,record=self.guarded(self.plan['budget']['mapping_seconds'],lambda:inner.mapping(
                self.reference,self.rt,self.reference,dict(price=q,population_scale=1.),prices,pensions,
                np.full(H,P.psi_child),self.plan['housing'],folder/f'map_{count}',
                self.plan['path']['cache_max_bytes']))
            valid=all(record['gates'].values())
            if count==1:
                check=inner.terminal_checks(self.reference,self.rt,self.reference,dict(price=q,population_scale=1.),
                    result,np.full(H,P.psi_child),dict(terminal_tolerances={k:1e-6 for k in inner.TERMINAL_KEYS},
                                                       raw_queue_relative_tolerance=1e-6))
                valid=valid and check['all_checks_pass'] and max(map(abs,record['market_residual']))<=2e-4 and max(map(abs,record['fiscal_residual']))<=1e-6
                inner.write(folder/'baseline_check.json',dict(valid=valid,terminal=check))
            return dict(mapping_valid=valid,market_residual=record['market_residual'],fiscal_residual=record['fiscal_residual'])
        try:
            matrix=inner.measure_jacobian(evaluate,np.full(H,q),np.full(H,P.pension),a['perturbed_date'],
                a['log_step'],folder/'measured',dict(reference_manifest_sha256=inner.MANIFEST_SHA,
                    source_pins=self.plan['source_pins'],housing=self.plan['housing'],closure='fixed_tax',
                    expectations=self.plan['expectations'],psi_child=P.psi_child,
                    reuse='approximate initialization only; nonlinear equilibrium and fresh replay remain mandatory'))
            self.seed_receipt=inner.read(folder/'measured/receipt.json')
            self.seed_receipt['scientific_validation']=dict(
                baseline_mapping_and_terminal_passed=True,
                unchanged_reference=True,baseline_mapping_count=count)
            self.seed_receipt['matrix_reconstruction']=dict(
                exact=True,finite=True,shape=list(matrix.shape))
            inner.write(folder/'measured/receipt.json',self.seed_receipt)
            inner.require(np.array_equal(matrix,self.initial_jacobian(H)),'Lag reconstruction changed measured matrix')
        finally:self.candidate_deadline=previous_deadline

    def _reconstruct_seed(self,receipt,horizon):
        from e5f_four_shock_acceleration import extend_measured_jacobian
        return extend_measured_jacobian(receipt,horizon)

    def _validate_seed_receipt(self,receipt):
        P=self.reference['parameters'];proof=receipt.get('scientific_validation',{})
        inner.require(receipt.get('mapping_count')==5 and
            receipt.get('reference_manifest_sha256')==inner.MANIFEST_SHA and
            receipt.get('source_pins')==self.plan['source_pins'] and
            receipt.get('housing')==self.plan['housing'] and
            receipt.get('closure')=='fixed_tax' and
            receipt.get('expectations')==self.plan['expectations'] and
            receipt.get('psi_child')==P.psi_child,
            'Pinned measured seed provenance or actual reference psi changed')
        inner.require(proof.get('baseline_mapping_and_terminal_passed') is True and
            proof.get('unchanged_reference') is True and proof.get('baseline_mapping_count')==5 and
            receipt.get('matrix_reconstruction',{}).get('exact') is True and
            receipt['matrix_reconstruction'].get('finite') is True and
            receipt['matrix_reconstruction'].get('shape')==[2*receipt.get('horizon',0)]*2,
            'Pinned seed lacks baseline scientific validation evidence')
        matrix=receipt.get('matrix',{})
        inner.require(set(matrix)=={'path','sha256'} and isinstance(matrix['sha256'],str) and len(matrix['sha256'])==64,
            'Pinned seed requires a matrix path and SHA-256')

    def initial_jacobian(self,horizon):
        inner.require(self.seed_receipt is not None,'Measure the current-reference cross-date Jacobian before fitting')
        self._validate_seed_receipt(self.seed_receipt)
        if horizon not in self.seed_matrices:
            self.seed_matrices[horizon]=self._reconstruct_seed(self.seed_receipt,horizon)
        return self.seed_matrices[horizon].copy()

    def endpoint(self,psi):
        import numpy as np
        from e5f_ssj_scaled_step_root import solve_price_path_scaled
        key=float(psi).hex()
        if key in self.endpoint_index:
            item=self.endpoint_index[key]
            with gzip.open(inner.pinned(item['checkpoint']),'rb') as stream:packet=pickle.load(stream)
            return packet,item
        self.endpoint_attempts+=1
        folder=self.out/'endpoints'/('endpoint_'+str(self.endpoint_attempts));folder.mkdir(parents=True,exist_ok=False)
        c=self.plan['endpoint'];q0=float(self.reference['solution'].p_eq[0]);latest={};count=0
        def evaluate(q):
            nonlocal count,latest
            count+=1;packet,record=self.stationary(psi,float(q[0]),folder/f'point_{count:03d}')
            latest=dict(packet=packet,record=record)
            inner.write(folder/'latest_completed.json',record)
            return dict(residual=np.array([record['renewal_residual']]),mapping_valid=all(record['gates'].values()),payload={'point':count})
        def progress(row):
            inner.write(folder/'root_progress.json',row)
            if row.get('new_best'):inner.write(folder/'best_so_far.json',latest['record'])
        root=solve_price_path_scaled(initial_prices=np.array([q0]),evaluate=evaluate,
            project=lambda q:np.clip(q,q0*c['price_bound_ratios'][0],q0*c['price_bound_ratios'][1]),
            slope=c['slope'],market_tolerance=c['renewal_tolerance'],max_log_step=c['max_log_step'],damping=c['damping'],
            max_evaluations=c['max_evaluations'],deadline_monotonic=min(self.candidate_deadline,time.monotonic()+self.plan['budget']['endpoint_seconds']),
            max_condition_number=self.plan['fit']['max_condition_number'],worsening_factor=self.plan['fit']['worsening_factor'],
            final_reproduction_tolerance=c['reproduction_tolerance'],callback=progress)
        inner.write(folder/'root.json',root)
        accepted(root['converged'],'Candidate terminal equilibrium did not converge')
        packet=latest['packet'];r=latest['record'];P=packet['parameters'];q=r['price']
        result,check=self.guarded(self.plan['budget']['mapping_seconds'],lambda:inner.mapping(packet,self.rt,packet,
            dict(price=q,population_scale=1.),np.array([q]),np.array([P.pension]),np.array([psi]),
            'elastic_reference',folder/'one_step',self.plan['path']['cache_max_bytes'],measure_fertility=True))
        terminal=inner.terminal_checks(packet,self.rt,packet,dict(price=q,population_scale=1.),result,[psi],
            dict(terminal_tolerances={k:1e-6 for k in inner.TERMINAL_KEYS},raw_queue_relative_tolerance=1e-6))
        accepted(all(check['gates'].values()) and terminal['all_checks_pass'] and
            max(map(abs,check['market_residual']))<=2e-4 and max(map(abs,check['fiscal_residual']))<=1e-6,
            'Terminal native one-step failed')
        path=folder/'checkpoint.pkl.gz';inner.dump_checkpoint(path,packet)
        item=dict(price=q,psi_child=psi,pension=P.pension,population_scale=r['population_scale'],
            checkpoint=dict(path=str(path),sha256=inner.sha(path)),repeat_verified=True,native_one_step_verified=True,
            source_manifest_sha256=self.m['source_manifest']['sha256'],reference_manifest_sha256=inner.MANIFEST_SHA,
            housing=self.plan['housing'],terminal=terminal)
        inner.write(folder/'receipt.json',item);self.endpoint_index[key]=item
        return packet,item

    def path(self,psi,terminal,endpoint,horizon,folder):
        import numpy as np
        from e5f_four_shock_acceleration import solve_joint_with_acceleration
        from e5f_social_security_root import CandidateDomainError
        p=self.plan['path'];P=self.reference['parameters'];q0=float(self.reference['solution'].p_eq[0]);qT=endpoint['price']
        spec=dict(kind=self.plan['kind'],start_year=self.year,psi=float(psi),
                  expectations='current_shock_permanent_until_next_surprise')
        path=np.full(horizon,psi);latest={};count=0
        warm=self.warm.get(horizon)
        prices=np.linspace(q0,qT,horizon) if warm is None else warm['prices']
        pensions=np.linspace(P.pension,terminal['parameters'].pension,horizon) if warm is None else warm['fiscal_values']
        def evaluate(q,b):
            nonlocal count,latest
            count+=1
            try:self.rt.rt['primitive'].pf.rents_from_asset_prices(q,qT,P)
            except ValueError as exc:raise CandidateDomainError(str(exc)) from exc
            result,record=self.guarded(self.plan['budget']['mapping_seconds'],lambda:inner.mapping(
                self.reference,self.rt,terminal,endpoint,q,b,path,self.plan['housing'],folder/f'map_{count:03d}',
                p['cache_max_bytes'],measure_fertility=True,initial_state=self.inherited,start_year=self.year))
            latest=dict(result=result,record=record,prices=q.copy(),pensions=b.copy())
            inner.write(folder/'latest_completed.json',record)
            inner.dump_checkpoint(folder/'latest_state.pkl.gz',dict(terminal_state=result.terminal_state,
                prices=q,pensions=b,psi_path=path,evaluation=count))
            return dict(mapping_valid=all(record['gates'].values()),market_residual=record['market_residual'],
                        fiscal_residual=record['fiscal_residual'],payload={'evaluation':count})
        def progress(row):
            inner.write(folder/'root_progress.json',row)
            if row.get('new_best'):
                inner.write(folder/'best_so_far.json',latest['record'])
                shutil.copyfile(folder/'latest_state.pkl.gz',folder/'best_state.pkl.gz')
        root=solve_joint_with_acceleration(closure='fixed_tax',initial_prices=prices,initial_fiscal_values=pensions,
            evaluate=evaluate,project_prices=lambda q:np.clip(q,q0*p['price_bound_ratios'][0],q0*p['price_bound_ratios'][1]),
            fiscal_bounds=[P.pension*v for v in p['pension_bound_ratios']],market_tolerance=p['market_tolerance'],
            fiscal_tolerance=p['fiscal_tolerance'],market_slope=p['market_slope'],fiscal_slope=p['fiscal_slope'],
            max_log_step=p['max_log_step'],damping=p['damping'],max_evaluations=p['max_evaluations'],
            deadline_monotonic=min(self.candidate_deadline,time.monotonic()+self.plan['budget']['path_seconds']),
            max_condition_number=self.plan['fit']['max_condition_number'],worsening_factor=self.plan['fit']['worsening_factor'],
            final_reproduction_tolerance=p['final_reproduction_tolerance'],callback=progress,
            initial_jacobian=self.initial_jacobian(horizon) if warm is None else warm['final_jacobian'])
        inner.write(folder/'root.json',root)
        accepted(root['converged'],'Candidate perfect-foresight equilibrium did not converge')
        check=inner.terminal_checks(self.reference,self.rt,terminal,endpoint,latest['result'],path,p)
        receipt=dict(reference_manifest_sha256=inner.MANIFEST_SHA,source_pins=self.plan['source_pins'],housing=self.plan['housing'],
            shock_contract=spec,horizon=horizon,root_and_terminal_pass=bool(root['converged'] and check['all_checks_pass']),
            terminal=check,rows=latest['record']['rows'],fertility=latest['record']['fertility'])
        inner.write(folder/'receipt.json',receipt)
        accepted(receipt['root_and_terminal_pass'],'Candidate path has not approached its terminal equilibrium')
        self.warm[horizon]=dict(prices=root['final']['prices'],fiscal_values=root['final']['fiscal_values'],
                                final_jacobian=np.asarray(root['final_jacobian'],dtype=float).copy())
        return receipt,latest

    def __call__(self,psi):
        psi=float(psi)
        self.count+=1;folder=self.out/f'candidate_{self.count:04d}';folder.mkdir(parents=True,exist_ok=False)
        self.candidate_deadline=min(self.deadline,time.monotonic()+self.plan['budget']['candidate_seconds'])
        inner.write(folder/'proposal.json',dict(psi=psi,kind=self.plan['kind'],start_year=self.year,
                    expectations='current_shock_permanent_until_next_surprise'))
        inner.write(self.out/'heartbeat.json',dict(phase='candidate',candidate=self.count,psi=psi,year=self.year,epoch=time.time()))
        try:
            terminal,endpoint=self.endpoint(psi);previous=None
            periods=self.plan['horizon_comparison_periods']
            for H in self.plan['horizons']:
                receipt,latest=self.path(psi,terminal,endpoint,H,folder/('horizon_'+str(H)))
                if previous is not None:
                    comparison=inner.compare_horizons(previous,receipt,periods,self.plan['horizon_relative_tolerance'])
                    a=[r['period_tfr_topcode_adjusted'] for r in previous['fertility'][:periods]]
                    b=[r['period_tfr_topcode_adjusted'] for r in receipt['fertility'][:periods]]
                    comparison['fertility_absolute_gaps']=[abs(x-y) for x,y in zip(a,b)]
                    comparison['passed']=comparison['passed'] and max(comparison['fertility_absolute_gaps'])<=self.plan['fit']['fertility_tolerance']/5
                    inner.write(folder/('horizon_'+str(H)+'_comparison.json'),comparison)
                    accepted(comparison['passed'],'Fertility/prices/population depend on the terminal horizon')
                previous=receipt
            models=[r['period_tfr_topcode_adjusted'] for r in receipt['fertility'][:periods]]
            target_index=3 if self.plan['kind']=='one_permanent' else self.stage
            model=models[-1];target=self.plan['target_contract']['rows'][target_index]['target']
            summary=dict(certified=True,model=model,payload=dict(candidate=self.count,path=str(folder),
                         models=models,start_year=self.year),endpoint=endpoint,psi=psi,horizon_verified=True,
                         target_index=target_index,gap=model-target,loss_contribution=(model-target)**2)
            inner.write(folder/'complete.json',summary);inner.write(self.out/'latest_completed.json',summary)
            bestpath=self.out/f'best_stage_{self.stage}.json'
            if not bestpath.exists() or summary['loss_contribution']<inner.read(bestpath)['loss_contribution']:
                inner.write(bestpath,summary)
                inner.write(self.out/'best_so_far.json',summary)
            self.latest=dict(terminal=terminal,endpoint=endpoint,latest=latest,psi=psi,summary=summary)
            return summary
        except (CandidateRejected,TimeoutError) as exc:
            inner.write(folder/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),certified=False,psi=psi))
            return dict(certified=False,model=None,payload=dict(candidate=self.count,path=str(folder),error=str(exc)))

    def start_stage(self,stage):
        inner.require(stage==self.stage and self.year==2007+4*stage,'Surprises must be fitted in order')
        return self

    def advance(self,stage,result,*,diagnostics=True):
        """Implement the accepted forecast prefix, before revealing any next shock.

        Its boundary is that vintage's next price and value function. Replacing
        them by a later surprise's price or the stationary endpoint changes
        today's rent, choices and estates, and is therefore forbidden.
        """
        import numpy as np
        from e5f_preference_shock_fit import fit_rows
        last=self.latest;final=result['root']['final'];latest=last['latest']
        inner.require(result['converged'] and stage==self.stage and final is not None and
            final['mapping_valid'] and float(final['prices'][0])==last['psi'] and
            final['payload']['candidate']==last['summary']['payload']['candidate'],
            'Only the freshly certified current-stage fit can advance')
        count=1 if self.plan['kind']=='four_successive' else 4
        inner.require(len(latest['prices'])>count,'Forecast must extend beyond the implemented prefix')
        boundary_q=float(latest['prices'][count]);boundary_V=latest['result'].values[count]
        boundary=dict(evaluation=SimpleNamespace(policy=SimpleNamespace(V=boundary_V)))
        folder=self.out/f'accepted_{self.year}';folder.mkdir(parents=True,exist_ok=False)
        self.candidate_deadline=self.deadline
        replay,record=self.guarded(self.plan['budget']['mapping_seconds'],lambda:inner.mapping(
            self.reference,self.rt,boundary,dict(price=boundary_q),latest['prices'][:count],
            latest['pensions'][:count],np.full(count,last['psi']),self.plan['housing'],folder/'implemented',
            self.plan['path']['cache_max_bytes'],capture=diagnostics,measure_fertility=True,
            initial_state=self.inherited,start_year=self.year))
        gaps=[]
        for actual,expected in zip(record['rows'],latest['record']['rows'][:count]):
            for key,value in expected.items():
                if isinstance(value,(int,float,np.number)):
                    gaps.append(abs(float(actual[key])-float(value))/max(1.,abs(float(value))))
        models=[row['period_tfr_topcode_adjusted'] for row in record['fertility']]
        wanted=last['summary']['payload']['models']
        gates=dict(scientific=all(record['gates'].values()),
            housing=max(map(abs,record['market_residual']))<=self.plan['path']['market_tolerance'],
            fiscal=max(map(abs,record['fiscal_residual']))<=self.plan['path']['fiscal_tolerance'],
            values=np.allclose(replay.values[0],latest['result'].values[0],rtol=0,atol=2e-10),
            forecast_rows=len(record['rows'])==count and bool(gaps) and max(gaps)<=2e-10,
            fertility=len(models)==count and np.max(np.abs(np.array(models)-wanted))<=self.plan['fit']['reproduction_tolerance'])
        inner.write(folder/'replay.json',dict(gates=gates,expected_next_price=boundary_q,
            source_candidate=last['summary']['payload'],start_year=self.year,end_year=self.year+4*count,
            expectations=self.plan['expectations']))
        inner.require(all(gates.values()),'Implemented prefix changed its accepted forecast')
        checkpoint=folder/'inherited_next.pkl.gz'
        inner.dump_checkpoint(checkpoint,dict(year=self.year+4*count,households=replay.terminal_state))
        with gzip.open(checkpoint,'rb') as stream:loaded=pickle.load(stream)
        pf=self.rt.rt['primitive'].pf
        inner.require(np.array_equal(loaded['households'].g_pre,replay.terminal_state.g_pre),'Inherited population changed on save')
        for name in ('scheduled_entries','scheduled_raw_entries'):
            inner.require(np.array_equal(pf.birth_queue_values(getattr(loaded['households'],name)),
                pf.birth_queue_values(getattr(replay.terminal_state,name))),'Inherited birth queue changed on save')
        if diagnostics:
            rendered=self.guarded(max(1.,self.deadline-time.monotonic()),lambda:inner.render_diagnostics(
                record['diagnostic_packets'],folder/'diagnostics',self.rt.rt['audit'],self.m['standard_diagnostic_names']))
            inner.write(folder/'diagnostics.json',rendered)
        inner.write(folder/'fit.json',result)
        inner.write(folder/'checkpoint.json',dict(path=str(checkpoint),sha256=inner.sha(checkpoint),
            year=loaded['year'],both_queues_preserved=True,population_rescaled=False))
        self.realized.extend(record['rows']);self.realized_fertility.extend(models)
        self.parameters.append(dict(parameter='psi_'+str(self.year),**result['parameter']))
        table(self.out/'estimated_path.csv',self.realized)
        table(self.out/'estimated_shocks.csv',self.parameters)
        table(self.out/'fertility_fit.csv',fit_rows(self.plan['target_contract']['rows'],
            self.realized_fertility+[None]*(4-len(self.realized_fertility)),self.plan['kind']))
        self.inherited=copy.deepcopy(replay.terminal_state);self.year+=4*count;self.stage+=1
        self.warm.clear()
        return record


def run(plan,output):
    validate_plan(plan,launching=True)
    inner.require(sys.platform=='linux' and os.environ.get('SLURM_JOB_ID','').isdigit(),'Torch Slurm only')
    inner.require(set(plan['source_pins'])==set(SOURCES),'Complete estimation source pins required')
    for name,digest in plan['source_pins'].items():inner.require(inner.sha(Path(__file__).parent/name)==digest,'Source changed: '+name)
    ready=inner.read(inner.pinned(plan['readiness_receipt']))
    inner.require(ready['status']=='PASS' and ready['estimator_sources']==plan['source_pins'] and
                  ready['native_endpoint_and_fertility_smoke_passed'],'Matching preparation proof required')
    contract=plan['target_contract'];fresh=target_contract(inner.pinned(contract['blocks']),inner.pinned(contract['annual']))
    inner.require(fresh==contract,'Complete target/measurement contract changed')
    output.mkdir(parents=True,exist_ok=False);inner.write(output/'plan.json',plan)
    stop=threading.Event()
    def heartbeat():
        while not stop.wait(60):
            inner.write(output/'process_heartbeat.json',dict(epoch=time.time(),phase='estimation_active'))
    worker=threading.Thread(target=heartbeat,daemon=True);worker.start()
    started=time.monotonic()
    try:
        import numpy as np
        from e5f_preference_shock_fit import fit_one,fit_sequence
        m,packet,evaluator=inner.load_reference(output/'reference')
        native=NativeEstimator(plan,output,m,packet,evaluator)
        native.deadline=started+plan['budget']['total_seconds'];native.candidate_deadline=native.deadline
        native.prepare_jacobian()
        initial=float(packet['parameters'].psi_child)
        bounds=initial*np.array(plan['search_bound_ratios'])
        controls=dict(plan['fit'],total_seconds=max(0.,native.deadline-time.monotonic()))
        if plan['kind']=='four_successive':
            result=fit_sequence(evaluate_factory=native.start_stage,advance=native.advance,
                targets=contract['rows'],initial_level=initial,bounds=bounds,controls=controls,
                callback=lambda stage,row:inner.write(output/'outer_progress.json',dict(stage=stage,**row)))
        else:
            result=fit_one(evaluate=native,target=contract['rows'][3]['target'],initial_level=initial,
                bounds=bounds,controls=controls,callback=lambda row:inner.write(output/'outer_progress.json',row))
            if result['converged']:native.advance(0,result)
        inner.write(output/'fit_result.json',result)
        for name in ('target_fit.csv','parameters.csv'):
            shutil.copyfile(Path(m['local_export'])/name,output/('fixed_reference_'+name))
        for name,digest in plan['source_pins'].items():inner.require(inner.sha(Path(__file__).parent/name)==digest,'Source changed during fit')
        inner.write(output/'complete.json',dict(status=result['status'],converged=result['converged'],
            kind=plan['kind'],reference_manifest_sha256=inner.MANIFEST_SHA,equilibrium_evaluations=result['equilibrium_evaluations'],
            full_fit='fertility_fit.csv',all_estimated_parameters='estimated_shocks.csv',
            standard_diagnostic_review_pending=True,production_eligible=False))
    except BaseException as exc:
        inner.write(output/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),completed_fit=False));raise
    finally:stop.set();worker.join(timeout=1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--plan-sha256');parser.add_argument('--output',type=Path);parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();plan=inner.read(args.plan)
    if not args.execute:
        print(json.dumps(dict(execution_enabled=False,unresolved=validate_plan(plan),estimates='shock values and their terminal equilibrium')));return
    inner.require(args.plan_sha256 and inner.sha(args.plan)==args.plan_sha256 and args.output is not None,'Pinned plan and new output required')
    run(plan,args.output)


if __name__=='__main__':main()
