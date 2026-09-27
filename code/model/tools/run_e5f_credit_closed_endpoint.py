#!/usr/bin/env python3
"""DRAFT: fixed-parameter closed natural-credit endpoint, no transition.

Matches the pre-September14 closed endpoint in
run_e5f_post2023_no_policy_continuations.solve_closed_stationary_endpoint:
root adjusted_births/(2.1*entry)-1 in price; population=supply/demand.
No H0, psi_child, entry-state distribution or survival changes. Lead historical
and numerical review required before launch; the CLI requires explicit consent.
"""
from __future__ import annotations
import argparse
import copy
import csv
import difflib
import gzip
import inspect
import json
import math
import os
from pathlib import Path
import pickle
import signal
import textwrap
import time
import types
import numpy as np
import run_e5f_credit_benchmark as io

REPLACEMENT = 2.1
ROOT_TOL = 2.5e-5


PRIMITIVE_NAMES = ('adult_entry_clock','entrant_conversion_factor','adult_entry_replacement_fertility',
    'use_age_survival','survival_probs','phi','q','R_gross','user_cost_rate','tau_H','tau_pay',
    'pension','property_tax_lump_sum_transfer','H0','r_bar','xi_supply','psi_child',
    'period_years','population_closure','N_target','entry_shares','z_grid','z_weights','Pi_z')


def primitive_snapshot(P):
    def plain(value):
        if isinstance(value,np.ndarray): return value.tolist()
        if isinstance(value,np.generic): return value.item()
        return value
    result={key:plain(getattr(P,key,None)) for key in PRIMITIVE_NAMES}
    result['annual_q']=(1+float(P.q))**(1/float(P.period_years))-1
    return result


def closed_accounting(entry, births, demand, supply, price):
    values = (entry, births, demand, supply, price)
    if not all(math.isfinite(float(x)) for x in values) or min(entry, demand, supply, price) <= 0 or births < 0:
        raise ValueError('Invalid closed endpoint flow/quantity')
    scale = supply / demand
    return dict(price=float(price), entry_per_normalized_household=float(entry),
        adjusted_birth_children_per_normalized_household=float(births),
        adjusted_births_per_entry=float(births/entry),
        renewal_residual=float(births/(REPLACEMENT*entry)-1),
        population_scale=float(scale), normalized_housing_demand=float(demand),
        absolute_housing_supply=float(supply), absolute_housing_demand=float(scale*demand),
        absolute_market_residual=float(scale*demand-supply),
        absolute_entry=float(scale*entry), absolute_adjusted_births=float(scale*births),
        outside_entry=0., retention=1., replacement=REPLACEMENT)


def bounded_root(evaluate, reference_price, *, low_ratio=.35, high_ratio=3., max_prices=16, tolerance=ROOT_TOL):
    """Small log-price schedule then safeguarded interpolation; no hidden solves."""
    if not 0 < low_ratio < 1 < high_ratio or not 5 <= max_prices <= 16 or not 0 < tolerance <= ROOT_TOL:
        raise ValueError('Invalid root budget/bounds/tolerance')
    cache = {}
    def ev(price):
        key = float(price)
        if key not in cache:
            if len(cache) >= max_prices: raise RuntimeError('Fixed-price solve budget exhausted')
            cache[key] = evaluate(key)
        return cache[key]
    for ratio in (1., low_ratio, math.sqrt(low_ratio), math.sqrt(high_ratio), high_ratio):
        ev(reference_price*ratio)
    def best(): return min(cache.values(), key=lambda r:abs(r['renewal_residual']))
    rows = sorted(cache.values(), key=lambda r:r['price'])
    exact = [r for r in rows if abs(r['renewal_residual']) <= tolerance]
    brackets = [(a,b) for a,b in zip(rows[:-1],rows[1:]) if a['renewal_residual']*b['renewal_residual'] < 0]
    if len(exact) > 1 or (exact and len(brackets)>1) or len(brackets)>1:
        return dict(status='multiple_candidate_roots_not_selected',selected=None,rows=rows)
    if exact:
        return dict(status='root_tolerance_met',selected=exact[0],rows=rows)
    if not brackets:
        return dict(status='no_sign_change_on_declared_schedule',selected=best(),rows=rows)
    a,b = brackets[0]
    while len(cache) < max_prices:
        la,lb = math.log(a['price']),math.log(b['price'])
        fa,fb = a['renewal_residual'],b['renewal_residual']
        fraction = min(.8,max(.2,-fa/(fb-fa)))
        row = ev(math.exp(la+fraction*(lb-la)))
        if abs(row['renewal_residual']) <= tolerance:
            return dict(status='root_tolerance_met',selected=row,rows=sorted(cache.values(),key=lambda r:r['price']))
        if fa*row['renewal_residual'] <= 0: b=row
        else: a=row
    return dict(status='budget_exhausted_no_certified_root',selected=best(),rows=sorted(cache.values(),key=lambda r:r['price']))


def normalized_report_method(runtime, scale, directory):
    """Clone only report units: inherited Habs becomes Habs/S for unit mass.

    Parameters and household solution are not changed. The absolute market
    calculation is saved independently. This source diff requires lead review.
    """
    original = inspect.getsource(type(runtime).evaluate_point)
    source = textwrap.dedent(original)
    old = 'float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]),float(P.xi_supply[0]))'
    new = 'float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0])/closed_population_scale,float(P.xi_supply[0]))'
    if source.count(old) != 1: raise ValueError('Report supply anchor changed')
    source = source.replace(old,new)
    path=directory/'normalized_report.generated.py';path.write_text(source)
    (directory/'normalized_report.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),source.splitlines(True))))
    ns=dict(type(runtime).evaluate_point.__globals__);ns['closed_population_scale']=float(scale)
    exec(compile(source,str(path),'exec'),ns)
    return types.MethodType(ns['evaluate_point'],runtime)


def compare_tables(first, second, tolerance):
    result={}
    for filename,key,value in (('target_fit.csv','moment','model'),('parameters.csv','parameter','estimate')):
        with (first/filename).open() as f:a={r[key]:float(r[value]) for r in csv.DictReader(f)}
        with (second/filename).open() as f:b={r[key]:float(r[value]) for r in csv.DictReader(f)}
        if a.keys()!=b.keys(): raise RuntimeError('Table row set differs')
        differences={k:b[k]-v for k,v in a.items()}
        if len(a)!=(14 if key=='moment' else 31): raise RuntimeError('Incomplete report table')
        limit=tolerance if key=='moment' else 0.
        if max(map(abs,differences.values()))>limit: raise RuntimeError('Reference/repeat table mismatch')
        result[filename]=differences
    return result


def run(plan_path, output, *, inspect_only=False):
    plan=io.read(plan_path)
    for name,pathkey,hashkey in [('runner','runner','runner_sha256'),('helper','helper','helper_sha256'),('contract','contract','contract_sha256')]:
        if io.sha(plan[pathkey])!=plan[hashkey]:raise RuntimeError(name+' pin changed')
    if Path(plan['runner']).resolve()!=Path(__file__).resolve():raise RuntimeError('Wrong runner')
    if io.sha(io.__file__)!=plan['io_sha256']:raise RuntimeError('Shared I/O helper pin changed')
    deadline=min(float(plan['deadline']),time.time()+1200)
    if deadline<=time.time():raise RuntimeError('Expired endpoint budget')
    output.mkdir(parents=True,exist_ok=False)
    def timeout(*_):raise TimeoutError('Closed endpoint hard walltime exhausted')
    signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,deadline-time.time())
    records=[]
    try:
        contract=Path(plan['contract']); c=io.read(contract)
        os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=io.sha(contract)
        os.environ['E5F_LOCAL_EXECUTION_AUTHORIZATION']=c['execution']['authorization_id']
        driver=io.module('closed_credit_frozen_driver',c['files']['driver']['path'])
        c,objective=driver.verify(contract);driver.verify_execution(c)
        case=Path(plan['reference_case']); receipt=io.read(case/'receipt.json')
        if io.sha(case/'receipt.json')!=plan['reference_receipt_sha256']:raise RuntimeError('Reference receipt changed')
        if io.sha(case/'initial_state.pkl.gz')!=receipt['case_checkpoint_sha256']:raise RuntimeError('Reference checkpoint changed')
        work=output/'runtime';work.mkdir()
        runtime,tax,selected,rt,*_=driver.setup(c,objective,receipt['point'],work)
        with gzip.open(case/'initial_state.pkl.gz','rb') as f:reference=pickle.load(f)
        helper=io.module('closed_credit_solvency',plan['helper'])
        model=rt['model']; grid=reference['b_grid']; P0=reference['parameters']
        if (P0.I!=1 or P0.population_closure!='normalized' or not np.isclose(model.housing_demand_normalizer(P0),1.,rtol=0,atol=1e-12)):
            raise ValueError('Expected inherited normalized one-household one-market units')
        psi=float(P0.psi_child); h0=np.asarray(P0.H0).copy()
        primitives=primitive_snapshot(P0);io.write(output/'reference_primitives.json',primitives)
        if inspect_only:
            io.write(output/'complete.json',dict(status='authenticated_reference_inspection_only',solves=0,reference_primitives=primitives));return
        original_estate=runtime.estate.audit
        def estate_report(*a,**kw):
            try:return original_estate(*a,**kw)
            except runtime.estate.EstateFundingShortfall as exc:return exc.audit
        runtime.estate.audit=estate_report
        def fixed(price,label):
            if time.time()>=deadline:raise TimeoutError('Deadline')
            io.write(output/'heartbeat.json',dict(phase='fixed_price',label=label,price=price,epoch=time.time(),completed=len(records),deadline=deadline))
            P=copy.deepcopy(P0);start=time.monotonic()
            sol=model.solve_markov_income_at_prices(np.array([price]),P,grid,verbose=False,fast_stats=False)
            if primitive_snapshot(P)!=primitives:raise RuntimeError('Fixed reference primitives changed')
            if not np.isclose(float(np.sum(sol.g)),1.,rtol=0,atol=1e-9):raise RuntimeError('Nonunit stationary population')
            rt['certify_initial_pension'](sol.g,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
            demand=float(np.asarray(sol.housing_demand).sum())
            supply=float(P.H0[0]*(P.user_cost_rate*price/P.r_bar[0])**P.xi_supply[0])
            row=closed_accounting(float(sol.entry_rate),float(sol.adult_entry_adjusted_birth_children),demand,supply,price)
            row.update(label=label,seconds=time.monotonic()-start,psi_child=psi)
            records.append(row);io.write(output/'price_schedule.json',records)
            io.write(output/'latest_completed.json',row)
            natural=[r for r in records if r['label'].startswith('natural')]
            if natural:io.write(output/'best_so_far.json',min(natural,key=lambda r:abs(r['renewal_residual'])))
            return sol,P,row
        def report(sol,P,row,destination,reference_control=False):
            scale=1. if reference_control else row['population_scale']
            completed=float(rt['chain'].extract_moments(sol,P)['tfr'])
            normalization=dict(status='fixed_psi_closed_endpoint',psi_child=psi,target=2.1,completed_fertility=completed,
                absolute_gap=abs(completed-2.1),stationary_solves=1,stationary_solve_seconds=row['seconds'])
            payroll,rule=runtime.adapter.pension_tax_from_demographics(P)
            def fixed_result(*_):return (sol,P,np.array([row['price']]),row['seconds'],normalization),rule,row
            runtime.normalize=fixed_result
            method=normalized_report_method(runtime,scale,work)
            result=method(tax=tax,objective=objective,selected=selected,runtime=rt,point=receipt['point'],
                output=destination,deadline_epoch=deadline,graphs=True)
            result.update(status='closed_endpoint_draft_verified_unit_mass_checks',closed_accounting=row,
                report_units='normalized household; supply divided by population scale only in reporter',
                absolute_market_audit={'demand':row['absolute_housing_demand'],'supply':row['absolute_housing_supply'],'residual':row['absolute_market_residual']},
                economic_H0_unchanged=True,psi_child_fixed=psi,reference_primitives=primitives,standard_plots_status='17_generated_pending_lead_visual_review')
            parameter_path=destination/'parameters.csv'
            with parameter_path.open() as f: parameter_rows=list(csv.DictReader(f))
            for parameter_row in parameter_rows:
                if parameter_row['parameter'] in receipt['point'] or parameter_row['parameter']=='psi_child':
                    parameter_row['status']='fixed reference estimate; no recalibration'
            with parameter_path.open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=parameter_rows[0]);writer.writeheader();writer.writerows(parameter_rows)
            io.write(destination/'receipt.json',result)
            return result
        baseprice=float(reference['solution'].p_eq[0])
        control= fixed(baseprice,'reference_smoke')
        report(*control,output/'reference_smoke',reference_control=True)
        io.write(output/'reference_comparison.json',compare_tables(case,output/'reference_smoke',1e-5))
        installation=helper.install(model,work,enabled=True);io.write(output/'installation.json',installation)
        rt['accounting'].audit_purchase_accounting=helper.audit_purchase_accounting
        best=[None]
        def evaluate(price):
            sol,P,row=fixed(price,'natural_'+str(sum(r['label'].startswith('natural') for r in records)))
            if best[0] is None or abs(row['renewal_residual'])<abs(best[0][2]['renewal_residual']):
                best[0]=(sol,P,row)
                with gzip.open(output/'best_fixed_price.pkl.gz','wb',compresslevel=1) as f:pickle.dump(dict(solution=sol,parameters=P,b_grid=grid,closed_accounting=row),f,protocol=5)
            return row
        root=bounded_root(evaluate,baseprice,low_ratio=float(plan.get('low_ratio',.35)),high_ratio=float(plan.get('high_ratio',3.)),max_prices=16)
        io.write(output/'root_search.json',root)
        if root['status']!='root_tolerance_met':
            if best[0] is not None:
                report(*best[0],output/'best_diagnostic_not_closed')
            io.write(output/'complete.json',dict(status=root['status'],usable_closed_root=False,root=root,solves=len(records),best_tables='best_diagnostic_not_closed',standard_plots_pending=False));return
        selected_tuple=best[0]
        if selected_tuple[2]['price']!=root['selected']['price']:raise RuntimeError('Selected payload mismatch')
        result=report(*selected_tuple,output/'selected')
        repeated=fixed(root['selected']['price'],'repeat')
        report(*repeated,output/'repeat')
        comparison=compare_tables(output/'selected',output/'repeat',0.)
        io.write(output/'repeat_comparison.json',comparison)
        funded=result['estate_funding']['status']=='funded'
        io.write(output/'complete.json',dict(status='verified_closed_root' if funded else 'demographic_root_estate_funding_failure',usable_closed_root=funded,
            endpoint=root['selected'],solves=len(records),repeat_verified=True,standard_plots_pending=False,not_a_transition=True))
    except BaseException as exc:
        io.write(output/'failure.json',dict(error_type=type(exc).__name__,error=str(exc),solves=len(records),epoch=time.time()))
        raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--lead-reviewed-draft',action='store_true')
    parser.add_argument('--inspect-only',action='store_true',help='Authenticate and export reference primitives; zero model solves')
    args=parser.parse_args()
    if not args.lead_reviewed_draft and not args.inspect_only:parser.error('Lead historical, mathematical and generated-source review is required before execution')
    run(args.plan,args.output,inspect_only=args.inspect_only)
