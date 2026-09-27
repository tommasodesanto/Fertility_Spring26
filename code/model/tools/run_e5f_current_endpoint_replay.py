#!/usr/bin/env python3
"""One native fixed-price replay of the authenticated closed credit endpoint.

No reroot, fertility normalization, or generated accounting/credit installation.
The reported stationary cross-section retains unit household mass; population
scales supply only in reporting and is separately checked against renewal.
"""
from __future__ import annotations
import argparse
import copy
import csv
import gzip
from pathlib import Path
import pickle
import signal
import time

import numpy as np
import e5f_current_transition_runtime as native

REFERENCE = native.ROOT/'output/model/daytime_calibration_20260927/credit_closed_endpoint/run_v1'
PRICE = .8523565697536034
POPULATION = 1.0352351249517946
PINS = {
    'complete.json':'2969db96931b6328721ba3ff44e432a2371483286a2d6bbcd179fbb35b2a3359',
    'selected/receipt.json':'2ab17fae4139c2d1273e34f7efc2f4189f27abf716740c67141badfd71c79bd6',
    'selected/initial_state.pkl.gz':'a12db0c90419b3565a82d4b0d539f334238000a9d2478dfe09f0f0a449d13131',
    'selected/parameters.csv':'a83a4a656fef6379e677718890a55100e4d287635d4d086ce4bafe9e182f3e9b',
    'selected/target_fit.csv':'d2dc7d7b8c8ce21aca54692b5efa4f6b893d28f3e0d428b622c6a1aaa2aaaed4',
}
HELPER_SHA='63fc1e1af117095c0db6277fdd141e7e23a5611aedff83bea3b6e7f817b8bd7f'


def table_rows(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream))


def authenticate_endpoint():
    for relative,expected in PINS.items():
        if native.sha(REFERENCE/relative)!=expected:
            raise RuntimeError('Closed endpoint reference changed: '+relative)
    complete=native.read(REFERENCE/'complete.json')
    if not complete['usable_closed_root'] or not complete['repeat_verified']:
        raise RuntimeError('Reference endpoint is not a repeated closed root')
    endpoint=complete['endpoint']
    if endpoint['price']!=PRICE or endpoint['population_scale']!=POPULATION:
        raise RuntimeError('Endpoint constants differ')
    if len(table_rows(REFERENCE/'selected/parameters.csv'))!=31:
        raise RuntimeError('Incomplete reference parameters')
    return complete


def closed_accounting(sol,P):
    entry=float(sol.entry_rate);births=float(sol.adult_entry_adjusted_birth_children)
    demand=float(np.asarray(sol.housing_demand).sum())
    supply=float(P.H0[0]*(P.user_cost_rate*PRICE/P.r_bar[0])**P.xi_supply[0])
    if not np.isfinite([entry,births,demand,supply]).all() or min(entry,demand,supply)<=0:
        raise RuntimeError('Invalid closed endpoint quantities')
    scale=supply/demand;renewal=births/(2.1*entry)-1
    if abs(renewal)>2.5e-5 or abs(scale-POPULATION)>1e-5:
        raise RuntimeError('Native endpoint does not reproduce closed renewal/population')
    return dict(price=PRICE,population_scale=POPULATION,recomputed_population_scale=scale,
        normalized_housing_demand=demand,absolute_housing_supply=supply,
        absolute_housing_demand=POPULATION*demand,absolute_market_residual=POPULATION*demand-supply,
        entry_per_normalized_household=entry,adjusted_birth_children_per_normalized_household=births,
        renewal_residual=renewal,outside_entry=0.,retention=1.,replacement=2.1)


def unchanged_parameters(old,new):
    changes=[]
    for key,value in vars(old).items():
        if key.startswith('_') or key.startswith('native_') or key=='eq_iter':continue
        other=getattr(new,key,None)
        try: same=np.array_equal(value,other,equal_nan=True)
        except (TypeError,ValueError):same=repr(value)==repr(other)
        if not same:changes.append(key)
    if changes:raise RuntimeError('Fixed endpoint parameters changed: '+repr(changes))


def compare_tables(output):
    result={}
    for filename,key,value,count,tol in [('target_fit.csv','moment','model',14,1e-5),
                                          ('parameters.csv','parameter','estimate',31,1e-12)]:
        a={r[key]:float(r[value]) for r in table_rows(REFERENCE/'selected'/filename)}
        b={r[key]:float(r[value]) for r in table_rows(output/filename)}
        if len(a)!=count or set(a)!=set(b):raise RuntimeError('Incomplete comparison: '+filename)
        gaps={name:b[name]-a[name] for name in a}
        result[filename]=dict(tolerance=tol,differences=gaps,passed=all(abs(x)<=tol for x in gaps.values()))
    return dict(tables=result,passed=all(table['passed'] for table in result.values()))


def run(output,*,seconds=600,inspect_only=False):
    if not 1<=seconds<=600:raise ValueError('At most 600 seconds for this single solve/reporter')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    reference_complete=authenticate_endpoint()
    prepared=native.setup(output/'preparation')
    rt=prepared['runtime'];model=rt['model'];P=copy.deepcopy(prepared['parameters'])
    P.native_solvency_credit=True
    model.validate_native_solvency_mode(P)
    helper_path=native.ROOT/'code/model/tools/e5f_solvency_credit_benchmark.py'
    if native.sha(helper_path)!=HELPER_SHA:raise RuntimeError('Reviewed independent credit audit changed')
    helper=native.load('native_endpoint_credit_audit_only',helper_path)
    rt['accounting'].audit_purchase_accounting=helper.audit_purchase_accounting
    with gzip.open(REFERENCE/'selected/initial_state.pkl.gz','rb') as stream:
        reference=pickle.load(stream)
    unchanged_parameters(reference['parameters'],P)
    evidence=dict(status='prepared_no_solve',endpoint_reference_files=PINS,
        price=PRICE,population_scale=POPULATION,credit_audit_sha256=HELPER_SHA,
        runner_sha256=native.sha(__file__),runtime_sha256=native.sha(native.__file__),
        economic_changes_relative_to_de0093=['Experimental solvency credit replaces artificial borrowing limits'],
        population_closure='Closed endogenous population; no outside entry or retention adjustment',
        economic_changes_relative_to_closed_endpoint=[],
        credit_approximation='Conservative feasible grid-node continuation and native value cutoff',
        all_reference_parameters_fixed=True,normalization=False,generated_installers=False)
    native.write(output/'preparation_endpoint.json',evidence)
    # Save exact driver bytes, independently of mutable working-source metadata.
    (output/'runner_snapshot.py').write_bytes(Path(__file__).read_bytes())
    if inspect_only:return evidence
    deadline=time.time()+seconds
    def timeout(*_):raise TimeoutError('Native endpoint replay hard budget exhausted')
    prior=signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,seconds)
    try:
        native.verify_current_sources(prepared['preparation'])
        start=time.monotonic()
        sol=model.solve_markov_income_at_prices(np.array([PRICE]),P,prepared['selected']['b_grid'],
            verbose=False,fast_stats=False)
        elapsed=time.monotonic()-start
        unchanged_parameters(reference['parameters'],P)
        rt['certify_initial_pension'](sol.g,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
        accounting=closed_accounting(sol,P)
        completed=float(rt['chain'].extract_moments(sol,P)['tfr'])
        normalization=dict(status='fixed_reference_closed_endpoint_no_normalization',psi_child=float(P.psi_child),
            target=2.1,completed_fertility=completed,absolute_gap=abs(completed-2.1),
            stationary_solves=1,stationary_solve_seconds=elapsed)
        objective=prepared['objective']
        _,fiscal_rule=objective.adapter.pension_tax_from_demographics(P)
        objective.normalize=lambda *_:((sol,P,np.array([PRICE]),elapsed,normalization),fiscal_rule,accounting)
        result=objective.evaluate_point(tax=prepared['tax'],objective=prepared['objective_definition'],
            selected=prepared['selected'],runtime=rt,point=prepared['reference_receipt']['point'],
            output=output/'case',deadline_epoch=deadline,graphs=True,report_population_scale=POPULATION)
        comparison=compare_tables(output/'case')
        native.write(output/'case/native_terminal_comparison.json',comparison)
        with gzip.open(output/'case/initial_state.pkl.gz','rb') as stream:candidate=pickle.load(stream)
        native.write(output/'case/native_terminal_arrays.json',native.compare_arrays(reference,candidate))
        native.verify_current_sources(prepared['preparation'])
        if native.sha(helper_path)!=HELPER_SHA or native.sha(__file__)!=evidence['runner_sha256']:
            raise RuntimeError('Endpoint driver/audit changed during replay')
        parameter_path=output/'case/parameters.csv'
        rows=table_rows(parameter_path)
        for row in rows:
            if row['parameter'] in prepared['reference_receipt']['point'] or row['parameter']=='psi_child':
                row['status']='fixed reference estimate; no recalibration or fertility normalization'
            if row['parameter']=='financed_share':
                row['status']='retained parameter; collateral rule replaced by experimental solvency credit'
        with parameter_path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
        result.update(status='native_terminal_tables_pass_array_review_required' if comparison['passed'] else 'native_terminal_tables_failed',
            closed_accounting=accounting,native_source_identity=prepared['preparation'],
            endpoint_replay_contract=evidence,economic_changes_relative_to_closed_endpoint=[],
            all_parameters_fixed=True,standard_plots_status='17 generated; lead visual review required')
        native.write(output/'case/receipt.json',result);native.write(output/'complete.json',result)
        if not comparison['passed']:raise RuntimeError('Native terminal table comparison failed')
        return result
    except BaseException as exc:
        native.write(output/'failure.json',dict(error_type=type(exc).__name__,error=str(exc)))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,prior)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--inspect-only',action='store_true')
    parser.add_argument('--seconds',type=int,default=600)
    args=parser.parse_args();run(args.output,seconds=args.seconds,inspect_only=args.inspect_only)
