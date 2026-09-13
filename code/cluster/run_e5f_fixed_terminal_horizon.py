"""Smoke-gated short transition using an already verified stationary endpoint.

The retained frozen solver supplies all numerical operations. This wrapper only
selects the horizon, records fertility and population diagnostics, and enforces
the new author-approved time budget. It never overwrites the carried population.
"""
import argparse
from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import threading
import time
from types import SimpleNamespace as NS
from unittest.mock import patch


def read(path):
    return json.loads(Path(path).read_text())


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', type=Path, required=True)
    args = p.parse_args()
    m = read(args.manifest)
    for path, digest in m['file_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('Changed short-horizon input: '+path)
    spec = read(m['spec'])
    sys.path.insert(0, str(Path(spec['batch'])/'source'))
    import run_e5f_original_queue_experiments as runner
    import numpy as np
    c = runner.load_context(m['spec'])
    c.spec_path = Path(m['spec'])
    smoke = read(spec['smoke_summary'])
    if smoke.get('status') != 'passed' or smoke['spec_sha256'] != c.driver.sha(m['spec']):
        raise ValueError('Original native smoke is missing')
    endpoint_receipt = read(m['endpoint_receipt'])
    with gzip.open(m['endpoint_pickle'], 'rb') as f:
        endpoint = pickle.load(f)
    if (not endpoint.verified or not endpoint_receipt['verified']
            or endpoint.receipt != endpoint_receipt
            or endpoint.parameters.psi_child != spec['permanent_psi']):
        raise ValueError('Verified endpoint does not match its receipt and permanent shock')
    count = m['count']
    if count != 10:
        raise ValueError('This authorized diagnostic is exactly ten periods')
    out = Path(m['output'])
    if (out/'experiment_contract.json').exists():
        raise ValueError('Refusing to overwrite a started experiment')
    end = min(time.time()+m['seconds'], spec['absolute_deadline_unix'])
    deadline = time.monotonic()+end-time.time()
    c.driver.save(out/'experiment_contract.json', dict(manifest_sha256=c.driver.sha(args.manifest),
        count=count, years=count*4, seconds=m['seconds'], absolute_deadline_unix=end,
        population_law=spec['population_law'], no_immigration=True,
        psi=spec['permanent_psi'], production_eligible=False))
    stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            c.driver.save(out/'controller_heartbeat.json', dict(remaining_seconds=end-time.time()))
            if time.time() >= end:
                c.driver.save(out/'controller_failure.json', dict(error='Short-test deadline'))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    import run_e5f_transition_calibration as fertility
    calendar = c.primitive.calendar
    def stationary_reference(end_point):
        P = end_point.parameters
        shared = calendar.model.precompute_shared(P, c.old.b_grid)
        e = calendar.evaluate_period(np.array([end_point.asset_price]), end_point.state.g_pre,
            P, c.old.b_grid, shared, calendar.SolveCounter(), supply_rule=c.old.supply_rule,
            supplied_policy=end_point.policy)
        return dict(fertility=fertility.period_fertility_diagnostics(e,P),
                    quantities=runner.reference(c,e,P))
    initial = NS(parameters=c.old.parameters, policy=c.old.policy,
        asset_price=float(c.old.policy.price[0]), state=c.old.initial_state,
        coordinates=runner.initial_coordinates(c,0))
    references = dict(initial=stationary_reference(initial), terminal=stationary_reference(endpoint),
        source_hashes=m['file_sha256'], stationary_endpoint_verified=True)

    def distance(path, target):
        actual = path.person_tail.terminal_state
        mass = float(target.g_pre.sum())
        return dict(distribution_relative_l1=float(np.abs(actual.g_pre-target.g_pre).sum())/mass,
            population_relative_gap=abs(float(actual.g_pre.sum())/mass-1),
            queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_entries)
                /np.asarray(target.scheduled_entries)-1))),
            raw_queue_relative_max=float(np.max(np.abs(np.asarray(actual.scheduled_raw_entries)
                /np.asarray(target.scheduled_raw_entries)-1))), production_eligible=False)

    @contextmanager
    def capture(folder, target):
        native = c.rebated.evaluate_forecast
        c.driver.save(folder/'stationary_reference.json', references)
        def evaluate(**kwargs):
            observations = []
            previous_observer = kwargs.get('observer')
            def observe(i,e,P,grid,shared):
                if previous_observer is not None:
                    previous_observer(i,e,P,grid,shared)
                observations.append(dict(period=i,calendar_year=2007+4*i,
                    **fertility.period_fertility_diagnostics(e,P)))
            kwargs['observer'] = observe
            result = native(**kwargs)
            case = folder/'last_sweep'
            c.driver.save(case/'rows.json', result.rows)
            c.driver.save(case/'fertility.json', observations)
            c.driver.save(case/'stationary_reference.json', references)
            c.driver.save(case/'terminal_distance.json', distance(result,target))
            c.driver.save(case/'root_receipt.json', dict(converged=False,
                status='Last full mapping only; equilibrium acceptance not established'))
            return result
        with patch.object(c.rebated,'evaluate_forecast',evaluate):
            yield

    try:
        with c.queue.original_queue_adapter(), runner.original_receipts(c), \
                c.cache.policy_cache(c.joined.pf,max_bytes=6*1024**3):
            with capture(out/'smoke', initial.state):
                test = runner.fixed_terminal_path(c,initial,out/'smoke',count,
                    c.old.parameters.psi_child,min(deadline,time.monotonic()+600))
            if test.next_state is None or not test.root_receipt['finite_horizon_market_fiscal_converged']:
                raise RuntimeError('Native ten-period no-shock root/replay failed')
            gaps = distance(test.path,initial.state)
            if any(gaps[k]>1e-5 for k in ['distribution_relative_l1','population_relative_gap',
                                          'queue_relative_max','raw_queue_relative_max']):
                raise RuntimeError('Ten-period no-shock stationary drift failed')
            c.driver.save(out/'smoke_summary.json',dict(passed=True,**gaps))
            with capture(out/'transition',endpoint.state):
                result = runner.fixed_terminal_path(c,endpoint,out/'transition',count,
                    spec['permanent_psi'],deadline)
            c.driver.save(out/'controller_complete.json',dict(
                finite_root_converged=bool(result.root_receipt['finite_horizon_market_fiscal_converged']),
                stationary_endpoint_verified=True, horizon_verified=False, production_eligible=False))
        c.driver.verify_pins(m['file_sha256'])
    except BaseException as exc:
        c.driver.save(out/'controller_failure.json',dict(error_type=type(exc).__name__,error=str(exc)))
        raise
    finally:
        stop.set()


if __name__ == '__main__':
    main()
