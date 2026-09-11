"""Bounded preference sensitivity with a read-only dated fertility observer.

Every original economic source and numerical gate remains pinned. Household-age
rates are diagnostic analogues, not female TFR: this does not estimate the shock.
"""
import argparse
import copy
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import run_candidate_path as pipeline
import resume_history

CASES = (('delta_m0025', -.025), ('delta_m010', -.10), ('delta_m005', -.05))
ACCEPTED = pipeline.HERE / 'recovery/delta_m005/accepted_history_6.json'


def jsonable(value):
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if hasattr(value, 'tolist'):
        return value.tolist()
    return value


def collect_observer(original, measure, save):
    mappings = []

    def observe(index, evaluation, parameters, grid, shared):
        if original is not None:
            original(index, evaluation, parameters, grid, shared)
        if index == 0:
            mappings.append([])
        if not mappings or index != len(mappings[-1]):
            raise RuntimeError('Dated fertility observer received an incomplete sequence')
        measured = measure(evaluation, parameters)
        mappings[-1].append(dict(calendar_year=2007+4*index, diagnostics=measured))
        save(mappings)
    return observe


def observed_solve(contract, output):
    sys.path[:0] = [str(pipeline.ROOT/'code/model/tools'), str(pipeline.ROOT/'code/model')]
    import e5f_balanced_history as adapter
    import run_e5f_candidate_history as driver
    import run_e5f_transition_calibration as fertility
    original = adapter.solve_balanced_history

    def solve(**kwargs):
        def save(mappings):
            pipeline.save(output/'dated_household_fertility.json', jsonable(dict(
                mappings=mappings,
                denominator='model adult households in each age cell, not female exposure',
                observer_source=pipeline.pin(fertility.__file__),
                reporting_wrapper=pipeline.pin(__file__),
                production_eligible=False, shock_estimated=False)))
        kwargs['observer'] = collect_observer(kwargs.get('observer'),
            fertility.period_fertility_diagnostics, save)
        return original(**kwargs)

    adapter.solve_balanced_history = solve
    try:
        return driver.run(SimpleNamespace(contract=contract,
            contract_sha256=pipeline.sha(contract), output=output))
    finally:
        adapter.solve_balanced_history = original


def run_stage(folder, name, c):
    cp = folder/'contracts'/f'{name}.json'
    pipeline.save(cp, c)
    output = folder/name
    pipeline.save(folder/'latest_stage.json', dict(stage=name, status='running'))
    with (folder/f'{name}.log').open('w') as log:
        run = subprocess.run([sys.executable, __file__, '--contract', str(cp),
            '--output', str(output)], stdout=log, stderr=subprocess.STDOUT,
            timeout=c['seconds']+45)
    summary = pipeline.read(output/'summary.json') if (output/'summary.json').exists() else {}
    pipeline.save(folder/'latest_stage.json', dict(stage=name, exit_code=run.returncode, summary=summary))
    # Continue only a fully checked but unconverged numerical root, never a crash.
    if not summary.get('mapping_replay_verified') or not summary.get('checkpoint_reload_verified'):
        raise RuntimeError(f'{name}: missing successful replay/checkpoint gates')
    if run.returncode not in (0, 1):
        raise RuntimeError(f'{name}: unexpected exit {run.returncode}')
    pipeline.compare(output, folder)
    return output, summary


def prepare(case):
    receipt = pipeline.read(ACCEPTED)
    for pin in receipt.values():
        if pipeline.sha(pin['path']) != pin['sha256']:
            raise ValueError('Accepted smoke receipt hash differs')
    summary = pipeline.read(receipt['summary']['path'])
    if not summary.get('finite_horizon_market_fiscal_converged'):
        raise ValueError('Accepted exact-loop smoke is not converged')
    c = pipeline.read(receipt['contract']['path'])
    c.pop('contract_sha256', None)
    c.pop('continuation')
    c.update(schema='e5f_candidate_history_root_v1')
    c['root_controls']['initial_jacobian'] = None
    root = pipeline.read(receipt['root_receipt']['path'])
    c['initial_prices'] = root['final']['prices']
    c['initial_pensions'] = root['final']['fiscal_values']
    label, delta = CASES[case]
    prior = pipeline.HERE/'results'/label
    terminal = prior/'terminal'
    c.update(psi_change_from_initial=delta,
        terminal_checkpoint=pipeline.pin(terminal/'terminal_state.pkl.gz'),
        terminal_summary=pipeline.pin(terminal/'summary.json'),
        terminal_contract=pipeline.pin(prior/'contracts/terminal.json'),
        terminal_root_receipt=pipeline.pin(terminal/'root_receipt.json'))
    return label, c


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--case', type=int, choices=range(3))
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--contract', type=Path)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    if args.contract:
        return observed_solve(args.contract, args.output)
    sys.path[:0] = [str(pipeline.ROOT/'code/model/tools'), str(pipeline.ROOT/'code/model')]
    import run_e5f_candidate_history as driver
    label, c = prepare(args.case)
    driver.validate_contract(c)
    driver.terminal_driver.verify_sources(c)
    for name in driver.PIN_NAMES:
        driver.verify(c[name]['path'], c[name]['sha256'])
    if args.preflight:
        print(f'{label}: complete source/input preflight passed', flush=True)
        return 0
    folder = pipeline.HERE/'observed_bracket_v2'/label
    folder.mkdir(parents=True, exist_ok=False)
    pipeline.save(folder/'plan.json', dict(delta=CASES[args.case][1],
        accepted_smoke=pipeline.pin(ACCEPTED), controller=pipeline.pin(__file__),
        maximum_mappings=8 if args.case==2 else 24,
        maximum_bellman_calls=96 if args.case==2 else 640,
        stages_seconds=[1800] if args.case==2 else [1800,1800,7200],
        observed_six_date_mapping_seconds=170,
        stop='One six-date solve, at most one audited continuation; 28 dates only after convergence. Baseline case only records six-date rates; its long run already exists.',
        production_eligible=False, shock_estimated=False,
        unresolved='Female exposure and maternal-age mapping; household-rate sensitivity is diagnostic.'))
    try:
        out, summary = run_stage(folder, 'history_6', c)
        if args.case != 2 and not summary.get('finite_horizon_market_fiscal_converged'):
            c = resume_history.continuation_contract(out)
            driver.validate_contract(c)
            out, summary = run_stage(folder, 'history_6_resumed', c)
        if not summary.get('finite_horizon_market_fiscal_converged'):
            raise RuntimeError('Bounded short root did not converge; no long stage')
        if args.case != 2:
            root = pipeline.read(out/'root_receipt.json')
            terminal = pipeline.read(c['terminal_root_receipt']['path'])
            # The common helper expects continuation metadata only to discard it.
            temp = copy.deepcopy(c); temp.setdefault('continuation', {})
            lc = resume_history.long_contract(temp, root, terminal)
            driver.validate_contract(lc)
            out, summary = run_stage(folder, 'history_28', lc)
        pipeline.save(folder/'completion.json', dict(summary=summary,
            production_eligible=False, shock_estimated=False, horizon_verified=False))
    except Exception as exc:
        pipeline.save(folder/'failure.json', dict(type=type(exc).__name__, message=str(exc)))
        raise
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
