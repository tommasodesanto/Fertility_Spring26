#!/usr/bin/env python3
"""Pinned block0506 replay and no-shock operators, exclusively on Torch.

This overlay changes no scientific source, preference, measurement or shock.
It never calls the calibration normalization routine. A passing suite is not
a shocked transition, credit experiment, terminal equilibrium or horizon test.
"""
from __future__ import annotations
import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import time
import traceback

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
HERE = Path(__file__).resolve().parent
REFERENCE = ROOT / 'output/model/fertility_identification_20260928'
MANIFEST = REFERENCE / 'fixed_reference_manifest.json'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
STAGES = (('replay', 360), ('no_shock_2', 180), ('no_shock_6', 360), ('fixed_stock_no_shock_2', 180))


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def table(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def load_setup(out):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    require(sha(MANIFEST) == MANIFEST_SHA, 'Frozen manifest identity')
    m = read(MANIFEST)
    c = read(m['contract']['path'])
    for pin in (m['contract'], m['objective'], m['source_manifest']):
        require(sha(pin['path']) == pin['sha256'], 'Pin differs: ' + pin['path'])
    export = Path(m['local_export'])
    require(sha(export / 'initial_state.pkl.gz') == m['checkpoint']['sha256'], 'block0506 checkpoint identity')
    require(len(table(export / 'target_fit.csv')) == 14 and len(table(export / 'parameters.csv')) == 31, 'Full reference tables')
    for name, digest in m['artifact_hashes'].items():
        require(sha(export / name) == digest, 'Reference artifact changed: ' + name)
    sys.path[:0] = [str(ROOT / 'code/model/tools'), str(ROOT / 'tmp/e5f_overnight_local_20260927/portable/tools_v4')]
    import e5f_evening_calibration_runtime as runtime
    objective = read(m['objective']['path'])
    evaluator = runtime.setup(dict(c, objective=m['objective']), objective, out / 'runtime')
    with gzip.open(export / 'initial_state.pkl.gz', 'rb') as stream:
        packet = pickle.load(stream)
    spec = importlib.util.spec_from_file_location('fixed_auth_helpers', REFERENCE / 'fixed_reference_authenticate.py')
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    require(helper.json_value(vars(packet['parameters'])) == m['actual_serialized_parameters'], 'All serialized parameters must match')
    require(float(packet['parameters'].psi_child) == m['economic_contract']['frozen_psi_child'], 'No fertility renormalization')
    write(out / 'identity.json', dict(status='PASS', label=m['label'], manifest_sha256=MANIFEST_SHA,
          checkpoint_sha256=m['checkpoint']['sha256'], parameter_fields=len(vars(packet['parameters'])),
          source_manifest_sha256=m['source_manifest']['sha256'], objective_sha256=m['objective']['sha256'],
          overlay_sha256=sha(__file__), model_source='unchanged frozen Torch snapshot'))
    return m, packet, evaluator, runtime


def replay(out, m, packet, evaluator, runtime, deadline):
    import numpy as np
    import e5f_current_transition_runtime as native
    rt = evaluator.rt
    old = read(Path(m['local_export']) / 'receipt.json')
    evaluator.selected = packet
    evaluator.prepared['selected'] = packet
    evaluator.prepared['reference_receipt'] = old
    evaluator.fixed_reference = True
    original = copy.deepcopy(packet['parameters'])
    point = {r['parameter']: float(r['estimate']) for r in table(Path(m['local_export']) / 'parameters.csv') if r['parameter'] in runtime.FREE}
    calls = []

    def fixed_solve(point_unused, output_unused, deadline_unused):
        require(not calls and time.time() < deadline, 'Exactly one fixed-price stationary solve allowed')
        P = copy.deepcopy(original)
        start = time.monotonic()
        sol = rt['model'].solve_markov_income_at_prices(packet['solution'].p_eq.copy(), P, packet['b_grid'], verbose=False, fast_stats=False)
        seconds = time.monotonic() - start
        calls.append(seconds)
        require(time.time() < deadline, 'Replay time budget')
        completed = float(rt['chain'].extract_moments(sol, P)['tfr'])
        norm = dict(status='fixed_saved_psi_no_renormalization', psi_child=float(P.psi_child),
                    target=2.1, completed_fertility=completed, absolute_gap=abs(completed - 2.1),
                    stationary_solves=1, stationary_solve_seconds=seconds)
        renewal = dict(status='reference_replay_only', relative_replacement_gap=float(sol.adult_entry_stationary_relative_gap))
        return (sol, P, packet['solution'].p_eq.copy(), seconds, norm), old['fiscal_rule'], renewal

    evaluator.normalize = fixed_solve
    result = evaluator.evaluate_point(tax=evaluator.tax, objective=evaluator.obj, selected=packet,
        runtime=rt, point=point, output=out / 'case', deadline_epoch=deadline, graphs=True)
    comparison = native.compare_arrays(packet, evaluator.last_packet)
    write(out / 'array_comparison.json', comparison)
    require(all(row.get('exact', False) for row in comparison['arrays'].values()), 'Replay numeric arrays must be exact')
    expected_fit = table(Path(m['local_export']) / 'target_fit.csv')
    expected_parameters = table(Path(m['local_export']) / 'parameters.csv')
    require(table(out / 'case/target_fit.csv') == expected_fit, 'All 14 target rows must be exact')
    require(table(out / 'case/parameters.csv') == expected_parameters, 'All 31 parameter rows must be exact')
    plot_hashes = {}
    for name in m['standard_diagnostic_names']:
        plot_hashes[name] = sha(out / 'case/standard_diagnostics' / name)
        require(plot_hashes[name] == m['artifact_hashes']['standard_diagnostics/' + name], 'Standard plot identity: ' + name)
    runtime.verify_sources(evaluator.c)
    write(out / 'result.json', dict(status='PASS', scope='exact fixed-price stationary replay',
        actual_stationary_solves=len(calls), seconds=calls[0], arrays_compared=comparison['array_count'],
        target_rows=14, parameter_rows=31, standard_plots=17, plot_hashes=plot_hashes,
        reference_label=m['label'], original_psi=float(original.psi_child), credit_changed=False,
        market_residual=result['market_residual'], economic_changes=[], production_transition_approved=False))


def mapping(stage, out, m, packet, evaluator, runtime, deadline):
    import numpy as np
    import e5f_overnight_estate_audit as estate
    rt = evaluator.rt
    pf = rt['primitive'].pf
    P = copy.deepcopy(packet['parameters'])
    grid = packet['b_grid']
    g0 = packet['stationary_g_pre']
    initial = pf.stationary_initial_state(g0, float(g0[:, :, :, 0].sum()), float(packet['evaluation'].births), P, 1/2.1)
    require(np.array_equal(initial.g_pre, g0), 'No initial population rescaling')
    horizon = 6 if stage == 'no_shock_6' else 2
    p0 = float(packet['solution'].p_eq[0])
    supply = packet['supply_rule']
    if stage == 'fixed_stock_no_shock_2':
        supply = pf.calendar.HousingSupplyRule('fixed-stock', p0, float(supply.quantity(np.array([p0]))[0]), 0.)
    audits = []

    def observe(period, e, parameters, bg, shared, next_entrant_cohort):
        require(time.time() < deadline, 'Dated mapping deadline')
        date = out / f'date_{period:02d}'
        date.mkdir()
        budget = rt['primitive'].dated_budget(e, parameters, shared, bg, float(P.user_cost_rate * p0))
        purchase = rt['accounting'].audit_purchase_accounting(e, parameters, shared, bg, rt['model'])
        funding = estate.audit(e, parameters, bg, next_entrant_cohort=next_entrant_cohort)
        arrays = rt['audit'].policy_array_audit(dict(parameters=parameters, b_grid=bg, evaluation=e), date)
        gates = dict(budget=float(budget['budget_excess_mass']) <= 2e-10,
            transaction_wealth=float(purchase['maximum_occupied_transaction_wealth_error']) <= 1e-9,
            funded=funding['status'] == 'funded', negative_estates=float(funding['estate']['totals']['net_negative']) <= 1e-10,
            next_entrant_cohort=funding['audit_id'] == 'estate_funded_dated_entry_provisional_net_v1',
            occupied_value=arrays['occupied_negative_steps'] == 0,
            probabilities=all(not a['nonfinite'] and 0 <= a['minimum'] <= a['maximum'] <= 1 for a in arrays['probabilities'].values()))
        for key, value in purchase.items():
            if key.endswith('violation_mass') or key in ('transaction_outside_grid_mass', 'negative_estate_exposure_mass', 'saving_outside_grid_mass'):
                gates[key] = abs(float(value)) <= 2e-10
        audits.append(dict(period=period, gates=gates, budget=budget, purchase=purchase, estate=funding, arrays=arrays))
        write(out / 'dated_audits.json', audits)
        write(out.parent / 'progress.json', dict(stage=stage, period=period, epoch=time.time()))
        require(all(gates.values()), 'Dated household/accounting gate')

    result = pf.evaluate_path_at_prices(prices=np.full(horizon, p0), psi_path=np.full(horizon, P.psi_child),
        terminal_price=p0, terminal_V=packet['evaluation'].policy.V, base_parameters=P, b_grid=grid,
        initial_state=initial, supply_rule=supply, birth_to_entry_conversion=1/2.1,
        transfer_path=np.full(horizon, P.property_tax_lump_sum_transfer), pension_path=np.full(horizon, P.pension),
        payroll_tax_path=np.full(horizon, P.tau_pay), dated_observer=observe)
    write(out / 'rows.json', result.rows)
    fiscal = max(abs(float(r['scaled_pension_budget_residual'])) for r in result.rows)
    gates = dict(market=result.maximum_market_residual <= 2e-4, fiscal=fiscal <= 1e-6,
        mass=result.maximum_mass_accounting_error <= 2e-8, projection=result.maximum_feasibility_projection_mass == 0,
        backward_forward=result.maximum_policy_reproduction_error <= 1e-10,
        calls=result.bellman_solves == 2*horizon)
    result_row = dict(status='PASS' if all(gates.values()) else 'FAIL', scope='prescribed constant-path operator; no shocked equilibrium',
        gates=gates, horizon=horizon, bellman_solves=result.bellman_solves, seconds=result.elapsed_seconds,
        housing_supply_mode=supply.mode, physical_stock=float(supply.initial_stock),
        market_residual=result.maximum_market_residual, fiscal_residual=fiscal,
        mass_accounting_error=result.maximum_mass_accounting_error,
        initial_mass=float(g0.sum()), final_mass=float(result.terminal_state.g_pre.sum()),
        population_l1=float(np.abs(result.terminal_state.g_pre-g0).sum()),
        adjusted_queue_max_change=float(np.max(np.abs(pf.birth_queue_values(result.terminal_state.scheduled_entries)-pf.birth_queue_values(initial.scheduled_entries)))),
        fixed_preferences=True, fertility_renormalized=False, production_transition_approved=False)
    write(out / 'result.json', result_row)
    runtime.verify_sources(evaluator.c)
    require(all(gates.values()), 'No-shock aggregate gate; diagnose without changing tolerance')


def suite():
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    out = HERE / 'runs' / os.environ['SLURM_JOB_ID']
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    write(out / 'plan.json', dict(stages=STAGES, workers=1, total_seconds=1200, stationary_solve_cap=1,
        dated_bellman_solve_cap=20, source_overlay_sha256=sha(__file__), manifest_sha256=MANIFEST_SHA,
        economic_changes='none in replay/no-shock; fixed-stock no-shock changes supply elasticity only at its exact reference quantity',
        reference_untouched=True, no_credit_experiment=True, no_historical_shock=True, no_retries=True))
    completed = []
    for stage, seconds in STAGES:
        remaining = min(seconds, 1200 - (time.time()-started))
        require(remaining > 0, 'Suite budget exhausted')
        case = out / stage
        write(out / 'progress.json', dict(stage=stage, epoch=time.time(), state='starting'))
        with (out / f'{stage}.log').open('w') as log:
            try:
                run = subprocess.run([sys.executable, __file__, '--stage', stage, '--out', str(case),
                    '--deadline', str(time.time()+remaining)], stdout=log, stderr=subprocess.STDOUT, timeout=remaining)
                require(run.returncode == 0, stage + ' failed; inspect preserved log')
            except BaseException as exc:
                write(out / 'suite_result.json', dict(status='FAIL', stage=stage, completed=completed, error=str(exc), elapsed=time.time()-started))
                raise
        completed.append(read(case / 'result.json'))
        write(out / 'latest_completed.json', dict(stage=stage, result=completed[-1]))
        write(out / 'best_so_far.json', dict(meaning='furthest completed gate, not an optimized candidate', completed_stages=[s for s,_ in STAGES[:len(completed)]]))
    write(out / 'suite_result.json', dict(status='PASS', completed=completed, elapsed=time.time()-started,
        production_transition_approved=False, no_shocked_cases_run=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', default='suite', choices=['suite']+[x[0] for x in STAGES])
    parser.add_argument('--out', type=Path)
    parser.add_argument('--deadline', type=float)
    args = parser.parse_args()
    if args.stage == 'suite':
        suite()
    else:
        args.out.mkdir(parents=True, exist_ok=False)
        try:
            m, packet, evaluator, runtime = load_setup(args.out)
            if args.stage == 'replay':
                replay(args.out, m, packet, evaluator, runtime, args.deadline)
            else:
                mapping(args.stage, args.out, m, packet, evaluator, runtime, args.deadline)
        except BaseException as exc:
            write(args.out / 'failure.json', dict(error_type=type(exc).__name__, error=str(exc), traceback=traceback.format_exc()))
            raise
