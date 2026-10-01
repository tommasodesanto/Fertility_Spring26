#!/usr/bin/env python3
"""Torch-only pension-accounting diagnostic at fixed house prices.

This is a bounded numerical check, not a shock estimator or a production
transition launcher.  It changes only the dated pension path used in a fresh
native mapping and stops at the first uncertified trial.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import os
import pickle
from pathlib import Path
import sys
import threading
import time

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
TOOLS = Path(__file__).resolve().parent.parent / 'model' / 'tools'
FROZEN_TOOLS = ROOT / 'tmp/e5f_overnight_local_20260927/portable/tools_v4'
ESTIMATOR_SOURCES = (
    'run_e5f_preference_transition.py', 'e5f_four_shock_acceleration.py',
    'e5f_exact_policy_cache.py', 'e5f_social_security_root.py',
    'e5f_matched_pf_path_root.py', 'e5f_ssj_scaled_step_root.py',
    'e5f_ssj_toeplitz_jacobian.py', 'e5f_preference_shock_fit.py',
    'run_e5f_preference_estimation.py')
SOURCE_NAMES = ESTIMATOR_SOURCES + (Path(__file__).name,)
GATES = dict(market=2e-4, fiscal=1e-6, replay=1e-10)
TERMINAL_TOLERANCES = ('population_relative_gap', 'normalized_distribution_l1',
                       'birth_queue_maximum_relative_gap', 'asset_price_relative_gap',
                       'renter_price_relative_gap', 'psi_absolute_gap')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def _json_default(value):
    if hasattr(value, 'tolist'):
        return value.tolist()
    raise TypeError(f'Object of type {type(value).__name__} is not JSON serializable')


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=_json_default) + '\n')
    temporary.replace(path)


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Pinned path/SHA-256 required')
    path = Path(item['path'])
    require(path.is_file() and len(item['sha256']) == 64 and sha(path) == item['sha256'],
            'Pinned artifact changed: ' + str(path))
    return path


def budget_update(rows, damping=1.0):
    """Balance each recorded old-distribution account without mutating rows."""
    require(isinstance(rows, list) and rows, 'At least one mapping row required')
    require(type(damping) in (int, float) and math.isfinite(damping) and 0 < damping <= 1,
            'Damping must be finite in (0,1]')
    updated = []
    for row in rows:
        try:
            revenue = float(row['payroll_tax_revenue']); outlays = float(row['pension_outlays'])
            pension = float(row['pension_period_units']); price = float(row['asset_price'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('Mapping row lacks pension-accounting fields') from exc
        require(all(math.isfinite(v) and v > 0 for v in (revenue, outlays, pension, price)),
                'Pension-accounting units must be finite and positive')
        updated.append(pension * math.exp(float(damping) * math.log(revenue / outlays)))
    return updated


def _load_inner():
    """Pin all code before importing a native dependency."""
    spec = importlib.util.spec_from_file_location('e5f_budget_inner', TOOLS / 'run_e5f_preference_transition.py')
    require(spec is not None and spec.loader is not None, 'Cannot load pinned transition helper')
    module = importlib.util.module_from_spec(spec)
    sys.path[:0] = [str(TOOLS), str(ROOT / 'code' / 'model' / 'tools'), str(FROZEN_TOOLS)]
    spec.loader.exec_module(module)
    return module


def validate_config(config):
    required = {'mode', 'reference_manifest_sha256', 'source_pins', 'endpoint_receipt', 'horizon', 'cache_max_bytes', 'mapping_seconds',
                'total_seconds', 'numerical_threads', 'market_tolerance', 'fiscal_tolerance',
                'replay_tolerance', 'smoke_pension_factor'}
    allowed = required | {'baseline_mapping'}
    require(set(config).issubset(allowed) and required.issubset(config), 'Diagnostic config schema differs')
    require(config['mode'] in ('smoke', 'full'), 'Mode must be smoke or full')
    horizon = config['horizon']; require(type(horizon) is int and horizon == (6 if config['mode'] == 'smoke' else 104), 'Mode horizon differs')
    require(config['reference_manifest_sha256'] == '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4', 'Wrong fixed reference')
    require(set(config['source_pins']) == set(SOURCE_NAMES), 'Exact nine estimator sources plus harness required')
    require(all(isinstance(v, str) and len(v) == 64 for v in config['source_pins'].values()), 'Complete source SHA-256 pins required')
    require(('baseline_mapping' not in config) if config['mode'] == 'smoke' else isinstance(config.get('baseline_mapping'), dict), 'Full mode requires baseline mapping only')
    if 'baseline_mapping' in config: pinned(config['baseline_mapping'])
    pinned(config['endpoint_receipt'])
    require(type(config['cache_max_bytes']) is int and 0 <= config['cache_max_bytes'] <= 64 * 1024**3, 'Invalid cache budget')
    require(all(type(config[k]) in (int, float) and math.isfinite(config[k]) and config[k] > 0 for k in ('mapping_seconds', 'total_seconds')), 'Finite time budgets required')
    require(config['mapping_seconds'] <= config['total_seconds'] and config['total_seconds'] <= 14400 and config['mapping_seconds'] <= 10800, 'Budget exceeds diagnostic cap')
    require(isinstance(config['numerical_threads'], str) and config['numerical_threads'].isdigit() and
            int(config['numerical_threads']) > 0 and config['numerical_threads'] == os.environ.get('NUMBA_NUM_THREADS'),
            'NUMBA thread count must be recorded and verified')
    require(all(config[k] == value for k, value in (('market_tolerance', GATES['market']), ('fiscal_tolerance', GATES['fiscal']), ('replay_tolerance', GATES['replay']))), 'Fixed acceptance gates changed')
    require(config['smoke_pension_factor'] == 1.00001, 'Smoke pension perturbation changed')


def check_sources(config):
    locations = {name: (Path(__file__) if name == Path(__file__).name else TOOLS / name) for name in SOURCE_NAMES}
    for name, path in locations.items():
        require(path.is_file() and sha(path) == config['source_pins'][name], 'Numerical source changed: ' + name)


def _all_true(value):
    if isinstance(value, dict): return bool(value) and all(_all_true(v) for v in value.values())
    return value is True


def _residuals(rows):
    market = []
    fiscal = []
    for row in rows:
        try:
            demand = float(row['housing_demand']); supply = float(row['housing_supply'])
            revenue = float(row['payroll_tax_revenue']); outlays = float(row['pension_outlays'])
            stored = float(row['scaled_pension_budget_residual'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('Mapping row lacks residual-accounting fields') from exc
        require(all(math.isfinite(x) for x in (demand, supply, revenue, outlays, stored)) and
                demand > 0 and supply > 0 and revenue > 0 and outlays > 0, 'Mapping residual accounts must be finite and positive')
        market.append((demand - supply) / supply)
        calculated = (revenue - outlays) / max(abs(revenue), abs(outlays))
        require(abs(stored - calculated) <= 1e-15, 'Recorded fiscal residual differs from accounts')
        fiscal.append(calculated)
    return market, fiscal


def _same_array(actual, expected):
    return len(actual) == len(expected) and all(abs(float(a) - float(b)) <= 1e-15 for a, b in zip(actual, expected))


def validate_baseline(record, horizon, reference):
    expected = {'backward_forward_policy_error', 'cache', 'diagnostic_packets', 'fertility', 'final_mass', 'fiscal_residual',
                'gates', 'initial_mass', 'market_residual', 'mass_error', 'population_l1', 'projection_mass', 'rows', 'seconds'}
    require(set(record) == expected, 'Baseline receipt schema differs')
    require(len(record.get('rows', [])) == horizon and len(record.get('fertility', [])) == horizon,
            'Baseline row/fertility horizon differs')
    require(isinstance(record['gates'], dict) and len(record['gates']) == 4 and _all_true(record['gates']),
            'Baseline lacks four passing recorded mapping gates')
    require(all(math.isfinite(float(record[key])) and float(record[key]) >= 0 for key in
                ('backward_forward_policy_error', 'final_mass', 'initial_mass', 'mass_error', 'population_l1', 'projection_mass', 'seconds')) and
            record['final_mass'] > 0 and record['initial_mass'] > 0 and record['seconds'] > 0,
            'Baseline receipt has nonfinite or invalid accounting values')
    rows = record['rows']; q = float(reference['q']); psi = float(reference['psi']); pension = float(reference['pension'])
    require(all(r.get('calendar_year') == 2007 + 4 * t and float(r['asset_price']) == q and
                float(r['psi_child']) == psi and float(r['pension_period_units']) == pension for t, r in enumerate(rows)),
            'Baseline changed fixed prices, preference, pension, or calendar')
    market, fiscal = _residuals(rows)
    require(_same_array(record['market_residual'], market) and _same_array(record['fiscal_residual'], fiscal),
            'Stored residual arrays differ from rows')
    require(max(map(abs, market)) <= GATES['market'], 'Baseline recorded housing gate fails')


def validate_record_residuals(record, horizon):
    rows = record.get('rows', [])
    require(len(rows) == horizon and len(record.get('fertility', [])) == horizon, 'Mapping row/fertility count differs')
    market, fiscal = _residuals(rows)
    require(_same_array(record.get('market_residual', []), market) and _same_array(record.get('fiscal_residual', []), fiscal),
            'Mapping residual fields/array differ')


def summarize(record, terminal):
    market, fiscal = record['market_residual'], record['fiscal_residual']
    return dict(max_market_error=max(map(abs, market)), max_fiscal_error=max(map(abs, fiscal)),
                residual_ratios=dict(market=max(map(abs, market)) / GATES['market'], fiscal=max(map(abs, fiscal)) / GATES['fiscal']),
                cache=record['cache'], gates=record['gates'], terminal=terminal)


class Diagnostic:
    def __init__(self, config, output, inner):
        self.config, self.output, self.inner = config, Path(output), inner
        self.deadline = time.monotonic() + config['total_seconds']
        self.reference = None; self.packet = self.evaluator = self.terminal = self.endpoint = None
        self.best_score = math.inf

    def guarded(self, call):
        import signal
        remaining = min(self.config['mapping_seconds'], self.deadline - time.monotonic())
        if remaining <= 0: raise TimeoutError('Diagnostic time budget exhausted')
        def timeout(*_): raise TimeoutError('Bounded mapping exhausted its time budget')
        previous = signal.signal(signal.SIGALRM, timeout); signal.setitimer(signal.ITIMER_REAL, remaining)
        try: return call()
        finally: signal.setitimer(signal.ITIMER_REAL, 0); signal.signal(signal.SIGALRM, previous)

    def setup(self):
        manifest, self.packet, self.evaluator = self.inner.load_reference(self.output / 'reference')
        sys.path[:0] = [str(TOOLS), str(ROOT / 'code' / 'model' / 'tools'), str(FROZEN_TOOLS)]
        psi = float(self.packet['parameters'].psi_child); pension = float(self.packet['parameters'].pension); q = float(self.packet['solution'].p_eq[0])
        self.reference = dict(manifest_sha256=self.config['reference_manifest_sha256'], source_manifest_sha256=manifest['source_manifest']['sha256'], psi=psi, pension=pension, q=q)
        receipt = self.inner.read(self.inner.pinned(self.config['endpoint_receipt']))
        require(set(receipt) == {'checkpoint', 'housing', 'native_one_step_verified', 'pension', 'population_scale', 'price',
                                 'psi_child', 'reference_manifest_sha256', 'repeat_verified', 'source_manifest_sha256', 'terminal'},
                'Endpoint receipt schema differs')
        checkpoint = self.inner.pinned(receipt['checkpoint'])
        with gzip.open(checkpoint, 'rb') as stream: self.terminal = pickle.load(stream)
        import numpy as np
        self.inner.check_endpoint_primitives(self.packet['parameters'], self.terminal['parameters'])
        require(np.array_equal(self.packet['b_grid'], self.terminal['b_grid']) and receipt['repeat_verified'] is True and
                np.array_equal(self.terminal['stationary_g_pre'], self.terminal['evaluation'].g_pre) and
                abs(float(np.sum(self.terminal['evaluation'].g_pre)) - 1.) <= 1e-9 and
                receipt['native_one_step_verified'] is True and receipt['terminal']['all_checks_pass'] is True and receipt['housing'] == 'fixed_stock' and
                receipt['psi_child'] == self.terminal['parameters'].psi_child == psi and receipt['price'] == self.terminal['solution'].p_eq[0] == q and
                receipt['pension'] == self.terminal['parameters'].pension == pension and receipt['reference_manifest_sha256'] == self.reference['manifest_sha256'] and
                receipt['source_manifest_sha256'] == self.reference['source_manifest_sha256'], 'Endpoint evidence failed authentication')
        scale = float(receipt['population_scale']); supply = float(self.packet['supply_rule'].quantity([q])[0]); demand = float(self.terminal['evaluation'].demand_by_loc[0])
        require(math.isfinite(scale) and scale > 0 and demand > 0 and abs(supply / demand / scale - 1) <= 1e-8,
                'Endpoint population scale inconsistent with fixed stock')
        self.endpoint = dict(receipt, checkpoint=dict(path=str(checkpoint), sha256=self.inner.sha(checkpoint)), population_scale=scale)

    def mapping(self, name, pensions):
        import numpy as np
        H = self.config['horizon']; q = self.reference['q']; psi = self.reference['psi']
        result, record = self.guarded(lambda: self.inner.mapping(self.packet, self.evaluator, self.terminal, self.endpoint,
            np.full(H, q), np.asarray(pensions), np.full(H, psi), 'fixed_stock', self.output / name,
            self.config['cache_max_bytes'], measure_fertility=True, initial_state=None, start_year=2007))
        validate_record_residuals(record, H)
        terminal = self.inner.terminal_checks(self.packet, self.evaluator, self.terminal, self.endpoint, result,
            [psi] * H, dict(terminal_tolerances={k: 1e-3 for k in TERMINAL_TOLERANCES}, raw_queue_relative_tolerance=1e-3))
        terminal = self.inner.plain(terminal)
        accepted = (all(record['gates'].values()) and terminal['all_checks_pass'] and
                    max(map(abs, record['market_residual'])) <= GATES['market'] and max(map(abs, record['fiscal_residual'])) <= GATES['fiscal'])
        receipt = dict(phase='mapping', last_completed=name, name=name, coordinates=dict(asset_price=q, psi_child=psi, pensions=list(map(float, pensions))),
                       accepted=accepted, **summarize(record, terminal))
        write(self.output / 'latest_completed.json', receipt)
        score = max(receipt['residual_ratios'].values())
        if score < self.best_score:
            self.best_score = score; write(self.output / 'best_so_far.json', dict(receipt, normalized_max_residual=score))
        return result, record, terminal, receipt

    def run(self):
        self.guarded(self.setup); H = self.config['horizon']; b0 = [self.reference['pension']] * H
        if self.config['mode'] == 'smoke':
            base_result, base, base_terminal, base_summary = self.mapping('baseline', [x * self.config['smoke_pension_factor'] for x in b0])
        else:
            base = read(pinned(self.config['baseline_mapping'])); validate_baseline(base, H, self.reference)
            base_terminal = base.get('terminal'); base_summary = summarize(base, base_terminal)
            base_result = None
            base_score = max(base_summary['residual_ratios'].values())
            self.best_score = base_score; write(self.output / 'best_so_far.json', dict(name='baseline', accepted=False, **base_summary, normalized_max_residual=base_score))
        baseline_fertility = [float(row['period_tfr_topcode_adjusted']) for row in base['fertility'][:4]]
        input_coordinates = dict(pensions=[float(row['pension_period_units']) for row in base['rows']])
        del base_result
        updated = budget_update(base['rows'])
        trial_result, trial, trial_terminal, trial_summary = self.mapping('trial', updated)
        trial_fertility = [float(row['period_tfr_topcode_adjusted']) for row in trial['fertility'][:4]]
        max_first_four_fertility_change = max(abs(a - b) for a, b in zip(baseline_fertility, trial_fertility))
        if not trial_summary['accepted']:
            return dict(numerical_certified=False, outcome='trial_not_certified', baseline=base_summary,
                        input=base_summary, input_coordinates=input_coordinates, updated=trial_summary, replay=None, baseline_fertility=baseline_fertility,
                        trial_fertility=trial_fertility, max_first_four_fertility_change=max_first_four_fertility_change,
                        comparison=('improved' if trial_summary['max_fiscal_error'] < base_summary['max_fiscal_error'] else 'worse_or_equal'))
        replay_result, replay, replay_terminal, replay_summary = self.mapping('fresh_replay', updated)
        differences = dict(market=max(abs(a-b) for a,b in zip(trial['market_residual'], replay['market_residual'])),
                           fiscal=max(abs(a-b) for a,b in zip(trial['fiscal_residual'], replay['fiscal_residual'])),
                           fertility=max(abs(float(a['period_tfr_topcode_adjusted'])-float(b['period_tfr_topcode_adjusted'])) for a,b in zip(trial['fertility'], replay['fertility'])))
        import numpy as np
        pf = self.evaluator.rt['primitive'].pf
        differences['g_pre'] = float(np.max(np.abs(trial_result.terminal_state.g_pre - replay_result.terminal_state.g_pre)))
        differences['queues'] = max(float(np.max(np.abs(pf.birth_queue_values(getattr(trial_result.terminal_state, name)) - pf.birth_queue_values(getattr(replay_result.terminal_state, name))))) for name in ('scheduled_entries', 'scheduled_raw_entries'))
        certified = (replay_summary['accepted'] and trial_summary['gates'] == replay_summary['gates'] and
                     trial_summary['terminal'] == replay_summary['terminal'] and max(differences.values()) <= GATES['replay'])
        return dict(numerical_certified=certified, outcome='certified' if certified else 'fresh_replay_not_exact',
                    baseline=base_summary, input=base_summary, input_coordinates=input_coordinates, updated=trial_summary, replay=replay_summary,
                    baseline_fertility=baseline_fertility, trial_fertility=trial_fertility,
                    max_first_four_fertility_change=max_first_four_fertility_change, replay_max_difference=differences)


def run(config, output):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    validate_config(config); check_sources(config)
    output = Path(output); require(not output.exists(), 'Output directory already exists'); output.mkdir(parents=True)
    write(output / 'config.json', config)
    stop = threading.Event()
    write(output / 'process_heartbeat.json', dict(epoch=time.time(), phase='created', last_completed='none'))
    def heartbeat():
        while not stop.wait(60): write(output / 'process_heartbeat.json', dict(epoch=time.time(), phase='budget_diagnostic_active', last_completed='see_latest_completed'))
    worker = threading.Thread(target=heartbeat, daemon=True); worker.start(); started = time.monotonic()
    try:
        inner = _load_inner(); result = Diagnostic(config, output, inner).run()
        result.update(runtime_seconds=time.monotonic()-started, unchanged_primitives='reference q, psi, grid, terminal, fixed-stock closure; pension path only',
                      production_ready=False, fitted_shocks=False)
        write(output / 'complete.json', result)
        require(result['numerical_certified'], 'Diagnostic completed but is not numerically certified')
    except BaseException as exc:
        write(output / 'failure.json', dict(error_type=type(exc).__name__, error=str(exc), numerical_certified=False))
        raise
    finally: stop.set(); worker.join(timeout=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True); parser.add_argument('--config-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True); args = parser.parse_args()
    require(sha(args.config) == args.config_sha256, 'Config SHA-256 differs')
    run(read(args.config), args.output)


if __name__ == '__main__': main()
