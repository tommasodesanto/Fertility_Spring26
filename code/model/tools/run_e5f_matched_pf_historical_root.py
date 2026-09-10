"""Bounded matched historical PF market root from a pinned price Jacobian.

Every trial uses the complete unchanged twelve-moment observer and a solved
stationary endpoint. Finite-horizon market convergence is not tail convergence.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import threading
import time
from types import SimpleNamespace

import run_e5f_matched_pf_baseline as baseline
import numpy as np
import e5f_matched_pf_path_root as root
import run_e5f_perfect_foresight_rebated_property_tax as rent_domain

primitive, joined, pf = baseline.primitive, baseline.joined, baseline.pf
MAXIMUM_ROOT_SECONDS = 21600


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', required=True, type=Path)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--arm', required=True, choices=('sequential', 'nested'))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    c, originating = joined.load_smoke_contract(args.contract, args.contract_sha256,
                                               args.arm, maximum_seconds=MAXIMUM_ROOT_SECONDS)
    if (c['experiment'] != 'normalized_historical_path_root' or c['probe_coordinate'] != -1
            or type(c['maximum_path_evaluations']) is not int
            or not 2 <= c['maximum_path_evaluations'] <= 6):
        raise ValueError('Explicit bounded six-evaluation historical root contract required')
    primitive.verify(c['jacobian_packet'], c['jacobian_packet_sha256'])
    packet = json.loads(Path(c['jacobian_packet']).read_text())
    # The remaining packet/scientific-source checks live below, before any solve.
    validate_jacobian_packet(packet, c, args.arm)
    initial_prices, initial_jacobian = load_restart(c, packet, args.arm)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    progress = dict(phase='load', arm=args.arm, evaluation=0)
    stop = threading.Event()
    def save(name, value):
        pf.write_json(out/name, value)
    def heartbeat():
        while not stop.wait(15):
            elapsed = time.monotonic() - started
            save('heartbeat.json', dict(progress, elapsed_seconds=elapsed))
            if elapsed > c['seconds']:
                save('failure.json', dict(progress, error='wall-time budget exhausted'))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        seed = baseline.load_seed(c, originating, args.arm)
        old, _ = baseline.load_normalized(c, args.arm)
        terminal_packet = baseline.load_terminal(c, args.arm, old)
        terminal = SimpleNamespace(parameters=terminal_packet['parameters'],
            asset_price=float(terminal_packet['policy'].price[0]))
        del old, terminal_packet
        def project(prices):
            return rent_domain.project_price_path_to_positive_rents(prices,
                terminal=terminal, minimum_rent_share=1e-6)[0]
        save('contract.json', dict(c, contract_sha256=args.contract_sha256, arm=args.arm,
            numerical_controls=dict(slope=1.63, market_tolerance=2e-4,
                max_log_step=.10, damping=1., max_condition_number=1e10,
                worsening_factor=1.5, final_reproduction_tolerance=2e-10),
            scope='finite-horizon dated housing-price root, not a re-estimation or production promotion'))
        def evaluate(prices):
            progress.update(evaluation=progress['evaluation']+1, phase='path')
            trial = out/f"evaluation_{progress['evaluation']:03d}"
            trial.mkdir()
            def trial_save(name, value):
                pf.write_json(trial/name, value)
            result = baseline.run_history_probe(seed, c, args, trial, progress, trial_save,
                time.monotonic(), prices_override=prices)
            return dict(residual=result['residual'], mapping_valid=result['mapping_valid'],
                payload=dict(directory=str(trial), loss=result['loss'],
                    terminal_distance=result['terminal_distance'],
                    summary_sha256=primitive.digest(trial/'summary.json')))
        def record(value):
            save('latest_completed.json', value)
            if value.get('new_best'):
                save('best_so_far.json', value)
        result = root.solve_price_path(initial_prices=initial_prices, evaluate=evaluate,
            project=project, slope=1.63, market_tolerance=2e-4, max_log_step=.10,
            damping=1., max_evaluations=c['maximum_path_evaluations'],
            deadline_monotonic=started+c['seconds'], max_condition_number=1e10,
            worsening_factor=1.5, final_reproduction_tolerance=2e-10,
            initial_jacobian=initial_jacobian, callback=record)
        best = result['final'] if result['converged'] else result['best']
        tail_pass = bool(best is not None and best['payload']['terminal_distance']['all_checks_pass'])
        save('root_history.json', result)
        save('summary.json', dict(status=result['status'], arm=args.arm,
            finite_horizon_market_converged=result['converged'],
            terminal_distance_passed=tail_pass, horizon_extension_verified=False,
            historical_equilibrium_certified=False, calibrated_history=False,
            production_promoted=False, evaluations=result['evaluations'],
            elapsed_seconds=time.monotonic()-started,
            final_reproduction_max_abs=result['final_reproduction_max_abs'],
            best=best, outstanding=['horizon extension and terminal convergence',
                'production empirical group/date alignment', 'matched recalibration']))
    except Exception as error:
        save('failure.json', dict(progress, error=repr(error), elapsed_seconds=time.monotonic()-started))
        raise
    finally:
        stop.set()


def load_restart(c, packet, arm):
    """Resume only a completed, reproduced root under the same economic inputs."""
    names = ('restart_contract', 'restart_history', 'restart_summary')
    if not any(name in c or name+'_sha256' in c for name in names):
        return packet['prices'], packet['jacobian']
    if not all(name in c and name+'_sha256' in c for name in names):
        raise ValueError('A restart requires all three pinned root receipts')
    parent = Path(c['restart_contract']).parent
    for name, filename in zip(names, ('contract.json', 'root_history.json', 'summary.json')):
        if Path(c[name]) != parent/filename:
            raise ValueError('Restart receipts must belong to one completed root directory')
    for name in names:
        if not Path(c[name]).is_absolute():
            raise ValueError('Restart receipt paths must be absolute')
        primitive.verify(c[name], c[name+'_sha256'])
    previous = json.loads(Path(c['restart_contract']).read_text())
    history = json.loads(Path(c['restart_history']).read_text())
    summary = json.loads(Path(c['restart_summary']).read_text())
    keys = ('arm', 'checkpoint_sha256', 'selected_summary_sha256',
        'normalized_checkpoint_sha256', 'normalized_summary_sha256', 'normalized_contract_sha256',
        'terminal_checkpoint_sha256', 'terminal_summary_sha256', 'terminal_contract_sha256',
        'demographic_sources', 'target_fingerprint', 'path_date_count',
        'initial_price_rule', 'terminal_preference_rule', 'probe_log_step', 'jacobian_packet_sha256')
    if any(previous[name] != c[name] for name in keys) or summary['arm'] != arm:
        raise ValueError('Restart scientific inputs differ')
    # These reviewed driver/test changes only add state saving, longer budgets
    # and receipt-based warm starts. All economic and numerical kernels agree.
    wrappers = {'code/model/tools/run_e5f_matched_pf_baseline.py',
        'code/model/tools/run_e5f_matched_pf_historical_root.py',
        'code/model/tools/test_run_e5f_matched_pf_historical_root.py'}
    for name, expected in previous['source_sha256'].items():
        if name not in wrappers and c['source_sha256'].get(name) != expected:
            raise ValueError(f'Restart economic source changed: {name}')
    best, final = history.get('best'), history.get('final')
    if (best is None or final is None or not best['mapping_valid'] or not final['mapping_valid']
            or summary.get('finite_horizon_market_converged') is not False
            or summary.get('final_reproduction_max_abs') is None
            or not np.isfinite(summary['final_reproduction_max_abs'])
            or summary['final_reproduction_max_abs'] > 2e-10
            or summary['best'] != best or history['evaluations'] != summary['evaluations']):
        raise ValueError('Restart must be a complete reproduced unfinished market root')
    n = c['path_date_count']
    prices = np.asarray(best['prices'], dtype=float)
    residual = np.asarray(best['residual'], dtype=float)
    matrix = np.asarray(history['final_jacobian'], dtype=float)
    final_residual = np.asarray(final['residual'], dtype=float)
    if (final_residual.shape != (n,) or not np.isfinite(final_residual).all()
            or prices.shape != (n,) or residual.shape != (n,) or matrix.shape != (n,n)
            or not np.isfinite(prices).all() or np.any(prices <= 0)
            or not np.isfinite(residual).all() or not np.isfinite(matrix).all()
            or not np.array_equal(prices, np.asarray(final['prices']))
            or not np.allclose(residual, final_residual, rtol=0, atol=2e-10)
            or not np.isclose(best['score'], np.max(np.abs(residual)), rtol=0, atol=1e-12)):
        raise ValueError('Restart price/residual/Jacobian or replay mismatch')
    # The prior approximate Jacobian is a preconditioner. The new root still
    # evaluates these prices afresh and reserves another final reproduction.
    return prices, matrix


def validate_jacobian_packet(packet, c, arm):
    """Fail closed until the collector's complete provenance is reconciled."""
    import collect_e5f_matched_pf_price_jacobian as collector
    if (packet.get('schema') != collector.SCHEMA
            or packet.get('status') != 'complete_validated_finite_difference_jacobian'
            or packet['arm'] != arm or c.get('arm') != arm
            or packet['target_fingerprint'] != c['target_fingerprint']):
        raise ValueError('Unexpected Jacobian receipt, arm or target contract')
    provenance = packet['provenance']
    repeated = collector.collect(provenance['anchor']['directory'],
        [row['directory'] for row in provenance['columns']])
    if repeated != packet:
        raise ValueError('Jacobian packet does not reproduce from its pinned complete probe panel')
    shared = packet['shared_contract']
    required = ('checkpoint_sha256', 'selected_summary_sha256', 'normalized_checkpoint_sha256',
        'normalized_summary_sha256', 'normalized_contract_sha256', 'terminal_checkpoint_sha256',
        'terminal_summary_sha256', 'terminal_contract_sha256', 'demographic_sources',
        'target_fingerprint', 'path_date_count', 'initial_price_rule',
        'terminal_preference_rule', 'probe_log_step')
    if any(shared[name] != c[name] for name in required):
        raise ValueError('Jacobian and root scientific input contracts differ')
    driver = 'code/model/tools/run_e5f_matched_pf_baseline.py'
    reviewed = c['reviewed_evaluator_driver_change']
    if reviewed != dict(path=driver, from_sha256=shared['source_sha256'][driver],
                        to_sha256=c['source_sha256'][driver],
                        scope='explicit supplied-price hook and optional dated-state checkpoint; unchanged economic evaluator'):
        raise ValueError('Explicit reviewed probe-to-root driver change required')
    runtime_paths = {
        'code/model/tools/run_e5f_matched_pf_historical_root.py',
        'code/model/tools/test_run_e5f_matched_pf_historical_root.py',
        'code/model/tools/run_e5f_matched_pf_history.py',
    }
    changed_runtime = {name for name in runtime_paths
        if name in shared['source_sha256']
        and c['source_sha256'].get(name) != shared['source_sha256'][name]}
    if changed_runtime or 'reviewed_root_runtime_changes' in c:
        expected_review = dict(
            scope='root runtime ceiling and explicit provenance checks only; unchanged economic evaluator and numerical root',
            files={name: dict(from_sha256=shared['source_sha256'][name],
                             to_sha256=c['source_sha256'].get(name))
                   for name in sorted(changed_runtime)})
        if (not changed_runtime or any(c['source_sha256'].get(name) is None for name in changed_runtime)
                or c.get('reviewed_root_runtime_changes') != expected_review):
            raise ValueError('Jacobian economic source changed: explicit reviewed root runtime source changes required')
        for name in changed_runtime:
            collector.require_hash(c['source_sha256'][name], name)
    for name, expected in shared['source_sha256'].items():
        if name != driver and name not in changed_runtime and c['source_sha256'].get(name) != expected:
            raise ValueError(f'Jacobian economic source changed before root: {name}')
    n = c['path_date_count']
    if (packet['years'] != list(2007 + 4*np.arange(n))
            or np.asarray(packet['jacobian']).shape != (n,n)
            or not np.isfinite(packet['jacobian']).all()):
        raise ValueError('Jacobian horizon or finite matrix invalid')


if __name__ == '__main__':
    main()
