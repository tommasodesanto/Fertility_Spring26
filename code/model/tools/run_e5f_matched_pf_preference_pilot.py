"""Bounded fixed-endpoint preference-shape diagnostics at inherited parameters.

The converged parent supplies a numerical price seed, never candidate convergence.
Every candidate is evaluated afresh; optional roots retain all existing gates.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import threading
import time
from types import SimpleNamespace

import numpy as np
import run_e5f_matched_pf_baseline as baseline
import e5f_matched_pf_path_root as price_root
import run_e5f_perfect_foresight_rebated_property_tax as rent_domain

primitive, joined, pf = baseline.primitive, baseline.joined, baseline.pf
DRIVER = 'code/model/tools/run_e5f_matched_pf_baseline.py'
CHANGE_SCOPE = 'strict optional fixed-endpoint preference-path hook; unchanged default household and population evaluator'
SHAPES = (-.5, 0., .5)


def preference_shape(old_path, count, coefficient):
    """psi0 + Delta * [x + a*x*(1-x)], with exact inherited endpoints."""
    if coefficient not in SHAPES or not np.isfinite(coefficient):
        raise ValueError('Only the three contracted diagnostic shape coefficients are allowed')
    reference = baseline.validated_preference_path(old_path, count)
    if coefficient == 0.:
        return reference  # Exact numerical reproduction of the inherited formula.
    x = np.minimum(np.arange(count, dtype=float) / 4., 1.)
    result = reference[0] + (reference[4] - reference[0]) * (x + coefficient*x*(1.-x))
    result[0], result[4:] = reference[0], reference[4]
    return baseline.validated_preference_path(old_path, count, result)


def validate_parent(parent, history, summary, c, arm):
    """Validate a pinned converged parent, with no restart/convergence waiver."""
    keys = ('checkpoint_sha256', 'selected_summary_sha256', 'normalized_checkpoint_sha256',
        'normalized_summary_sha256', 'normalized_contract_sha256', 'terminal_checkpoint_sha256',
        'terminal_summary_sha256', 'terminal_contract_sha256', 'demographic_sources',
        'target_fingerprint', 'initial_price_rule', 'terminal_preference_rule', 'probe_log_step')
    if (parent['arm'] != arm or summary['arm'] != arm or c['arm'] != arm
            or any(parent[k] != c[k] for k in keys)):
        raise ValueError('Pilot parent economic input/arm contract mismatch')
    if (summary.get('finite_horizon_market_converged') is not True
            or summary.get('status') != 'converged' or not history.get('converged')
            or summary.get('final_reproduction_max_abs') != 0.
            or history.get('final_reproduction_max_abs') != 0.):
        raise ValueError('Pilot seed must be a completed exactly reproduced finite-horizon root')
    final = history.get('final')
    if final is None or not final.get('mapping_valid') or summary['best'] != final:
        raise ValueError('Pilot parent final mapping receipt mismatch')
    n = parent['path_date_count']
    prices = np.asarray(final['prices'], dtype=float)
    residual = np.asarray(final['residual'], dtype=float)
    if (prices.shape != (n,) or residual.shape != (n,) or not np.isfinite(prices).all()
            or np.any(prices <= 0) or not np.isfinite(residual).all()
            or np.max(np.abs(residual)) > 2e-4
            or not np.isclose(final['score'], np.max(np.abs(residual)), rtol=0, atol=1e-15)):
        raise ValueError('Pilot parent prices/residuals fail the unchanged market gate')
    source = c['source_sha256']
    required = (DRIVER, 'code/model/tools/run_e5f_matched_pf_preference_pilot.py',
        'code/model/tools/test_run_e5f_matched_pf_preference_pilot.py',
        'code/model/tools/e5f_matched_pf_path_root.py', 'code/model/tools/e5f_matched_pf_moments.py')
    if any(name not in source for name in required):
        raise ValueError('Missing pilot/observer/root source pins')
    for name, expected in parent['source_sha256'].items():
        if name != DRIVER and source.get(name) != expected:
            raise ValueError(f'Pilot scientific source changed: {name}')
    if c.get('reviewed_preference_driver_change') != dict(path=DRIVER,
            from_sha256=parent['source_sha256'][DRIVER], to_sha256=source[DRIVER], scope=CHANGE_SCOPE):
        raise ValueError('Exact reviewed optional-driver source change required')
    count = c['path_date_count']
    if count != n and not (count == 6 and c.get('scope') == 'six_date_plumbing_smoke_only'):
        raise ValueError('Only the parent horizon or explicit six-date plumbing smoke is allowed')
    if c['preconditioner'] == 'verified_parent_broyden':
        J = np.asarray(history['final_jacobian'], dtype=float)
        if count != n or J.shape != (n, n) or not np.isfinite(J).all():
            raise ValueError('Parent Broyden preconditioner requires its complete matching horizon')
    elif c['preconditioner'] == 'diagonal':
        J = None
    else:
        raise ValueError('Explicit parent or diagonal preconditioner required')
    return prices[:count].copy(), J, residual[:count].copy()


def load_parent(c, arm):
    names = {'parent_contract': 'contract.json', 'parent_history': 'root_history.json',
             'parent_summary': 'summary.json'}
    directory = Path(c['parent_contract']).parent
    values = []
    for key, filename in names.items():
        path = Path(c[key])
        if not path.is_absolute() or path != directory / filename:
            raise ValueError('Parent receipts must share one absolute completed-root directory')
        primitive.verify(path, c[key + '_sha256'])
        values.append(json.loads(path.read_text()))
    validated = validate_parent(*values, c, arm)
    trial = Path(values[2]['best']['payload']['directory'])
    if not trial.is_absolute() or trial.parent != directory or not trial.name.startswith('evaluation_'):
        raise ValueError('Parent final evaluation must belong to the pinned completed-root directory')
    primitive.verify(trial/'summary.json', values[2]['best']['payload']['summary_sha256'])
    final_summary = json.loads((trial/'summary.json').read_text())
    if (final_summary.get('mapping_valid') is not True
            or final_summary.get('target_fingerprint') != c['target_fingerprint']
            or not all(row['passed'] for row in final_summary['gates'].values())):
        raise ValueError('Parent final evaluation mapping/target receipt mismatch')
    return (*validated, final_summary['artifact_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', required=True, type=Path)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--arm', required=True, choices=('sequential', 'nested'))
    parser.add_argument('--shape', required=True, type=float, choices=SHAPES)
    parser.add_argument('--mode', required=True, choices=('replay', 'conditional', 'root'))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    c, originating = joined.load_smoke_contract(args.contract, args.contract_sha256,
        args.arm, maximum_seconds=10800)
    if (c.get('experiment') != 'fixed_endpoint_preference_shape_pilot'
            or c.get('shape_coefficient') != args.shape or c.get('pilot_mode') != args.mode
            or c.get('probe_coordinate') != -1 or c.get('probe_log_step') != .01
            or c.get('terminal_preference_rule') != 'hold_normalized_2023_intercept'
            or type(c.get('maximum_path_evaluations')) is not int
            or c['maximum_path_evaluations'] != (3 if args.mode == 'root' else 1)
            or (args.mode == 'replay' and args.shape != 0.)):
        raise ValueError('Explicit fixed-endpoint, fixed-budget pilot contract required')
    prices, jacobian, parent_residual, parent_artifacts = load_parent(c, args.arm)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(out)
    out.mkdir(parents=True, exist_ok=True)
    started, stop = time.monotonic(), threading.Event()
    progress = dict(phase='load', evaluation=0, arm=args.arm, shape=args.shape)
    def save(name, value):
        pf.write_json(out/name, value)
    def heartbeat():
        while not stop.wait(15):
            elapsed = time.monotonic() - started
            save('heartbeat.json', dict(progress, elapsed_seconds=elapsed))
            if elapsed > c['seconds']:
                save('failure.json', dict(progress, error='pilot wall-time budget exhausted'))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        seed = baseline.load_seed(c, originating, args.arm)
        old, _ = baseline.load_normalized(c, args.arm)
        terminal_packet = baseline.load_terminal(c, args.arm, old)
        terminal = SimpleNamespace(parameters=terminal_packet['parameters'],
            asset_price=float(terminal_packet['policy'].price[0]))
        psi = preference_shape(old.psi_path, c['path_date_count'], args.shape)
        save('contract.json', dict(c, contract_sha256=args.contract_sha256, psi_path=psi,
            parent_price_role='numerical starting vector only',
            preconditioner_role='approximate matrix; not a derivative measured at this candidate',
            shape_status='fixed diagnostic coordinate; no parameter estimation',
            scope=c.get('scope', 'fixed-endpoint inherited-parameter PF diagnostic')))
        pf.write_csv(out/'diagnostic_coordinates.csv', [dict(parameter='preference_shape_coefficient',
            value=args.shape, is_free_parameter=False,
            status='fixed_diagnostic_scenario_coordinate_not_estimated')])
        del old, terminal_packet
        def project(p):
            return rent_domain.project_price_path_to_positive_rents(p, terminal=terminal,
                minimum_rent_share=1e-6)[0]
        def evaluate(p):
            progress.update(phase='path', evaluation=progress['evaluation'] + 1)
            trial = out / f"evaluation_{progress['evaluation']:03d}"
            trial.mkdir()
            result = baseline.run_history_probe(seed, c, args, trial, progress,
                lambda name, value: pf.write_json(trial/name, value), time.monotonic(),
                prices_override=p, psi_path_override=psi)
            row = dict(evaluation=progress['evaluation'], score=result['maximum_market_residual'],
                directory=str(trial), shape=args.shape, elapsed_seconds=time.monotonic()-started)
            save('latest_completed.json', row)
            return dict(residual=result['residual'], mapping_valid=result['mapping_valid'],
                payload=dict(directory=str(trial), summary_sha256=primitive.digest(trial/'summary.json'),
                    terminal_distance=result['terminal_distance']))
        seed_difference = None
        if args.mode == 'root':
            def record(row):
                save('latest_completed.json', row)
                if row.get('new_best'):
                    save('best_so_far.json', row)
            result = price_root.solve_price_path(initial_prices=prices, evaluate=evaluate,
                project=project, slope=1.63, market_tolerance=2e-4, max_log_step=.10,
                damping=1., max_evaluations=3, deadline_monotonic=started+c['seconds'],
                max_condition_number=1e10, worsening_factor=1.5,
                final_reproduction_tolerance=2e-10, initial_jacobian=jacobian, callback=record)
            save('root_history.json', result)
            converged = bool(result['converged'])
            status = result['status']
            replay_gap = result['final_reproduction_max_abs']
        else:
            result = evaluate(prices)
            seed_difference = float(np.max(np.abs(np.asarray(result['residual']) - parent_residual)))
            full_parent_horizon = c.get('scope') != 'six_date_plumbing_smoke_only'
            replay_gap = seed_difference if args.mode == 'replay' and full_parent_horizon else None
            if replay_gap is not None and replay_gap > 2e-10:
                raise RuntimeError('Inherited price/preference replay differs from pinned parent residuals')
            if args.mode == 'replay' and full_parent_horizon:
                for name in ('target_fit.csv', 'parameters.csv', 'measurement.json', 'transition_path.csv'):
                    primitive.verify(Path(result['payload']['directory'])/name, parent_artifacts[name])
            save('best_so_far.json', result)
            converged = False
            status = 'passed_conditional_preference_pilot'
        save('summary.json', dict(status=status, arm=args.arm, shape=args.shape,
            mode=args.mode, evaluations=progress['evaluation'], elapsed_seconds=time.monotonic()-started,
            scope=c.get('scope', 'fixed-endpoint inherited-parameter PF diagnostic'),
            finite_horizon_market_converged=converged, residual_replay_gap=replay_gap,
            residual_difference_from_parent_at_seed=seed_difference,
            historical_equilibrium_certified=False, horizon_extension_verified=False,
            calibrated_history=False, production_promoted=False,
            outstanding=['terminal and horizon convergence', 'empirical observation mapping',
                         'full target and checkpoint replay verification']))
    except Exception as error:
        save('failure.json', dict(progress, error=repr(error), elapsed_seconds=time.monotonic()-started))
        raise
    finally:
        stop.set()


if __name__ == '__main__':
    main()
