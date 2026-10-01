#!/usr/bin/env python3
"""Bounded historical numerical diagnostic; never a current-floor transition.

Fresh unchanged endpoint root, one permanent illustrative 2007 psi change,
then at most three six-date fixed-price pension-accounting mappings. No fit,
long horizon, credit repair, gate change or reference promotion is permitted.
"""
import argparse
import importlib
import json
import os
from pathlib import Path
import signal
import sys
import threading
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REFERENCE_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
PSI = .1355551166583114


def reserve_mapping(completed):
    from selected_adapter import require
    require(type(completed) is int and 0 <= completed < 3, 'Three-map diagnostic cap reached')


def import_pinned():
    from selected_adapter import require
    sys.path.insert(0, str(HERE / 'pinned_tools'))
    pins = json.loads((HERE / 'legacy_source_pins.json').read_text())
    modules = {}
    for name in pins:
        module = importlib.import_module(name[:-3])
        require(Path(module.__file__).resolve() == (HERE / 'pinned_tools' / name).resolve(),
                'Imported numerical module outside snapshot: ' + name)
        modules[name] = module
    plan = modules['run_e5f_preference_estimation.py'].draft_plan('one_permanent')
    require(plan['endpoint']['max_evaluations'] == 24, 'Original endpoint call cap changed')
    require(callable(modules['run_e5f_preference_estimation.py'].NativeEstimator.endpoint), 'Native endpoint callable missing')
    return modules


def preflight():
    # The adapter helpers have no native imports or numerical side effects.
    from selected_adapter import sha, require
    pins = json.loads((HERE / 'legacy_source_pins.json').read_text())
    for name, digest in pins.items():
        require(sha(HERE / 'pinned_tools' / name) == digest, 'Legacy snapshot changed: ' + name)
    execution = json.loads((HERE / 'execution_pins.json').read_text())
    for name, digest in execution.items():
        require(sha(HERE / name) == digest, 'Readiness execution source changed: ' + name)
    require(sha(ROOT / 'output/model/fertility_identification_20260928/fixed_reference_manifest.json') == REFERENCE_SHA,
            'Reference manifest changed')
    return dict(status='PASS', model_calls=0, economic_change={'psi_child': {'from': PSI, 'to': PSI * 1.001,
                'classification': 'experimental', 'dates': [2007], 'expectations': 'permanent until next surprise'}},
                all_other_primitives='unchanged authenticated block0506', horizon=6, maximum_six_date_maps=3,
                maximum_endpoint_stationary_calls=24, maximum_endpoint_one_step_maps=1,
                maximum_native_policy_calls_conservative=62, total_seconds=3550, numerical_threads=1,
                cache_max_bytes=64 * 1024**3, historical_repayment_limitation=True,
                current_floor_validation=False, production_ready=False, joint_price_pension_root_verified=False,
                budget_estimate={'previous_three_six_date_maps_seconds': 1170.57,
                    'endpoint_setup_render_reserve_seconds': 2379.43,
                    'changed_endpoint_observed_timing': None, 'completion_guaranteed': False})


def run(output):
    started = time.monotonic()
    from selected_adapter import require
    proof = preflight()
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    for name in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS'):
        require(os.environ.get(name) == '1', 'One numerical thread required: ' + name)
    require(not output.exists(), 'Output already exists')
    output.mkdir(parents=True)
    modules = import_pinned()
    inner = modules['run_e5f_preference_transition.py']
    estimator = modules['run_e5f_preference_estimation.py']
    diagnostic = modules['run_e5f_preference_budget_diagnostic.py']
    deadline = started + 3550; stop = threading.Event(); native = None; last_record = None
    inner.write(output / 'preflight.json', proof)
    inner.write(output / 'latest_completed.json', dict(status='not_started'))
    inner.write(output / 'best_so_far.json', dict(status='no_completed_mapping'))
    def heartbeat():
        while not stop.wait(30):
            inner.write(output / 'heartbeat.json', dict(epoch=time.time(), elapsed=time.monotonic()-started, phase='active'))
    inner.write(output / 'heartbeat.json', dict(epoch=time.time(), elapsed=0, phase='setup'))
    worker = threading.Thread(target=heartbeat, daemon=True); worker.start()
    def timeout(*_):
        raise TimeoutError('Hard 3550-second historical diagnostic budget exhausted')
    previous = signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, max(.001, deadline-time.monotonic()))
    # Independent absolute watchdog survives the unchanged inner SIGALRM guards.
    previous_watchdog = signal.signal(signal.SIGUSR1, timeout)
    def watchdog():
        if not stop.wait(max(0, deadline-time.monotonic())):
            os.kill(os.getpid(), signal.SIGUSR1)
    deadline_worker = threading.Thread(target=watchdog, daemon=True); deadline_worker.start()
    try:
        import numpy as np
        manifest, packet, evaluator = inner.load_reference(output / 'reference')
        require(float(packet['parameters'].psi_child) == PSI, 'Reference psi differs')
        plan = estimator.draft_plan('one_permanent')
        plan['housing'] = 'fixed_stock'
        plan['budget'].update(total_seconds=3550, candidate_seconds=3550, endpoint_seconds=2000,
                              mapping_seconds=1200, path_seconds=1600, jacobian_seconds=1)
        plan['path']['cache_max_bytes'] = 64 * 1024**3
        native = estimator.NativeEstimator(plan, output, manifest, packet, evaluator)
        native.deadline = started + 3550; native.candidate_deadline = native.deadline
        terminal, endpoint = native.endpoint(PSI * 1.001)
        inner.write(output / 'latest_completed.json', dict(phase='endpoint', endpoint=endpoint))
        path = inner.shock_path(dict(kind='one_permanent', years=[2007], period_years=4,
            levels=[PSI * 1.001], expectations='perfect_foresight_announced_at_start',
            interpretation='explicit_illustrative_transport'), 6)
        q = float(packet['solution'].p_eq[0]); b = float(packet['parameters'].pension)
        numerics = plan['path']; best = float('inf'); maps = []
        def mapping(name, pensions):
            nonlocal best, last_record
            reserve_mapping(len(maps))
            result, record = native.guarded(1200, lambda: inner.mapping(packet, evaluator, terminal, endpoint,
                np.full(6, q), np.asarray(pensions), path, 'fixed_stock', output / name,
                64 * 1024**3, capture=True, measure_fertility=True, start_year=2007))
            diagnostic.validate_record_residuals(record, 6)
            check = inner.terminal_checks(packet, evaluator, terminal, endpoint, result, path, numerics)
            accepted = all(record['gates'].values()) and check['all_checks_pass'] and max(map(abs, record['market_residual'])) <= 2e-4 and max(map(abs, record['fiscal_residual'])) <= 1e-6
            receipt = dict(name=name, accepted=accepted, **diagnostic.summarize(record, inner.plain(check)))
            maps.append(receipt); last_record = record
            inner.write(output / 'latest_completed.json', receipt)
            score = max(receipt['residual_ratios'].values())
            if score < best:
                best = score; inner.write(output / 'best_so_far.json', receipt)
            inner.dump_checkpoint(output / (name + '_checkpoint.pkl.gz'), dict(terminal_state=result.terminal_state,
                prices=np.full(6, q), pensions=pensions, psi_path=path, rows=result.rows))
            return result, record, check, accepted
        base_result, base, _, _ = mapping('changed_psi_input', [b] * 6)
        del base_result
        updated = diagnostic.budget_update(base['rows'])
        trial_result, trial, trial_terminal, accepted = mapping('pension_trial', updated)
        certified = False; differences = None
        if accepted:
            replay_result, replay, replay_terminal, replay_ok = mapping('fresh_replay', updated)
            pf = evaluator.rt['primitive'].pf
            differences = dict(market=max(abs(a-b) for a,b in zip(trial['market_residual'], replay['market_residual'])),
                fiscal=max(abs(a-b) for a,b in zip(trial['fiscal_residual'], replay['fiscal_residual'])),
                fertility=max(abs(float(a['period_tfr_topcode_adjusted'])-float(b['period_tfr_topcode_adjusted'])) for a,b in zip(trial['fertility'], replay['fertility'])),
                g_pre=float(np.max(np.abs(trial_result.terminal_state.g_pre-replay_result.terminal_state.g_pre))),
                queues=max(float(np.max(np.abs(pf.birth_queue_values(getattr(trial_result.terminal_state, name))-pf.birth_queue_values(getattr(replay_result.terminal_state, name))))) for name in ('scheduled_entries', 'scheduled_raw_entries')))
            certified = bool(replay_ok and trial['gates'] == replay['gates'] and inner.plain(trial_terminal) == inner.plain(replay_terminal) and max(differences.values()) <= 1e-10)
        visual = native.guarded(max(1, native.deadline-time.monotonic()), lambda: inner.render_diagnostics(
            trial['diagnostic_packets'], output / 'diagnostics', evaluator.rt['audit'], manifest['standard_diagnostic_names']))
        inner.write(output / 'complete.json', dict(numerical_certified=certified, maps=maps, diagnostics=visual,
            replay_max_difference=differences, elapsed_seconds=time.monotonic()-started, production_ready=False,
            current_floor_validation=False, full_horizon_verified=False, fitted_shocks=False,
            joint_price_pension_root_verified=False, interpretation='fixed-price historical diagnostic; does not validate joint path solver'))
        require(certified, 'Changed-psi six-date diagnostic was not certified; stop without repair')
    except BaseException as exc:
        # A completed input/trial retains its capture packet even if a later map fails.
        if last_record is not None and native is not None and time.monotonic() < deadline and not (output / 'diagnostics').exists():
            try:
                visual = native.guarded(deadline-time.monotonic(), lambda: inner.render_diagnostics(
                    last_record['diagnostic_packets'], output / 'diagnostics', evaluator.rt['audit'], manifest['standard_diagnostic_names']))
                inner.write(output / 'diagnostics_receipt.json', visual)
            except BaseException as render_exc:
                inner.write(output / 'diagnostics_failure.json', dict(error=str(render_exc), captured_packets_retained=True))
        inner.write(output / 'failure.json', dict(error_type=type(exc).__name__, error=str(exc),
            elapsed_seconds=time.monotonic()-started, production_ready=False, current_floor_validation=False))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0); signal.signal(signal.SIGALRM, previous)
        stop.set(); worker.join(timeout=1); deadline_worker.join(timeout=1)
        signal.signal(signal.SIGUSR1, previous_watchdog)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--native-preflight', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.preflight or args.native_preflight:
        proof = preflight()
        if args.native_preflight:
            import_pinned()
            proof['native_source_imports_verified'] = True
        print(json.dumps(proof, indent=2, sort_keys=True))
    else:
        parser.error('--output is required') if args.output is None else run(args.output)


if __name__ == '__main__':
    main()
