#!/usr/bin/env python3
"""One explicit 24-date native map from a pinned frozen v5 reference.

Diagnostic only: no endpoint solve, scalar fit, root, replay, or adoption.
The supplied unverified solve contributes only its genuine terminal value V;
its separate completed stationary receipts authenticate the endpoint boundary.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import pickle
import sys
import threading
import time
import traceback
from types import SimpleNamespace

HORIZON = 24
PSI = .14736308634876963
QT = .6482672107267026
PENSION = .917784047463731
GROWTH = 1.1802793938887526
RENT = .13689881249028354
REFERENCE_SHA = 'a2c6b2b266ef524fa8d020da1623b5d38e55caa4acfd364a07cf2260762244fe'
TERMINAL_SHA = '295a586d2dbd52e931b5224f0233623db3f65f52f63e485b2413c66704190923'
MANIFEST_SHA = '14c54b9d4fbb4a80ca8c71121268ac0fff8cf79868316cd348c01528782a0f46'
INITIAL_SHA = '6a7cdad9dc8bd32a22b832756b746c106acb070cd62b19770a0093be743c10b3'
TERMINAL_V_SHA = 'fbe3eacff8a14d71c663a766fa1a238f350320df215086d874e1c1c6de687d20'
PIN_NAMES = ('checkpoint', 'raw_receipt', 'stationary', 'root', 'latest_completed', 'one_step_record')
ENDOGENOUS_FIELDS = {'birth_count_action_probs', 'birth_count_realized_probs',
                     'birth_count_pre_distribution', 'birth_count_post_distribution',
                     'birth_count_first_birth_tagged_distribution'}


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def pinned(item, root=None):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Exact path/SHA pin required')
    path = Path(item['path']).resolve(strict=True)
    require(path.is_file() and path.is_absolute(), 'Pinned file missing')
    if root is not None:
        require(path.is_relative_to(root), 'Pin escapes frozen package: ' + str(path))
    require(digest(path) == item['sha256'], 'Pinned hash differs: ' + str(path))
    return path


def read_pin(item, root=None):
    return json.loads(pinned(item, root).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False, default=str) + '\n')
    temporary.replace(path)


def expected_paths():
    base = RENT / (GROWTH - 1.)
    return [base + (QT - base) * GROWTH ** -(HORIZON - t) for t in range(HORIZON)]


def validate_config(config):
    require(set(config) == {'package_root', 'fit_manifest', 'driver', 'reference_receipt',
                           'endpoint_pins', 'q_path', 'pension_path', 'psi_path'}, 'Config keys differ')
    root = Path(config['package_root']).resolve(strict=True)
    require(root.is_dir() and root.name == 'execution_smoke_v5', 'Frozen v5 package required')
    require(pinned(config['driver']).resolve() == Path(__file__).resolve(), 'Executing wrapper differs from self pin')
    manifest_path = pinned(config['fit_manifest'], root)
    require(manifest_path == root / 'inputs/fit_manifest.json' and digest(manifest_path) == MANIFEST_SHA,
            'Frozen fit manifest identity differs')
    source = root / 'frozen/source/code/model/experiments/birth_count_choice'
    require(source.is_dir(), 'Frozen model source missing')
    receipt_path = pinned(config['reference_receipt'], root)
    require(receipt_path == root / 'results/fit_fit_v1/run/reference/reference_reconstruction.json',
            'Wrong reconstructed reference receipt')
    receipt = json.loads(receipt_path.read_text())
    require(receipt['checkpoint_sha256'] == REFERENCE_SHA, 'Wrong reference native checkpoint')
    pins = config['endpoint_pins']
    require(isinstance(pins, dict) and set(pins) == set(PIN_NAMES), 'Complete endpoint pin set required')
    paths = {name: pinned(pin, root) for name, pin in pins.items()}
    require(paths['checkpoint'].name == 'native_solve_unverified.pkl.gz' and
            pins['checkpoint']['sha256'] == TERMINAL_SHA, 'Wrong terminal native checkpoint')
    require(paths['raw_receipt'] == paths['checkpoint'].with_suffix('').with_suffix('.json'),
            'Raw solve receipt must accompany checkpoint')
    require(paths['stationary'].name == 'stationary.json' and
            paths['root'].name == 'root.json' and
            paths['latest_completed'].name == 'latest_completed.json' and
            paths['one_step_record'].name == 'native_record.json', 'Endpoint proof names differ')
    require(paths['root'].parent == paths['latest_completed'].parent and
            paths['one_step_record'].parent == paths['root'].parent / 'one_step' and
            paths['stationary'].parent == paths['checkpoint'].parent,
            'Endpoint proof locations differ')
    for key, target in [('q_path', expected_paths()), ('pension_path', [PENSION] * HORIZON),
                        ('psi_path', [PSI] * HORIZON)]:
        values = config[key]
        require(isinstance(values, list) and len(values) == HORIZON, 'Path length differs: ' + key)
        require(all(type(x) in (float, int) and math.isfinite(x) and x > 0 for x in values),
                'Nonfinite/nonpositive path: ' + key)
        require(all(abs(a-b) <= 2e-14 * max(1., abs(b)) for a, b in zip(values, target)),
                'Explicit path formula or level differs: ' + key)
    return root, source, paths


def validate_endpoint(config, paths, identity):
    raw = json.loads(paths['raw_receipt'].read_text())
    root = json.loads(paths['root'].read_text())
    latest = json.loads(paths['latest_completed'].read_text())
    one = json.loads(paths['one_step_record'].read_text())
    stationary = json.loads(paths['stationary'].read_text())
    require(stationary == latest, 'Stationary point differs from completed root record')
    raw_checkpoint = raw.get('checkpoint', {})
    require(raw.get('schema') == 'current_floor_unverified_native_solve_v1' and
            raw_checkpoint.get('sha256') == TERMINAL_SHA and
            Path(raw_checkpoint.get('path', '')).name == paths['checkpoint'].name and
            raw.get('identity') == identity and raw.get('reference_verified') is False,
            'Raw terminal solve receipt differs')
    require(raw.get('price') == QT and raw.get('psi_child') == PSI, 'Raw terminal levels differ')
    require(root.get('converged') is True and root.get('status') == 'converged' and
            root.get('final_reproduction_max_abs') == 0. and
            root['final']['prices'] == [QT], 'Stationary endpoint root proof failed')
    require(latest.get('price') == QT and latest.get('psi_child') == PSI and
            latest.get('pension') == PENSION and latest.get('accounting_valid') is True,
            'Completed endpoint levels/accounting differ')
    require(stationary.get('price') == QT and stationary.get('psi_child') == PSI and
            abs(stationary.get('renewal_residual', math.inf)) <= 1e-6,
            'Stationary endpoint proof differs')
    require(one.get('accounting_valid') is True and all(one.get('gates', {}).values()) and
            all(math.isfinite(float(x)) for k in ('market_residual', 'fiscal_residual') for x in one[k]) and
            max(map(abs, one['market_residual'])) <= 2e-4 and
            max(map(abs, one['fiscal_residual'])) <= 1e-6,
            'One-step native gates/residuals failed')
    with gzip.open(paths['checkpoint'], 'rb') as stream:
        saved = pickle.load(stream)
    require(saved['identity'] == identity and saved['stage'] == 'native_solve_completed_reconstruction_pending' and
            saved['price'] == QT and
            float(saved['solution'].p_eq[0]) == QT and
            float(saved['parameters'].psi_child) == PSI,
            'Native terminal checkpoint identity/levels differ')
    return saved, dict(root=root, latest_completed=latest, stationary=stationary,
                       one_step_record=one, raw_receipt=raw)


def public_parameters(parameters):
    import numpy as np
    return {key: value for key, value in vars(parameters).items()
            if not key.startswith('_') and key != 'native_inherited_distribution_evidence_dir'
            and key not in ENDOGENOUS_FIELDS}


def same_public_parameters(endpoint, reference):
    import numpy as np
    a, b = public_parameters(endpoint), public_parameters(reference)
    require(set(a) == set(b), 'Terminal public parameter keys differ')
    for key in a:
        if key == 'psi_child':
            require(float(a[key]) == PSI, 'Terminal psi differs')
        else:
            require(np.array_equal(np.asarray(a[key]), np.asarray(b[key])),
                    'Terminal public parameter differs: ' + key)


def bind_reference_cache_outputs(rt, receipt_pin, package_root, folder):
    """Restore five authenticated solve outputs before the unchanged reference gate.

    The constructor can recreate the economic primitives but not these arrays,
    which the frozen household and distribution routines populate during a solve.
    No solve, state, queue, primitive, or acceptance check is changed here.
    """
    import numpy as np
    receipt = read_pin(receipt_pin, package_root)
    require(receipt.get('schema') == 'current_floor_reference_reconstruction_v1' and
            receipt.get('status') == 'passed' and receipt.get('identity') == rt.identity(),
            'Reference reconstruction receipt identity differs')
    checkpoint_pin = receipt.get('checkpoint')
    require(isinstance(checkpoint_pin, dict) and checkpoint_pin.get('sha256') == REFERENCE_SHA and
            receipt.get('checkpoint_sha256') == REFERENCE_SHA,
            'Authenticated reference native checkpoint differs')
    checkpoint = pinned(checkpoint_pin, package_root)
    with gzip.open(checkpoint, 'rb') as stream:
        packet = pickle.load(stream)
    require(np.array_equal(packet['b_grid'], rt.grid), 'Restored reference grid differs')
    saved, current = vars(packet['parameters']), vars(rt.P)
    require(ENDOGENOUS_FIELDS <= set(saved) and ENDOGENOUS_FIELDS <= set(current),
            'Expected five native cached output arrays missing')
    a, b = public_parameters(packet['parameters']), public_parameters(rt.P)
    require(set(a) == set(b), 'Reference public primitive keys differ')
    require(float(a['psi_child']) == float(b['psi_child']) == .17892072066041628,
            'Baseline reference psi differs')
    for key in a:
        require(np.array_equal(np.asarray(a[key]), np.asarray(b[key])),
                'Reference public primitive differs: ' + key)
    def array_receipt(value):
        array = np.asarray(value)
        require(array.dtype.kind in 'biuf' and np.isfinite(array).all(),
                'Nonfinite native cached output array')
        return dict(dtype=str(array.dtype), shape=list(array.shape),
                    sha256=hashlib.sha256(array.tobytes()).hexdigest())
    before = {name: array_receipt(current[name]) for name in sorted(ENDOGENOUS_FIELDS)}
    source = {name: array_receipt(saved[name]) for name in sorted(ENDOGENOUS_FIELDS)}
    for name in sorted(ENDOGENOUS_FIELDS):
        setattr(rt.P, name, copy.deepcopy(saved[name]))
    after = {name: array_receipt(getattr(rt.P, name)) for name in sorted(ENDOGENOUS_FIELDS)}
    require(after == source, 'Cached native output binding differs from authenticated checkpoint')
    receipt = dict(kind='diagnostic_reference_cached_native_outputs_v1',
                   reference_receipt=receipt_pin, checkpoint=checkpoint_pin,
                   source_assignments=['frozen household.py:1086-1087',
                                       'frozen distribution.py:1339-1341'],
                   fields={name:dict(before=before[name], source=source[name], after=after[name])
                           for name in sorted(ENDOGENOUS_FIELDS)},
                   primitive_comparison='all other public fields exact',
                   reference_restore_and_gates_still_required=True,
                   native_calls=getattr(rt, 'total_native_calls', 0))
    write(Path(folder) / 'reference_cache_binding.json', receipt)
    return receipt


def restore_reference_with_bridge_redirect(rt, receipt_pin, package_root, folder):
    """Preserve frozen comparisons while writing two bridge receipts into output."""
    receipt = read_pin(receipt_pin, package_root)
    require(receipt.get('schema') == 'current_floor_reference_reconstruction_v1' and
            receipt.get('status') == 'passed' and receipt.get('identity') == rt.identity(),
            'Reconstruction receipt identity differs before bridge redirect')
    reports = [Path(path).resolve(strict=True) for path in receipt['reports']]
    require(len(reports) == 2 and len(set(reports)) == 2 and
            all(path.is_relative_to(package_root) for path in reports),
            'Exactly two immutable reference reports required')
    by_target = {}
    for index, report in enumerate(reports):
        require(f'repeat_{index}' in report.parts and report.parts[-2:] == ('phase_b_ge','selected_root'),
                'Original repeat report order/location differs')
        by_target[report / 'saved_platform_bridge.json'] = Path(folder) / f'repeat_{index}_bridge.json'
    comparison = rt.runner.compare_repeated
    namespace = comparison.__func__.__globals__
    require(Path(namespace['__file__']).resolve() ==
            package_root / 'frozen/source/code/model/experiments/birth_count_choice/transition_runtime.py',
            'Frozen comparison implementation differs')
    original_write = namespace['write']
    seen = {}
    def redirected_write(path, payload):
        target = Path(path).resolve()
        if target not in by_target:
            return original_write(path, payload)
        require(target not in seen, 'Duplicate saved-platform bridge write')
        require(isinstance(payload, dict) and
                payload.get('schema') == 'current_estate_a_saved_platform_bridge_v1' and
                payload.get('status') == 'passed' and payload.get('identity') == rt.identity(),
                'Saved-platform bridge payload differs')
        destination = by_target[target]
        original_write(destination, payload)
        seen[target] = dict(original_path=str(target), redirected_path=str(destination),
                            redirected_sha256=digest(destination))
    try:
        namespace['write'] = redirected_write
        restored = rt.load_reconstructed_reference(receipt_pin, folder)
        require(set(seen) == set(by_target), 'Both authenticated bridge writes required')
        write(Path(folder) / 'bridge_redirect.json',
              dict(kind='diagnostic_reference_bridge_redirect_v1',
                   writes=[seen[path] for path in by_target],
                   original_comparisons_unchanged=True, original_reports_unchanged=True,
                   native_calls=getattr(rt, 'total_native_calls', 0)))
        return restored
    finally:
        namespace['write'] = original_write


def run(config, output):
    root, source, paths = validate_config(config)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    deadline = started + 1200.
    phase = 'preflight'
    runtime = None
    stop = threading.Event()
    progress_lock = threading.Lock()
    def progress(name, **extra):
        with progress_lock:
            write(output / 'progress.json', dict(phase=name, epoch=time.time(), elapsed_seconds=time.monotonic()-started,
                                                 actual_native_calls=getattr(getattr(runtime, 'rt', None), 'total_native_calls', 0), **extra))
    def pulse():
        while not stop.wait(180.):
            progress('native_in_progress' if phase == 'mapping' else phase)
    thread = threading.Thread(target=pulse, name='one-map-heartbeat', daemon=True)
    thread.start()
    try:
        progress(phase, diagnostic_only=True, root_converged=False, empirical_fitted=False,
                 scientific_validation=False, production_ready=False)
        sys.path.insert(0, str(source))
        import two_shock as d
        import two_shock_runtime as native_module
        require(Path(d.__file__).resolve() == source / 'two_shock.py' and
                Path(native_module.__file__).resolve() == source / 'two_shock_runtime.py',
                'Mutable or cached frozen driver import')
        plan = read_pin(config['fit_manifest'], root)
        with native_module.retained.watchdog(deadline):
            preflight = d.preflight(plan)
            require(preflight['native_calls'] == 0, 'Preflight made native calls')
            phase = 'constructor'; progress(phase)
            runtime = native_module.NativeRuntime(plan, output / 'runtime')
            rt = runtime.rt
            require(rt.total_native_calls == 0, 'Constructor made native calls')
            phase = 'bind_reference_cache_outputs'; progress(phase)
            bind_reference_cache_outputs(rt, config['reference_receipt'], root, output)
            require(rt.total_native_calls == 0, 'Reference cache binding made native calls')
            phase = 'restore_reference'; progress(phase)
            restore_reference_with_bridge_redirect(rt, config['reference_receipt'], root,
                                                   output / 'reference_restored')
            require(rt.total_native_calls == 0, 'Reference restore made native calls')
            initial = copy.deepcopy(rt.initial_state)
            import numpy as np
            population_hash = hashlib.sha256(np.asarray(initial.g_pre).tobytes()).hexdigest()
            require(population_hash == INITIAL_SHA, 'Initial population hash differs')
            initial_hash = d.state_hash(initial, rt.pf.birth_queue_values)
            phase = 'verify_endpoint'; progress(phase, initial_population_sha256=population_hash,
                                                initial_state_sha256=initial_hash)
            saved, proofs = validate_endpoint(config, paths, rt.identity())
            require(native_module.retained.stationary_mapping_valid(proofs['stationary']) and
                    native_module.retained.stationary_mapping_valid(proofs['latest_completed']),
                    'Original stationary native gates failed')
            require(np.array_equal(saved['b_grid'], rt.grid), 'Terminal native grid differs')
            same_public_parameters(saved['parameters'], rt.packet['parameters'])
            terminal_v_sha = hashlib.sha256(np.asarray(saved['solution'].V).tobytes()).hexdigest()
            require(terminal_v_sha == TERMINAL_V_SHA, 'Pre-reviewed terminal V bytes differ')
            terminal = dict(evaluation=SimpleNamespace(policy=SimpleNamespace(V=saved['solution'].V)))
            endpoint = dict(price=QT, population_scale=proofs['latest_completed']['population_scale'])
            write(output / 'terminal_extraction.json', dict(kind='zero_advance_raw_native_solution_V_only',
                  raw_solve=config['endpoint_pins']['checkpoint'], raw_receipt=config['endpoint_pins']['raw_receipt'],
                  completed_proofs={k: config['endpoint_pins'][k] for k in ('stationary','root','latest_completed','one_step_record')},
                  terminal_V_sha256=terminal_v_sha, terminal_certification_performed=False, native_calls=0))
            phase = 'mapping'; progress(phase, q0=config['q_path'][0], q_terminal=QT)
            original_guard = rt._guard_native_call
            def counted_guard():
                original_guard()
                progress('native_call', actual_call_started=True)
            rt._guard_native_call = counted_guard
            try:
                with rt.native_budget(deadline, 64):
                    native, record = rt.mapping(terminal, endpoint, config['q_path'], config['pension_path'],
                                                config['psi_path'], output / 'map_001',
                                                initial_state=copy.deepcopy(initial), start_year=2007)
            finally:
                rt._guard_native_call = original_guard
            require(d.state_hash(initial, rt.pf.birth_queue_values) == initial_hash and
                    d.state_hash(rt.initial_state, rt.pf.birth_queue_values) == initial_hash,
                    'Reference initial state or queues mutated')
            require(record['accounting_valid'] is True and all(record['gates'].values()),
                    'Dated native accounting gates failed')
            require(record['policy_calls'] == rt.total_native_calls and rt.total_native_calls <= 64,
                    'Actual native call count differs')
            require(len(record['rows']) == HORIZON and len(native.dated_states) == HORIZON,
                    'Dated native map length differs')
            finite = all(math.isfinite(float(x)) for key in ('market_residual','fiscal_residual') for x in record[key])
            require(finite, 'Nonfinite physical residuals')
            result = dict(status='mapping_completed', diagnostic_only=True, empirical_fitted=False,
                          root_converged=False, scientific_validation=False, production_ready=False,
                          terminal_certification_performed=False, identity=rt.identity(),
                          initial_state_sha256=initial_hash, reference_state_unchanged=True,
                          actual_native_calls=rt.total_native_calls, policy_calls=record['policy_calls'],
                          accounting_valid=record['accounting_valid'], gates=record['gates'],
                          market_maximum_residual=max(map(abs,record['market_residual'])),
                          fiscal_maximum_residual=max(map(abs,record['fiscal_residual'])),
                          original_market_tolerance=plan['gates']['market_tolerance'],
                          original_fiscal_tolerance=plan['gates']['fiscal_tolerance'],
                          source_manifest=config['fit_manifest'], endpoint_pins=config['endpoint_pins'],
                          reference_receipt=config['reference_receipt'], q_path=config['q_path'],
                          pension_path=config['pension_path'], psi_path=config['psi_path'])
            write(output / 'result.json', result)
            phase = 'complete'; progress(phase)
            return result
    except BaseException as exc:
        write(output / 'failure.json', dict(status='failed', phase=phase, error=repr(exc),
              traceback=traceback.format_exc(), actual_native_calls=getattr(getattr(runtime,'rt',None),'total_native_calls',0),
              scientific_validation=False, production_ready=False))
        progress('failed')
        raise
    finally:
        stop.set(); thread.join(timeout=1.)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)
    config = json.loads(args.config.read_text())
    run(config, args.output)


if __name__ == '__main__':
    main()
