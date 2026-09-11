"""Pinned checkpoint observation only. No submissions or equilibrium solves."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import threading
import time


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def write(path, value):
    path = Path(path)
    temp = path.with_name(path.name + f'.{threading.get_ident()}.tmp')
    temp.write_bytes(json.dumps(value, indent=2, sort_keys=True, allow_nan=False).encode() + b'\n')
    os.replace(temp, path)


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def verify(path, expected):
    require(digest(path) == expected, f'Fingerprint mismatch: {path}')


def bounded_path(root, relative):
    path = (root / relative).resolve()
    require(path.is_relative_to(root), f'Path outside bundle: {relative}')
    return path


def parent_gates(case, bundle, top, gate):
    """Validate the native two-repetition joint smoke; never impersonate a panel row."""
    import csv
    import importlib.util
    parent_path = bounded_path(bundle, case['parent_contract'])
    parent = json.loads(parent_path.read_text())
    require(gate['status'] == 'verified' and gate['repetitions'] == 2, 'Original joint gate missing')
    require(gate['exact_early_equality'] is True and gate['every_GE_market_verified_both_repetitions'] is True, 'Original joint replay/market gate failed')
    require(gate['source_pins_verified'] == 634 and gate['graph_count'] == 17, 'Original joint source/graph gate failed')
    require(gate['quantity_report_repetition'] == 2, 'Wrong final report repetition')
    require(case['case_id'] == parent['case_id'] == 'joint_smoke_analysis_47', 'Original joint case identity changed')
    verify(parent_path, gate['contract_sha256'])
    require(parent['repetitions'] == 2 and parent['maximum_GE_solves'] == 16, 'Original two-loop contract changed')
    require(parent['run_input_fingerprint'] == gate['input_fingerprint'], 'Joint input fingerprint mismatch')
    actual = json.loads(bounded_path(bundle, case['actual_output_contract']).read_text())
    require(actual == dict(parent, contract_sha256=gate['contract_sha256'], case='new_balanced'), 'Original joint output contract mismatch')
    require(top['status'] == 'passed_initial_candidate_loop' and top['normalized'] is True and top['repetitions'] == 2, 'Original joint smoke failed')
    require(top['stationary_solves'] == gate['stationary_solves'] <= 16 and top['elapsed_seconds'] <= 1800, 'Original joint budget failed')
    smoke = bounded_path(bundle, 'inputs/smoke')
    for name, pin in gate['artifact_sha256'].items():
        if name.endswith(('.json', '.csv')):
            verify(bounded_path(smoke, name), pin)
    spec = importlib.util.spec_from_file_location('original_panel_validator', bounded_path(bundle, 'inputs/panel_validator.py'))
    validator = importlib.util.module_from_spec(spec); spec.loader.exec_module(validator)
    reps = []; early = []; solve_count = 0
    for number in (1, 2):
        p = smoke / f'repetition_{number:02d}'
        final = json.loads((p/'summary.json').read_text())
        observation = json.loads((p/'early_measurement.json').read_text())
        ges = json.loads((p/'stationary_solves.json').read_text())
        wrapper = dict(top, repetitions=1, stationary_solves=len(ges), final=final)
        with (p/'parameters.csv').open() as stream:
            validator.validate_summary(wrapper, parent, observation, list(csv.DictReader(stream)), ges)
        require(final['checkpoint_sha256'] == gate['artifact_sha256'][f'repetition_{number:02d}/initial_state.pkl.gz'], 'Original joint checkpoint claim mismatch')
        reps.append(final); early.append(observation); solve_count += len(ges)
    validator.validate_market(json.loads((smoke/'repetition_02/market_quantity_units.json').read_text()))
    require(early[0] == early[1], 'Original joint early measurements differ')
    for key in ('price', 'legacy_stationary_moments'):
        require(reps[0][key] == reps[1][key], 'Original joint numerical replay differs: '+key)
    require({k:v for k,v in reps[0]['normalization'].items() if k != 'stationary_solve_seconds'} == {k:v for k,v in reps[1]['normalization'].items() if k != 'stationary_solve_seconds'}, 'Original normalization differs')
    require(top['final'] == reps[1] and solve_count == top['stationary_solves'], 'Original aggregate joint summary differs')
    require(case['checkpoint'] == gate['proposal_reuse']['output_path'] + '/initial_state.pkl.gz', 'Joint checkpoint path changed')
    require(case['checkpoint_sha256'] == reps[1]['checkpoint_sha256'], 'Selected joint checkpoint pin differs')
    return parent


# Explicit compiled transport override is confined to the startup-test child.
# The unchanged 13 tests ordinarily choose the uncompiled scatter path.
COMPILED_TESTS = r'''
import sys, unittest
from unittest.mock import patch
import test_e5f_recent_parent_flow_observer as tests
original = tests.model.realize_current_cross_section
def compiled(*args, **kwargs):
    kwargs['use_compiled_scatter'] = True
    return original(*args, **kwargs)
suite = unittest.defaultTestLoader.loadTestsFromModule(tests)
assert suite.countTestCases() == 13
with patch.object(tests.model, 'realize_current_cross_section', compiled):
    result = unittest.TextTestRunner(verbosity=2).run(suite)
assert result.testsRun == 13 and result.wasSuccessful()
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify-only', action='store_true', help='No imports, checkpoint load, tests or observation')
    args = parser.parse_args()
    started = time.monotonic()
    verify(args.contract, args.contract_sha256)
    contract = json.loads(args.contract.read_text())
    bundle = args.contract.resolve().parent
    root = args.source_root.resolve()
    require(contract['schema'] == 'e5f_recent_parent_joint_smoke_probe_v1', 'Wrong schema')
    require(contract['calibrated_smm'] is False and contract['maximum_model_solves'] == 0, 'Read-only contract required')
    require(0 < contract['work_seconds'] <= 220 and contract['hard_seconds'] == 240, 'Wall limit changed')
    verify(Path(__file__), contract['driver_sha256'])
    for name, pin in contract['bundle_input_sha256'].items():
        verify(bounded_path(bundle, name), pin)
    sources = contract['source_sha256']
    actual_sources = {str(p.relative_to(root)) for p in (root / 'code/model').rglob('*')
                      if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc', '.nbc', '.nbi')}
    require(actual_sources == set(sources), 'Archive source inventory is incomplete or has additional files')
    for name, pin in sources.items():
        verify(bounded_path(root, name), pin)
    require(hashlib.sha256(canonical(sources)).hexdigest() == contract['source_fingerprint'], 'Source fingerprint mismatch')
    raw = json.loads(bounded_path(bundle, contract['joint_parent_summary']).read_text())
    panel = json.loads(bounded_path(bundle, contract['joint_smoke_gate_receipt']).read_text())
    empirical = json.loads(bounded_path(bundle, contract['empirical_target_file']).read_text())
    target = next(x for x in empirical['records'] if x['id'] == contract['empirical_target_id'])
    require(target['estimate'] == contract['empirical_target_value'] == 0.16289550916123285, 'Empirical target changed')
    ids = [c['case_id'] for c in contract['cases']]
    require(len(ids) == len(set(ids)) and len(ids) > 0, 'Duplicate or empty cases')
    for case in contract['cases']:
        require(type(case['passes']) is int and case['passes'] in (1, 2), 'Invalid pass count')
        parent = parent_gates(case, bundle, raw, panel)
        require(parent['source_commit'] == '7e872053', 'Wrong inherited source')
        require(all(sources.get(k) == v for k, v in parent['source_sha256'].items()), 'Inherited source changed')
        require(len(parent['source_sha256']) == contract['inherited_source_count'], 'Incomplete inherited manifest')
        if not args.verify_only:
            verify(case['checkpoint'], case['checkpoint_sha256'])
    if args.verify_only:
        print(json.dumps(dict(status='local_preflight_passed', source_files=len(sources),
                             inherited_source_files=contract['inherited_source_count'],
                             checkpoint_bytes_not_checked=True, model_solves=0)))
        return
    args.output.mkdir(parents=True, exist_ok=False)
    out = args.output
    stopped = threading.Event()
    phase = {'phase': 'startup', 'completed_passes': 0}
    def heartbeat():
        while not stopped.wait(15):
            write(out / 'heartbeat.json', {**phase, 'elapsed_seconds': time.monotonic() - started})
    def hard_stop():
        write(out / 'failure.json', dict(status='timeout', elapsed_seconds=time.monotonic() - started, **phase))
        os._exit(124)
    timer = threading.Timer(max(1., 240 - (time.monotonic() - started)), hard_stop)
    timer.daemon = True; timer.start()
    thread = threading.Thread(target=heartbeat, daemon=True); thread.start()
    write(out / 'contract.json', contract)
    write(out / 'best_so_far.json', dict(status='not_applicable', reason='Observation has no optimization or best candidate'))
    write(out / 'heartbeat.json', phase)
    def remaining():
        seconds = contract['work_seconds'] - (time.monotonic() - started)
        require(seconds > 0, 'Work deadline reached; reporting reserve retained')
        return seconds
    try:
        phase['phase'] = 'compiled_startup_tests'
        env = dict(os.environ, NUMBA_DISABLE_JIT='0', PYTHONPATH=f'{root}/code/model/tools:{root}/code/model')
        with (out / 'compiled_startup_tests.log').open('w') as stream:
            subprocess.run([sys.executable, '-c', COMPILED_TESTS], env=env, cwd=root,
                           stdout=stream, stderr=subprocess.STDOUT, timeout=min(110., remaining()), check=True)
        for path in (root / 'code/model', root / 'code/model/tools'):
            sys.path.insert(0, str(path))
        import run_e5f_open_population_transition as transition
        from intergen_eqscale_seq_optimized import solver as model
        from e5f_recent_parent_flow_observer import observe_recent_parent_flow, SNAPSHOT, AGE_PROJECTION
        # Explicit driver binding; no chain loader, solver or hidden observer setup.
        transition.calendar.model = model
        require(os.environ.get('NUMBA_DISABLE_JIT') == '0', 'Compiled probe environment required')
        records = []
        for case in contract['cases']:
            remaining(); phase.update(phase='loading_checkpoint', case_id=case['case_id'])
            with gzip.open(case['checkpoint'], 'rb') as stream:
                checkpoint = pickle.load(stream)
            P, evaluation = checkpoint['parameters'], checkpoint['evaluation']
            require(P.use_numba_scatter is True, 'Checkpoint must already request compiled scatter; no parameter mutation')
            first_bytes = None
            for repetition in range(1, case['passes'] + 1):
                remaining(); phase.update(phase='observing', repetition=repetition)
                write(out / 'heartbeat.json', {**phase, 'elapsed_seconds': time.monotonic() - started})
                result = observe_recent_parent_flow(evaluation, P, diagnostic_enabled=True,
                    snapshot=SNAPSHOT, age_projection=AGE_PROJECTION, diagnostic_allow_residence_proxy=True,
                    input_provenance=dict(case_id=case['case_id'], checkpoint_sha256=case['checkpoint_sha256'],
                                          parent_source_commit='7e872053', observer_source_commit=contract['source_commit']))
                encoded = canonical(result)
                if first_bytes is not None:
                    require(encoded == first_bytes, 'Observer passes do not reproduce exactly')
                first_bytes = encoded
                name = f"{case['case_id']}_pass_{repetition:02d}.json"
                write(out / name, result)
                record = dict(case_id=case['case_id'], repetition=repetition, file=name,
                              sha256=digest(out / name), canonical_result_sha256=hashlib.sha256(encoded).hexdigest(),
                              diagnostic_moment=result['moment'], diagnostic_value=result['model_value'])
                records.append(record); phase['completed_passes'] += 1
                write(out / 'latest_completed.json', record)
            del checkpoint, P, evaluation
        summary = dict(status='passed_read_only_observation', calibrated_smm=False,
            household_or_equilibrium_solves=0, compiled_startup_tests=13,
            compiled_transport_override_tested=True, source_files_verified=len(sources),
            source_fingerprint=contract['source_fingerprint'], contract_sha256=args.contract_sha256,
            empirical_reference=dict(id=contract['empirical_target_id'], value=target['estimate'],
                                     exact_model_row_filled=False, weight=None, loss_contribution=None),
            observations=records, exact_reproduction_passed=all(c['passes'] == 2 for c in contract['cases']),
            elapsed_seconds=time.monotonic() - started)
        write(out / 'summary.json', summary)
        receipt = dict(status='passed', summary_sha256=digest(out / 'summary.json'),
                       contract_sha256=args.contract_sha256, startup_log_sha256=digest(out / 'compiled_startup_tests.log'))
        write(out / 'receipt.json', receipt)
        print(json.dumps({**receipt, 'receipt_sha256': digest(out / 'receipt.json')}), flush=True)
        phase['phase'] = 'complete'
        write(out / 'heartbeat.json', {**phase, 'elapsed_seconds': time.monotonic() - started})
    except BaseException as exc:
        write(out / 'failure.json', dict(status='failed', error=repr(exc), **phase,
                                         elapsed_seconds=time.monotonic() - started))
        raise
    finally:
        stopped.set(); timer.cancel()


if __name__ == '__main__':
    main()
