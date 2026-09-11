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


def parent_gates(case, bundle, raw, collection, *, verify_remote_files):
    """Validate a native successful fresh joint case and its original receipts."""
    import csv
    import importlib.util
    require(raw['case_count'] == 23 and collection['case_count_verified'] == 23,
            'Require final 23-success collection, including the reused smoke')
    require(raw['status'] == collection['status'] == 'partial',
            'Preserve the original failed-case status; do not relabel the 24-case round')
    require(collection['array_job_id'] == '17362324' and collection['source_pins_verified'] == 634,
            'Wrong original joint collection')
    failed = [s for s in collection['states'] if s['status'] != 'verified']
    require(len(failed) == 1 and failed[0]['index'] == 1 and failed[0]['status'] == 'failed',
            'Expected exactly the original failed index 1')
    verify(bounded_path(bundle, 'inputs/raw_case_summary.json'), collection['raw_case_summary_sha256'])
    verify(bounded_path(bundle, 'inputs/original_collector.py'), collection['collector_sha256'])
    verify(bounded_path(bundle, 'inputs/original_cases.json'), collection['cases_file_sha256'])
    plan_path = bounded_path(bundle, 'inputs/original_run_plan.json')
    verify(plan_path, collection['plan_sha256'])
    plan = json.loads(plan_path.read_text())
    meta = json.loads(bounded_path(bundle, 'inputs/original_preparation_metadata.json').read_text())
    verify(bounded_path(bundle, 'inputs/panel_validator.py'), plan['artifact_sha256']['inputs/panel_validator.py'])
    require(collection['smoke_gate_receipt']['exact_early_equality'] is True,
            'Original prerequisite smoke replay did not pass')
    entries = [r for r in raw['cases'] if r['case_id'] == case['case_id']]
    require(len(entries) == 1, 'Missing or duplicate native joint record')
    record = entries[0]
    require(record['status'] == 'verified' and record['numeric_gates_verified'] is True,
            'Original numerical gates did not pass')
    require(record['smoke_reuse'] is False and record['index'] not in (1, 21),
            'Only successful fresh cases are observed in this batch')
    require(case['original_index'] == record['index'] and record['source_commit'] == '7e872053',
            'Original joint source or index changed')
    states = [s for s in collection['states'] if s['case_id'] == case['case_id']]
    require(len(states) == 1 and states[0]['status'] == 'verified', 'Case absent from final verified states')
    original_cases = json.loads(bounded_path(bundle, 'inputs/original_cases.json').read_text())
    proposal = next(x for x in original_cases if x['case_id'] == case['case_id'])
    parent_path = bounded_path(bundle, case['parent_contract'])
    parent = json.loads(parent_path.read_text())
    verify(parent_path, proposal['contract_sha256'])
    require(record['actual_executed_contract_sha256'] == record['proposal_contract_sha256']
            == proposal['contract_sha256'], 'Fresh case must retain its own executed contract')
    require(record['input_fingerprint'] == parent['run_input_fingerprint']
            == meta['run_input_fingerprint'] == plan['input_fingerprint'], 'Joint input fingerprints differ')
    require(record['run_plan_sha256'] == collection['plan_sha256'], 'Run-plan fingerprint differs')
    require(parent['case_id'] == case['case_id'] and parent['repetitions'] == 1
            and parent['maximum_GE_solves'] == 8, 'Original fresh-case contract changed')
    require(parent['source_sha256'] == meta['source_pins'] and len(parent['source_sha256']) == 634,
            'Original economic source manifest differs')
    dest = bounded_path(bundle, case['collected_directory'])
    require(json.loads((dest/'collection_receipt.json').read_text()) == record,
            'Case receipt differs from original raw collection')
    for name, pin in record['file_sha256'].items():
        verify(bounded_path(dest, name), pin)
    actual = json.loads((dest/'contract.json').read_text())
    require(actual == dict(parent, contract_sha256=proposal['contract_sha256'], case='new_balanced'),
            'Original output contract mismatch')
    seed = json.loads((dest/'seed_mapping.json').read_text())
    require(seed['source_files_verified'] == 634 and seed['initial_psi'] == parent['initial_psi'],
            'Original seed/source mapping changed')
    top = json.loads((dest/'summary.json').read_text())
    rep = dest/'repetition_01'
    final = json.loads((rep/'summary.json').read_text())
    early = json.loads((rep/'early_measurement.json').read_text())
    ges = json.loads((rep/'stationary_solves.json').read_text())
    require(top == record['summary'] and final == record['repetition_summary'] and top['final'] == final,
            'Original summary objects differ')
    spec = importlib.util.spec_from_file_location('original_panel_validator', bounded_path(bundle, 'inputs/panel_validator.py'))
    validator = importlib.util.module_from_spec(spec); spec.loader.exec_module(validator)
    with (rep/'parameters.csv').open() as stream:
        parameters = validator.validate_summary(top, parent, early, list(csv.DictReader(stream)), ges)
    require(parameters == record['parameters'] and len(ges) == record['stationary_solves'],
            'Original parameter or solve accounting differs')
    market = validator.validate_market(json.loads((rep/'market_quantity_units.json').read_text()))
    require(market == record['market_residual'], 'Original market receipt differs')
    checkpoint = record['checkpoint']
    require(case['checkpoint'] == checkpoint['path'] == record['remote_output'] + '/repetition_01/initial_state.pkl.gz',
            'Original checkpoint path changed')
    require(case['checkpoint_sha256'] == checkpoint['sha256'] == final['checkpoint_sha256']
            and case['checkpoint_bytes'] == checkpoint['bytes'], 'Original checkpoint fingerprint changed')
    graphs = record['original_graphs']
    require(len(graphs) == 17 and sorted(x['filename'] for x in graphs)
            == sorted(meta['expected_graph_filenames']), 'Original 17-graph receipt incomplete')
    require(len({x['path'] for x in graphs}) == 17, 'Duplicate original graph paths')
    if verify_remote_files:
        require(Path(case['checkpoint']).stat().st_size == case['checkpoint_bytes'], 'Checkpoint byte count differs')
        for graph in graphs:
            expected = record['remote_output'] + '/repetition_01/standard_diagnostics/' + graph['filename']
            require(graph['path'] == expected, 'Original graph path changed')
            path = Path(graph['path'])
            require(path.stat().st_size == graph['bytes'], 'Original graph byte count differs')
            verify(path, graph['sha256'])
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
    require(contract['schema'] == 'e5f_recent_parent_joint_cases_probe_v1', 'Wrong schema')
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
    raw = json.loads(bounded_path(bundle, contract['raw_case_summary']).read_text())
    panel = json.loads(bounded_path(bundle, contract['joint_collection_receipt']).read_text())
    empirical = json.loads(bounded_path(bundle, contract['empirical_target_file']).read_text())
    target = next(x for x in empirical['records'] if x['id'] == contract['empirical_target_id'])
    require(target['estimate'] == contract['empirical_target_value'] == 0.16289550916123285, 'Empirical target changed')
    ids = [c['case_id'] for c in contract['cases']]
    require(len(ids) == len(set(ids)) and len(ids) > 0, 'Duplicate or empty cases')
    require(len(ids) == 22 and {c['original_index'] for c in contract['cases']} == set(range(24)) - {1, 21},
            'Require exactly 22 successful fresh cases; failed and reused indices excluded')
    for case in contract['cases']:
        require(type(case['passes']) is int and case['passes'] in (1, 2), 'Invalid pass count')
        parent = parent_gates(case, bundle, raw, panel, verify_remote_files=not args.verify_only)
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
