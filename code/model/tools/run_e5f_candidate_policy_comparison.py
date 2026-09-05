#!/usr/bin/env python3
"""Observe unchanged policy drivers in a pinned, verified source snapshot.

The observer copies solved states for the established diagnostic packet. It
does not replace the household solver, transition, target, or policy equations.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import os
import pickle
import sys
import threading
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code/model'), str(ROOT / 'code/model/tools')]
CASES = ('baseline', 'supply-plus-20', 'dependent-child-ltv95', 'combined',
         'property-tax-2pct-no-rebate')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def verify_plan(path, expected):
    if sha(path) != expected:
        raise RuntimeError('Plan hash changed')
    plan = read(path)
    if plan['schema'] != 'e5f_candidate_policy_comparison_v1':
        raise RuntimeError('Unknown comparison plan')
    if time.time() > plan['launch_deadline_epoch']:
        raise RuntimeError('Declared launch deadline exceeded')
    for name, digest in plan['source_hashes'].items():
        if sha(ROOT / name) != digest:
            raise RuntimeError(f'Source hash changed: {name}')
    for name, digest in plan['input_hashes'].items():
        if sha(name) != digest:
            raise RuntimeError(f'Input hash changed: {name}')
    import run_e5f_transition_calibration as calibration
    from intergen_eqscale_seq_optimized import solver
    if calibration.code_fingerprint_contract(solver)['bundle_sha256'] != plan['code_bundle_sha256']:
        raise RuntimeError('Scientific source bundle changed')
    return plan


def validate_smoke(root, plan):
    import collect_e5f_post2023_policy_mechanisms as collector
    summary = read(root / 'smoke/report/summary.json')
    if summary['status'] != 'complete_post2023_policy_mechanism_comparison':
        raise RuntimeError('Five-case smoke is incomplete')
    if summary['cases'] != list(CASES) or summary['post_2023_periods'] != 1:
        raise RuntimeError('Smoke case set or horizon differs')
    for case in CASES:
        s, _ = collector.validate_case(root / 'smoke' / case, case)
        if s['input_hashes']['selected_report_sha256'] != plan['selected_report_sha256']:
            raise RuntimeError('Smoke selected calibration differs')
        receipt = read(root / 'smoke_diagnostics' / case / 'observer_receipt.json')
        if receipt['status'] != 'complete' or receipt['completed_dates'] != 2:
            raise RuntimeError('Smoke observer loop is incomplete')
    return sha(root / 'smoke/report/summary.json')


def call_main(module, arguments):
    previous = sys.argv
    try:
        sys.argv = [str(module.__file__), *arguments]
        module.main()
    finally:
        sys.argv = previous


def run(args):
    plan = verify_plan(args.plan, args.plan_sha256)
    root = Path(plan['output_root'])
    if args.stage != 'smoke':
        smoke_sha = validate_smoke(root, plan)
    else:
        smoke_sha = None
    import run_e5f_post2023_policy_mechanisms as policy
    import run_e5f_independent_numerical_audit as audit
    import numpy as np
    baseline = policy.baseline
    diagnostic_root = root / f'{args.stage}_diagnostics' / args.case
    if diagnostic_root.exists():
        raise FileExistsError(diagnostic_root)
    diagnostic_root.mkdir(parents=True)
    started = time.time()
    progress = {'phase': 'source_and_input_gates_passed', 'completed_dates': 0}
    stopped = threading.Event()
    def heartbeat():
        while not stopped.is_set():
            write(diagnostic_root / 'heartbeat.json', {
                **progress, 'epoch': time.time(), 'elapsed_seconds': time.time()-started,
                'job_id': os.environ.get('SLURM_JOB_ID'), 'case': args.case})
            stopped.wait(60)
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    receipts = []
    active = False
    original_advance = baseline.advance_from_evaluation
    original_path = policy.run_policy_path

    def observed_advance(**kwargs):
        result = original_advance(**kwargs)
        if not active:
            return result
        year = 2007 + 4 * int(kwargs['period_from_2007'])
        out = diagnostic_root / str(year)
        out.mkdir()
        progress['phase'] = f'diagnostics_{year}'
        packet = copy.deepcopy(dict(parameters=kwargs['P'], b_grid=kwargs['b_grid'],
            evaluation=kwargs['evaluation'], shared=kwargs['shared'],
            supply_rule=kwargs['supply_rule'], state=kwargs['state'],
            calendar_year=year, policy_case=args.case))
        checkpoint = out / 'dated_state.pkl.gz'
        with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
            pickle.dump(packet, stream, protocol=5)
        audit.standard_diagnostics(packet, out, validate_production_young=False)
        budget = audit.budget_audit(packet, out)
        arrays = audit.policy_array_audit(packet, out)
        for bounds in arrays['probabilities'].values():
            if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
                raise RuntimeError('Invalid probability in solved policy')
        # Read-only copying and graph generation must leave the live state intact.
        e, original_e = packet['evaluation'], kwargs['evaluation']
        for field in ('g_pre', 'g_current', 'g_post_fertility'):
            if not np.array_equal(getattr(e, field), getattr(original_e, field)):
                raise RuntimeError(f'Diagnostic mutation of {field}')
        for field in ('V', 'bp_pol', 'hR_pol', 'tenure_probs', 'fert_probs', 'fert2_probs'):
            if not np.array_equal(getattr(e.policy, field), getattr(original_e.policy, field), equal_nan=True):
                raise RuntimeError(f'Diagnostic mutation of policy {field}')
        graphs = sorted((out / 'standard_diagnostics').glob('*.png'))
        if len(graphs) != 17:
            raise RuntimeError('Standard seventeen-graph set is incomplete')
        record = dict(calendar_year=year, status='complete', budget=budget,
            policy_arrays=arrays, grid_fallback=bool(kwargs['grid_fallback']),
            artifact_hashes={str(p.relative_to(out)): sha(p) for p in [checkpoint, *graphs]})
        write(out / 'receipt.json', record)
        receipts.append(record)
        progress.update(phase=f'completed_{year}', completed_dates=len(receipts))
        write(diagnostic_root / 'latest_completed_case.json', record)
        write(diagnostic_root / 'best_so_far.json', {
            'scope': 'fixed selected candidate; no search or promotion',
            'calibration_loss': plan['candidate_loss'], 'completed_dates': len(receipts)})
        return result

    def observed_path(*positional, **keyword):
        nonlocal active
        active = True
        try:
            return original_path(*positional, **keyword)
        finally:
            active = False

    try:
        common = plan['common_arguments']
        if args.stage == 'rebate':
            import run_e5f_post2023_rebated_property_tax_smoke as impact
            import run_e5f_rebated_tax_shapley_diagnosis as shapley
            progress['phase'] = 'rebated_impact_three_equilibria'
            call_main(impact, [*common, '--output-dir', str(root / 'rebated_impact')])
            impact_summary = read(root / 'rebated_impact/summary.json')
            if impact_summary['status'] != 'complete_exact_rebated_property_tax_smoke':
                raise RuntimeError('Exact rebated-tax impact failed')
            progress['phase'] = 'shapley_eight_cells'
            call_main(shapley, [*common, '--rebated-path-csv', str(root / 'rebated_impact/impact_results.csv'),
                '--output-dir', str(root / 'shapley')])
            if read(root / 'shapley/summary.json')['status'] != 'complete':
                raise RuntimeError('Shapley decomposition failed')
        else:
            baseline.advance_from_evaluation = observed_advance
            policy.run_policy_path = observed_path
            periods = 1 if args.stage == 'smoke' else 10
            call_main(policy, [*common, '--policy-case', args.case, '--post-2023-periods', str(periods),
                '--output-dir', str(root / args.stage / args.case)])
            if len(receipts) != periods + 1:
                raise RuntimeError('Observer date count differs from policy loop')
            if args.stage == 'full':
                # Same first two dates must exactly nest the short-loop smoke.
                before = baseline.read_csv(root / 'smoke' / args.case / 'policy_path.csv')
                after = baseline.read_csv(root / 'full' / args.case / 'policy_path.csv')[:2]
                for left, right in zip(before, after):
                    for key in left:
                        if not any(x in key for x in ('elapsed', 'seconds')) and left[key] != right[key]:
                            raise RuntimeError(f'Smoke/full prefix differs: {key}')
        progress['phase'] = 'complete'
        write(diagnostic_root / 'observer_receipt.json', dict(status='complete',
            completed_dates=len(receipts), stage=args.stage, case=args.case,
            plan_sha256=args.plan_sha256, smoke_summary_sha256=smoke_sha,
            observer_sha256=sha(__file__), elapsed_seconds=time.time()-started,
            production_promoted=False, source_bundle=plan['code_bundle_sha256']))
    except BaseException as error:
        progress['phase'] = 'failed'
        write(diagnostic_root / 'failure.json', dict(error=str(error), traceback=traceback.format_exc()))
        raise
    finally:
        baseline.advance_from_evaluation = original_advance
        policy.run_policy_path = original_path
        stopped.set()
        thread.join(timeout=2)
        write(diagnostic_root / 'heartbeat.json', {**progress, 'elapsed_seconds': time.time()-started})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--plan-sha256', required=True)
    parser.add_argument('--stage', choices=('smoke', 'full', 'rebate'), required=True)
    parser.add_argument('--case', choices=CASES, default='baseline')
    run(parser.parse_args())
