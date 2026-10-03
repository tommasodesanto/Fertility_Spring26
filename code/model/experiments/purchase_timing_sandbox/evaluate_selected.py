"""Evaluate the pinned transaction-timing experiment at a saved soft candidate.

This is one fixed-coordinate normalized equilibrium, not a calibration search.
Run in a fresh interpreter through the authenticated local source overlay.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import difflib
import json
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
PACKETS = ROOT / 'output/model/fixed_reference_economics_20260928'
TIMING = PACKETS / 'purchase_timing_sandbox_v1'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selection', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise RuntimeError('Refusing existing output directory')
    args.out.mkdir(parents=True)
    started = time.time()
    deadline = started + 1200

    def alarm(*_):
        raise TimeoutError('20-minute fixed-coordinate evaluation limit')

    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, 1200)
    spec = importlib.util.spec_from_file_location('pinned_timing_driver', TIMING / 'run.py')
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    v2 = driver.v2
    try:
        _, manifest = driver.verify()
        saved = json.loads(args.selection.read_text())
        source = ROOT / saved['source']
        assert hashlib.sha256(source.read_bytes()).hexdigest() == saved['source_sha256']
        selected = saved['selected']
        assert json.loads(source.read_text())[saved['source_key']]['best'] == selected
        assert selected['weight_fingerprint'] == manifest['weight_fingerprint']
        assert v2.native.target_identity(selected['target_fit']) == v2.CONFIG['base_target_contract']
        point = selected['parameters']
        lane = 'floor_s0'
        _, bounds, _ = v2.inputs.seed_and_bounds(lane)
        bounds = {k: tuple(v) for k, v in bounds.items()}
        assert bounds['h_P'] == (.1, 2.3)
        bounds.update(h_P=(.1, 2.6), psi_child=tuple(v2.CONFIG['psi_bounds']))
        assert set(point) == set(v2.inputs.parameters(lane)) | {'psi_child'}
        assert len(point) == 10 and all(bounds[k][0] <= point[k] <= bounds[k][1] for k in point)
        v2.inputs.LANES[lane].update(seed=dict(point), bounds=bounds, free_coordinates=list(point))
        P, grid = v2.inputs.proposal(lane)
        P, _ = v2.inputs.entry(P, grid, 'nonnegative_mean')
        assert P.native_purchase_income and P.native_due_stayer_credit and not P.joint_nested_choice
        assert P.N_target == 1. and P.R_gross > 1.
        change = ('Experimental transaction timing: b_next=R*b+net_sale-purchase+income-consumption-costs; '
                  'seller solvency tested at R*b+net_sale>=0. Ending debt floors, owner stayers, '
                  'earnings, entrants, utility, ten coordinates, targets and weights unchanged. '
                  'Price clears renewal and H0 is derived at N0=1, as in the soft reference.')
        v2.write(args.out / 'start.json', dict(status='running', started_epoch=started,
            deadline_epoch=deadline, lifecycle_cap=32, fixed_parameters=point,
            selection_sha256=driver.sha(args.selection), source_manifest_sha256=driver.sha(TIMING / 'manifest.json'),
            experimental_change=change, target_fingerprint=manifest['target_fingerprint'],
            weight_fingerprint=manifest['weight_fingerprint']))
        Q = v2.native.utility_checks(P, grid, lane, args.out)
        # The frozen calendar observer has its own solver reference. It must
        # use the experiment's transaction map when reconstructing households.
        # Its independent transaction audit must test the new dated budget.
        install_metadata = v2.native.install_observer_metadata

        def timing_metadata(ge, arm, source_bounds):
            install_metadata(ge, arm, source_bounds)
            observe = ge.observe_price

            def timing_observe(ctx, live, label, *, final=False):
                rt = ctx['prepared'].rt
                if not getattr(rt['accounting'], '_transaction_timing_installed', False):
                    cal = rt['primitive'].pf.calendar
                    old_map = cal.model.build_forward_tenure_transition_maps
                    new_map = driver.sandbox_household.build_forward_tenure_transition_maps
                    old_audit = rt['accounting']._inherited_audit()
                    before = inspect.getsource(old_audit)
                    replacements = {
                        'x = grid if old == new else grid + sale[old] - costs[new]':
                        'x = grid if old == new else grid + (sale[old] - costs[new]) / float(P.R_gross)',
                        'invalid = x + y / float(P.R_gross) < floor - 1e-10':
                        'invalid = float(P.R_gross) * x + y < floor - 1e-10',
                    }
                    after = before
                    for old, new in replacements.items():
                        assert after.count(old) == 1, old
                        after = after.replace(old, new)
                    namespace = dict(old_audit.__globals__)
                    path = args.out / 'timing_purchase_audit.py'
                    path.write_text(after)
                    (args.out / 'timing_purchase_audit.diff').write_text(''.join(
                        difflib.unified_diff(before.splitlines(True), after.splitlines(True),
                            fromfile='original', tofile='transaction_timing')))
                    exec(compile(after, str(path), 'exec'), namespace)
                    rt['accounting']._INHERITED = namespace['audit_purchase_accounting']
                    cal.model.build_forward_tenure_transition_maps = new_map
                    rt['model'].build_forward_tenure_transition_maps = new_map
                    rt['accounting']._transaction_timing_installed = True
                    v2.write(args.out / 'timing_observer_receipt.json', dict(
                        original_map_file=inspect.getsourcefile(old_map),
                        experiment_map_file=inspect.getsourcefile(new_map),
                        original_audit_sha256=hashlib.sha256(before.encode()).hexdigest(),
                        timing_audit_sha256=hashlib.sha256(after.encode()).hexdigest(),
                        acceptance_tolerances_unchanged=True,
                        audit_budget='R*b + net_sale - purchase + income - consumption - costs',
                        unrelated_audit_statements_unchanged=True))
                return observe(ctx, live, label, final=final)

            ge.observe_price = timing_observe

        v2.native.install_observer_metadata = timing_metadata
        evaluate = v2.normalized_objective.make_evaluator(args.out, lane, Q, grid, deadline,
            float(selected['price']), native_runner=v2.native)
        review = args.out / 'normalization_source_review/receipt.json'
        receipt = json.loads(review.read_text())
        receipt.update(household_solver_unchanged=False,
            isolated_purchase_timing_source_manifest_sha256=driver.sha(TIMING / 'manifest.json'))
        v2.write(review, receipt)
        result = evaluate('alternative_timing', point, deadline)
        assert result['status'] == 'passed', result
        report = Path(result['report'])
        fits = v2.native.readtable(report / 'target_fit.csv')
        parameters = v2.native.readtable(report / 'parameters.csv')
        assert len(fits) == 14 and len(parameters) == 31
        assert v2.native.target_identity(fits) == v2.CONFIG['base_target_contract']
        repeat = v2.native.compare_repeated(report, report.parent / 'selected_repeat_final')
        result.update(status='fixed_coordinate_timing_experiment_passed', target_fit=fits,
            parameters=parameters, repeat=repeat, loss=sum(float(r['loss_contribution'] or 0) for r in fits),
            fixed_parameters=point, experimental_change=change, experimental_not_adopted=True,
            elapsed_seconds=time.time()-started)
        v2.write(args.out / 'completed.json', result)
    except BaseException as exc:
        v2.write(args.out / 'failure.json', dict(type=type(exc).__name__, message=str(exc),
            elapsed_seconds=time.time()-started, no_auto_retry=True))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
