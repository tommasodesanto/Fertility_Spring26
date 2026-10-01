"""Prepare, or execute once, a fixed-old-parameters floor full GE.

Default mode authenticates sources, binds inputs, and constructs the reviewed
native evaluator with zero lifecycle calls. --run-once calls it exactly once.
"""
from __future__ import annotations
import argparse, builtins, csv, hashlib, importlib.machinery, io, json, os, runpy, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
OUT = Path(__file__).resolve().parent
ROUND2 = ROOT / 'output/model/fixed_reference_economics_20260928/utility_floor_round2_v1'
REFERENCE = ROOT / 'output/model/fixed_reference_economics_20260928/utility_calibration_round1_v1/deployment/attempt2/comparison_requested'
CANDIDATE = ROOT / 'output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/chain_11/results/0026_nm/phase_b_ge/selected_root'

OLD_PARAMETERS = REFERENCE / 'with_A_parameters.csv'
OLD_FIT = REFERENCE / 'with_A_target_fit.csv'
HP = 1.9105885246262313
PRICE_START = 0.678660685

# Activate the already reviewed, read-only local source overlay used by the
# running utility-floor chains before loading the native runner.
BOOTSTRAP = ROUND2 / 'local_run/bootstrap.py'
_overlay_source = BOOTSTRAP.read_text()
_overlay_tail = "runpy.run_path(str(HERE/'runner_local.py'),run_name='__main__')"
assert _overlay_source.count(_overlay_tail) == 1
exec(compile(_overlay_source.replace(_overlay_tail, 'pass'), str(BOOTSTRAP), 'exec'),
     {'__file__': str(BOOTSTRAP), '__name__': '_old_parameters_floor_overlay'})
import inputs
import runner as native


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_parameters():
    with OLD_PARAMETERS.open(newline='') as f:
        rows = list(csv.DictReader(f))
    return {r['parameter']: float(r['estimate']) for r in rows}, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-once', action='store_true', help='Call exactly one full native GE evaluator')
    ap.add_argument('--go-reviewed', action='store_true', help='Required explicit lead GO after reviewing preflight artifacts')
    ap.add_argument('--deadline-epoch', type=float)
    ap.add_argument('--out', type=Path, default=OUT / 'results')
    args = ap.parse_args()
    run_started = time.time()
    deadline = min(run_started + 1200, args.deadline_epoch or float('inf'))
    out = args.out.resolve()
    inputs.require(not out.exists(), 'Refusing existing results directory')
    out.mkdir(parents=True)
    inputs.require(not args.run_once or args.go_reviewed, 'Numerical execution requires explicit lead GO after review')
    inputs.require(args.run_once or not args.go_reviewed, '--go-reviewed only applies to --run-once')
    native.verify_sources()
    old, rows = read_parameters()
    inputs.require(old['h_P'] == 0.0, 'With-A reference must retain its original h_P=0')
    # Preserve the reference's original h_P=0; the physical floor is separately pinned.
    lane = 'floor_s0'
    seed, bounds, _ = inputs.seed_and_bounds(lane)
    seed.update({k: old[k] for k in inputs.parameters(lane) if k != 'h_P'})
    seed['h_P'] = HP
    inputs.check_point(seed, bounds)
    inputs.LANES[lane]['seed'] = seed
    native.PLAN['reference_parameter_table'] = rows
    P, grid = inputs.proposal(lane)
    P, entry = inputs.entry(P, grid, 'nonnegative_mean')
    # Carry the old point's fixed common values from its authenticated 31-row table.
    import numpy as np
    np.testing.assert_array_equal(np.asarray(P.H0), np.asarray(old['H0']))
    np.testing.assert_array_equal(np.asarray(P.psi_child), np.asarray(old['psi_child']))
    # Confirm the requested, strictly scoped economic changes through the existing
    # reviewed binding: floor arm sets A(m)=0, child-share loading=0, and h_P=HP.
    Q = native.utility_checks(P, grid, lane, out)
    np.testing.assert_array_equal(np.asarray(Q.H0), np.asarray(old['H0']))
    np.testing.assert_array_equal(np.asarray(Q.psi_child), np.asarray(old['psi_child']))
    inputs.require(Q.child_room_floor and Q.delta_alpha_jump == 0. and Q.delta_alpha == 0., 'Floor utility mapping differs')
    inputs.require(Q.hbar_first_child_jump == HP, 'Pinned floor value differs')
    target_sha = inputs.canonical(native.PLAN['target_contract'])
    inputs.require(entry['dimensions'] == [120, 9] and entry['negative_wealth_share'] == 0., 'Entry/grid/no-borrowing contract differs')
    # Construct the native evaluator only: source auth + mappings, zero lifecycle calls.
    evaluator = native.native_evaluator(out, lane, Q, grid, deadline, PRICE_START)
    inputs.require(callable(evaluator), 'Native evaluator initialization failed')
    closure = dict(zip(evaluator.__code__.co_freevars, (c.cell_contents for c in evaluator.__closure__)))
    inputs.require('base' in closure and 'ge' in closure, 'Native evaluator cap hooks not exposed')
    base, ge = closure['base'], closure['ge']
    OriginalBudget = base.ArmBudget
    class Budget12(OriginalBudget):
        @property
        def max_lifecycle(self): return self._old_floor_max_lifecycle
        @max_lifecycle.setter
        def max_lifecycle(self, value): self._old_floor_max_lifecycle = min(int(value), 12)
    base.ArmBudget = Budget12
    original_phase_b = ge.run_phase_b
    def run_phase_b_12(context, *a, **kw):
        context['phase_b_max_new_lifecycle'] = 12
        return original_phase_b(context, *a, **kw)
    ge.run_phase_b = run_phase_b_12
    pins = {
        'old_reference_parameters_csv': {'path': str(OLD_PARAMETERS.relative_to(ROOT)), 'sha256': sha(OLD_PARAMETERS)},
        'old_reference_target_fit_csv': {'path': str(OLD_FIT.relative_to(ROOT)), 'sha256': sha(OLD_FIT)},
        'candidate_floor_parameters_csv': {'path': str(CANDIDATE.relative_to(ROOT) / 'parameters.csv'), 'sha256': sha(CANDIDATE / 'parameters.csv')},
        'candidate_floor_target_fit_csv': {'path': str(CANDIDATE.relative_to(ROOT) / 'target_fit.csv'), 'sha256': sha(CANDIDATE / 'target_fit.csv')},
        'native_floor_source_pins': str((ROUND2 / 'source_pins.json').relative_to(ROOT)),
        'native_floor_source_pins_sha256': sha(ROUND2 / 'source_pins.json'),
        'coarse_bundle_sha256': inputs.sha(ROOT / native.PLAN['bundle_sources']['coarse']['folder'] / 'bundle.json'),
        'coarse_arrays_sha256': inputs.sha(ROOT / native.PLAN['bundle_sources']['coarse']['folder'] / 'arrays.npz'),
        'old_reference_parameter_values': old,
        'effective_expected_parameters': native.expected_parameters(seed, (120, 9), 'floor'),
        'economic_differences': ['A(m) disabled', f'delta_alpha_jump {old["delta_alpha_jump"]:.17g} -> 0', f'h_P 0 -> {HP:.17g} physical rooms'],
        'h_P_pinned_from_chain_11_0026_nm': HP,
        'selected_price_start': PRICE_START,
        'target_contract_sha256': target_sha,
        'entry_contract': entry,
    }
    inputs.require(old['psi_child'] == 0.1355551166583114 and old['delta_alpha_jump'] == 0.1303383207216736, 'Historical reference value mismatch')
    inputs.require(old['child_benefit_curvature'] == 0.07222157735996547, 'Historical curvature mismatch')
    # The primitive benefit coefficient is derived from the fixed benefit and curvature.
    derived = (1. - old['child_benefit_curvature']) * old['psi_child']
    inputs.require(abs(derived - old['child_benefit_CRRA_coefficient']) < 1e-15, 'Derived benefit coefficient mismatch')
    native.write(out / 'source_input_pins.json', pins)
    native.write(out / 'initialized.json', {
        'status': 'zero_solve_native_initializer_passed', 'lifecycle_solves': 0,
        'native_evaluator_constructed': True, 'evaluator_called': False,
        'enforced_lifecycle_cap': 12,
        'cap_binding': 'Budget12.max_lifecycle property clamps runner assignment to 12; run_phase_b_12 binds phase_b_max_new_lifecycle=12',
        'old_reference': 'with-A nonnegative_mean_120x9, 037_r2_gn_1.0, loss 18.128926; provisional historical snapshot',
        'old_parameter_csv_sha256': sha(OLD_PARAMETERS), 'h_P': HP,
        'shared_parameters': 'exact old 31-row values; psi_child and derived benefit coefficient retained',
        'economic_differences': [
            'A(m) disabled',
            f'delta_alpha_jump changed {old["delta_alpha_jump"]:.17g} to 0 (child-share loading disabled)',
            f'h_P changed 0 to {HP:.17g} physical rooms (pinned from chain_11/results/0026_nm)',
        ],
        'same_grid_entry_no_borrowing_rate_targets_closure': True,
        'dimensions': [120, 9], 'annual_rate': 0.02,
        'target_contract_sha256': target_sha,
        'selected_price_start': PRICE_START,
        'time_budget_seconds': 1200, 'lifecycle_call_cap': 12,
        'cpu_threads': 2, 'rss_limit_gib': 12,
        'run_command_after_lead_GO': 'python prepare.py --run-once --go-reviewed --out <fresh-results-directory>',
        'provenance': pins,
    })
    if args.run_once:
        # launch.py sets thread environment variables before imports; preserve the
        # native compatibility overlay's existing Numba limit without resetting it.
        result = evaluator('000_fixed_old_parameters', seed, deadline)
        inputs.require(result.get('status') == 'passed', 'One full native GE did not pass: ' + str(result.get('status')))
        inputs.require(int(result.get('lifecycle_solves', 0)) <= 12, 'Native evaluator exceeded the 12-call cap')
        report = Path(result['report'])
        rows = native.readtable(report / 'target_fit.csv')
        params = native.readtable(report / 'parameters.csv')
        inputs.require(len(rows) == 14 and len(params) == 31, 'Full GE report tables incomplete')
        plot_hashes = native.report_hashes(report)
        inputs.require(len(plot_hashes) == 20, '17 plot + 3 report-file hashes incomplete')
        native.write(out / 'completed.json', dict(status='one_full_ge_passed', result=result, full_target_rows=14, full_parameter_rows=31, standard_plot_count=17, lifecycle_cap=12, elapsed_seconds=time.time()-run_started, target_fit=rows, parameters=params, report_hashes=plot_hashes))

if __name__ == '__main__':
    main()
