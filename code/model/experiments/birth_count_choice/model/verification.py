"""One fresh-process execution or strict reference/production artifact comparison."""
from __future__ import annotations
import argparse
import csv
import functools
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import socket
import sys
import time
from pathlib import Path
from typing import Any
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PACKETS = ROOT / 'output/model/fixed_reference_economics_20260928'
SELECTION = PACKETS / 'soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
NATIVE_RECEIPT = SELECTION.parent / 'native_postcheck/completed.json'
REQUIRED = {'V', 'bp_pol', 'c_pol', 'hR_pol', 'g', 'tenure_probs', 'fert_probs', 'fert2_probs', 'b_grid', 'shared.c_bar', 'shared.h_bar', 'shared.alpha_flat', 'shared.escale_flat'}
POLICIES = {'V', 'bp_pol', 'c_pol', 'hR_pol', 'tenure_probs', 'fert_probs', 'fert2_probs'}


def _json(value: Any):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, Path): return str(value)
    raise TypeError(type(value).__name__)


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _selected():
    saved = json.loads(SELECTION.read_text())
    native = json.loads(NATIVE_RECEIPT.read_text())
    selected = dict(saved['selected'])
    if saved['arm'] != 'alternative' or saved['chain'] != 13 or saved['status'] != 'selected_numerically_verified':
        raise RuntimeError('Expected the authenticated post-interest chain13 continuation anchor')
    receipt = native['selected_postcheck']
    if native['status'] != 'full_native_postcheck_passed' or receipt['status'] != 'passed':
        raise RuntimeError('Missing passed chain13 native receipt')
    if not receipt['price'] == selected['price'] == .7760569760205563:
        raise RuntimeError('Chain13 native price drift')
    if native['native_loss'] != saved['native_loss']:
        raise RuntimeError('Chain13 native loss drift')
    selected['price'] = receipt['price']
    return selected


def _point(case):
    point = {key: float(value) for key, value in _selected()['parameters'].items()}
    external = None
    if case == 'beta': point['beta_annual'] -= .001
    elif case == 'external': external = {'sigma': 2.01}
    elif case != 'unchanged': raise ValueError(case)
    return point, external


def _install_reference_overlay():
    if any(name.startswith('production.engine') for name in sys.modules):
        raise RuntimeError('Reference process already imported the production engine')
    path = PACKETS / 'purchase_rules_overnight_v1/local_runtime/bootstrap.py'
    prefix = path.read_text().split("if '--preflight-context' in sys.argv:", 1)[0]
    if 'DIGESTS=' not in prefix or 'Frozen overlay write forbidden' not in prefix:
        raise RuntimeError('Authenticated local overlay structure drift')
    exec(compile(prefix, str(path), 'exec'), {'__file__': str(path), '__name__': 'authenticated_reference_overlay'})


def _reference(case, out, budget, max_lifecycle, preflight=False):
    point, external = _point(case)
    if external is not None:
        raise RuntimeError('Reference external edits require a separate expected-parameter observer contract; use production external propagation check')
    if max_lifecycle != 32: raise ValueError('Historical oracle contract requires max-lifecycle=32')
    _install_reference_overlay()
    calibrate = _load(ROOT / 'code/model/experiments/purchase_timing_sandbox/calibrate.py', 'timing_reference_calibrate')
    v2, timing, manifest, _ = calibrate.checked_inputs('alternative')
    selected = _selected()
    if manifest['target_fingerprint'] != json.loads(SELECTION.read_text())['target_fingerprint']:
        raise RuntimeError('Reference target contract does not match chain13')
    lane = 'floor_s0'
    _, bounds, _ = v2.inputs.seed_and_bounds(lane)
    bounds = {key: tuple(value) for key, value in bounds.items()}
    bounds.update(h_P=(.1, 2.6), psi_child=tuple(v2.CONFIG['psi_bounds']))
    v2.inputs.check_point(point, bounds)
    v2.inputs.LANES[lane].update(seed=point, bounds=bounds, free_coordinates=list(point))
    P, grid = v2.inputs.proposal(lane)
    P, _ = v2.inputs.entry(P, grid, 'nonnegative_mean')
    calibrate.install_timing_observer(v2, timing, out, P)
    Q = v2.native.utility_checks(P, grid, lane, out)
    deadline = time.time() + budget
    evaluator = v2.normalized_objective.make_evaluator(out, lane, Q, grid, deadline,
        selected['price'], native_runner=v2.native, exploratory=False)
    if preflight:
        return dict(status='preflight_passed_zero_solves', lifecycle_solves=0, point=point, price_start=selected['price'])
    result = evaluator('reference', point, deadline)
    if result['status'] != 'passed': raise RuntimeError(f'Reference failed: {result}')
    if any(name.startswith('production.engine') for name in sys.modules):
        raise RuntimeError('Reference execution imported production engine')
    return result


def _production(case, out, budget, max_lifecycle, preflight=False):
    point, external = _point(case)
    from .inputs import load_inputs
    from .equilibrium import solve_stationary_ge
    P, grid = load_inputs(parameters=point, external_inputs=external)
    if external and float(P.sigma) != external['sigma']: raise RuntimeError('External sigma failed input propagation')
    if preflight:
        return dict(status='preflight_passed_zero_solves', lifecycle_solves=0, point=point, price_start=_selected()['price'])
    result = solve_stationary_ge(P, grid, out=out / 'production', price_start=_selected()['price'],
        budget_seconds=budget, max_lifecycle=max_lifecycle, closure='population_one')
    if result.get('status', 'passed') != 'passed': raise RuntimeError(f'Production failed: {result}')
    return dict(status='passed', report=str(result['report_directory']), lifecycle_solves=result['lifecycle_solves'])


def _rows(path):
    with Path(path).open(newline='') as stream: return list(csv.DictReader(stream))


def _inventory(report):
    # The native normalization adapter writes this final exact-price repeat.
    # Backend directory prefixes must never become comparison keys.
    archive = report.parent / 'selected_repeat/stage/solution_arrays.npz'
    if not archive.is_file(): raise RuntimeError(f'Missing final canonical array archive: {archive}')
    with np.load(archive, allow_pickle=False) as saved:
        arrays = {name: saved[name].copy() for name in saved.files}
    _validate_arrays(arrays)
    inventory = {name: dict(shape=list(value.shape), dtype=str(value.dtype),
                    sha256=hashlib.sha256(value.tobytes()).hexdigest()) for name, value in arrays.items()}
    return arrays, inventory


@functools.lru_cache(maxsize=1)
def _expected_array_names():
    archive = SELECTION.parent / "native_postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz"
    with np.load(archive, allow_pickle=False) as saved:
        names = set(saved.files)
    if not REQUIRED <= names: raise RuntimeError("Anchor lacks required full policy/shared arrays")
    return names


def _validate_arrays(arrays):
    missing = (REQUIRED | _expected_array_names()) - set(arrays)
    empty = [name for name, value in arrays.items() if value.size == 0]
    if missing or empty or not arrays: raise RuntimeError(f'Incomplete policy/shared array inventory: missing={sorted(missing)} empty={empty}')
    for name, value in arrays.items():
        if value.dtype.kind not in 'biufc': raise RuntimeError('Non-numeric canonical array: ' + name)
        if not np.isfinite(value).all(): raise RuntimeError('Nonfinite canonical array: ' + name)


def _runtime():
    import importlib.metadata
    packages = {}
    for name in ('numpy', 'scipy', 'numba', 'matplotlib'):
        try: packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: packages[name] = None
    return dict(hostname=socket.gethostname(), platform=platform.platform(), machine=platform.machine(),
                python=platform.python_version(), executable=str(Path(sys.executable).resolve()), packages=packages,
                threads={key: os.environ.get(key) for key in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')})


def _plots(report):
    result = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((report / 'standard_diagnostics').glob('*.png'))}
    if len(result) != 17: raise RuntimeError('Expected exactly 17 standard diagnostic PNGs')
    return result


def _closure(report):
    value = json.loads((report / 'closure.json').read_text())
    required = {'H0_derived', 'price', 'population_scale', 'renewal_residual', 'absolute_housing_residual',
                'actual_paygo_residual', 'native_population_step', 'occupied_renter_floor'}
    if required - set(value): raise RuntimeError('Missing required closure/residual fields: ' + str(sorted(required - set(value))))
    return value


def execute(backend, case, out, budget, max_lifecycle=32, preflight=False):
    if out.exists(): raise RuntimeError(f'Refusing existing output path: {out}')
    out.mkdir(parents=True)
    result = (_reference if backend == 'reference' else _production)(case, out, budget, max_lifecycle, preflight)
    if preflight:
        (out / 'preflight.json').write_text(json.dumps(dict(result, runtime=_runtime()), indent=2))
        return
    report = Path(result['report'])
    fits, parameters = _rows(report / 'target_fit.csv'), _rows(report / 'parameters.csv')
    if len(fits) != 14 or len(parameters) != 31: raise RuntimeError('Full 14/31 report row-count contract failed')
    arrays, inventory = _inventory(report)
    plots, closure = _plots(report), _closure(report)
    point, external = _point(case)
    if external:
        sigma = [row for row in parameters if row['parameter'] == 'sigma']
        if len(sigma) != 1 or float(sigma[0]['estimate']) != external['sigma']:
            raise RuntimeError('External sigma did not reach the executed parameter report')
    shutil.copytree(report, out / 'report')
    np.savez_compressed(out / 'arrays.npz', **arrays)
    receipt = dict(backend=backend, case=case, point=point, external_inputs=external, report='report',
        array_inventory=inventory, shared_arrays_observed=True, plot_hashes=plots, closure=closure,
        status=result['status'], lifecycle_solves=int(result['lifecycle_solves']), created_epoch=time.time(), runtime=_runtime(),
        closure_contract='population_one', selection_source=str(SELECTION.relative_to(ROOT)),
        selection_sha256=hashlib.sha256(SELECTION.read_bytes()).hexdigest(),
        native_receipt_sha256=hashlib.sha256(NATIVE_RECEIPT.read_bytes()).hexdigest())
    (out / 'receipt.json').write_text(json.dumps(receipt, default=_json, indent=2, allow_nan=False))


# Explicit descriptive transitions approved for ordinary production solves.
# Numeric bounds/restrictions and all target-table fields remain exact.
SUPPLIED_ROLE = 'supplied primitive; reference calibration bounds advisory'
H0_ROLE = 'derived housing supply coefficient at N0=1; reference bounds advisory'
PARAMETER_STATUS_TRANSITIONS = {
    **{name: ('free in experimental utility calibration', SUPPLIED_ROLE) for name in
       ('beta_annual', 'chi', 'first_birth_fixed_cost', 'kappa_fert',
        'kappa_fert_continuation', 'theta0', 'child_benefit_curvature',
        'tenure_choice_kappa', 'psi_child', 'h_P')},
    'H0': ('derived calibrated housing supply coefficient at N0=1', H0_ROLE),
    'child_benefit_CRRA_coefficient': ('derived from fixed benefit and proposed curvature', 'derived from supplied benefit and curvature'),
}
CLOSURE_ROLE_TRANSITION = ('derived calibrated coefficient at N0=1', 'derived coefficient at N0=1')


def _table_comparison(one, two, count, *, parameter_status_transitions=None):
    if len(one) != count or len(two) != count: return dict(equal=False, reason='missing table rows', allowed_differences=[])
    differences, allowed = [], []
    for index, (xrow, yrow) in enumerate(zip(one, two)):
        if set(xrow) != set(yrow): differences.append(dict(row=index, column='schema', left=sorted(xrow), right=sorted(yrow)))
        for key in set(xrow) | set(yrow):
            x, y = xrow.get(key), yrow.get(key)
            if x == y: continue
            item = dict(row=index, parameter=xrow.get('parameter'), column=key, before=x, after=y)
            transition = (parameter_status_transitions or {}).get(xrow.get('parameter'))
            if key == 'status' and xrow.get('parameter') == yrow.get('parameter') and transition == (x, y):
                allowed.append(dict(item, reason='descriptive role of supplied production primitive; reference search bounds are advisory'))
                continue
            try:
                xn, yn = float(x), float(y)
                equal = xn == yn
            except (TypeError, ValueError): equal = False
            if not equal: differences.append(item)
    return dict(equal=not differences, all_fields_equal=not differences and not allowed,
                differences=differences, allowed_differences=allowed)


def _closure_comparison(one, two, *, reference_first):
    reference, production = (one, two) if reference_first else (two, one)
    allowed, differences = [], []
    extras = {'closure_mode', 'implied_H0_at_population_one', 'fixed_h0_population_scale'}
    for key in sorted(set(reference) | set(production)):
        x, y = reference.get(key), production.get(key)
        if key in reference and key in production and x == y: continue
        item = dict(field=key, before=one.get(key), after=two.get(key))
        if key == 'housing_supply_coefficient_role' and (x, y) == CLOSURE_ROLE_TRANSITION:
            allowed.append(dict(item, reason='descriptive supply-coefficient role, unchanged normalization'))
        elif key not in reference and key in production and key in extras:
            if key == 'closure_mode': valid = y == 'population_one'
            elif key == 'implied_H0_at_population_one': valid = y == production['H0_derived']
            else: valid = y == _selected()['H0_derived'] / production['H0_derived']
            if valid: allowed.append(dict(item, reason='additional production scale metadata, validated against retained H0 and population-one coefficient'))
            else: differences.append(item)
        else: differences.append(item)
    return dict(equal=not differences, all_fields_equal=not differences and not allowed,
                differences=differences, allowed_differences=allowed)


def compare(left, right, out):
    if out.exists(): raise RuntimeError(f'Refusing existing output path: {out}')
    out.mkdir(parents=True)
    try:
        a, b = (json.loads((path / 'receipt.json').read_text()) for path in (left, right))
        same_host_runtime = a['runtime'] == b['runtime']
        same_inputs = all(a[key] == b[key] for key in ('case', 'point', 'external_inputs', 'closure_contract', 'selection_sha256', 'native_receipt_sha256'))
        with np.load(left / 'arrays.npz', allow_pickle=False) as x, np.load(right / 'arrays.npz', allow_pickle=False) as y:
            xa, ya = ({name: source[name] for name in source.files} for source in (x, y))
            _validate_arrays(xa); _validate_arrays(ya)
            comparisons = {name: name in xa and name in ya and xa[name].dtype == ya[name].dtype and np.array_equal(xa[name], ya[name]) for name in sorted(set(xa) | set(ya))}
        fits = _table_comparison(_rows(left / 'report/target_fit.csv'), _rows(right / 'report/target_fit.csv'), 14)
        reference_first = a['backend'] == 'reference' and b['backend'] == 'production'
        production_first = a['backend'] == 'production' and b['backend'] == 'reference'
        transitions = PARAMETER_STATUS_TRANSITIONS if reference_first else ({key: tuple(reversed(value)) for key, value in PARAMETER_STATUS_TRANSITIONS.items()} if production_first else {})
        params = _table_comparison(_rows(left / 'report/parameters.csv'), _rows(right / 'report/parameters.csv'), 31, parameter_status_transitions=transitions)
        plots_a, plots_b = _plots(left / 'report'), _plots(right / 'report')
        closure_a, closure_b = _closure(left / 'report'), _closure(right / 'report')
        closure_check = _closure_comparison(closure_a, closure_b, reference_first=not production_first) if reference_first or production_first else dict(equal=closure_a == closure_b, all_fields_equal=closure_a == closure_b, differences=[], allowed_differences=[])
        source_integrity = closure_a == a['closure'] and closure_b == b['closure']
        receipt = dict(same_host_runtime=same_host_runtime, same_inputs=same_inputs,
            policy_arrays_equal=all(comparisons[name] for name in POLICIES), all_arrays_equal=all(comparisons.values()),
            array_comparison=comparisons, target_fit=fits, parameters=params,
            plot_hashes_equal=plots_a == plots_b == a['plot_hashes'] == b['plot_hashes'],
            closure_equal=closure_check['equal'] and source_integrity, closure_comparison=closure_check,
            parameter_csv_bytes_equal=(left / 'report/parameters.csv').read_bytes() == (right / 'report/parameters.csv').read_bytes(),
            target_csv_bytes_equal=(left / 'report/target_fit.csv').read_bytes() == (right / 'report/target_fit.csv').read_bytes(),
            allowed_metadata_differences=params['allowed_differences'] + closure_check['allowed_differences'],
            expected_scale_only_differences=[item for item in closure_check['allowed_differences'] if item['field'] != 'housing_supply_coefficient_role'],
            comparison_policy='Exact same-host/runtime arrays, targets, numeric parameter restrictions and common numeric closure fields; only enumerated descriptive roles and validated additional scale metadata may differ')
        passed = all(receipt[key] for key in ('same_host_runtime', 'same_inputs', 'policy_arrays_equal', 'all_arrays_equal', 'plot_hashes_equal', 'closure_equal')) and fits['equal'] and params['equal']
        receipt['status'] = 'passed' if passed else 'failed'
        (out / 'comparison.json').write_text(json.dumps(receipt, indent=2, default=_json, allow_nan=False))
        if not passed: raise RuntimeError('Artifact parity failed; inspect ' + str(out / 'comparison.json'))
    except Exception as error:
        if not (out / 'comparison.json').exists():
            (out / 'comparison.json').write_text(json.dumps(dict(status='failed', reason=str(error)), indent=2))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    run = sub.add_parser('execute')
    run.add_argument('--backend', choices=('reference', 'production'), required=True)
    run.add_argument('--case', choices=('unchanged', 'beta', 'external'), required=True)
    run.add_argument('--out', type=Path, required=True)
    run.add_argument('--budget-seconds', type=float, default=1800.)
    run.add_argument('--max-lifecycle', type=int, default=32)
    run.add_argument('--preflight', action='store_true', help='Authenticate and initialize only, zero lifecycle solves')
    diff = sub.add_parser('compare')
    diff.add_argument('--left', type=Path, required=True); diff.add_argument('--right', type=Path, required=True); diff.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'execute': execute(args.backend, args.case, args.out, args.budget_seconds, args.max_lifecycle, args.preflight)
    else: compare(args.left, args.right, args.out)

if __name__ == '__main__': main()
