"""Calibration adapter for the common production stationary GE solver.

Search bounds and the derived-H0 acceptance interval belong here, never in GE.
The full target/weight contract is authenticated before invoking a solve.
"""
from __future__ import annotations
import copy
import csv
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any, Callable, Mapping
import numpy as np
from .inputs import DEFAULT_PRICE, load_inputs

ROOT = Path(__file__).resolve().parents[5]
ANCHOR = ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
SNAPSHOT = ROOT / 'code/model/production/reference_inputs/bundle.json'
H0_BOUNDS = (.2, 80.)

class H0BoundError(RuntimeError):
    """Derived supply coefficient is outside the calibration acceptance interval."""

def _canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def _contract():
    anchor = json.loads(ANCHOR.read_text())
    snapshot = json.loads(SNAPSHOT.read_text())
    target = [{k: row[k] for k in ('moment', 'target', 'weight', 'role')} for row in anchor['target_fit']]
    fingerprint = _canonical(target)
    weights = _canonical(dict(base_contract=target, multipliers={}))
    if not fingerprint == snapshot['target_fingerprint'] == anchor['target_fingerprint']:
        raise RuntimeError('Production complete target-and-weight fingerprint drift')
    if not weights == snapshot['weight_fingerprint'] == anchor['weight_fingerprint']:
        raise RuntimeError('Production weight fingerprint drift')
    bounds = {row['parameter']: (float(row['lower']), float(row['upper'])) for row in anchor['parameters'] if row['parameter'] in anchor['selected']['parameters']}
    return target, fingerprint, weights, bounds

def _rows(path):
    with Path(path).open(newline='') as stream: return list(csv.DictReader(stream))

def residual_from_report(report_directory, expected_contract=None):
    report = Path(report_directory)
    fits, parameters = _rows(report / 'target_fit.csv'), _rows(report / 'parameters.csv')
    if len(fits) != 14 or len(parameters) != 31:
        raise RuntimeError('production report requires 14 target rows and 31 parameter rows')
    target = [{k: row[k] for k in ('moment', 'target', 'weight', 'role')} for row in fits]
    if expected_contract is None: expected_contract = _contract()[0]
    if target != expected_contract: raise RuntimeError('Solved target-and-weight contract drift')
    scored = [row for row in fits if row['role'] == 'scored']
    values = np.asarray([math.sqrt(float(row['weight'])) * float(row['gap']) for row in scored])
    loss = sum(float(row['loss_contribution']) for row in scored)
    if values.shape != (10,) or not np.isfinite(values).all() or abs(float(values @ values) - loss) > 1e-8:
        raise RuntimeError('ten finite scored residuals with exact loss arithmetic required')
    return values, fits, parameters

def _serial(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, Path): return str(value)
    if isinstance(value, Mapping): return {str(k): _serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [_serial(v) for v in value]
    return value

def effective_input_fingerprint(P, grid):
    """Hash exact candidate-bound caller inputs before any solver mutations.

    Only the native diagnostics-output directory is normalized. Array identity
    includes dtype, shape and contents; optional scalar infinities/NaNs have
    deterministic binary representations rather than invalid JSON numbers.
    """
    import struct

    def encode(value):
        if isinstance(value, np.ndarray):
            content = ([encode(item) for item in value.flat] if value.dtype.hasobject
                       else hashlib.sha256(value.tobytes(order='C')).hexdigest())
            return ['ndarray', value.dtype.str, list(value.shape), content]
        if isinstance(value, np.generic):
            return ['numpy_scalar', value.dtype.str, value.tobytes().hex()]
        if value is None or isinstance(value, (str, bool, int)):
            return [type(value).__name__, value]
        if isinstance(value, float):
            return ['float', struct.pack('>d', value).hex()]
        if isinstance(value, Path):
            return ['Path', str(value)]
        if isinstance(value, Mapping):
            if any(not isinstance(key, str) for key in value):
                raise TypeError('input identity requires string mapping keys')
            return ['mapping', {key: encode(item) for key, item in value.items()}]
        if isinstance(value, (list, tuple)):
            return [type(value).__name__, [encode(item) for item in value]]
        raise TypeError('unsupported input identity type: ' + type(value).__name__)

    fields = dict(vars(P))
    if 'native_inherited_distribution_evidence_dir' in fields:
        fields['native_inherited_distribution_evidence_dir'] = ''
    return _canonical(dict(schema='canonical_caller_inputs_v1',
                           parameters=encode(fields), wealth_grid=encode(np.asarray(grid))))


def make_evaluator(out, lane, P, grid, deadline, price_start=None, *, native_runner=None,
                   exploratory=False, solver: Callable[..., Mapping[str, Any]] | None = None,
                   target_fingerprint=None, weight_fingerprint=None, bounds_override=None):
    """Preserve caller primitives/grid and bind only the ten search coordinates."""
    del lane, native_runner, exploratory
    target, target_pin, weight_pin, bounds = _contract()
    if bounds_override is not None:
        revised = {k: tuple(map(float, v)) for k, v in bounds_override.items()}
        expected = dict(bounds, beta_annual=(.93, .99))
        if revised != expected:
            raise ValueError("Only estate-A beta bound .93-.99 is supported")
        bounds = revised
    if target_fingerprint is not None and target_fingerprint != target_pin:
        raise RuntimeError('Caller target fingerprint differs from production contract')
    if weight_fingerprint is not None and weight_fingerprint != weight_pin:
        raise RuntimeError('Caller weight fingerprint differs from production contract')
    if P is None:
        if grid is not None: raise ValueError('grid requires caller P')
        base_P, base_grid = load_inputs()
    else:
        if grid is None: raise ValueError('caller P requires its grid')
        base_P, base_grid = copy.deepcopy(P), np.asarray(grid).copy()
    if solver is None:
        from .equilibrium import solve_stationary_ge
        solver = solve_stationary_ge
    destination = Path(out)
    start_price = DEFAULT_PRICE if price_start is None else float(price_start)

    def evaluate(label, point, end):
        effective_end = min(float(end), float(deadline))
        if time.time() >= effective_end:
            return dict(status='budget_exhausted', reason='evaluation deadline reached', lifecycle_solves=0)
        if set(point) != set(bounds): raise ValueError('Wrong ten calibration coordinates')
        for key, value in point.items():
            if not math.isfinite(float(value)) or not bounds[key][0] <= float(value) <= bounds[key][1]:
                raise ValueError('Out of calibration bounds: ' + key)
        from .inputs import bind_parameters
        b_grid = base_grid.copy()
        model_P = bind_parameters(base_P, b_grid, point)
        input_fingerprint = effective_input_fingerprint(model_P, b_grid)
        result = solver(model_P, b_grid, out=destination / str(label), price_start=start_price,
                        budget_seconds=effective_end - time.time(), max_lifecycle=32, closure='population_one')
        if result.get('status', 'passed') != 'passed':
            return _serial({k: result[k] for k in ('status', 'reason', 'lifecycle_solves', 'price_search') if k in result})
        report = Path(result['report_directory'])
        residual, fits, parameters = residual_from_report(report, target)
        closure = result['closure']
        h0 = float(closure['H0_derived'])
        if not H0_BOUNDS[0] <= h0 <= H0_BOUNDS[1]:
            return dict(status='inadmissible_numerical', reason=f'Derived normalized H0={h0:.16g} outside {H0_BOUNDS}',
                        lifecycle_solves=int(result['lifecycle_solves']), derived_H0_bound_rejection=True,
                        rejection_kind='derived_H0_constraint')
        loss = float(residual @ residual)
        receipt = dict(status='passed', report=str(report), residual=residual.tolist(), loss=loss, objective=loss,
            population=1., normalization='N0=1', H0_derived=h0, price=float(result['price']), starting_price=start_price,
            target_fit=fits, parameter_table=parameters, closure=closure, lifecycle_solves=int(result['lifecycle_solves']),
            target_fingerprint=target_pin, weight_fingerprint=weight_pin,
            input_fingerprint=input_fingerprint, input_fingerprint_schema='canonical_caller_inputs_v1',
            parameter_report_bound_semantics='reference calibration search bounds are advisory in GE; enforced by the calibration adapter only')
        receipt = _serial(receipt)
        json.dumps(receipt, allow_nan=False)
        return receipt
    return evaluate


def check_adapter_without_solves():
    """Exercise input propagation and optimizer receipt compatibility with a mock."""
    anchor = json.loads(ANCHOR.read_text())
    report = ANCHOR.parent / 'native_postcheck/selected_postcheck/phase_b_ge/selected_root'
    P, grid = load_inputs()
    P.sigma = 2.01
    P.adapter_check_array = np.asarray([123.])
    expected_grid = grid.copy()
    calls = []
    def mock_solver(Q, supplied_grid, **kwargs):
        if Q is P or Q.sigma != 2.01 or Q.adapter_check_array is P.adapter_check_array:
            raise AssertionError('caller primitives were reset or aliased')
        np.testing.assert_array_equal(supplied_grid, expected_grid)
        if kwargs['closure'] != 'population_one': raise AssertionError('wrong calibration closure')
        calls.append(kwargs)
        return dict(report_directory=str(report), closure=json.loads((report / 'closure.json').read_text()),
                    price=anchor['selected']['price'], lifecycle_solves=2)
    end = time.time() + 120.
    evaluate = make_evaluator('unused_mock_destination', 'floor_s0', P, grid, end, solver=mock_solver)
    receipt = evaluate('mock', anchor['selected']['parameters'], end)
    if receipt['lifecycle_solves'] != 2 or receipt['price'] != anchor['selected']['price'] or 'solution' in receipt:
        raise AssertionError('optimizer receipt compatibility failed')
    if abs(receipt['loss'] - anchor['native_loss']) > 1e-12 or len(calls) != 1:
        raise AssertionError('residual/loss arithmetic or common solver call failed')
    json.dumps(dict(label='mock', parameters=anchor['selected']['parameters'], **receipt), allow_nan=False)
    try:
        make_evaluator('unused', 'floor_s0', P, grid, end, solver=mock_solver, target_fingerprint='wrong')
    except RuntimeError: pass
    else: raise AssertionError('unmatched target fingerprint did not fail before solving')
    return dict(status='passed_zero_solves', actual_lifecycle_solves=0, mocked_lifecycle_count=receipt['lifecycle_solves'],
                supplied_P_grid_preserved=True, common_solver_called=True, receipt_serializable=True,
                full_target_pin_rejected=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(check_adapter_without_solves(), indent=2))


def export_verified_parameters(destination, completed_receipt, P, grid, *, budget_seconds=1800):
    """Export a run-local input file after the final fresh native acceptance.

    Only the adopted post-interest/old-target production contract is supported.
    Preserve caller fixed primitives, the verified price and derived H0. Reload
    the file and compare every effective input/grid before publishing; unsupported
    contracts or primitive mappings fail rather than reset to reference values.
    This never promotes or overwrites code/model/parameters/best_params.py.
    """
    import pprint
    from .inputs import bind_parameters
    from .parameter_files import BEST_FILE, PARAMETER_ROOT, load_parameter_file

    destination = Path(destination).resolve()
    if destination.name != 'best_params.py' or destination.is_relative_to(PARAMETER_ROOT.resolve()):
        raise ValueError('calibration export must be a run-local best_params.py')
    if destination == BEST_FILE.resolve() or destination.exists():
        raise ValueError('refusing to overwrite an existing parameter file')
    source = Path(completed_receipt).resolve(strict=True)
    receipt = json.loads(source.read_text())
    if receipt.get('status') not in ('selected_numerically_verified', 'selected_native_passed_search_loss_differs'):
        raise ValueError('parameter export requires a completed fresh native selected-point check')
    if receipt.get('arm') != 'alternative' or receipt.get('selected_evaluator') != 'full_native':
        raise ValueError('only the adopted post-interest production selected point can be exported')
    target, target_pin, weight_pin, bounds = _contract()
    if receipt.get('target_fingerprint') != target_pin or receipt.get('weight_fingerprint') != weight_pin:
        raise RuntimeError('calibration export target/weight contract differs from production')
    selected = receipt.get('selected', {})
    verification = receipt.get('selected_postcheck', {})
    repeat = receipt.get('repeat', {})
    if selected.get('status') != 'passed' or verification.get('status') != 'passed' or repeat.get('status') != 'exact_full_ge_repeat_passed':
        raise ValueError('accepted selected point and exact native repeat are required')
    if repeat.get('target_rows') != 14 or repeat.get('parameter_rows') != 31 or len(repeat.get('standard_plot_hashes', {})) != 17:
        raise ValueError('native repeat table/plot verification is incomplete')
    point = selected.get('parameters', {})
    if set(point) != set(bounds):
        raise ValueError('export requires the complete ten-coordinate selected point')
    for key, value in point.items():
        if not math.isfinite(float(value)) or not bounds[key][0] <= float(value) <= bounds[key][1]:
            raise ValueError('selected export coordinate out of bounds: ' + key)
    if verification.get('input_fingerprint_schema') != 'canonical_caller_inputs_v1' or not verification.get('input_fingerprint'):
        raise ValueError('native receipt lacks authenticated caller-input identity; historical receipts cannot export')
    effective = bind_parameters(P, np.asarray(grid), point)
    candidate_fingerprint = effective_input_fingerprint(effective, grid)
    if candidate_fingerprint != verification['input_fingerprint']:
        raise ValueError('export caller inputs/grid differ from the fresh native verified candidate')
    report = Path(verification['report'])
    residual, fits, rows = residual_from_report(report, target)
    loss = float(residual @ residual)
    if abs(loss - float(receipt['native_loss'])) > 1e-8:
        raise RuntimeError('completed native export loss differs from verified report')
    reported = {row['parameter']: float(row['estimate']) for row in rows}
    for key, value in point.items():
        if reported.get(key) != float(value):
            raise RuntimeError('native report selected-coordinate drift: ' + key)
    h0 = float(verification['H0_derived'])
    price = float(verification['price'])
    if not H0_BOUNDS[0] <= h0 <= H0_BOUNDS[1] or not math.isfinite(price) or price <= 0:
        raise ValueError('verified native H0/price is inadmissible')
    if reported.get('H0') != h0:
        raise RuntimeError('native report derived H0 differs from selected verification')
    effective.H0 = np.full_like(effective.H0, h0)
    # This is a diagnostics output directory, not an economic primitive. Native
    # solves always replace it with the new run's stage directory.
    effective.native_inherited_distribution_evidence_dir = ''
    _, reference_grid = load_inputs(point)
    if not np.array_equal(np.asarray(grid), reference_grid):
        raise ValueError('cannot export a changed wealth grid using the authenticated input snapshot')

    def same(a, b):
        if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
            return isinstance(a, np.ndarray) and isinstance(b, np.ndarray) and a.dtype == b.dtype and np.array_equal(a, b, equal_nan=True)
        if isinstance(a, Mapping) and isinstance(b, Mapping):
            return set(a) == set(b) and all(same(a[k], b[k]) for k in a)
        if isinstance(a, (list, tuple)) and isinstance(b, type(a)):
            return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
        if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
            return True
        return type(a) == type(b) and a == b

    editable = load_parameter_file(BEST_FILE)
    external = {key: _serial(getattr(effective, key)) for key in editable['external_inputs']}
    reconstructed, _ = load_inputs(point, external)
    derived = {'q', 'user_cost_rate', 'rho', 'rho_hat', 'eps_fert', 'income', 'pension', 'pension_by_loc'}
    mapped = {'beta', 'hbar_first_child_jump'} | (set(point) - {'beta_annual', 'h_P'})
    native = {}
    differences = []
    for key, value in vars(effective).items():
        if not hasattr(reconstructed, key):
            differences.append(key + ' (absent from reference snapshot)')
        elif not same(value, getattr(reconstructed, key)):
            if key in derived | mapped:
                differences.append(key + ' (derived input cannot be exported independently)')
            else:
                native[key] = _serial(value)
    if differences:
        raise ValueError('cannot faithfully export effective inputs: ' + ', '.join(differences))
    provenance = dict(role='run_local_verified_calibration', completed_receipt=str(source),
        completed_receipt_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        target_fingerprint=target_pin, weight_fingerprint=weight_pin, native_loss=loss,
        arm='alternative', chain=receipt.get('chain'), global_default_promoted=False,
        normalized_population=1., native_repeat='exact_full_ge_repeat_passed',
        input_snapshot_sha256=hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest(),
        native_candidate_input_fingerprint=candidate_fingerprint,
        native_candidate_input_fingerprint_schema='canonical_caller_inputs_v1',
        exported_input_fingerprint=effective_input_fingerprint(effective, grid),
        wealth_grid_sha256=hashlib.sha256(np.asarray(grid).tobytes()).hexdigest(),
        effective_field_count=len(vars(effective)),
        output_only_evidence_directory_reset=True)
    config = dict(PARAMETERS=dict(point), EXTERNAL_INPUTS=external, NATIVE_OVERRIDES=native,
                  PRICE_GUESS=price, CLOSURE='fixed_h0', BUDGET_SECONDS=budget_seconds, PROVENANCE=provenance)
    text = ('"""Run-local verified calibration inputs; no global promotion.\n'
            'The stationary solver uses these inputs without searching parameters.\n"""\n\n')
    text += '\n\n'.join(key + ' = ' + pprint.pformat(value, sort_dicts=False, width=88) for key, value in config.items()) + '\n'
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name('.best_params_pending.py')
    if temporary.exists():
        raise ValueError('pending parameter export already exists')
    try:
        temporary.write_text(text)
        loaded = load_parameter_file(temporary)
        restored, restored_grid = load_inputs(loaded['parameters'], loaded['external_inputs'], loaded['native_overrides'])
        mismatches = sorted(set(vars(effective)) ^ set(vars(restored)))
        mismatches += [key for key in vars(effective) if hasattr(restored, key) and not same(getattr(effective, key), getattr(restored, key))]
        if mismatches or not np.array_equal(restored_grid, np.asarray(grid)):
            raise ValueError('export roundtrip differs from effective inputs: ' + ', '.join(mismatches))
        # Exclusive creation also rejects a file created while validation ran.
        with destination.open('x') as stream:
            stream.write(text)
    finally:
        temporary.unlink(missing_ok=True)
    return dict(path=str(destination), sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
                effective_field_count=len(vars(effective)), full_input_roundtrip=True,
                global_default_promoted=False, target_fingerprint=target_pin, weight_fingerprint=weight_pin)
