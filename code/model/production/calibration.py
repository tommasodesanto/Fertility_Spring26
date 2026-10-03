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

ROOT = Path(__file__).resolve().parents[3]
ANCHOR = ROOT / 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
SNAPSHOT = Path(__file__).with_name('reference_inputs') / 'bundle.json'
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

def make_evaluator(out, lane, P, grid, deadline, price_start=None, *, native_runner=None,
                   exploratory=False, solver: Callable[..., Mapping[str, Any]] | None = None,
                   target_fingerprint=None, weight_fingerprint=None):
    """Preserve caller primitives/grid and bind only the ten search coordinates."""
    del lane, native_runner, exploratory
    target, target_pin, weight_pin, bounds = _contract()
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
