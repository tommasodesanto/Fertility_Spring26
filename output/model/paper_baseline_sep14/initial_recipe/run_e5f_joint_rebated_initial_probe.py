"""Joint initial root for price, fertility normalization, and equal rebate.

This is an adapter for the pinned initial raw driver.  It replaces the nested
price GE -> rebate root -> fertility root with fixed-price household and KFE
evaluations of three simultaneous residuals.  The analytic stationary PAYGO
pension remains separately bound and certified by the saved source routines.
"""
from __future__ import annotations

import argparse
import copy
import importlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

MAX_EVALUATIONS = 20
HOUSING_TOLERANCE = 2.5e-5
FERTILITY_TOLERANCE = 5e-4
REBATE_TOLERANCE = 1e-6


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def load_module(name, path):
    path = Path(path).resolve(); sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def relative_gap(residual, *flows):
    return abs(float(residual)) / max(*(abs(float(value)) for value in flows), 1e-12)


def residual_passes(residual, *, housing_tolerance=HOUSING_TOLERANCE):
    return (abs(float(residual[0])) <= housing_tolerance
            and abs(float(residual[1])) <= FERTILITY_TOLERANCE
            and abs(float(residual[2])) <= REBATE_TOLERANCE)


def solve_three_residual_root(evaluate, start, *, maximum_evaluations=MAX_EVALUATIONS,
                              housing_tolerance=HOUSING_TOLERANCE, progress=None):
    """Bounded damped Newton root with explicit economically sized differences."""
    if not 4 <= maximum_evaluations <= MAX_EVALUATIONS:
        raise ValueError("Joint root evaluation budget must lie in [4,20]")
    x = np.asarray(start, dtype=float).copy()
    if x.shape != (3,) or not np.isfinite(x).all():
        raise ValueError("Joint root requires three finite transformed starts")
    count = 0; records = []

    def call(point, label):
        nonlocal count
        if count >= maximum_evaluations:
            raise TimeoutError("Joint initial root exhausted its fixed-price evaluation budget")
        result = evaluate(np.asarray(point, dtype=float))
        residual = np.asarray(result["residual"], dtype=float)
        if residual.shape != (3,) or not np.isfinite(residual).all():
            raise RuntimeError("Joint initial residual is incomplete or nonfinite")
        count += 1
        record = dict(evaluation=count, label=label, transformed=np.asarray(point).tolist(),
                      residual=residual.tolist(), score=float(residual @ residual))
        records.append(record)
        if progress is not None: progress(record, result)
        return result

    current = call(x, "start")
    steps = np.array([0.02, 0.02, 0.05])
    caps = np.array([0.30, 0.25, 0.70])
    while not residual_passes(current["residual"], housing_tolerance=housing_tolerance):
        base = np.asarray(current["residual"], dtype=float)
        jacobian = np.empty((3, 3))
        for column in range(3):
            probe = x.copy(); probe[column] += steps[column]
            trial = call(probe, f"jacobian_{column}")
            jacobian[:, column] = (np.asarray(trial["residual"]) - base) / steps[column]
        if not np.isfinite(jacobian).all() or np.linalg.matrix_rank(jacobian) < 3:
            raise RuntimeError("Joint initial residual Jacobian is rank deficient")
        direction = np.linalg.solve(jacobian.T @ jacobian + 1e-8 * np.eye(3), -jacobian.T @ base)
        direction = np.clip(direction, -caps, caps)
        accepted = None
        for scale in (1.0, 0.5, 0.25):
            candidate = x + scale * direction
            candidate[0] = np.clip(candidate[0], math.log(1e-4), math.log(100.0))
            candidate[2] = np.clip(candidate[2], math.log(1e-10), math.log(8.0))
            trial = call(candidate, f"line_{scale:g}")
            if float(np.asarray(trial["residual"]) @ np.asarray(trial["residual"])) < float(base @ base):
                accepted = (candidate, trial); break
        if accepted is None:
            raise RuntimeError("Joint initial root failed its residual-decrease gate")
        x, current = accepted
    return dict(status="converged", transformed=x, result=current,
                evaluations=count, records=records)


def joint_initial_solution(*, model, parameters, b_grid, initial_prices,
                           payroll_tax, marginal_tolerance, fiscal_tolerance,
                           bind_pension, certify_pension, progress=None):
    """Solve the initial stationary closure with one joint three-variable root."""
    if int(parameters.I) != 1 or float(parameters.tau_H) != 0.04:
        raise ValueError("Joint initial probe requires the saved one-market 1% annual tax economy")
    if fiscal_tolerance > REBATE_TOLERANCE or float(parameters.tol_eq) > HOUSING_TOLERANCE:
        raise ValueError("Saved fiscal or housing gate was relaxed")
    calibration = importlib.import_module(model.__package__ + ".calibration")

    def fixed(price, psi, transfer):
        P = copy.deepcopy(parameters)
        P.psi_child = float(psi)
        P.property_tax_lump_sum_transfer = float(transfer)
        P, predicted = bind_pension(P, payroll_tax=payroll_tax)
        solution = model.solve_markov_income_at_prices(
            np.array([price]), P, b_grid, verbose=False, fast_stats=False)
        solution = model.attach_markov_market_accounting(solution, P, b_grid)
        demand, _ = model.markov_market_housing_demand(solution, P, b_grid)
        supply = float(np.asarray(solution.housing_supply).reshape(-1)[0])
        housing = (float(demand[0]) - supply) / max(abs(supply), 1e-12)
        fertility = float(calibration.extract_moments(solution, P)["tfr"])
        revenue = float(solution.property_tax_revenue)
        outlays = float(solution.property_tax_transfer_outlays)
        rebate = float(solution.property_tax_budget_residual)
        return dict(solution=solution, parameters=P, price=np.array([price]),
                    predicted_pension=predicted, completed_fertility=fertility,
                    property_tax_revenue=revenue, rebate_outlays=outlays,
                    residual=np.array([housing, fertility - 2.1,
                                       rebate / max(abs(revenue), abs(outlays), 1e-12)]))

    initial_price = float(np.asarray(initial_prices).reshape(-1)[0])
    pilot = fixed(initial_price, float(parameters.psi_child), 0.0)
    pilot_mass = float(np.asarray(pilot["solution"].g).sum())
    transfer_start = pilot["property_tax_revenue"] / max(pilot_mass, 1e-12)

    def evaluate(transformed):
        return fixed(math.exp(float(transformed[0])), float(transformed[1]),
                     math.exp(float(transformed[2])))

    root = solve_three_residual_root(evaluate,
        np.array([math.log(initial_price), float(parameters.psi_child), math.log(transfer_start)]),
        housing_tolerance=min(float(parameters.tol_eq), HOUSING_TOLERANCE), progress=progress)
    selected = root["result"]
    solution, P, price = selected["solution"], selected["parameters"], selected["price"]
    pension = certify_pension(solution.g, P, marginal_tolerance=marginal_tolerance,
                              fiscal_tolerance=fiscal_tolerance)
    residual = np.asarray(selected["residual"], dtype=float)
    actual_housing_tolerance = min(float(P.tol_eq), HOUSING_TOLERANCE)
    if not residual_passes(residual, housing_tolerance=actual_housing_tolerance):
        raise RuntimeError("Joint initial root returned without all three unchanged gates")
    solution.converged = True
    solution.timings = {**getattr(solution, "timings", {}),
        "strict_converged": True, "accepted": True,
        "best_eq_error": abs(float(residual[0])),
        "joint_initial_root_evaluations": int(root["evaluations"]),
        "convergence_reason": "joint_price_fertility_rebate_root"}
    receipt = dict(pension)
    receipt["joint_initial_root"] = dict(
        status="verified", evaluations=int(root["evaluations"]),
        fixed_price_pilot_evaluations=1, price=float(price[0]),
        psi_child=float(P.psi_child), completed_fertility=float(selected["completed_fertility"]),
        transfer_period_units=float(P.property_tax_lump_sum_transfer),
        housing_tolerance=actual_housing_tolerance,
        housing_relative_residual=float(residual[0]),
        fertility_absolute_gap=float(residual[1]),
        rebate_relative_residual=float(residual[2]),
        property_tax_revenue=float(selected["property_tax_revenue"]),
        rebate_outlays=float(selected["rebate_outlays"]),
        method="simultaneous fixed-price household/KFE residual root")
    return solution, P, price, receipt


def install_on_initial_driver(driver, progress_path=None):
    """Patch only the raw driver's initial closure and normalization dispatcher."""
    bind_pension = driver.bind_initial_balanced_pension
    certify_pension = driver.certify_initial_pension

    def progress(record, _):
        if progress_path is not None: write(progress_path, record)

    def joint(**kwargs):
        return joint_initial_solution(**kwargs, bind_pension=bind_pension,
                                      certify_pension=certify_pension, progress=progress)

    def one_call(chain, base_overrides, *, initial_psi, completed_fertility_target,
                 completed_fertility_tolerance, normalize):
        if (base_overrides or not normalize or completed_fertility_target != 2.1
                or completed_fertility_tolerance != FERTILITY_TOLERANCE):
            raise ValueError("Joint adapter requires the unchanged normalized initial contract")
        started = time.monotonic()
        solution, P, price = chain.run_model_cp_dt({"psi_child": float(initial_psi)})
        moments = chain.extract_moments(solution, P)
        gap = abs(float(moments["tfr"]) - 2.1)
        if gap > FERTILITY_TOLERANCE:
            raise RuntimeError("Joint adapter failed the final completed-fertility gate")
        normalization = dict(status="derived_intercept", psi_child=float(P.psi_child),
            completed_fertility=float(moments["tfr"]), target=2.1, absolute_gap=gap,
            stationary_solves=1, stationary_solve_seconds=time.monotonic()-started,
            fixed_price_root=True, normalization_method="joint_price_fertility_rebate_root")
        return solution, P, price, normalization["stationary_solve_seconds"], normalization

    driver.solve_balanced_initial_equilibrium = joint
    driver.calibration.solve_old_steady_state = one_call
    original_writer = driver.primitive.pf.calendar.write_json_atomic

    def receipt_writer(path, value):
        # The saved driver creates the row before calling the joint closure.
        # Preserve that requested intercept and record the actual solved one.
        if Path(path).name == "stationary_solves.json" and isinstance(value, list):
            value = copy.deepcopy(value)
            for row in value:
                joint_receipt = row.get("fiscal", {}).get("joint_initial_root")
                if joint_receipt is not None:
                    row["requested_psi_child"] = row.get("psi_child")
                    row["solved_psi_child"] = joint_receipt["psi_child"]
                    row["psi_child"] = joint_receipt["psi_child"]
        return original_writer(path, value)

    driver.primitive.pf.calendar.write_json_atomic = receipt_writer
    return driver


def raw_mode(source_root, driver_path, driver_args, progress):
    source_root = Path(source_root).resolve(); driver_path = Path(driver_path).resolve()
    if not driver_path.is_relative_to(source_root):
        raise ValueError("Pinned raw driver escapes source root")
    sys.path[:0] = [str(source_root / "code/model/tools"), str(source_root / "code/model")]
    driver = load_module("joint_rebated_initial_raw_driver", driver_path)
    install_on_initial_driver(driver, progress)
    sys.argv = [str(driver_path), *driver_args]
    driver.main()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--progress", type=Path)
    parser.add_argument("driver_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    values = args.driver_args[1:] if args.driver_args[:1] == ["--"] else args.driver_args
    raw_mode(args.source_root, args.driver, values, args.progress)


if __name__ == "__main__": main()
