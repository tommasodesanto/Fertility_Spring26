#!/usr/bin/env python3
"""Run one checkpointed task from the optimized-package promotion battery."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any

for _thread_variable in (
    "NUMBA_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
):
    os.environ.setdefault(_thread_variable, "1")

import numpy as np

from .m5_profile import M5_PRICE, m5_overrides, m5_target_system
from .promotion_contract import (
    ARRAY_NAMES,
    CALIBRATION_CASES,
    DEMAND_PRICES,
    PARITY_TOLERANCES,
    ROOT_CASES,
    ROOT_TOLERANCES,
    bound_cases,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "task",
        choices=("bound-parity", "demand-point", "root-case", "strict-repeat", "calibration"),
    )
    parser.add_argument("--index", type=int, default=1, help="One-based task/case index.")
    parser.add_argument("--package", choices=("production", "optimized"), default="optimized")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="Use Nb=30 but preserve the exact loop shape.")
    return parser.parse_args()


def package_api(package: str):
    if package == "production":
        from intergen_housing_fertility.calibration import extract_moments
        from intergen_housing_fertility.solver import run_model_cp_dt
    else:
        from intergen_housing_fertility_optimized.calibration import extract_moments
        from intergen_housing_fertility_optimized.solver import run_model_cp_dt
    return run_model_cp_dt, extract_moments


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def exception_record(exc: Exception) -> dict[str, Any]:
    return {
        "status": "infeasible" if type(exc).__name__ == "InfeasibleThetaError" else "error",
        "exception_type": type(exc).__name__,
        "message": str(exc),
        "stage": getattr(exc, "stage", None),
        "dead_mass": getattr(exc, "dead_mass", None),
    }


def run_model(
    package: str,
    changes: dict[str, Any],
    *,
    mode: str,
    price: float = M5_PRICE,
    smoke: bool = False,
) -> tuple[dict[str, Any], Any | None]:
    run_model_cp_dt, extract_moments = package_api(package)
    optimized = package == "optimized"
    overrides = m5_overrides(tight=mode == "ge", optimized=optimized)
    overrides.update(changes)
    if smoke:
        overrides["Nb"] = 30
    if mode == "pe":
        overrides.update(
            solve_mode="pe",
            p_fixed=np.array([price]),
            w_fixed=np.array([1.0]),
            entry_shares_fixed=np.array([1.0]),
        )
    started = time.perf_counter()
    try:
        solution, parameters, equilibrium_price = run_model_cp_dt(overrides, verbose=False)
    except Exception as exc:  # Case-level checkpoints must survive invalid bound points.
        record = exception_record(exc)
        record["elapsed_seconds"] = time.perf_counter() - started
        return record, None
    elapsed = time.perf_counter() - started
    moments = extract_moments(solution, parameters)
    target_system = m5_target_system()
    selected_moments = {name: float(moments[name]) for name in target_system.moment_names}
    loss = target_system.loss(selected_moments)
    record = {
        "status": "ok",
        "elapsed_seconds": elapsed,
        "price": np.asarray(equilibrium_price, dtype=float).reshape(-1),
        "market_residual": float(getattr(solution, "best_max_abs_rel_excess", np.nan)),
        "mass": float(np.sum(solution.g)),
        "loss": float(loss),
        "moments": selected_moments,
        "timings": getattr(solution, "timings", {}),
        "housing_demand": float(np.asarray(solution.housing_demand).reshape(-1)[0]),
        "housing_supply": float(np.asarray(solution.housing_supply).reshape(-1)[0]),
    }
    record["excess_demand"] = record["housing_demand"] - record["housing_supply"]
    return record, solution


def optional_array_difference(left: Any, right: Any) -> float:
    if left is None and right is None:
        return 0.0
    if left is None or right is None:
        return float("inf")
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def compare_failures(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    same_status = left["status"] == right["status"]
    matching_infeasibility = same_status and left["status"] == "infeasible"
    same_type = left.get("exception_type") == right.get("exception_type")
    same_stage = left.get("stage") == right.get("stage")
    lm = left.get("dead_mass")
    rm = right.get("dead_mass")
    if lm is None or rm is None:
        mass_close = lm is rm
        mass_difference = None
    else:
        mass_difference = abs(float(lm) - float(rm))
        mass_close = mass_difference <= max(
            PARITY_TOLERANCES["exception_dead_mass_abs"],
            PARITY_TOLERANCES["exception_dead_mass_rel"] * max(abs(float(lm)), abs(float(rm))),
        )
    return {
        "same_status": same_status,
        "same_exception_type": same_type,
        "same_stage": same_stage,
        "dead_mass_abs_difference": mass_difference,
        "pass": bool(matching_infeasibility and same_type and same_stage and mass_close),
    }


def run_bound_parity(index: int, output: Path, smoke: bool) -> None:
    case = bound_cases()[index - 1]
    production, production_solution = run_model("production", case["changes"], mode="pe", smoke=smoke)
    optimized, optimized_solution = run_model("optimized", case["changes"], mode="pe", smoke=smoke)
    if production["status"] != "ok" or optimized["status"] != "ok":
        comparison = compare_failures(production, optimized)
    else:
        array_differences = {
            name: optional_array_difference(
                getattr(production_solution, name, None), getattr(optimized_solution, name, None)
            )
            for name in ARRAY_NAMES
        }
        moment_difference = max(
            abs(production["moments"][name] - optimized["moments"][name])
            for name in production["moments"]
        )
        mass_difference = abs(production["mass"] - optimized["mass"])
        comparison = {
            "array_max_abs_differences": array_differences,
            "target_moment_max_abs_difference": moment_difference,
            "mass_abs_difference": mass_difference,
            "pass": bool(
                max(array_differences.values()) <= PARITY_TOLERANCES["array_max_abs"]
                and moment_difference <= PARITY_TOLERANCES["target_moment_max_abs"]
                and mass_difference <= PARITY_TOLERANCES["mass_abs"]
            ),
        }
    write_json(
        output,
        {
            "task": "bound-parity",
            "index": index,
            "case": case,
            "smoke": smoke,
            "production": production,
            "optimized": optimized,
            "comparison": comparison,
        },
    )


def run_demand_point(index: int, output: Path, smoke: bool) -> None:
    price = DEMAND_PRICES[index - 1]
    production, _ = run_model("production", {}, mode="pe", price=price, smoke=smoke)
    optimized, _ = run_model("optimized", {}, mode="pe", price=price, smoke=smoke)
    if production["status"] != "ok" or optimized["status"] != "ok":
        comparison = compare_failures(production, optimized)
    else:
        difference = abs(production["excess_demand"] - optimized["excess_demand"])
        comparison = {
            "excess_demand_abs_difference": difference,
            "sign_equal": bool(np.signbit(production["excess_demand"]) == np.signbit(optimized["excess_demand"])),
            "pass": bool(difference <= ROOT_TOLERANCES["fixed_map_excess_abs"]),
        }
    write_json(
        output,
        {
            "task": "demand-point",
            "index": index,
            "price": price,
            "smoke": smoke,
            "production": production,
            "optimized": optimized,
            "comparison": comparison,
        },
    )


def run_root_case(index: int, output: Path, smoke: bool) -> None:
    case = ROOT_CASES[index - 1]
    production, _ = run_model("production", case["changes"], mode="ge", smoke=smoke)
    optimized, _ = run_model("optimized", case["changes"], mode="ge", smoke=smoke)
    both_ok = production["status"] == optimized["status"] == "ok"
    comparison = (
        {
            "price_abs_difference": abs(production["price"][0] - optimized["price"][0]),
            "loss_abs_difference": abs(production["loss"] - optimized["loss"]),
            "optimized_speedup": production["elapsed_seconds"] / optimized["elapsed_seconds"],
            "production_strict": production["market_residual"] <= ROOT_TOLERANCES["strict_residual"],
            "optimized_strict": optimized["market_residual"] <= ROOT_TOLERANCES["strict_residual"],
            "pass": bool(
                production["market_residual"] <= ROOT_TOLERANCES["strict_residual"]
                and optimized["market_residual"] <= ROOT_TOLERANCES["strict_residual"]
            ),
        }
        if both_ok
        else compare_failures(production, optimized)
    )
    if not both_ok:
        comparison["pass"] = False
        comparison["reason"] = "root_case_must_be_feasible_in_both_packages"
    write_json(
        output,
        {
            "task": "root-case",
            "index": index,
            "case": case,
            "smoke": smoke,
            "production": production,
            "optimized": optimized,
            "comparison": comparison,
        },
    )


def run_strict_repeat(index: int, output: Path, smoke: bool) -> None:
    optimized, _ = run_model("optimized", {}, mode="ge", smoke=smoke)
    write_json(
        output,
        {"task": "strict-repeat", "index": index, "smoke": smoke, "optimized": optimized},
    )


def run_calibration(package: str, output: Path, smoke: bool) -> None:
    records: list[dict[str, Any]] = []
    started = time.perf_counter()
    for index, case in enumerate(CALIBRATION_CASES, start=1):
        result, _ = run_model(package, case["changes"], mode="ge", smoke=smoke)
        record = {"index": index, "case": case, "result": result}
        records.append(record)
        write_json(
            output,
            {
                "task": "calibration",
                "package": package,
                "smoke": smoke,
                "status": "running",
                "elapsed_seconds": time.perf_counter() - started,
                "records": records,
            },
        )
    feasible = [record for record in records if record["result"]["status"] == "ok"]
    best = min(feasible, key=lambda record: record["result"]["loss"]) if feasible else None
    write_json(
        output,
        {
            "task": "calibration",
            "package": package,
            "smoke": smoke,
            "status": "complete",
            "elapsed_seconds": time.perf_counter() - started,
            "records": records,
            "best": best,
        },
    )


def main() -> None:
    args = parse_args()
    if args.task == "bound-parity":
        run_bound_parity(args.index, args.output, args.smoke)
    elif args.task == "demand-point":
        run_demand_point(args.index, args.output, args.smoke)
    elif args.task == "root-case":
        run_root_case(args.index, args.output, args.smoke)
    elif args.task == "strict-repeat":
        run_strict_repeat(args.index, args.output, args.smoke)
    else:
        run_calibration(args.package, args.output, args.smoke)


if __name__ == "__main__":
    main()
