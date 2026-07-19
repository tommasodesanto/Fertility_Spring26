#!/usr/bin/env python3
"""Fixed-price production/optimized parity panel around the M5 winner."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from .m5_profile import M5_PRICE, m5_overrides, m5_target_system


CASES: tuple[tuple[str, dict[str, float]], ...] = (
    ("baseline", {}),
    ("beta_lower", {"beta": 0.955}),
    ("consumption_floor_lower", {"c_bar_0": 1.15}),
    ("child_housing_lower", {"psi_child": 0.10}),
    ("deterministic_tenure", {"tenure_choice_kappa": 0.0}),
    ("larger_additional_child_space", {"h_bar_n": 0.35}),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def package_api(package: str):
    if package == "production":
        from intergen_housing_fertility.calibration import extract_moments
        from intergen_housing_fertility.solver import run_model_cp_dt
    else:
        from intergen_housing_fertility_optimized.calibration import extract_moments
        from intergen_housing_fertility_optimized.solver import run_model_cp_dt
    return run_model_cp_dt, extract_moments


def run_case(package: str, changes: dict[str, float]) -> tuple[Any, dict[str, float], float]:
    run_model_cp_dt, extract_moments = package_api(package)
    overrides = m5_overrides(tight=True, optimized=package == "optimized")
    overrides.update(changes)
    overrides.update(
        solve_mode="pe",
        p_fixed=np.array([M5_PRICE]),
        w_fixed=np.array([1.0]),
        entry_shares_fixed=np.array([1.0]),
    )
    start = time.perf_counter()
    solution, parameters, _ = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - start
    moments = extract_moments(solution, parameters)
    system = m5_target_system()
    return solution, {name: float(moments[name]) for name in system.moment_names}, elapsed


def max_abs_array_difference(left: Any, right: Any) -> float:
    """Compare optional policy arrays without hiding representation mismatches."""
    if left is None and right is None:
        return 0.0
    if left is None or right is None:
        return float("inf")
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def main() -> None:
    args = parse_args()
    records: list[dict[str, Any]] = []
    arrays = ("V", "c_pol", "hR_pol", "bp_pol", "tenure_probs", "loc_probs", "fert_probs", "g")
    for name, changes in CASES:
        production, production_moments, production_seconds = run_case("production", changes)
        optimized, optimized_moments, optimized_seconds = run_case("optimized", changes)
        records.append(
            {
                "case": name,
                "changes": changes,
                "production_seconds": production_seconds,
                "optimized_seconds": optimized_seconds,
                "speedup": production_seconds / optimized_seconds,
                "max_abs_target_moment_difference": max(
                    abs(production_moments[key] - optimized_moments[key])
                    for key in production_moments
                ),
                "array_max_abs_differences": {
                    key: max_abs_array_difference(
                        getattr(production, key), getattr(optimized, key)
                    )
                    for key in arrays
                },
                "production_mass": float(np.sum(production.g)),
                "optimized_mass": float(np.sum(optimized.g)),
            }
        )
    payload = {
        "cases": records,
        "max_abs_target_moment_difference": max(
            record["max_abs_target_moment_difference"] for record in records
        ),
        "minimum_speedup": min(record["speedup"] for record in records),
        "median_speedup": float(np.median([record["speedup"] for record in records])),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
