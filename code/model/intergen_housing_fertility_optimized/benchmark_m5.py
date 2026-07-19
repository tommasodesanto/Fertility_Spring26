#!/usr/bin/env python3
"""Benchmark and parity-check production against the isolated optimized M5."""

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

from .m5_profile import M5_LOSS, M5_PRICE, m5_overrides, m5_target_system


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", choices=("production", "optimized", "both"), default="both")
    parser.add_argument("--mode", choices=("fixed-price", "tight-ge"), default="fixed-price")
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def _package_api(package: str):
    if package == "production":
        from intergen_housing_fertility.calibration import diagnostic_loss, extract_moments
        from intergen_housing_fertility.solver import run_model_cp_dt
    else:
        from intergen_housing_fertility_optimized.calibration import diagnostic_loss, extract_moments
        from intergen_housing_fertility_optimized.solver import run_model_cp_dt
    return run_model_cp_dt, extract_moments, diagnostic_loss


def run_case(package: str, mode: str) -> tuple[dict[str, Any], Any]:
    run_model_cp_dt, extract_moments, diagnostic_loss = _package_api(package)
    optimized = package == "optimized"
    overrides = m5_overrides(tight=True, optimized=optimized)
    if mode == "fixed-price":
        overrides.update(
            solve_mode="pe",
            p_fixed=np.array([M5_PRICE]),
            w_fixed=np.array([1.0]),
            entry_shares_fixed=np.array([1.0]),
        )
    start = time.perf_counter()
    solution, parameters, price = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - start
    moments = extract_moments(solution, parameters)
    target_system = m5_target_system()
    loss = diagnostic_loss(
        moments,
        targets=target_system.targets_dict(),
        weights=target_system.weights_dict(),
    )
    result = {
        "package": package,
        "mode": mode,
        "elapsed_seconds": elapsed,
        "price": np.asarray(price, dtype=float).tolist(),
        "loss": float(loss),
        "canonical_loss": M5_LOSS,
        "market_residual": float(solution.best_max_abs_rel_excess),
        "mass": float(np.sum(solution.g)),
        "timings": getattr(solution, "timings", {}),
        "moments": {name: float(moments[name]) for name in target_system.moment_names},
    }
    return result, solution


def array_parity(production: Any, optimized: Any) -> dict[str, Any]:
    comparisons: dict[str, Any] = {}
    for name in (
        "V",
        "c_pol",
        "hR_pol",
        "bp_pol",
        "tenure_choice",
        "tenure_probs",
        "loc_probs",
        "fert_probs",
        "fert_value",
        "g",
    ):
        left = np.asarray(getattr(production, name))
        right = np.asarray(getattr(optimized, name))
        comparisons[name] = {
            "shape_equal": left.shape == right.shape,
            "max_abs_difference": float(np.max(np.abs(left - right))),
            "exactly_equal": bool(np.array_equal(left, right)),
        }
    return comparisons


def main() -> None:
    args = parse_args()
    if args.package == "both" and args.mode != "fixed-price":
        raise ValueError("--package both is restricted to fixed-price parity; benchmark GE separately")
    packages = ("production", "optimized") if args.package == "both" else (args.package,)
    records: dict[str, Any] = {}
    solutions: dict[str, Any] = {}
    for package in packages:
        records[package], solutions[package] = run_case(package, args.mode)
    payload: dict[str, Any] = {"runs": records}
    if args.package == "both":
        payload["array_parity"] = array_parity(solutions["production"], solutions["optimized"])
        payload["speedup"] = (
            records["production"]["elapsed_seconds"]
            / records["optimized"]["elapsed_seconds"]
        )
        payload["max_abs_target_moment_difference"] = max(
            abs(records["production"]["moments"][name] - records["optimized"]["moments"][name])
            for name in records["production"]["moments"]
        )
    encoded = json.dumps(payload, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded + "\n")
    print(encoded)


if __name__ == "__main__":
    main()
