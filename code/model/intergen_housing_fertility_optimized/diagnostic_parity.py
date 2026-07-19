#!/usr/bin/env python3
"""Generate and compare production/optimized M5 diagnostic packets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
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

import matplotlib.image as mpimg
import numpy as np

from .m5_profile import M5_PRICE, M5_PROFILE_NAME, M5_THETA, m5_overrides, m5_target_system


NUMERIC_TOLERANCE = 5.0e-10
PIXEL_TOLERANCE = 1.0 / 255.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true", help="Use Nb=30 for exact-loop preflight.")
    return parser.parse_args()


def package_api(package: str):
    if package == "production":
        from intergen_housing_fertility.calibration import extract_moments
        from intergen_housing_fertility.diagnostics import write_diagnostics
        from intergen_housing_fertility.solver import run_model_cp_dt
    else:
        from intergen_housing_fertility_optimized.calibration import extract_moments
        from intergen_housing_fertility_optimized.diagnostics import write_diagnostics
        from intergen_housing_fertility_optimized.solver import run_model_cp_dt
    return run_model_cp_dt, extract_moments, write_diagnostics


def solve_and_write(package: str, outdir: Path, smoke: bool) -> dict[str, Any]:
    from tools import build_intergen_mechanics_packet as mechanics

    run_model_cp_dt, extract_moments, write_diagnostics = package_api(package)
    overrides = m5_overrides(tight=True, optimized=package == "optimized")
    if smoke:
        overrides["Nb"] = 30
    overrides.update(
        solve_mode="pe",
        p_fixed=np.array([M5_PRICE]),
        w_fixed=np.array([1.0]),
        entry_shares_fixed=np.array([1.0]),
        diagnostic_policy_ages=np.array([30.0, 42.0]),
    )
    started = time.perf_counter()
    solution, parameters, price = run_model_cp_dt(overrides, verbose=False)
    elapsed = time.perf_counter() - started
    moments = extract_moments(solution, parameters)
    target_system = m5_target_system()
    selected_moments = {name: float(moments[name]) for name in target_system.moment_names}
    baseline = {
        "label": "canonical_fixed_price",
        "sol": solution,
        "P": parameters,
        "p_eq": np.asarray(price, dtype=float).reshape(-1),
        "moments": moments,
        "rank_loss": target_system.loss(selected_moments),
        "market_residual": float(getattr(solution, "best_max_abs_rel_excess", np.nan)),
        "elapsed_sec": elapsed,
        "overrides": overrides,
    }
    source = {
        "theta": dict(M5_THETA),
        "source_meta": {
            "profile": M5_PROFILE_NAME,
            "price": M5_PRICE,
            "mode": "canonical_fixed_price",
        },
    }
    outdir.mkdir(parents=True, exist_ok=True)
    write_diagnostics(solution, parameters, outdir / "diagnostics")
    mechanics.write_core_outputs(
        outdir,
        baseline,
        source,
        target_system.targets_dict(),
        target_system.weights_dict(),
        include_contact_sheet=True,
        quick_first_look_only=False,
        write_csv_sidecars=True,
    )
    mechanics.write_classic_draft_outputs(outdir / "classic_draft", baseline)
    return {
        "package": package,
        "elapsed_seconds": elapsed,
        "loss": baseline["rank_loss"],
        "market_residual": baseline["market_residual"],
        "mass": float(np.sum(solution.g)),
    }


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_number(text: str) -> float | None:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def compare_csv(left: Path, right: Path) -> dict[str, Any]:
    with left.open(newline="") as stream:
        left_rows = list(csv.reader(stream))
    with right.open(newline="") as stream:
        right_rows = list(csv.reader(stream))
    if len(left_rows) != len(right_rows) or any(
        len(left_row) != len(right_row)
        for left_row, right_row in zip(left_rows, right_rows)
    ):
        return {"pass": False, "reason": "shape_mismatch"}
    maximum = 0.0
    text_equal = True
    for left_row, right_row in zip(left_rows, right_rows):
        for left_value, right_value in zip(left_row, right_row):
            left_number = parse_number(left_value)
            right_number = parse_number(right_value)
            if left_number is not None and right_number is not None:
                maximum = max(maximum, abs(left_number - right_number))
            elif left_value != right_value:
                text_equal = False
    return {
        "max_abs_numeric_difference": maximum,
        "nonnumeric_text_equal": text_equal,
        "pass": bool(maximum <= NUMERIC_TOLERANCE and text_equal),
    }


IGNORED_JSON_KEYS = {
    "elapsed_sec",
    "elapsed_seconds",
    "elapsed_sec_excluding_model_solve",
    "solver_timings",
    "timings",
    "package",
}


def compare_json_values(left: Any, right: Any) -> tuple[float, bool]:
    if isinstance(left, dict) and isinstance(right, dict):
        left_keys = set(left) - IGNORED_JSON_KEYS
        right_keys = set(right) - IGNORED_JSON_KEYS
        if left_keys != right_keys:
            return float("inf"), False
        maximum = 0.0
        for key in left_keys:
            difference, equal = compare_json_values(left[key], right[key])
            maximum = max(maximum, difference)
            if not equal:
                return maximum, False
        return maximum, True
    if isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            return float("inf"), False
        maximum = 0.0
        for left_item, right_item in zip(left, right):
            difference, equal = compare_json_values(left_item, right_item)
            maximum = max(maximum, difference)
            if not equal:
                return maximum, False
        return maximum, True
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if not math.isfinite(float(left)) or not math.isfinite(float(right)):
            return 0.0, str(left) == str(right)
        difference = abs(float(left) - float(right))
        return difference, difference <= NUMERIC_TOLERANCE
    return 0.0, left == right


def compare_json(left: Path, right: Path) -> dict[str, Any]:
    difference, equal = compare_json_values(
        json.loads(left.read_text()), json.loads(right.read_text())
    )
    return {"max_abs_numeric_difference": difference, "pass": bool(equal)}


def compare_png(left: Path, right: Path) -> dict[str, Any]:
    if file_digest(left) == file_digest(right):
        return {"sha256_equal": True, "max_abs_pixel_difference": 0.0, "pass": True}
    left_image = np.asarray(mpimg.imread(left), dtype=float)
    right_image = np.asarray(mpimg.imread(right), dtype=float)
    if left_image.shape != right_image.shape:
        return {"sha256_equal": False, "shape_equal": False, "pass": False}
    difference = float(np.max(np.abs(left_image - right_image)))
    return {
        "sha256_equal": False,
        "shape_equal": True,
        "max_abs_pixel_difference": difference,
        "pass": bool(difference <= PIXEL_TOLERANCE),
    }


def compare_packets(left_root: Path, right_root: Path) -> dict[str, Any]:
    suffixes = {".csv", ".json", ".md", ".png"}
    left_files = {
        path.relative_to(left_root): path
        for path in left_root.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    }
    right_files = {
        path.relative_to(right_root): path
        for path in right_root.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    }
    common = sorted(set(left_files) & set(right_files))
    records: dict[str, Any] = {}
    for relative in common:
        left = left_files[relative]
        right = right_files[relative]
        if relative.suffix == ".csv":
            comparison = compare_csv(left, right)
        elif relative.suffix == ".json":
            comparison = compare_json(left, right)
        elif relative.suffix == ".png":
            comparison = compare_png(left, right)
        else:
            comparison = {"text_equal": left.read_text() == right.read_text()}
            comparison["pass"] = comparison["text_equal"]
        records[str(relative)] = comparison
    missing_left = sorted(str(path) for path in set(right_files) - set(left_files))
    missing_right = sorted(str(path) for path in set(left_files) - set(right_files))
    return {
        "production_file_count": len(left_files),
        "optimized_file_count": len(right_files),
        "missing_from_production": missing_left,
        "missing_from_optimized": missing_right,
        "files": records,
        "pass": bool(
            not missing_left
            and not missing_right
            and records
            and all(record["pass"] for record in records.values())
        ),
    }


def main() -> None:
    args = parse_args()
    if args.outdir.exists():
        shutil.rmtree(args.outdir)
    production_root = args.outdir / "production"
    optimized_root = args.outdir / "optimized"
    production = solve_and_write("production", production_root, args.smoke)
    optimized = solve_and_write("optimized", optimized_root, args.smoke)
    comparison = compare_packets(production_root, optimized_root)
    payload = {
        "smoke": args.smoke,
        "production": production,
        "optimized": optimized,
        "comparison": comparison,
    }
    (args.outdir / "comparison.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not comparison["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
