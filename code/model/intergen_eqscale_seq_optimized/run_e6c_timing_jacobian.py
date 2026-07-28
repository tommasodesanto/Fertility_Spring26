#!/usr/bin/env python3
"""Strict local E6c timing Jacobian at the combined certified winner."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.e6c_profile import E6C_DOMAIN, E6C_SEED  # noqa: E402
from intergen_eqscale_seq_optimized.local_panel import jsonable  # noqa: E402


MOMENTS = ("mean_age_first_birth", "share_first_births_age30plus")
PARAMETERS = ("readiness_location_age", "readiness_spread_years")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--location-step", type=float, default=0.25)
    parser.add_argument("--spread-step", type=float, default=0.10)
    return parser.parse_args()


def configure_environment() -> None:
    os.environ.update(
        {
            "E3_L4": "1",
            "E5": "1",
            "E6A": "1",
            "E6B": "1",
            "E6C": "1",
            "E3_TFR_TOP_BIN_WEIGHT": "3.602359422009",
        }
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.location_step <= 0.0 or args.spread_step <= 0.0:
        raise ValueError("finite-difference steps must be positive")
    payload = json.loads(args.source.read_text())
    winner = ((payload.get("winners") or {}).get("E1"))
    if not isinstance(winner, dict) or not isinstance(winner.get("theta"), dict):
        raise ValueError(f"{args.source} does not contain winners.E1.theta")
    metadata = payload.get("winner_metadata") or {}
    if metadata.get("arm", winner.get("arm")) != "E6AB":
        raise ValueError("E6c Jacobian source must be the certified E6AB winner")
    repeat = payload.get("winner_repeat_check") or {}
    if not (
        repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", math.inf)) == 0.0
        and float(repeat.get("max_abs_moment_difference", math.inf)) == 0.0
    ):
        raise ValueError("E6c Jacobian source must preserve exact strict repeats")

    configure_environment()
    from intergen_eqscale_seq_optimized import run_e1_chain
    from intergen_eqscale_seq_optimized.e5_profile import e5_target_system

    chain = importlib.reload(run_e1_chain)
    chain.load_runtime()
    runtime_args = argparse.Namespace(
        J=int(metadata.get("J", 17)),
        Nb=int(metadata.get("Nb", 120)),
        max_iter_eq=40,
        tol_eq=2.5e-5,
    )
    overrides = chain.common_overrides(runtime_args)
    baseline = {name: float(value) for name, value in winner["theta"].items()}
    baseline.update(E6C_SEED)
    steps = {
        "readiness_location_age": float(args.location_step),
        "readiness_spread_years": float(args.spread_step),
    }
    cases: list[tuple[str, dict[str, float]]] = [("baseline", baseline)]
    for parameter in PARAMETERS:
        for direction, sign in (("minus", -1.0), ("plus", 1.0)):
            theta = dict(baseline)
            theta[parameter] += sign * steps[parameter]
            cases.append((f"{parameter}_{direction}", theta))

    system = e5_target_system()
    targets, weights = system.targets_dict(), system.weights_dict()
    records: list[dict[str, Any]] = []
    by_label: dict[str, dict[str, Any]] = {}
    for label, theta in cases:
        started = time.perf_counter()
        solution, parameters, prices = chain.run_model_cp_dt(
            {**overrides, **theta},
            verbose=False,
        )
        elapsed = time.perf_counter() - started
        moments = chain.extract_moments(solution, parameters)
        residual = float(getattr(solution, "best_max_abs_rel_excess", math.inf))
        timings = dict(getattr(solution, "timings", {}))
        strict = bool(
            timings.get(
                "strict_converged",
                getattr(solution, "converged", False),
            )
            and residual <= 2.5e-5
        )
        if not strict:
            raise RuntimeError(
                f"{label} is not a strict solve: residual={residual:.6g}"
            )
        record = {
            "label": label,
            "strict_converged": strict,
            "market_residual": residual,
            "elapsed_sec": elapsed,
            "price": float(np.asarray(prices, dtype=float).reshape(-1)[0]),
            "loss": float(
                chain.diagnostic_loss(
                    moments,
                    targets=targets,
                    weights=weights,
                )
            ),
            **{name: float(theta[name]) for name in PARAMETERS},
            **{name: float(moments[name]) for name in MOMENTS},
        }
        records.append(record)
        by_label[label] = record

    jacobian = np.zeros((len(MOMENTS), len(PARAMETERS)))
    for column, parameter in enumerate(PARAMETERS):
        minus = by_label[f"{parameter}_minus"]
        plus = by_label[f"{parameter}_plus"]
        for row, moment in enumerate(MOMENTS):
            jacobian[row, column] = (
                float(plus[moment]) - float(minus[moment])
            ) / (2.0 * steps[parameter])
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    rank = int(np.linalg.matrix_rank(jacobian))
    condition_number = float(
        singular_values[0] / singular_values[-1]
        if singular_values[-1] > 0.0
        else math.inf
    )
    timing_standard_errors = np.array([0.25, 0.008], dtype=float)
    parameter_ranges = np.array(
        [upper - lower for _name, lower, upper, _kind in E6C_DOMAIN],
        dtype=float,
    )
    standardized_jacobian = (
        np.diag(1.0 / timing_standard_errors)
        @ jacobian
        @ np.diag(parameter_ranges)
    )
    standardized_singular_values = np.linalg.svd(
        standardized_jacobian,
        compute_uv=False,
    )
    standardized_condition_number = float(
        standardized_singular_values[0] / standardized_singular_values[-1]
        if standardized_singular_values[-1] > 0.0
        else math.inf
    )
    summary = {
        "source": str(args.source),
        "arm": "E6ABC",
        "baseline": {name: baseline[name] for name in PARAMETERS},
        "steps": steps,
        "moments": list(MOMENTS),
        "parameters": list(PARAMETERS),
        "jacobian": jacobian,
        "singular_values": singular_values,
        "rank": rank,
        "condition_number": condition_number,
        "timing_standard_errors": timing_standard_errors,
        "parameter_ranges": parameter_ranges,
        "standardized_jacobian": standardized_jacobian,
        "standardized_singular_values": standardized_singular_values,
        "standardized_condition_number": standardized_condition_number,
        "strict_case_count": len(records),
        "all_cases_strict": all(record["strict_converged"] for record in records),
    }
    if rank != 2:
        raise RuntimeError(
            f"E6c timing Jacobian is not rank two: rank={rank}, matrix={jacobian}"
        )
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "strict_cases.csv", records)
    write_csv(
        args.outdir / "jacobian.csv",
        [
            {
                "moment": moment,
                **{
                    parameter: float(jacobian[row, column])
                    for column, parameter in enumerate(PARAMETERS)
                },
            }
            for row, moment in enumerate(MOMENTS)
        ],
    )
    (args.outdir / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )


if __name__ == "__main__":
    main()
