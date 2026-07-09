#!/usr/bin/env python3
"""Compare the exact frozen source specification with repaired timing."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from intergen_housing_fertility.calibration import (
    PERIOD_YEARS,
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.local_panel import income_process_overrides
from intergen_housing_fertility.production_profile import (
    PRODUCTION_INCOME_STATES,
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    PRODUCTION_VERIFY_NB,
    comparison_arm_switches,
    production_profile_metadata,
    production_profile_overrides,
    validate_frozen_source_theta,
    validate_production_profile,
)
from intergen_housing_fertility.solver import run_model_cp_dt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta-json", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--grid", choices=["search", "verify"], default="search")
    parser.add_argument("--max-iter-eq", type=int, default=PRODUCTION_MAX_ITER_EQ)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def load_theta(path: Path) -> tuple[dict[str, float], dict[str, Any]]:
    payload = json.loads(path.read_text())
    raw = payload.get("theta", payload)
    if not isinstance(raw, dict):
        raise ValueError(f"{path} does not contain a theta object")
    return {str(key): float(value) for key, value in raw.items()}, payload


def run_case(
    label: str,
    theta: dict[str, float],
    *,
    nb: int,
    max_iter_eq: int,
    targets: dict[str, float],
    weights: dict[str, float],
    verbose: bool,
) -> dict[str, Any]:
    timing_switches = comparison_arm_switches(label)
    overrides = {
        **base_overrides(J=PRODUCTION_J, Nb=nb, n_house=5, max_iter_eq=max_iter_eq),
        **income_process_overrides(PRODUCTION_INCOME_STATES),
        **production_profile_overrides(),
        **theta,
        **timing_switches,
    }
    started = time.perf_counter()
    sol, p, prices = run_model_cp_dt(overrides, verbose=verbose)
    elapsed = time.perf_counter() - started
    moments = extract_moments(sol, p)
    residual = float(moments.get("market_residual", math.inf))
    timings = dict(getattr(sol, "timings", {}))
    strict = bool(
        timings.get("strict_converged", getattr(sol, "converged", False))
        and np.isfinite(residual)
        and residual <= float(getattr(p, "tol_eq", 1e-4))
    )
    rows = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights.get(name, 1.0))
        rows.append(
            {
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap**2,
            }
        )
    return {
        "label": label,
        "theta": theta,
        "timing_switches": timing_switches,
        "prices": jsonable(prices),
        "market_residual": residual,
        "strict_converged": strict,
        "loss": diagnostic_loss(moments, targets=targets, weights=weights),
        "elapsed_sec": elapsed,
        "moments": moments,
        "target_fit": rows,
        "timings": jsonable(timings),
        "has_reproducible_beginning_asset_distribution": hasattr(
            sol, "g_beginning_assets_by_current_choice"
        ),
    }


def parameter_table(theta: dict[str, float]) -> list[dict[str, Any]]:
    rows = []
    for bound_name, lower, upper in PRODUCTION_SEARCH_BOUNDS:
        source_name = "beta" if bound_name == "beta_annual" else bound_name
        estimate = float(theta[source_name])
        reported_estimate = estimate ** (1.0 / PERIOD_YEARS) if bound_name == "beta_annual" else estimate
        width = max(float(upper - lower), 1e-12)
        distance = min(reported_estimate - lower, upper - reported_estimate) / width
        rows.append(
            {
                "parameter": bound_name,
                "model_key": source_name,
                "estimate": reported_estimate,
                "lower": float(lower),
                "upper": float(upper),
                "near_bound": bool(distance <= 0.02),
                "comparison_status": "held at supplied seed; free in the planned polish",
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    nb = PRODUCTION_SEARCH_NB if args.grid == "search" else PRODUCTION_VERIFY_NB
    validate_production_profile(
        PRODUCTION_PROFILE_NAME,
        J=PRODUCTION_J,
        Nb=nb,
        n_house=5,
        income_states=PRODUCTION_INCOME_STATES,
        target_set=PRODUCTION_TARGET_SET,
        max_iter_eq=args.max_iter_eq,
        stage=args.grid,
    )
    theta, theta_source = load_theta(args.theta_json)
    validate_frozen_source_theta(theta)
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    args.outdir.mkdir(parents=True, exist_ok=True)
    frozen_source = run_case(
        "frozen_source_repro",
        theta,
        nb=nb,
        max_iter_eq=args.max_iter_eq,
        targets=targets,
        weights=weights,
        verbose=not args.quiet,
    )
    repaired = run_case(
        "repaired_timing",
        theta,
        nb=nb,
        max_iter_eq=args.max_iter_eq,
        targets=targets,
        weights=weights,
        verbose=not args.quiet,
    )
    differences = [
        {
            "moment": name,
            "frozen_source_repro": float(frozen_source["moments"].get(name, math.nan)),
            "repaired_timing": float(repaired["moments"].get(name, math.nan)),
            "repaired_minus_frozen": float(repaired["moments"].get(name, math.nan))
            - float(frozen_source["moments"].get(name, math.nan)),
        }
        for name in targets
    ]
    payload = {
        "status": "frozen_source_vs_repaired_timing_comparison",
        "profile": production_profile_metadata(),
        "grid_stage": args.grid,
        "theta_source_path": str(args.theta_json.resolve()),
        "theta_source_metadata": theta_source,
        "frozen_source_repro": frozen_source,
        "repaired_timing": repaired,
        "differences": differences,
        "parameters": parameter_table(theta),
    }
    (args.outdir / "comparison.json").write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True))
    write_csv(args.outdir / "target_fit_frozen_source_repro.csv", frozen_source["target_fit"])
    write_csv(args.outdir / "target_fit_repaired_timing.csv", repaired["target_fit"])
    write_csv(args.outdir / "moment_differences.csv", differences)
    write_csv(args.outdir / "parameter_ledger.csv", payload["parameters"])
    print(json.dumps(jsonable(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
