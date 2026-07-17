#!/usr/bin/env python3
"""Fixed-parameter A/B and bounded polish for bequest normalization."""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Any

from intergen_housing_fertility.calibration import (
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.local_panel import (
    income_process_overrides,
    run_local_polish,
)
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_BOUNDS,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    REPAIRED_TIMING_SWITCHES,
    production_profile_overrides,
    validate_production_profile,
)
from intergen_housing_fertility.solver import run_model_cp_dt


MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
DEFAULT_SEED = (
    Path(__file__).resolve().parents[3]
    / "output/model/intergen_bequest_normalization_ab_20260710/fixed_theta_seed.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--theta-json", type=Path, default=DEFAULT_SEED)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("../../output/model/intergen_bequest_normalization_ab_20260710"),
    )
    parser.add_argument("--run-polish", action="store_true")
    parser.add_argument("--polish-max-evals", type=int, default=50)
    parser.add_argument("--polish-minutes", type=float, default=30.0)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def load_theta(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text())
    raw = payload.get("theta", payload)
    if not isinstance(raw, dict):
        raise ValueError(f"{path} does not contain a theta object")
    theta = {str(key): float(value) for key, value in raw.items()}
    expected = {"beta" if name == "beta_annual" else name for name, _, _ in PRODUCTION_SEARCH_BOUNDS}
    if set(theta) != expected:
        raise ValueError(f"theta keys differ from production bounds: {sorted(set(theta) ^ expected)}")
    return theta


def target_rows(
    moments: dict[str, float], targets: dict[str, float], weights: dict[str, float]
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights[name])
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
    return rows


def parameter_rows(theta: dict[str, float]) -> list[dict[str, Any]]:
    rows = []
    for name, lower, upper in PRODUCTION_SEARCH_BOUNDS:
        key = "beta" if name == "beta_annual" else name
        estimate = float(theta[key]) ** 0.25 if name == "beta_annual" else float(theta[key])
        span = float(upper - lower)
        near = min(estimate - float(lower), float(upper) - estimate) <= 0.02 * span
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower": float(lower),
                "upper": float(upper),
                "near_bound": bool(near),
            }
        )
    return rows


def run_case(
    label: str,
    normalized: bool,
    theta: dict[str, float],
    targets: dict[str, float],
    weights: dict[str, float],
    quiet: bool,
) -> dict[str, Any]:
    income = income_process_overrides(
        5,
        "rouwenhorst",
        MATCHED_ANNUAL_INNOVATION_SD,
        MATCHED_ANNUAL_RHO,
    )
    overrides = {
        **base_overrides(
            J=PRODUCTION_J,
            Nb=PRODUCTION_SEARCH_NB,
            n_house=5,
            max_iter_eq=PRODUCTION_MAX_ITER_EQ,
        ),
        **production_profile_overrides(),
        **REPAIRED_TIMING_SWITCHES,
        **income,
        **theta,
        "normalize_bequest_utility": bool(normalized),
    }
    started = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=not quiet)
    moments = extract_moments(sol, P)
    residual = float(moments["market_residual"])
    timings = dict(getattr(sol, "timings", {}))
    strict = bool(
        timings.get("strict_converged", getattr(sol, "converged", False))
        and math.isfinite(residual)
        and residual <= float(P.tol_eq)
    )
    rows = target_rows(moments, targets, weights)
    return {
        "label": label,
        "normalize_bequest_utility": bool(normalized),
        "strict_converged": strict,
        "market_residual": residual,
        "rank_loss": float(diagnostic_loss(moments, targets=targets, weights=weights)),
        "elapsed_sec": float(time.perf_counter() - started),
        "p_eq": jsonable(p_eq),
        "theta": dict(theta),
        "moments": jsonable(moments),
        "target_fit": rows,
        "timings": jsonable(timings),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    validate_production_profile(
        PRODUCTION_PROFILE_NAME,
        J=PRODUCTION_J,
        Nb=PRODUCTION_SEARCH_NB,
        n_house=5,
        income_states=5,
        target_set=PRODUCTION_TARGET_SET,
        max_iter_eq=PRODUCTION_MAX_ITER_EQ,
        stage="search",
    )
    theta = load_theta(args.theta_json)
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    cases = [
        run_case("unnormalized", False, theta, targets, weights, args.quiet),
        run_case("normalized", True, theta, targets, weights, args.quiet),
    ]
    if not all(case["strict_converged"] for case in cases):
        raise RuntimeError("fixed-parameter A/B did not strictly converge")
    if cases[0]["theta"] != cases[1]["theta"]:
        raise RuntimeError("A/B parameter vectors differ")

    fit_rows = [
        {"case": case["label"], **row}
        for case in cases
        for row in case["target_fit"]
    ]
    params = parameter_rows(theta)
    payload: dict[str, Any] = {
        "status": "fixed_parameter_bequest_normalization_ab",
        "specification_adopted": False,
        "theta_source": str(args.theta_json.resolve()),
        "parameters_identical_across_cases": True,
        "production_profile": PRODUCTION_PROFILE_NAME,
        "target_set": PRODUCTION_TARGET_SET,
        "Nb": PRODUCTION_SEARCH_NB,
        "income_process": {
            "method": "rouwenhorst",
            "states": 5,
            "annual_rho": MATCHED_ANNUAL_RHO,
            "annual_innovation_sd": MATCHED_ANNUAL_INNOVATION_SD,
        },
        "parameters": params,
        "cases": cases,
    }
    write_csv(args.outdir / "fixed_ab_target_fit.csv", fit_rows)
    write_csv(args.outdir / "fixed_ab_parameters.csv", params)
    (args.outdir / "fixed_ab_comparison.json").write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True)
    )

    if args.run_polish:
        polish = run_local_polish(
            args.outdir / "normalized_local_polish",
            max_evals=min(int(args.polish_max_evals), 50),
            seed=20260710,
            J=PRODUCTION_J,
            Nb=PRODUCTION_SEARCH_NB,
            n_house=5,
            max_iter_eq=PRODUCTION_MAX_ITER_EQ,
            minutes=min(float(args.polish_minutes), 30.0),
            income_states=5,
            income_process="rouwenhorst",
            income_annual_rho=MATCHED_ANNUAL_RHO,
            income_innovation_sd=MATCHED_ANNUAL_INNOVATION_SD,
            target_set=PRODUCTION_TARGET_SET,
            seed_theta=theta,
            method="pattern",
            initial_step=0.04,
            min_step=0.003,
            shrink=0.5,
            profile_name=PRODUCTION_PROFILE_NAME,
            normalize_bequest_utility=True,
            progress=not args.quiet,
        )
        payload["normalized_local_polish"] = polish
        (args.outdir / "fixed_ab_comparison.json").write_text(
            json.dumps(jsonable(payload), indent=2, sort_keys=True)
        )

    print(
        json.dumps(
            {
                "outdir": str(args.outdir),
                "cases": [
                    {
                        "label": case["label"],
                        "rank_loss": case["rank_loss"],
                        "market_residual": case["market_residual"],
                    }
                    for case in cases
                ],
                "polish_ran": bool(args.run_polish),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
