#!/usr/bin/env python3
"""Compare the production income process with annual-AR(1) Rouwenhorst alternatives.

The Rouwenhorst cases start from the Sommer--Sullivan annual specification
``rho=0.90`` and annual innovation standard deviation ``0.20``.  They are
sampled every four model years before discretization; ``0.20`` is therefore
never treated as an unconditional standard deviation.

Examples
--------
Construction checks only (no solve)::

    PYTHONPATH=$PWD .venv/bin/python tools/compare_intergen_income_processes.py \
      --construction-only

Fixed-theta comparison (the production default is ``Nb=120``)::

    PYTHONPATH=$PWD .venv/bin/python tools/compare_intergen_income_processes.py \
      --outdir ../../output/model/income_process_comparison
"""

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
    base_overrides,
    diagnostic_loss,
    extract_moments,
    get_target_set,
    jsonable,
)
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    REPAIRED_TIMING_SWITCHES,
    production_profile_overrides,
    validate_production_profile,
)
from intergen_housing_fertility.solver import run_model_cp_dt


ANNUAL_RHO = 0.90
ANNUAL_INNOVATION_SD = 0.20
PERIOD_YEARS = 4
DEFAULT_THETA_PATH = (
    Path(__file__).resolve().parents[3]
    / "output/model/intergen_repair_full13_nb120_best6942_policy_battery_20260710/policy_summary.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construction-only", action="store_true", help="Print process diagnostics and do not solve.")
    parser.add_argument("--outdir", type=Path, help="Directory for comparison JSON and CSV files (required unless construction-only).")
    parser.add_argument("--theta-json", type=Path, default=DEFAULT_THETA_PATH, help="Policy summary containing benchmark.theta.")
    parser.add_argument("--nb", type=int, default=PRODUCTION_SEARCH_NB, help="Asset-grid nodes for fixed-theta solves (default: 120).")
    parser.add_argument("--max-iter-eq", type=int, default=PRODUCTION_MAX_ITER_EQ)
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=("current", "rouwenhorst5", "rouwenhorst7"),
        default=("current", "rouwenhorst5"),
        help="Cases to solve; the default keeps the income-state count fixed at five.",
    )
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args()


def stationary_distribution(Pi: np.ndarray) -> np.ndarray:
    """Return the invariant row-vector distribution, with numerical checks."""
    Pi = np.asarray(Pi, dtype=float)
    if Pi.ndim != 2 or Pi.shape[0] != Pi.shape[1]:
        raise ValueError("transition matrix must be square")
    if np.max(np.abs(Pi.sum(axis=1) - 1.0)) > 1e-12:
        raise ValueError("transition matrix rows must sum to one")
    system = np.vstack((Pi.T - np.eye(Pi.shape[0]), np.ones(Pi.shape[0])))
    rhs = np.zeros(Pi.shape[0] + 1)
    rhs[-1] = 1.0
    weights = np.linalg.lstsq(system, rhs, rcond=None)[0]
    weights = np.maximum(weights, 0.0)
    return weights / weights.sum()


def rouwenhorst_log_process(nz: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Four-year Rouwenhorst process implied by the stated annual AR(1)."""
    rho4 = ANNUAL_RHO**PERIOD_YEARS
    sigma_eps4 = ANNUAL_INNOVATION_SD * math.sqrt(sum(ANNUAL_RHO ** (2 * k) for k in range(PERIOD_YEARS)))
    sigma_x = sigma_eps4 / math.sqrt(1.0 - rho4**2)
    p = (1.0 + rho4) / 2.0
    Pi = np.array([[p, 1.0 - p], [1.0 - p, p]], dtype=float)
    for size in range(3, int(nz) + 1):
        old = Pi
        Pi = np.zeros((size, size), dtype=float)
        Pi[:-1, :-1] += p * old
        Pi[:-1, 1:] += (1.0 - p) * old
        Pi[1:, :-1] += (1.0 - p) * old
        Pi[1:, 1:] += p * old
        Pi[1:-1] *= 0.5
    Pi /= Pi.sum(axis=1, keepdims=True)
    log_grid = np.linspace(-sigma_x * math.sqrt(nz - 1), sigma_x * math.sqrt(nz - 1), nz)
    weights = stationary_distribution(Pi)
    z_grid = np.exp(log_grid)
    z_grid /= float(weights @ z_grid)  # Level productivity has stationary mean exactly one.
    return z_grid, weights, Pi, {
        "annual_rho": ANNUAL_RHO,
        "annual_innovation_sd": ANNUAL_INNOVATION_SD,
        "period_years": PERIOD_YEARS,
        "rho_four_year": rho4,
        "innovation_sd_four_year": sigma_eps4,
        "unconditional_log_sd": sigma_x,
        "unconditional_log_variance": sigma_x**2,
    }


def process_diagnostics(z_grid: np.ndarray, weights: np.ndarray, Pi: np.ndarray) -> dict[str, Any]:
    z_grid, weights, Pi = map(lambda x: np.asarray(x, dtype=float), (z_grid, weights, Pi))
    weights = weights / weights.sum()
    log_z = np.log(z_grid)
    log_mean = float(weights @ log_z)
    log_var = float(weights @ (log_z - log_mean) ** 2)
    log_cov1 = float(np.sum(weights[:, None] * Pi * (log_z[:, None] - log_mean) * (log_z[None, :] - log_mean)))
    return {
        "z_grid": z_grid.tolist(),
        "stationary_weights": weights.tolist(),
        "transition_matrix": Pi.tolist(),
        "stationary_weight_error_max": float(np.max(np.abs(weights @ Pi - weights))),
        "transition_row_sum_error_max": float(np.max(np.abs(Pi.sum(axis=1) - 1.0))),
        "weighted_level_mean": float(weights @ z_grid),
        "weighted_level_variance": float(weights @ (z_grid - weights @ z_grid) ** 2),
        "weighted_log_mean": log_mean,
        "weighted_log_variance": log_var,
        "first_order_log_autocorrelation": float(log_cov1 / log_var) if log_var > 0 else math.nan,
    }


def all_processes() -> dict[str, dict[str, Any]]:
    current = production_profile_overrides()
    result: dict[str, dict[str, Any]] = {
        "current": {
            "description": "Current five-state production stay/redraw process.",
            "overrides": {key: current[key] for key in ("use_income_types", "income_type_transition", "z_grid", "z_weights", "Pi_z")},
        }
    }
    for nz in (5, 7):
        z_grid, weights, Pi, conversion = rouwenhorst_log_process(nz)
        result[f"rouwenhorst{nz}"] = {
            "description": f"{nz}-state Rouwenhorst approximation to the four-year sampled annual AR(1).",
            "conversion": conversion,
            "overrides": {
                "use_income_types": True,
                "income_type_transition": "markov",
                "z_grid": z_grid,
                "z_weights": weights,
                "Pi_z": Pi,
                "income_shock_persistence": conversion["rho_four_year"],
            },
        }
    for spec in result.values():
        override = spec["overrides"]
        spec["diagnostics"] = process_diagnostics(override["z_grid"], override["z_weights"], override["Pi_z"])
        spec["overrides"] = jsonable(override)
    return result


def load_benchmark_theta(path: Path) -> dict[str, float]:
    payload = json.loads(path.read_text())
    # The policy battery records the common fixed-theta source under
    # ``source_record`` (its cases contain policy perturbations, not the
    # benchmark parameter vector).
    raw = payload.get("source_record", {}).get("theta")
    if not isinstance(raw, dict):
        raise ValueError(f"Could not find source_record.theta in {path}")
    theta = {str(key): float(value) for key, value in raw.items()}
    if len(theta) != 13:
        raise ValueError(f"Benchmark theta is ambiguous or incomplete: expected 13 values, found {len(theta)}")
    return theta


def target_rows(moments: dict[str, float], targets: dict[str, float], weights: dict[str, float]) -> list[dict[str, Any]]:
    rows = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights.get(name, 1.0))
        rows.append({"moment": name, "target": float(target), "model": model, "gap": gap, "weight": weight, "loss_contribution": weight * gap**2})
    return rows


def run_case(label: str, spec: dict[str, Any], theta: dict[str, float], args: argparse.Namespace, targets: dict[str, float], weights: dict[str, float]) -> dict[str, Any]:
    income = spec["overrides"]
    overrides = {
        **base_overrides(J=PRODUCTION_J, Nb=args.nb, n_house=5, max_iter_eq=args.max_iter_eq),
        **production_profile_overrides(),
        **REPAIRED_TIMING_SWITCHES,
        **theta,
        **income,
    }
    started = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=not args.quiet)
    moments = extract_moments(sol, P)
    timings = dict(getattr(sol, "timings", {}))
    residual = float(moments["market_residual"])
    target_fit = target_rows(moments, targets, weights)
    headline_keys = ("tfr", "childless_rate", "own_rate", "aggregate_own_rate", "own_rate_2534", "own_rate_3544", "old_age_own_rate", "own_family_gap", "young_liquid_wealth_to_income", "mean_age_first_birth")
    return {
        "label": label,
        "status": "ok",
        "rank_loss": float(diagnostic_loss(moments, targets=targets, weights=weights)),
        "market_residual": residual,
        "elapsed_sec": float(time.perf_counter() - started),
        "strict_converged": bool(timings.get("strict_converged", getattr(sol, "converged", False)) and residual <= float(P.tol_eq)),
        "p_eq": jsonable(p_eq),
        "income_process": spec["diagnostics"],
        "headline_statistics": {key: float(moments.get(key, math.nan)) for key in headline_keys},
        "moments": jsonable(moments),
        "target_fit": target_fit,
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
    processes = all_processes()
    construction = {label: {"description": spec["description"], "conversion": spec.get("conversion"), "diagnostics": spec["diagnostics"]} for label, spec in processes.items()}
    if args.construction_only:
        print(json.dumps(construction, indent=2, sort_keys=True))
        return
    if args.outdir is None:
        raise SystemExit("--outdir is required unless --construction-only is used")
    if args.nb < 2:
        raise SystemExit("--nb must be at least 2")
    # The profile validator deliberately enforces the exact published
    # ``Nb=120, max_iter_eq=10`` configuration.  Retain that check for the
    # production default, but allow explicitly requested small-grid smoke and
    # sensitivity comparisons without changing any production source file.
    if args.nb == PRODUCTION_SEARCH_NB and args.max_iter_eq == PRODUCTION_MAX_ITER_EQ:
        validate_production_profile(
            PRODUCTION_PROFILE_NAME,
            J=PRODUCTION_J,
            Nb=args.nb,
            n_house=5,
            income_states=5,
            target_set=PRODUCTION_TARGET_SET,
            max_iter_eq=args.max_iter_eq,
            stage="search",
        )
    theta = load_benchmark_theta(args.theta_json)
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    args.outdir.mkdir(parents=True, exist_ok=True)
    cases = [run_case(label, processes[label], theta, args, targets, weights) for label in args.cases]
    summary_rows = [{"case": case["label"], "rank_loss": case["rank_loss"], "market_residual": case["market_residual"], "elapsed_sec": case["elapsed_sec"], **case["headline_statistics"]} for case in cases]
    fit_rows = [{"case": case["label"], **row} for case in cases for row in case["target_fit"]]
    payload = {
        "status": "fixed_theta_income_process_comparison",
        "theta_source_path": str(args.theta_json.resolve()),
        "theta": theta,
        "nb": args.nb,
        "max_iter_eq": args.max_iter_eq,
        "production_target_set": PRODUCTION_TARGET_SET,
        "annual_ar1_specification": {"rho": ANNUAL_RHO, "innovation_sd": ANNUAL_INNOVATION_SD},
        "construction": construction,
        "cases": cases,
    }
    (args.outdir / "income_process_comparison.json").write_text(json.dumps(jsonable(payload), indent=2, sort_keys=True))
    write_csv(args.outdir / "income_process_summary.csv", summary_rows)
    write_csv(args.outdir / "income_process_target_fit.csv", fit_rows)
    print(json.dumps(jsonable({"outdir": str(args.outdir), "summary": summary_rows}), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
