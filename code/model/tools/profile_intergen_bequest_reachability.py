#!/usr/bin/env python3
"""Profile the attainable late-life estate dispersion over a wide theta1 grid.

The 11 non-bequest coordinates are fixed at a strict M3 winner.  For each
fixed theta1 cell, theta0 and theta_n are re-optimized against the late-life
median estate ratio and the 2+-minus-1-child median gap.  The p90/median ratio
is an outcome, not part of the cell objective, so the completed grid is a
direct reachability frontier for the dispersion target.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from intergen_housing_fertility.calibration import extract_moments
from intergen_housing_fertility.local_panel import jsonable
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    validate_production_profile,
)
from intergen_housing_fertility.solver import InfeasibleThetaError, run_model_cp_dt
from tools.run_intergen_bequest_exit_chain import (
    THETA0_DOMAIN,
    THETA1_DOMAIN,
    THETA_N_DOMAIN,
    arm_contract,
    common_overrides,
    load_theta,
    target_system,
)
from tools.run_intergen_combined_wide_discovery import inverse, transform


THETA1_GRID = (0.02, 0.05, 0.10, 0.20, 0.35, 0.55, 0.80, 1.20, 2.0, 4.0, 8.0, 16.0)
RETIREMENT_INCOME_SCALE_GRID = (0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0)
MEDIAN = "old_total_estate_wealth_to_annual_income_median_7684"
DISPERSION = "old_total_estate_wealth_to_annual_income_p90_p50_7684"
FAMILY_GAP = "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575"
ACTIVE = (THETA0_DOMAIN, THETA_N_DOMAIN)
JOINT_ACTIVE = (THETA0_DOMAIN, THETA1_DOMAIN, THETA_N_DOMAIN)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--winner-json", type=Path, required=True)
    parser.add_argument("--winner-arm", type=str, default="M3")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cell", type=int, required=True)
    parser.add_argument(
        "--experiment",
        choices=("theta1_frontier", "retirement_income_dispersion"),
        default="theta1_frontier",
    )
    parser.add_argument("--seed", type=int, default=2026071600)
    parser.add_argument("--max-evals", type=int, default=150)
    parser.add_argument("--minutes", type=float, default=75.0)
    parser.add_argument("--initial-step", type=float, default=0.12)
    parser.add_argument("--min-step", type=float, default=0.002)
    parser.add_argument("--shrink", type=float, default=0.5)
    parser.add_argument("--J", type=int, default=PRODUCTION_J)
    parser.add_argument("--Nb", type=int, default=PRODUCTION_SEARCH_NB)
    parser.add_argument("--max-iter-eq", type=int, default=PRODUCTION_MAX_ITER_EQ)
    parser.add_argument("--tol-eq", type=float, default=1e-4)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def unit_from_theta(
    theta: dict[str, float],
    active: tuple[tuple[str, float, float, str], ...] = ACTIVE,
) -> np.ndarray:
    return np.asarray(
        [inverse(float(theta[name]), lo, hi, kind) for name, lo, hi, kind in active],
        dtype=float,
    )


def bequest_from_unit(
    unit: np.ndarray,
    active: tuple[tuple[str, float, float, str], ...] = ACTIVE,
) -> dict[str, float]:
    return {
        name: transform(float(u), lo, hi, kind)
        for u, (name, lo, hi, kind) in zip(np.asarray(unit, dtype=float), active)
    }


def main() -> None:
    args = parse_args()
    grid = THETA1_GRID if args.experiment == "theta1_frontier" else RETIREMENT_INCOME_SCALE_GRID
    if not 1 <= int(args.cell) <= len(grid):
        raise ValueError(f"cell must lie in [1,{len(grid)}] for {args.experiment}")
    if args.smoke:
        args.Nb = 60
        args.max_iter_eq = 2
        args.tol_eq = 0.25
        args.max_evals = min(int(args.max_evals), 6)
        args.minutes = min(float(args.minutes), 8.0)
    else:
        validate_production_profile(
            PRODUCTION_PROFILE_NAME,
            J=args.J,
            Nb=args.Nb,
            n_house=5,
            income_states=5,
            target_set=PRODUCTION_TARGET_SET,
            max_iter_eq=args.max_iter_eq,
            stage="search",
        )
        if not math.isclose(float(args.tol_eq), 1e-4):
            raise ValueError("production reachability cells require tol_eq=1e-4")

    is_retirement_diagnostic = args.experiment == "retirement_income_dispersion"
    cell_value = float(grid[int(args.cell) - 1])
    seed_theta = load_theta(args.winner_json, args.winner_arm)
    base_theta = dict(seed_theta)
    theta1 = float(seed_theta["theta1"]) if is_retirement_diagnostic else cell_value
    retirement_income_scale = cell_value if is_retirement_diagnostic else 0.0
    if not is_retirement_diagnostic:
        base_theta["theta1"] = theta1
    base_theta["tenure_choice_kappa"] = 0.0
    active = JOINT_ACTIVE if is_retirement_diagnostic else ACTIVE
    objective_names = (MEDIAN, DISPERSION, FAMILY_GAP) if is_retirement_diagnostic else (MEDIAN, FAMILY_GAP)
    x_best = unit_from_theta(seed_theta, active)

    contract = SimpleNamespace(
        arm="M3",
        ltv_terminal=0.0,
        theta1=theta1,
        fixed_theta0=None,
        J=int(args.J),
        Nb=int(args.Nb),
        max_iter_eq=int(args.max_iter_eq),
        tol_eq=float(args.tol_eq),
    )
    _, _, mechanism = arm_contract(contract)
    overrides = common_overrides(contract, mechanism)
    overrides["retirement_income_z_scale"] = retirement_income_scale
    targets, weights = target_system("candidate_replacement_bequest_internal_v1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path = args.outdir / "cases.jsonl"
    cases_path.write_text("")
    started = time.perf_counter()
    tight_reserve = 0.0 if args.smoke else 180.0
    deadline = started + max(1.0, float(args.minutes) * 60.0 - tight_reserve)
    max_evals = max(1, int(args.max_evals))
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None

    metadata = {
        "status": "smoke" if args.smoke else f"bequest_{args.experiment}_cell",
        "experiment": args.experiment,
        "cell": int(args.cell),
        "theta1": theta1,
        "theta1_grid": list(THETA1_GRID) if not is_retirement_diagnostic else None,
        "retirement_income_z_scale": retirement_income_scale,
        "retirement_income_scale_grid": (
            list(RETIREMENT_INCOME_SCALE_GRID) if is_retirement_diagnostic else None
        ),
        "mechanism_interpretation": (
            "persistent pension/retirement-income heterogeneity linked to the existing income state; "
            "the weighted mean pension is unchanged"
            if is_retirement_diagnostic
            else None
        ),
        "fixed_nonbequest_theta": {
            key: value
            for key, value in base_theta.items()
            if key not in {"theta0", "theta1", "theta_n"}
        },
        "active_bounds": [
            {"name": name, "lower": lo, "upper": hi, "transform": kind}
            for name, lo, hi, kind in active
        ],
        "cell_objective": list(objective_names),
        "reported_frontier_outcome": None if is_retirement_diagnostic else DISPERSION,
        "targets": {name: targets[name] for name in (MEDIAN, DISPERSION, FAMILY_GAP)},
        "weights": {name: weights[name] for name in objective_names},
        "winner_json": str(args.winner_json),
        "winner_arm": args.winner_arm,
        "J": int(args.J),
        "Nb": int(args.Nb),
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "max_evals": max_evals,
        "minutes": float(args.minutes),
        "tight_evaluator": {"max_iter_eq": 40, "tol_eq": 2.5e-5, "repeats": 2},
    }
    (args.outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2, sort_keys=True))

    def score(record: dict[str, Any] | None) -> float:
        if record is None or not record.get("strict_converged", False):
            return math.inf
        return float(record["profile_loss"])

    def can_eval() -> bool:
        return len(records) < max_evals and time.perf_counter() < deadline

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any], *, tight: bool = False) -> dict[str, Any]:
        nonlocal best
        clipped = np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)
        theta = {**base_theta, **bequest_from_unit(clipped, active)}
        if not is_retirement_diagnostic:
            theta["theta1"] = theta1
        run_overrides = dict(overrides)
        if tight:
            run_overrides.update(max_iter_eq=40, tol_eq=2.5e-5)
        t0 = time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt({**run_overrides, **theta}, verbose=False)
            moments = extract_moments(sol, P)
            residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
            timings = dict(getattr(sol, "timings", {}))
            strict = bool(
                timings.get("strict_converged", getattr(sol, "converged", False))
                and math.isfinite(residual)
                and residual <= float(P.tol_eq)
            )
            gaps = {name: float(moments[name]) - float(targets[name]) for name in (MEDIAN, DISPERSION, FAMILY_GAP)}
            loss = sum(float(weights[name]) * gaps[name] ** 2 for name in objective_names)
            status, error = "ok", ""
            price = float(np.asarray(p_eq).reshape(-1)[0])
        except InfeasibleThetaError as exc:
            moments, gaps, timings = {}, {}, {}
            residual, price, strict, loss = math.inf, math.nan, False, math.inf
            status, error = "infeasible_theta", str(exc)
        except Exception as exc:  # noqa: BLE001 - persist diagnostic failures.
            moments, gaps, timings = {}, {}, {}
            residual, price, strict, loss = math.inf, math.nan, False, math.inf
            status, error = f"failed:{type(exc).__name__}", str(exc)
        record = {
            "case": len(records),
            "label": label,
            "status": status,
            "strict_converged": strict,
            "profile_loss": loss,
            "market_residual": residual,
            "price": price,
            "theta": theta,
            "retirement_income_z_scale": retirement_income_scale,
            "moments": {name: moments.get(name) for name in (MEDIAN, DISPERSION, FAMILY_GAP)},
            "all_moments": moments,
            "gaps": gaps,
            "unit_vector": clipped,
            "origin": origin,
            "timings": timings,
            "error": error,
            "evaluator": {
                "max_iter_eq": int(run_overrides.get("max_iter_eq", args.max_iter_eq)),
                "tol_eq": float(run_overrides.get("tol_eq", args.tol_eq)),
            },
            "elapsed_sec": time.perf_counter() - t0,
        }
        if not tight:
            with cases_path.open("a") as handle:
                handle.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
            records.append(record)
            (args.outdir / "latest_completed_case.json").write_text(
                json.dumps(jsonable(record), indent=2, sort_keys=True)
            )
            if score(record) < score(best):
                best = record
                (args.outdir / "best.json").write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
            print(
                f"cell={args.cell} theta1={theta1:g} ret_scale={retirement_income_scale:g} "
                f"eval={len(records)}/{max_evals} "
                f"{label} loss={score(record):.6g} resid={residual:.2e} "
                f"best={score(best):.6g} elapsed={(time.perf_counter()-started)/60:.1f}m",
                flush=True,
            )
        return record

    best_record = evaluate(x_best, "seed", {"phase": "seed"})
    step = min(max(float(args.initial_step), 1e-5), 0.5)
    min_step = min(max(float(args.min_step), 1e-6), step)
    shrink = min(max(float(args.shrink), 0.1), 0.9)
    iteration = 0
    while can_eval() and step >= min_step:
        improved = False
        for dim in range(len(active)):
            if not can_eval():
                break
            incumbent = (score(best_record), x_best.copy(), best_record)
            candidates = [incumbent]
            for direction in (-1.0, 1.0):
                if not can_eval():
                    break
                trial = x_best.copy()
                trial[dim] = np.clip(trial[dim] + direction * step, 0.0, 1.0)
                rec = evaluate(
                    trial,
                    f"pattern_i{iteration:03d}_d{dim}_{'m' if direction < 0 else 'p'}",
                    {"phase": "coordinate", "iteration": iteration, "dim": dim, "direction": direction, "step": step},
                )
                candidates.append((score(rec), trial, rec))
            candidates.sort(key=lambda item: item[0])
            if candidates[0][0] + 1e-10 < score(best_record):
                _, x_best, best_record = candidates[0]
                x_best = x_best.copy()
                improved = True
        if not improved:
            step *= shrink
        iteration += 1

    tight_records: list[dict[str, Any]] = []
    if not args.smoke and best is not None:
        x_tight = np.asarray(best["unit_vector"], dtype=float)
        for repeat in range(2):
            tight = evaluate(
                x_tight,
                f"tight_repeat_{repeat+1}",
                {"phase": "tight_repeat", "repeat": repeat + 1},
                tight=True,
            )
            tight_records.append(tight)
            print(
                f"cell={args.cell} tight={repeat+1}/2 strict={tight['strict_converged']} "
                f"loss={tight['profile_loss']:.6g} resid={tight['market_residual']:.2e}",
                flush=True,
            )
        (args.outdir / "tight_cases.jsonl").write_text(
            "".join(json.dumps(jsonable(row), sort_keys=True) + "\n" for row in tight_records)
        )

    repeat_check = None
    best_tight = None
    if len(tight_records) == 2:
        names = (MEDIAN, DISPERSION, FAMILY_GAP)
        repeat_check = {
            "both_strict": all(row["strict_converged"] for row in tight_records),
            "loss_abs_difference": abs(tight_records[0]["profile_loss"] - tight_records[1]["profile_loss"]),
            "max_abs_moment_difference": max(
                abs(float(tight_records[0]["moments"][name]) - float(tight_records[1]["moments"][name]))
                for name in names
            ),
        }
        if repeat_check["both_strict"]:
            best_tight = min(tight_records, key=lambda row: row["profile_loss"])
            (args.outdir / "best_tight.json").write_text(
                json.dumps(jsonable(best_tight), indent=2, sort_keys=True)
            )

    summary = {
        "status": "complete" if len(records) >= max_evals else "time_or_step_budget_reached",
        "metadata": metadata,
        "best_search": best,
        "best_tight": best_tight,
        "tight_repeat_check": repeat_check,
        "n_cases_completed": len(records),
        "elapsed_sec": time.perf_counter() - started,
        "final_step": step,
        "iterations": iteration,
    }
    (args.outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
