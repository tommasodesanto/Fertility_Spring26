#!/usr/bin/env python3
"""Run one recoverable six-hour M5 continuation chain.

This driver belongs to the isolated optimized package.  It preserves the live
M5 model, 15-moment target system, 14-parameter search domain, and external
restriction ``theta_n=0``.  Search evaluations use the established loose
equilibrium tolerance; only two fresh strict repeats can promote a winner.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any, Callable

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from .calibration import extract_moments
from .local_panel import jsonable
from .m5_profile import (
    M5_LOSS,
    M5_OPTIMIZED_STRICT_ROOT_LOSS,
    M5_PROFILE_NAME,
    M5_THETA,
    m5_overrides,
    m5_target_system,
)
from .new_moment_profile import (
    NEW_MOMENT_PROFILE_NAME,
    NEW_MOMENT_PROFILE_RUNNABLE,
    NEW_MOMENT_SEED,
    new_moment_overrides,
    new_moment_target_system,
)
from .promotion_contract import FREE_PARAMETER_BOUNDS
from .solver import InfeasibleThetaError, run_model_cp_dt


PERIOD_YEARS = 4.0
STRICT_RESIDUAL = 2.5e-5
SEARCH_RESIDUAL = 1.0e-4
DOMAIN: tuple[tuple[str, float, float, str], ...] = (
    ("beta_annual", 0.80, 0.9995, "discount"),
    ("alpha_cons", 0.02, 0.98, "logit"),
    ("c_bar_0", 0.0, 2.0, "softzero"),
    ("c_bar_n", 0.0, 3.0, "softzero"),
    ("h_bar_0", 0.05, 5.80, "log"),
    ("h_bar_jump", 0.0, 8.0, "softzero"),
    ("h_bar_n", 0.0, 5.0, "softzero"),
    ("psi_child", -3.0, 3.0, "asinh"),
    ("kappa_fert", 0.02, 50.0, "log"),
    ("chi", 0.10, 5.0, "log"),
    ("H0", 0.20, 80.0, "log"),
    ("theta0", 0.0, 8.0, "softzero"),
    ("theta1", 0.02, 16.0, "log"),
    ("tenure_choice_kappa", 0.0, 0.12, "softzero"),
)


def validate_contract(profile: str = "m5") -> None:
    if profile == "new-moments" and not NEW_MOMENT_PROFILE_RUNNABLE:
        raise RuntimeError(
            "new-moments is disabled until its empirical auxiliaries and "
            "identification gates pass"
        )
    declared = tuple((name, lower, upper) for name, lower, upper, _ in DOMAIN)
    if declared != FREE_PARAMETER_BOUNDS:
        raise RuntimeError("calibration search domain differs from the promotion contract")
    target_system = m5_target_system() if profile == "m5" else new_moment_target_system()
    target_system.require_identified(free_parameter_count=len(DOMAIN))
    expected_count = 15 if profile == "m5" else 14
    if target_system.count != expected_count:
        raise RuntimeError(f"{profile} requires its declared {expected_count}-moment target system")
    overrides = m5_overrides if profile == "m5" else new_moment_overrides
    if not math.isclose(float(overrides(tight=False)["tol_eq"]), SEARCH_RESIDUAL):
        raise RuntimeError(f"{profile} loose evaluator tolerance drifted")
    if not math.isclose(float(overrides(tight=True)["tol_eq"]), STRICT_RESIDUAL):
        raise RuntimeError(f"{profile} strict evaluator tolerance drifted")


def transform(unit: float, lower: float, upper: float, kind: str) -> float:
    unit = float(np.clip(unit, 0.0, 1.0))
    if kind == "log":
        return lower * (upper / lower) ** unit
    if kind == "discount":
        return lower + (upper - lower) * (1.0 - (1.0 - unit) ** 2)
    if kind == "logit":
        return lower + (upper - lower) / (1.0 + math.exp(-8.0 * (unit - 0.5)))
    if kind == "asinh":
        return math.sinh((1.0 - unit) * math.asinh(lower) + unit * math.asinh(upper))
    if kind == "softzero":
        return lower + (upper - lower) * unit * unit
    raise ValueError(f"unknown transform {kind!r}")


def inverse(value: float, lower: float, upper: float, kind: str) -> float:
    value = float(np.clip(value, lower, upper))
    if kind == "log":
        return math.log(value / lower) / math.log(upper / lower)
    if kind == "discount":
        return 1.0 - math.sqrt(max(0.0, 1.0 - (value - lower) / (upper - lower)))
    if kind == "logit":
        share = float(np.clip((value - lower) / (upper - lower), 1e-12, 1.0 - 1e-12))
        return 0.5 + math.log(share / (1.0 - share)) / 8.0
    if kind == "asinh":
        return (math.asinh(value) - math.asinh(lower)) / (
            math.asinh(upper) - math.asinh(lower)
        )
    if kind == "softzero":
        return math.sqrt(max(0.0, (value - lower) / (upper - lower)))
    raise ValueError(f"unknown transform {kind!r}")


def theta_from_unit(unit: np.ndarray) -> dict[str, float]:
    theta = {"theta_n": 0.0}
    for coordinate, (name, lower, upper, kind) in zip(np.asarray(unit), DOMAIN):
        value = transform(float(coordinate), lower, upper, kind)
        theta["beta" if name == "beta_annual" else name] = (
            value**PERIOD_YEARS if name == "beta_annual" else value
        )
    return theta


def unit_from_theta(theta: dict[str, float]) -> np.ndarray:
    unit = []
    for name, lower, upper, kind in DOMAIN:
        value = float(theta["beta" if name == "beta_annual" else name])
        if name == "beta_annual":
            value **= 1.0 / PERIOD_YEARS
        unit.append(inverse(value, lower, upper, kind))
    return np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)


def parameter_rows(theta: dict[str, float]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, lower, upper, kind in DOMAIN:
        key = "beta" if name == "beta_annual" else name
        estimate = float(theta[key])
        if name == "beta_annual":
            estimate **= 1.0 / PERIOD_YEARS
        rows.append(
            {
                "parameter": name,
                "estimate": estimate,
                "lower": lower,
                "upper": upper,
                "transform": kind,
                "role": "estimated",
                "near_bound": min(estimate - lower, upper - estimate)
                <= 0.02 * (upper - lower),
            }
        )
    rows.append(
        {
            "parameter": "theta_n",
            "estimate": 0.0,
            "lower": None,
            "upper": None,
            "transform": "external restriction",
            "role": "fixed",
            "near_bound": False,
        }
    )
    return rows


def target_fit_rows(
    moments: dict[str, float], system=None
) -> list[dict[str, float | str]]:
    system = m5_target_system() if system is None else system
    return [
        {
            "moment": name,
            "target": target,
            "model": float(moments[name]),
            "gap": float(moments[name]) - target,
            "weight": weight,
            "loss_contribution": weight * (float(moments[name]) - target) ** 2,
        }
        for name, target, weight in zip(
            system.moment_names,
            system.target_values,
            system.weights,
        )
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--profile", choices=("m5", "new-moments"), default="m5")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--method", choices=("nelder-mead", "pattern"), required=True)
    parser.add_argument("--minutes", type=float, default=350.0)
    parser.add_argument("--strict-reserve-minutes", type=float, default=10.0)
    parser.add_argument("--max-evals", type=int, default=2000)
    parser.add_argument("--initial-step", type=float, default=0.015)
    parser.add_argument("--minimum-step", type=float, default=2.5e-4)
    parser.add_argument("--start-mix", type=float, default=0.0)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validate_contract(args.profile)
    if args.smoke:
        args.minutes = min(float(args.minutes), 8.0)
        args.strict_reserve_minutes = 0.0
        args.max_evals = min(int(args.max_evals), len(DOMAIN) + 1)
    if args.minutes <= args.strict_reserve_minutes:
        raise ValueError("total minutes must exceed the strict reserve")

    if args.profile == "m5":
        target_system = m5_target_system()
        profile_name = M5_PROFILE_NAME
        seed_theta = M5_THETA
        override_factory = m5_overrides
        run_status = "proper_joint_m5_continuation"
    else:
        target_system = new_moment_target_system()
        profile_name = NEW_MOMENT_PROFILE_NAME
        seed_theta = NEW_MOMENT_SEED
        override_factory = new_moment_overrides
        run_status = "provisional_joint_new_moment_calibration"
    rng = np.random.default_rng(int(args.seed))
    seed_unit = unit_from_theta(seed_theta)
    start_mix = float(np.clip(args.start_mix, 0.0, 0.25))
    start_unit = np.clip(
        (1.0 - start_mix) * seed_unit + start_mix * rng.random(len(DOMAIN)),
        0.0,
        1.0,
    )
    search_overrides = override_factory(tight=False, optimized=True)
    if args.smoke:
        search_overrides.update(Nb=60, max_iter_eq=2, tol_eq=0.25)

    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path = args.outdir / "cases.jsonl"
    cases_path.write_text("")
    metadata = {
        "status": "exact_loop_smoke" if args.smoke else run_status,
        "profile": profile_name,
        "target_system": target_system.name,
        "target_fingerprint": target_system.fingerprint,
        "target_count": target_system.count,
        "free_parameter_count": len(DOMAIN),
        "external_restrictions": {"theta_n": 0.0},
        "domain": [
            {"name": name, "lower": lower, "upper": upper, "transform": kind}
            for name, lower, upper, kind in DOMAIN
        ],
        "canonical_m5_loss": M5_LOSS,
        "optimized_strict_root_m5_loss": M5_OPTIMIZED_STRICT_ROOT_LOSS,
        "objective_note": (
            "legacy M5 weighting"
            if args.profile == "m5"
            else "sum of squared relative gaps; provisional until joint auxiliary covariance is available"
        ),
        "method": args.method,
        "seed": int(args.seed),
        "start_mix": start_mix,
        "initial_step": float(args.initial_step),
        "minimum_step": float(args.minimum_step),
        "minutes": float(args.minutes),
        "strict_reserve_minutes": float(args.strict_reserve_minutes),
        "max_evals": int(args.max_evals),
        "search_evaluator": {
            "J": int(search_overrides["J"]),
            "Nb": int(search_overrides["Nb"]),
            "max_iter_eq": int(search_overrides["max_iter_eq"]),
            "tol_eq": float(search_overrides["tol_eq"]),
            "optimized": True,
        },
        "strict_evaluator": {"J": 17, "Nb": 120, "max_iter_eq": 40, "tol_eq": 2.5e-5},
    }
    (args.outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n"
    )

    started = time.perf_counter()
    search_seconds = 60.0 * (float(args.minutes) - float(args.strict_reserve_minutes))
    deadline = started + max(1.0, search_seconds)
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    eval_index = 0

    def score(record: dict[str, Any] | None) -> float:
        if record is None or not bool(record.get("strict_converged", False)):
            return math.inf
        return float(record.get("rank_loss", math.inf))

    def can_evaluate() -> bool:
        return eval_index < int(args.max_evals) and time.perf_counter() < deadline

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any]) -> dict[str, Any]:
        nonlocal best, eval_index
        clipped = np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)
        theta = theta_from_unit(clipped)
        began = time.perf_counter()
        try:
            solution, parameters, price = run_model_cp_dt(
                {**search_overrides, **theta}, verbose=False
            )
            moments = extract_moments(solution, parameters)
            loss = float(target_system.loss(moments))
            residual = float(getattr(solution, "best_max_abs_rel_excess", math.inf))
            timings = dict(getattr(solution, "timings", {}))
            tolerance = float(search_overrides["tol_eq"])
            converged = bool(
                timings.get("strict_converged", getattr(solution, "converged", False))
                and math.isfinite(residual)
                and residual <= tolerance
            )
            status, error, census = "ok", "", []
            price_value = float(np.asarray(price).reshape(-1)[0])
        except InfeasibleThetaError as exc:
            moments, timings = {}, {}
            loss, residual, price_value, converged = math.inf, math.inf, math.nan, False
            status, error, census = "infeasible_theta", str(exc), list(exc.census)
        except Exception as exc:  # Persist unexpected failures instead of losing a chain.
            moments, timings = {}, {}
            loss, residual, price_value, converged = math.inf, math.inf, math.nan, False
            status, error, census = f"failed:{type(exc).__name__}", str(exc), []
        record = {
            "case": eval_index,
            "label": label,
            "status": status,
            "strict_converged": converged,
            "rank_loss": loss,
            "market_residual": residual,
            "price": price_value,
            "theta": theta,
            "moments": moments,
            "unit_vector": clipped,
            "origin": origin,
            "timings": timings,
            "feasibility_census": census,
            "error": error,
            "elapsed_sec": time.perf_counter() - began,
        }
        with cases_path.open("a") as handle:
            handle.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
        records.append(record)
        (args.outdir / "latest_completed_case.json").write_text(
            json.dumps(jsonable(record), indent=2, sort_keys=True) + "\n"
        )
        if score(record) < score(best):
            best = record
            (args.outdir / "best_so_far.json").write_text(
                json.dumps(jsonable(best), indent=2, sort_keys=True) + "\n"
            )
        eval_index += 1
        print(
            f"eval={eval_index}/{args.max_evals} method={args.method} label={label} "
            f"status={status} loss={score(record):.9g} residual={residual:.3e} "
            f"best={score(best):.9g} elapsed_min={(time.perf_counter()-started)/60.0:.2f}",
            flush=True,
        )
        return record

    if args.method == "nelder-mead":
        run_nelder_mead(start_unit, evaluate, can_evaluate, score, rng, float(args.initial_step))
    else:
        run_pattern_search(
            start_unit,
            evaluate,
            can_evaluate,
            score,
            rng,
            float(args.initial_step),
            float(args.minimum_step),
        )

    tight_records: list[dict[str, Any]] = []
    if not args.smoke and best is not None:
        tight_overrides = override_factory(tight=True, optimized=True)
        winning_unit = np.asarray(best["unit_vector"], dtype=float)
        winning_theta = theta_from_unit(winning_unit)
        tight_path = args.outdir / "tight_cases.jsonl"
        tight_path.write_text("")
        for repeat in range(2):
            began = time.perf_counter()
            try:
                solution, parameters, price = run_model_cp_dt(
                    {**tight_overrides, **winning_theta}, verbose=False
                )
                moments = extract_moments(solution, parameters)
                loss = float(target_system.loss(moments))
                residual = float(getattr(solution, "best_max_abs_rel_excess", math.inf))
                timings = dict(getattr(solution, "timings", {}))
                converged = bool(
                    timings.get("strict_converged", getattr(solution, "converged", False))
                    and math.isfinite(residual)
                    and residual <= STRICT_RESIDUAL
                )
                status, error, census = "ok", "", []
                price_value = float(np.asarray(price).reshape(-1)[0])
            except InfeasibleThetaError as exc:
                moments, timings = {}, {}
                loss, residual, price_value, converged = math.inf, math.inf, math.nan, False
                status, error, census = "infeasible_theta", str(exc), list(exc.census)
            except Exception as exc:
                moments, timings = {}, {}
                loss, residual, price_value, converged = math.inf, math.inf, math.nan, False
                status, error, census = f"failed:{type(exc).__name__}", str(exc), []
            tight = {
                "case": repeat,
                "label": f"strict_winner_repeat_{repeat + 1}",
                "status": status,
                "strict_converged": converged,
                "rank_loss": loss,
                "market_residual": residual,
                "price": price_value,
                "theta": winning_theta,
                "moments": moments,
                "target_fit": target_fit_rows(moments, target_system) if moments else [],
                "parameters": parameter_rows(winning_theta),
                "unit_vector": winning_unit,
                "origin": {"search_case": int(best["case"])},
                "timings": timings,
                "feasibility_census": census,
                "error": error,
                "evaluator": {"max_iter_eq": 40, "tol_eq": STRICT_RESIDUAL},
                "elapsed_sec": time.perf_counter() - began,
            }
            with tight_path.open("a") as handle:
                handle.write(json.dumps(jsonable(tight), sort_keys=True) + "\n")
            tight_records.append(tight)
            (args.outdir / "latest_tight_case.json").write_text(
                json.dumps(jsonable(tight), indent=2, sort_keys=True) + "\n"
            )
            print(
                f"strict_repeat={repeat + 1}/2 status={status} strict={converged} "
                f"loss={loss:.9g} residual={residual:.3e}",
                flush=True,
            )

    eligible = bool(
        len(tight_records) == 2
        and all(record.get("strict_converged") for record in tight_records)
        and tight_records[0].get("rank_loss") == tight_records[1].get("rank_loss")
        and tight_records[0].get("price") == tight_records[1].get("price")
        and tight_records[0].get("moments") == tight_records[1].get("moments")
    )
    if eligible:
        (args.outdir / "best_tight.json").write_text(
            json.dumps(jsonable(tight_records[0]), indent=2, sort_keys=True) + "\n"
        )
    summary = {
        "status": (
            "smoke_complete"
            if args.smoke
            else "evaluation_cap_reached"
            if eval_index >= int(args.max_evals)
            else "time_budget_or_step_stop"
        ),
        "eligible": eligible,
        "best_search": best,
        "best_tight": tight_records[0] if eligible else None,
        "tight_repeats": tight_records,
        "n_cases_completed": len(records),
        "n_converged": sum(bool(record.get("strict_converged")) for record in records),
        "n_infeasible": sum(record.get("status") == "infeasible_theta" for record in records),
        "n_failed": sum(str(record.get("status", "")).startswith("failed:") for record in records),
        "elapsed_sec": time.perf_counter() - started,
        "metadata": metadata,
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True) + "\n"
    )
    if args.smoke and (
        len(records) != len(DOMAIN) + 1
        or summary["n_converged"] != len(records)
        or summary["n_infeasible"] != 0
        or summary["n_failed"] != 0
    ):
        raise RuntimeError(
            "exact-loop smoke requires 15/15 converged cases with no infeasibility or program errors"
        )


def run_nelder_mead(
    start: np.ndarray,
    evaluate: Callable[[np.ndarray, str, dict[str, Any]], dict[str, Any]],
    can_evaluate: Callable[[], bool],
    score: Callable[[dict[str, Any] | None], float],
    rng: np.random.Generator,
    initial_step: float,
) -> None:
    simplex = [(start.copy(), evaluate(start, "seed", {"phase": "seed"}))]
    for dimension in range(len(DOMAIN)):
        if not can_evaluate():
            break
        step = initial_step * (0.75 + 0.5 * rng.random())
        trial = start.copy()
        trial[dimension] = (
            trial[dimension] + step
            if trial[dimension] + step <= 1.0
            else max(0.0, trial[dimension] - step)
        )
        simplex.append(
            (
                trial,
                evaluate(
                    trial,
                    f"nm_init_{dimension:02d}",
                    {"phase": "initial_simplex", "dimension": dimension, "step": step},
                ),
            )
        )
    iteration = 0
    while can_evaluate() and len(simplex) == len(DOMAIN) + 1:
        simplex.sort(key=lambda item: score(item[1]))
        centroid = np.mean([unit for unit, _ in simplex[:-1]], axis=0)
        worst_unit, worst_record = simplex[-1]
        reflected = np.clip(centroid + (centroid - worst_unit), 0.0, 1.0)
        reflected_record = evaluate(
            reflected, f"nm_reflect_{iteration:05d}", {"phase": "reflect", "iteration": iteration}
        )
        reflected_loss = score(reflected_record)
        best_loss = score(simplex[0][1])
        second_worst_loss = score(simplex[-2][1])
        if reflected_loss < best_loss and can_evaluate():
            expanded = np.clip(centroid + 2.0 * (reflected - centroid), 0.0, 1.0)
            expanded_record = evaluate(
                expanded, f"nm_expand_{iteration:05d}", {"phase": "expand", "iteration": iteration}
            )
            simplex[-1] = (
                (expanded, expanded_record)
                if score(expanded_record) < reflected_loss
                else (reflected, reflected_record)
            )
        elif reflected_loss < second_worst_loss:
            simplex[-1] = (reflected, reflected_record)
        else:
            if reflected_loss < score(worst_record):
                contracted = np.clip(centroid + 0.5 * (reflected - centroid), 0.0, 1.0)
                threshold, phase = reflected_loss, "outside_contract"
            else:
                contracted = np.clip(centroid - 0.5 * (centroid - worst_unit), 0.0, 1.0)
                threshold, phase = score(worst_record), "inside_contract"
            contracted_record = (
                evaluate(
                    contracted,
                    f"nm_contract_{iteration:05d}",
                    {"phase": phase, "iteration": iteration},
                )
                if can_evaluate()
                else worst_record
            )
            if score(contracted_record) < threshold:
                simplex[-1] = (contracted, contracted_record)
            else:
                best_unit = simplex[0][0].copy()
                reduced = [simplex[0]]
                for slot, (unit, old_record) in enumerate(simplex[1:], start=1):
                    if not can_evaluate():
                        reduced.append((unit, old_record))
                        continue
                    shrunk = np.clip(best_unit + 0.5 * (unit - best_unit), 0.0, 1.0)
                    reduced.append(
                        (
                            shrunk,
                            evaluate(
                                shrunk,
                                f"nm_shrink_{iteration:05d}_{slot:02d}",
                                {"phase": "shrink", "iteration": iteration, "slot": slot},
                            ),
                        )
                    )
                simplex = reduced
        iteration += 1


def run_pattern_search(
    start: np.ndarray,
    evaluate: Callable[[np.ndarray, str, dict[str, Any]], dict[str, Any]],
    can_evaluate: Callable[[], bool],
    score: Callable[[dict[str, Any] | None], float],
    rng: np.random.Generator,
    initial_step: float,
    minimum_step: float,
) -> None:
    incumbent_unit = start.copy()
    incumbent_record = evaluate(start, "seed", {"phase": "seed"})
    step = float(initial_step)
    sweep = 0
    while can_evaluate() and step >= minimum_step:
        improved = False
        for dimension in rng.permutation(len(DOMAIN)):
            for direction in rng.permutation(np.asarray((-1.0, 1.0))):
                if not can_evaluate():
                    break
                trial = incumbent_unit.copy()
                trial[dimension] = np.clip(trial[dimension] + direction * step, 0.0, 1.0)
                if trial[dimension] == incumbent_unit[dimension]:
                    continue
                record = evaluate(
                    trial,
                    f"pattern_{sweep:04d}_{dimension:02d}_{int(direction):+d}",
                    {
                        "phase": "coordinate_poll",
                        "sweep": sweep,
                        "dimension": int(dimension),
                        "direction": float(direction),
                        "step": step,
                    },
                )
                if score(record) < score(incumbent_record):
                    incumbent_unit, incumbent_record = trial, record
                    improved = True
                    break
        if improved:
            step = min(initial_step, step * 1.15)
        else:
            step *= 0.5
        sweep += 1


if __name__ == "__main__":
    main()
