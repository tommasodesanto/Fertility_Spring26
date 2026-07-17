#!/usr/bin/env python3
"""Transformed Nelder--Mead polish for relaxed repaired-objective basins."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any

import numpy as np

from intergen_housing_fertility.calibration import get_target_set
from intergen_housing_fertility.local_panel import (
    append_jsonl,
    income_process_overrides,
    is_better_record,
    jsonable,
    record_selection_loss,
    record_sort_key,
    run_local_panel_case,
    strict_records,
)
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_MAX_ITER_EQ,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
    production_profile_overrides,
    validate_production_profile,
)
from tools.run_intergen_combined_wide_discovery import (
    MATCHED_ANNUAL_INNOVATION_SD,
    MATCHED_ANNUAL_RHO,
    ROOMS_TARGET,
    ROOMS_WEIGHT,
    arm_spec,
    theta_from_unit,
    unit_from_theta,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-theta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--arm",
        choices=("positive_all", "deterministic_tenure", "deterministic_tenure_no_bequest"),
        default=os.getenv("WIDE_POLISH_ARM", "deterministic_tenure"),
    )
    parser.add_argument("--fixed-h-bar-0", type=float, default=None)
    parser.add_argument("--fixed-chi", type=float, default=None)
    parser.add_argument("--seed", type=int, default=int(os.getenv("WIDE_POLISH_SEED", "20260713")))
    parser.add_argument("--max-evals", type=int, default=int(os.getenv("WIDE_POLISH_MAX_EVALS", "2000")))
    parser.add_argument("--minutes", type=float, default=float(os.getenv("WIDE_POLISH_MINUTES", "230")))
    parser.add_argument("--initial-step", type=float, default=float(os.getenv("WIDE_POLISH_INITIAL_STEP", "0.015")))
    parser.add_argument("--start-mix", type=float, default=float(os.getenv("WIDE_POLISH_START_MIX", "0")))
    parser.add_argument("--shrink", type=float, default=float(os.getenv("WIDE_POLISH_SHRINK", "0.5")))
    args = parser.parse_args()

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
    payload = json.loads(args.seed_theta.read_text())
    seed_theta = dict(payload.get("theta", payload))
    active, fixed = arm_spec(args.arm)
    if args.fixed_h_bar_0 is not None:
        if not (0.05 <= args.fixed_h_bar_0 < 6.0):
            raise ValueError("fixed h_bar_0 must be in [0.05, 6)")
        fixed["h_bar_0"] = float(args.fixed_h_bar_0)
    if args.fixed_chi is not None:
        if args.fixed_chi <= 0.0:
            raise ValueError("fixed chi must be positive")
        fixed["chi"] = float(args.fixed_chi)
    fixed_names = set(fixed)
    active = [row for row in active if ("beta" if row[0] == "beta_annual" else row[0]) not in fixed_names]
    x0 = unit_from_theta(seed_theta, active)
    rng = np.random.default_rng(args.seed)
    start_mix = float(np.clip(args.start_mix, 0.0, 0.25))
    if start_mix > 0.0:
        x0 = np.clip((1.0 - start_mix) * x0 + start_mix * rng.random(len(active)), 0.0, 1.0)

    max_evals = max(1, int(args.max_evals))
    initial_step = float(np.clip(args.initial_step, 1e-5, 0.25))
    shrink = float(np.clip(args.shrink, 0.1, 0.9))
    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path = args.outdir / "cases.jsonl"
    best_path = args.outdir / "best.json"
    cases_path.write_text("")

    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = ROOMS_WEIGHT
    income = income_process_overrides(
        5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO
    )
    overrides = production_profile_overrides()
    overrides.update(
        {
            "normalize_bequest_utility": True,
            "q": (1.0 + 0.02) ** 4 - 1.0,
            "delta": 1.0 - (1.0 - 0.011) ** 4,
            "eta_supply": np.array([1.75]),
            "lambda_d": 0.0,
            "debt_taper_start_age": 42.0,
            "debt_taper_end_age": 62.0,
        }
    )
    metadata = {
        "status": "relaxed_transformed_local_polish",
        "algorithm": "bounded_nelder_mead_in_transformed_wide_coordinates",
        "arm": args.arm,
        "seed": args.seed,
        "seed_theta_path": str(args.seed_theta),
        "seed_theta": seed_theta,
        "fixed_theta": fixed,
        "active_dimension": len(active),
        "active_domain": [
            {"name": name, "lower": lo, "upper": hi, "transform": kind}
            for name, lo, hi, kind in active
        ],
        "initial_unit_vector": x0,
        "initial_step": initial_step,
        "start_mix": start_mix,
        "shrink": shrink,
        "max_evals": max_evals,
        "minutes": args.minutes,
        "J": PRODUCTION_J,
        "Nb": PRODUCTION_SEARCH_NB,
        "target_set": PRODUCTION_TARGET_SET,
        "target_count": len(targets),
    }
    (args.outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True)
    )

    start = time.perf_counter()
    deadline = start + max(1.0, float(args.minutes) * 60.0)
    records: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    eval_idx = 0

    def loss(record: dict[str, Any] | None) -> float:
        return record_selection_loss(record) if record is not None else math.inf

    def can_eval() -> bool:
        return eval_idx < max_evals and time.perf_counter() < deadline

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any]) -> dict[str, Any]:
        nonlocal eval_idx, best
        clipped = np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)
        theta = theta_from_unit(clipped, active, fixed)
        record = run_local_panel_case(
            eval_idx,
            {"label": label, "theta": theta},
            PRODUCTION_J,
            PRODUCTION_SEARCH_NB,
            5,
            PRODUCTION_MAX_ITER_EQ,
            income,
            targets,
            weights,
            overrides,
        )
        record["algorithm"] = "wide_transformed_nelder_mead"
        record["polish_arm"] = args.arm
        record["origin"] = jsonable(origin)
        record["unit_vector"] = jsonable(clipped)
        append_jsonl(cases_path, record)
        records.append(record)
        if is_better_record(record, best):
            best = record
            best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
        elapsed = time.perf_counter() - start
        print(
            f"eval {eval_idx + 1}/{max_evals} {label}: status={record.get('status')} "
            f"rank={loss(record):.6g} resid={float(record.get('market_residual', math.inf)):.2e} "
            f"best={loss(best):.6g} elapsed={elapsed / 60:.1f}m",
            flush=True,
        )
        eval_idx += 1
        return record

    dim = len(active)
    simplex: list[tuple[np.ndarray, dict[str, Any]]] = []
    simplex.append((x0.copy(), evaluate(x0, "seed", {"phase": "seed"})))
    for j in range(dim):
        if not can_eval():
            break
        step_j = initial_step * (0.75 + 0.5 * rng.random())
        trial = x0.copy()
        trial[j] = trial[j] + step_j if trial[j] + step_j <= 1.0 else max(0.0, trial[j] - step_j)
        rec = evaluate(
            trial,
            f"nm_init_{j:02d}",
            {"phase": "initial_simplex", "dim": j, "step": step_j},
        )
        simplex.append((trial, rec))

    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, shrink
    iteration = 0
    while can_eval() and len(simplex) >= 2:
        simplex.sort(key=lambda item: loss(item[1]))
        centroid = np.mean([unit for unit, _ in simplex[:-1]], axis=0)
        worst_unit, worst_record = simplex[-1]
        second_worst_loss = loss(simplex[-2][1])
        best_loss = loss(simplex[0][1])
        reflected = np.clip(centroid + alpha * (centroid - worst_unit), 0.0, 1.0)
        reflected_record = evaluate(
            reflected, f"nm_reflect_{iteration:04d}", {"phase": "reflect", "iteration": iteration}
        )
        reflected_loss = loss(reflected_record)
        if reflected_loss < best_loss and can_eval():
            expanded = np.clip(centroid + gamma * (reflected - centroid), 0.0, 1.0)
            expanded_record = evaluate(
                expanded, f"nm_expand_{iteration:04d}", {"phase": "expand", "iteration": iteration}
            )
            simplex[-1] = (expanded, expanded_record) if loss(expanded_record) < reflected_loss else (reflected, reflected_record)
        elif reflected_loss < second_worst_loss:
            simplex[-1] = (reflected, reflected_record)
        else:
            if reflected_loss < loss(worst_record):
                contracted = np.clip(centroid + rho * (reflected - centroid), 0.0, 1.0)
                phase, threshold = "outside_contract", reflected_loss
            else:
                contracted = np.clip(centroid - rho * (centroid - worst_unit), 0.0, 1.0)
                phase, threshold = "inside_contract", loss(worst_record)
            contract_record = (
                evaluate(
                    contracted,
                    f"nm_contract_{iteration:04d}",
                    {"phase": phase, "iteration": iteration},
                )
                if can_eval()
                else worst_record
            )
            if loss(contract_record) < threshold:
                simplex[-1] = (contracted, contract_record)
            else:
                best_unit = simplex[0][0].copy()
                new_simplex = [simplex[0]]
                for slot, (unit, old_record) in enumerate(simplex[1:], start=1):
                    if not can_eval():
                        new_simplex.append((unit, old_record))
                        continue
                    shrunk = np.clip(best_unit + sigma * (unit - best_unit), 0.0, 1.0)
                    shrunk_record = evaluate(
                        shrunk,
                        f"nm_shrink_{iteration:04d}_{slot:02d}",
                        {"phase": "shrink", "iteration": iteration, "slot": slot},
                    )
                    new_simplex.append((shrunk, shrunk_record))
                simplex = new_simplex
        iteration += 1

    ordered = sorted(records, key=record_sort_key)
    strict = strict_records(ordered)
    summary = {
        "best": strict[0] if strict else None,
        "top_records": strict[:10],
        "n_cases_completed": len(records),
        "n_strict": len(strict),
        "n_infeasible": sum(r.get("status") == "infeasible_theta" for r in records),
        "elapsed_sec": time.perf_counter() - start,
        "stopped_by_time_budget": eval_idx < max_evals,
        "iterations_completed": iteration,
        "metadata": metadata,
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(jsonable(summary), indent=2, sort_keys=True)
    )


if __name__ == "__main__":
    main()
