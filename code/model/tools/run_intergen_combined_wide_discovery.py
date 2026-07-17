#!/usr/bin/env python3
"""Wide-domain transformed DE discovery for the repaired combined objective.

This is deliberately a discovery outer loop, not a change to the production
parameter restrictions.  Every proposal is evaluated by the same repaired
Nb=120 model and 15-moment objective used by the production calibration.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any

import numpy as np

from intergen_housing_fertility.local_panel import (
    append_jsonl,
    income_process_overrides,
    is_better_record,
    jsonable,
    latin_hypercube,
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
from intergen_housing_fertility.calibration import get_target_set


MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768
ROOMS_TARGET = 5.779970481941968
ROOMS_WEIGHT = 6.0

# Intentionally much wider than the production restrictions.  The transforms
# allocate resolution near economically important boundaries and avoid invalid
# values (for example beta >= 1, alpha outside (0,1), or nonpositive H0).
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
    ("tenure_choice_kappa", 0.0001, 2.0, "log"),
    ("chi", 0.10, 5.0, "log"),
    ("theta0", 0.00001, 5.0, "log"),
    ("theta_n", -0.49, 8.0, "asinh"),
    ("H0", 0.20, 80.0, "log"),
)


def transform(u: float, lo: float, hi: float, kind: str) -> float:
    u = float(np.clip(u, 0.0, 1.0))
    if kind == "log":
        return float(math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo))))
    if kind == "logit":
        a, b = math.log(lo / (1.0 - lo)), math.log(hi / (1.0 - hi))
        z = a + u * (b - a)
        return float(1.0 / (1.0 + math.exp(-z)))
    if kind == "softzero":
        scale = 0.25
        return float(scale * math.expm1(u * math.log1p(hi / scale)))
    if kind == "asinh":
        a, b = math.asinh(lo), math.asinh(hi)
        return float(math.sinh(a + u * (b - a)))
    if kind == "discount":
        # Log-uniform in the annual discount rate rho=-log(beta_annual).
        rlo, rhi = -math.log(hi), -math.log(lo)
        rho = math.exp(math.log(rlo) + u * (math.log(rhi) - math.log(rlo)))
        return float(math.exp(-rho))
    raise ValueError(f"unknown transform: {kind}")


def inverse(value: float, lo: float, hi: float, kind: str) -> float:
    value = float(np.clip(value, lo, hi))
    if kind == "log":
        return float((math.log(value) - math.log(lo)) / (math.log(hi) - math.log(lo)))
    if kind == "logit":
        z = math.log(value / (1.0 - value))
        a, b = math.log(lo / (1.0 - lo)), math.log(hi / (1.0 - hi))
        return float((z - a) / (b - a))
    if kind == "softzero":
        scale = 0.25
        return float(math.log1p(value / scale) / math.log1p(hi / scale))
    if kind == "asinh":
        a, b = math.asinh(lo), math.asinh(hi)
        return float((math.asinh(value) - a) / (b - a))
    if kind == "discount":
        rho, rlo, rhi = -math.log(value), -math.log(hi), -math.log(lo)
        return float((math.log(rho) - math.log(rlo)) / (math.log(rhi) - math.log(rlo)))
    raise ValueError(f"unknown transform: {kind}")


def arm_spec(arm: str) -> tuple[list[tuple[str, float, float, str]], dict[str, float]]:
    fixed: dict[str, float] = {}
    if arm == "positive_all":
        pass
    elif arm == "deterministic_tenure":
        fixed["tenure_choice_kappa"] = 0.0
    elif arm == "deterministic_tenure_no_bequest":
        fixed.update(tenure_choice_kappa=0.0, theta0=0.0, theta_n=0.0)
    else:
        raise ValueError(f"unknown discovery arm: {arm}")
    active = [row for row in DOMAIN if row[0] not in fixed]
    return active, fixed


def theta_from_unit(unit: np.ndarray, active: list[tuple[str, float, float, str]], fixed: dict[str, float]) -> dict[str, float]:
    theta = dict(fixed)
    for u, (name, lo, hi, kind) in zip(unit, active):
        value = transform(float(u), lo, hi, kind)
        theta["beta" if name == "beta_annual" else name] = value**4 if name == "beta_annual" else value
    return theta


def unit_from_theta(theta: dict[str, Any], active: list[tuple[str, float, float, str]]) -> np.ndarray:
    out = []
    for name, lo, hi, kind in active:
        key = "beta" if name == "beta_annual" else name
        value = float(theta[key]) ** 0.25 if name == "beta_annual" else float(theta[key])
        out.append(inverse(value, lo, hi, kind))
    return np.clip(np.asarray(out, dtype=float), 0.0, 1.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-theta", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--arm", choices=("positive_all", "deterministic_tenure", "deterministic_tenure_no_bequest"), default=os.getenv("WIDE_ARM", "positive_all"))
    parser.add_argument("--seed", type=int, default=int(os.getenv("WIDE_SEED", "20260713")))
    parser.add_argument("--max-evals", type=int, default=int(os.getenv("WIDE_MAX_EVALS", "1800")))
    parser.add_argument("--minutes", type=float, default=float(os.getenv("WIDE_MINUTES", "460")))
    parser.add_argument("--pop-size", type=int, default=int(os.getenv("WIDE_POP_SIZE", "28")))
    parser.add_argument("--mutation", type=float, default=float(os.getenv("WIDE_MUTATION", "0.80")))
    parser.add_argument("--crossover", type=float, default=float(os.getenv("WIDE_CROSSOVER", "0.75")))
    parser.add_argument("--local-seed-count", type=int, default=int(os.getenv("WIDE_LOCAL_SEED_COUNT", "7")))
    parser.add_argument("--local-seed-mix", type=float, default=float(os.getenv("WIDE_LOCAL_SEED_MIX", "0.10")))
    args = parser.parse_args()

    validate_production_profile(PRODUCTION_PROFILE_NAME, J=PRODUCTION_J, Nb=PRODUCTION_SEARCH_NB, n_house=5, income_states=5, target_set=PRODUCTION_TARGET_SET, max_iter_eq=PRODUCTION_MAX_ITER_EQ, stage="search")
    payload = json.loads(args.seed_theta.read_text())
    seed_theta = dict(payload.get("theta", payload))
    active, fixed = arm_spec(args.arm)
    seed_unit = unit_from_theta(seed_theta, active)
    dim = len(active)
    pop_size = max(dim + 1, int(args.pop_size))
    max_evals = max(1, int(args.max_evals))
    rng = np.random.default_rng(args.seed)

    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path, best_path = args.outdir / "cases.jsonl", args.outdir / "best.json"
    cases_path.write_text("")
    targets, weights = get_target_set(PRODUCTION_TARGET_SET)
    targets["aggregate_mean_occupied_rooms_18_85"] = ROOMS_TARGET
    weights["aggregate_mean_occupied_rooms_18_85"] = ROOMS_WEIGHT
    income = income_process_overrides(5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO)
    overrides = production_profile_overrides()
    overrides.update({
        "normalize_bequest_utility": True,
        "q": (1.0 + 0.02) ** 4 - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** 4,
        "eta_supply": np.array([1.75]),
        "lambda_d": 0.0,
        "debt_taper_start_age": 42.0,
        "debt_taper_end_age": 62.0,
    })
    metadata = {
        "status": "wide_transformed_discovery_not_production_restrictions",
        "algorithm": "hybrid_seeded_latin_hypercube_de_rand_1_bin",
        "arm": args.arm,
        "seed": args.seed,
        "J": PRODUCTION_J,
        "Nb": PRODUCTION_SEARCH_NB,
        "target_set": PRODUCTION_TARGET_SET,
        "target_count": len(targets),
        "active_dimension": dim,
        "active_domain": [{"name": n, "lower": lo, "upper": hi, "transform": kind} for n, lo, hi, kind in active],
        "fixed_theta": fixed,
        "seed_theta": seed_theta,
        "seed_unit": seed_unit,
        "max_evals": max_evals,
        "minutes": args.minutes,
        "pop_size": pop_size,
        "mutation": args.mutation,
        "crossover": args.crossover,
        "local_seed_count": args.local_seed_count,
        "local_seed_mix": args.local_seed_mix,
    }
    (args.outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2, sort_keys=True))

    start = time.perf_counter()
    deadline = start + max(1.0, args.minutes * 60.0)
    pop = latin_hypercube(rng, pop_size, dim)
    pop[0] = seed_unit
    local_n = min(max(0, args.local_seed_count), pop_size - 1)
    for i in range(1, local_n + 1):
        mix = float(np.clip(args.local_seed_mix * (0.5 + rng.random()), 0.0, 0.5))
        pop[i] = np.clip((1.0 - mix) * seed_unit + mix * pop[i], 0.0, 1.0)
    pop_loss = np.full(pop_size, math.inf)
    best: dict[str, Any] | None = None
    records: list[dict[str, Any]] = []
    eval_idx = 0
    generation = 0

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any]) -> dict[str, Any]:
        nonlocal eval_idx, best
        clipped = np.clip(np.asarray(unit, dtype=float), 0.0, 1.0)
        theta = theta_from_unit(clipped, active, fixed)
        rec = run_local_panel_case(eval_idx, {"label": label, "theta": theta}, PRODUCTION_J, PRODUCTION_SEARCH_NB, 5, PRODUCTION_MAX_ITER_EQ, income, targets, weights, overrides)
        rec["algorithm"] = "wide_transformed_de"
        rec["discovery_arm"] = args.arm
        rec["origin"] = jsonable(origin)
        rec["unit_vector"] = jsonable(clipped)
        append_jsonl(cases_path, rec)
        records.append(rec)
        if is_better_record(rec, best):
            best = rec
            best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
        elapsed = time.perf_counter() - start
        print(f"eval {eval_idx + 1}/{max_evals} {label}: status={rec.get('status')} rank={record_selection_loss(rec):.5g} resid={float(rec.get('market_residual', math.inf)):.2e} best={record_selection_loss(best) if best else math.inf:.5g} elapsed={elapsed/60:.1f}m", flush=True)
        eval_idx += 1
        return rec

    for i in range(pop_size):
        if eval_idx >= max_evals or time.perf_counter() >= deadline:
            break
        phase = "production_seed" if i == 0 else ("local_seed" if i <= local_n else "latin_hypercube")
        rec = evaluate(pop[i], f"init_{i:03d}", {"phase": phase, "slot": i})
        pop_loss[i] = record_selection_loss(rec)

    while eval_idx < max_evals and time.perf_counter() < deadline:
        order = rng.permutation(pop_size)
        for raw_i in order:
            if eval_idx >= max_evals or time.perf_counter() >= deadline:
                break
            i = int(raw_i)
            choices = [j for j in range(pop_size) if j != i]
            r1, r2, r3 = rng.choice(choices, 3, replace=False)
            mutant = pop[r1] + args.mutation * (pop[r2] - pop[r3])
            bad = (mutant < 0.0) | (mutant > 1.0) | ~np.isfinite(mutant)
            mutant[bad] = rng.random(int(np.sum(bad)))
            mask = rng.random(dim) < args.crossover
            mask[rng.integers(0, dim)] = True
            trial = pop[i].copy()
            trial[mask] = mutant[mask]
            rec = evaluate(trial, f"de_g{generation:04d}_i{i:03d}", {"phase": "de_rand_1_bin", "generation": generation, "slot": i})
            trial_loss = record_selection_loss(rec)
            if trial_loss < pop_loss[i]:
                pop[i], pop_loss[i] = trial, trial_loss
        generation += 1

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
        "generations_completed": generation,
        "metadata": metadata,
    }
    (args.outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))
    if not strict:
        raise RuntimeError("wide discovery produced no strict-converged candidate")


if __name__ == "__main__":
    main()
