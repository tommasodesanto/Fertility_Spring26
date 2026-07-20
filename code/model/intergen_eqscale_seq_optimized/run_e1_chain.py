#!/usr/bin/env python3
"""Run one recoverable E1 eqscale/sequential-fertility calibration chain."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[3]
M5_RESULTS = ROOT / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json"
PERIOD_YEARS = 4.0
POSTRETIREMENT_SURVIVAL = {
    66.0: 0.9391263063710125,
    70.0: 0.9184976343249724,
    74.0: 0.8849521927812863,
    78.0: 0.8300468061015381,
}
TARGET_SET = "candidate_replacement_income_disciplined_v1"
DOMAIN: tuple[tuple[str, float, float, str], ...] = (
    ("beta_annual", 0.80, 0.9995, "discount"),
    ("alpha_cons", 0.02, 0.98, "logit"),
    ("delta_alpha", 0.0, 0.25, "softzero"),
    ("delta_alpha_jump", 0.0, 0.25, "softzero"),
    ("gamma_e", 0.0, 3.0, "softzero"),
    ("psi_child", -3.0, 3.0, "asinh"),
    ("kappa_fert", 0.02, 50.0, "log"),
    ("chi", 0.10, 5.0, "log"),
    ("H0", 0.20, 80.0, "log"),
    ("theta0", 0.0, 8.0, "softzero"),
    ("theta1", 0.02, 16.0, "log"),
    ("tenure_choice_kappa", 0.0, 0.12, "softzero"),
)
FIXED = {"theta_n": 0.0}
DEFAULT_J = 17
DEFAULT_NB = 120
DEFAULT_MAX_ITER_EQ = 10


def load_runtime() -> None:
    """Delay model imports so the standard-library ``--help`` remains usable."""
    global np, base_overrides, diagnostic_loss, external_entry_wealth_overrides_1824
    global extract_moments, get_target_set, income_process_overrides, jsonable
    global production_profile_overrides, InfeasibleThetaError, run_model_cp_dt
    import numpy as np  # type: ignore[no-redef]
    from intergen_eqscale_seq_optimized.calibration import (
        base_overrides, diagnostic_loss, external_entry_wealth_overrides_1824,
        extract_moments, get_target_set,
    )
    from intergen_eqscale_seq_optimized.local_panel import income_process_overrides, jsonable
    from intergen_eqscale_seq_optimized.production_profile import production_profile_overrides
    from intergen_eqscale_seq_optimized.solver import InfeasibleThetaError, run_model_cp_dt


def transform(u: float, lo: float, hi: float, kind: str) -> float:
    u = float(np.clip(u, 0.0, 1.0))
    if kind == "log":
        return lo * (hi / lo) ** u
    if kind == "discount":
        return lo + (hi - lo) * (1.0 - (1.0 - u) ** 2)
    if kind == "logit":
        return lo + (hi - lo) / (1.0 + math.exp(-8.0 * (u - 0.5)))
    if kind == "asinh":
        return math.sinh((1.0 - u) * math.asinh(lo) + u * math.asinh(hi))
    if kind == "softzero":
        return lo + (hi - lo) * u * u
    raise ValueError(f"unknown transform {kind!r}")


def inverse(value: float, lo: float, hi: float, kind: str) -> float:
    value = float(np.clip(value, lo, hi))
    if kind == "log":
        return math.log(value / lo) / math.log(hi / lo)
    if kind == "discount":
        return 1.0 - math.sqrt(max(0.0, 1.0 - (value - lo) / (hi - lo)))
    if kind == "logit":
        q = np.clip((value - lo) / (hi - lo), 1e-12, 1.0 - 1e-12)
        return 0.5 + math.log(q / (1.0 - q)) / 8.0
    if kind == "asinh":
        return (math.asinh(value) - math.asinh(lo)) / (math.asinh(hi) - math.asinh(lo))
    if kind == "softzero":
        return math.sqrt(max(0.0, (value - lo) / (hi - lo)))
    raise ValueError(f"unknown transform {kind!r}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=2026071801)
    parser.add_argument("--start-mix", type=float, default=0.0)
    parser.add_argument("--max-evals", type=int, default=1000)
    parser.add_argument("--minutes", type=float, default=225.0)
    parser.add_argument("--initial-step", type=float, default=0.015)
    parser.add_argument("--shrink", type=float, default=0.5)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--J", type=int, default=DEFAULT_J)
    parser.add_argument("--Nb", type=int, default=DEFAULT_NB)
    parser.add_argument("--max-iter-eq", type=int, default=DEFAULT_MAX_ITER_EQ)
    parser.add_argument("--tol-eq", type=float, default=1e-4)
    return parser.parse_args()


def build_seed_theta(seed_path: Path = M5_RESULTS) -> dict[str, float]:
    """Return the E1 seed: M5 shared coordinates plus E1-only primitives.

    E2 continuation: when E2_SEED_RECORD is set, seed from that record's
    winners.E1 theta, which is already in the 12-parameter contract shape.
    """
    e2_record = os.environ.get("E2_SEED_RECORD", "")
    if e2_record:
        allowed = {"beta" if name == "beta_annual" else name for name, *_ in DOMAIN}
        required = allowed | {"theta_n"}
        payload_e2: Any = json.loads(Path(e2_record).read_text())
        winner_e2 = (payload_e2.get("winners") or {}).get("E1")
        if not isinstance(winner_e2, dict) or not isinstance(winner_e2.get("theta"), dict):
            raise ValueError(f"{e2_record} does not contain winners.E1.theta")
        theta = {k: float(v) for k, v in winner_e2["theta"].items() if k in required}
        if set(theta) != required:
            raise RuntimeError(f"E2 seed keys differ from contract: {sorted(theta)}")
        return theta
    payload: Any = json.loads(seed_path.read_text())
    winner = (payload.get("winners") or {}).get("M5")
    if not isinstance(winner, dict) or not isinstance(winner.get("theta"), dict):
        raise ValueError(f"{seed_path} does not contain winners.M5.theta")
    source = winner["theta"]
    allowed = {"beta" if name == "beta_annual" else name for name, *_ in DOMAIN}
    theta = {name: float(source[name]) for name in allowed - {"delta_alpha", "delta_alpha_jump", "gamma_e"}}
    theta.update(delta_alpha=0.05, delta_alpha_jump=0.10, gamma_e=0.5, theta_n=0.0)
    required = allowed | {"theta_n"}
    if set(theta) != required:
        raise RuntimeError(f"E1 seed keys differ from contract: {sorted(theta)}")
    return theta


def theta_from_unit(unit: np.ndarray) -> dict[str, float]:
    theta = dict(FIXED)
    for u, (name, lo, hi, kind) in zip(np.asarray(unit, dtype=float), DOMAIN):
        value = transform(float(u), lo, hi, kind)
        theta["beta" if name == "beta_annual" else name] = value ** PERIOD_YEARS if name == "beta_annual" else value
    return theta


def unit_from_theta(theta: dict[str, float]) -> np.ndarray:
    values = []
    for name, lo, hi, kind in DOMAIN:
        value = float(theta["beta" if name == "beta_annual" else name])
        if name == "beta_annual":
            value **= 1.0 / PERIOD_YEARS
        values.append(inverse(value, lo, hi, kind))
    return np.clip(np.asarray(values, dtype=float), 0.0, 1.0)


def survival_schedule(J: int) -> np.ndarray:
    schedule = np.ones(int(J) - 1, dtype=float)
    for j in range(int(J) - 1):
        age = 18.0 + j * PERIOD_YEARS
        if age in POSTRETIREMENT_SURVIVAL:
            schedule[j] = POSTRETIREMENT_SURVIVAL[age]
    return schedule


def common_overrides(args: argparse.Namespace) -> dict[str, Any]:
    return {
        **base_overrides(J=args.J, Nb=args.Nb, n_house=5, max_iter_eq=args.max_iter_eq),
        **production_profile_overrides(),
        **income_process_overrides(5, "rouwenhorst", 0.20, 0.9601845894041878),
        **external_entry_wealth_overrides_1824(),
        "preference_spec": "eqscale", "sequential_births": True,
        "fecundity_omega1": 0.02, "fecundity_omega2": 0.134,
        "fecundity_terminal_age": 45.0,
        "transfer_floor_G0": 0.0, "transfer_floor_Gn": 0.0,
        "use_age_survival": True, "entry_wealth_censor_to_frontier": True,
        "bequest_spec": "linear_child_scale", "normalize_bequest_utility": True,
        "owner_ltv_taper": False, "owner_ltv_taper_start_age": 66.0,
        "owner_ltv_taper_end_age": 82.0, "owner_ltv_terminal_share": 0.0,
        "max_iter_eq": int(args.max_iter_eq), "tol_eq": float(args.tol_eq),
        "q": (1.0 + 0.02) ** PERIOD_YEARS - 1.0,
        "delta": 1.0 - (1.0 - 0.011) ** PERIOD_YEARS,
        "eta_supply": np.array([1.75]), "lambda_d": 0.0,
        "debt_taper_start_age": 42.0, "debt_taper_end_age": 62.0,
        "survival_probs": survival_schedule(args.J),
    }


def target_fit(moments: dict[str, Any], targets: dict[str, float], weights: dict[str, float]) -> list[dict[str, float | str]]:
    return [{"moment": name, "target": float(target), "model": float(moments.get(name, math.nan)),
             "gap": float(moments.get(name, math.nan)) - float(target), "weight": float(weights[name]),
             "loss_contribution": float(weights[name]) * (float(moments.get(name, math.nan)) - float(target)) ** 2}
            for name, target in targets.items()]


def main() -> None:
    args = parse_args()
    load_runtime()
    if args.smoke:
        args.Nb, args.max_iter_eq, args.tol_eq = 60, 2, 0.25
        args.max_evals, args.minutes = min(args.max_evals, 13), min(args.minutes, 8.0)
    if not math.isclose(float(args.tol_eq), 1e-4) and not args.smoke:
        raise ValueError("E1 search requires tol_eq=1e-4; winners are repeated at 40/2.5e-5")
    targets, weights = get_target_set(TARGET_SET)
    # Atomic in the optimized target registry; assignment remains idempotent
    # for compatibility with the original runner contract.
    targets["aggregate_mean_occupied_rooms_18_85"] = 5.779970481941968
    weights["aggregate_mean_occupied_rooms_18_85"] = 6.0
    if len(targets) != 15 or set(targets) != set(weights):
        raise ValueError("E1 requires the unchanged 15-moment income-disciplined target system")
    seed_theta = build_seed_theta()
    x0 = unit_from_theta(seed_theta)
    rng = np.random.default_rng(args.seed)
    start_mix = float(np.clip(args.start_mix, 0.0, 0.25))
    if start_mix:
        x0 = np.clip((1.0 - start_mix) * x0 + start_mix * rng.random(len(DOMAIN)), 0.0, 1.0)
    overrides = common_overrides(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    cases_path, best_path = args.outdir / "cases.jsonl", args.outdir / "best.json"
    cases_path.write_text("")
    metadata = {"status": "smoke" if args.smoke else "proper_joint_smm_chain", "arm": "E1",
                "free_parameter_count": len(DOMAIN), "target_count": len(targets),
                "active_domain": [{"name": n, "lower": lo, "upper": hi, "transform": k} for n, lo, hi, k in DOMAIN],
                "fixed_parameters": FIXED, "target_set": TARGET_SET, "targets": targets, "weights": weights,
                "seed": args.seed, "start_mix": start_mix, "initial_unit_vector": x0,
                "J": args.J, "Nb": args.Nb, "max_iter_eq": args.max_iter_eq, "tol_eq": args.tol_eq,
                "income_process": {"states": 5, "process": "rouwenhorst", "annual_rho": 0.9601845894041878, "annual_innovation_sd": 0.20},
                "tight_winner_evaluator": {"max_iter_eq": 40, "tol_eq": 2.5e-5, "repeats": 2}}
    (args.outdir / "metadata.json").write_text(json.dumps(jsonable(metadata), indent=2, sort_keys=True))
    started, records, best, eval_idx = time.perf_counter(), [], None, 0
    reserve = 0.0 if args.smoke else min(300.0, max(0.0, args.minutes * 60.0 - 1.0))
    deadline = started + max(1.0, args.minutes * 60.0 - reserve)

    def score(record: dict[str, Any] | None) -> float:
        return float(record.get("rank_loss", math.inf)) if record and record.get("strict_converged") else math.inf

    def can_eval() -> bool:
        return eval_idx < args.max_evals and time.perf_counter() < deadline

    def evaluate(unit: np.ndarray, label: str, origin: dict[str, Any], tight: bool = False, tight_case: int | None = None) -> dict[str, Any]:
        nonlocal eval_idx, best
        theta, began = theta_from_unit(np.clip(unit, 0.0, 1.0)), time.perf_counter()
        try:
            sol, P, p_eq = run_model_cp_dt({**(overrides if not tight else {**overrides, "max_iter_eq": 40, "tol_eq": 2.5e-5}), **theta}, verbose=False)
            moments = extract_moments(sol, P); loss = float(diagnostic_loss(moments, targets=targets, weights=weights)
            )
            residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf)); timings = dict(getattr(sol, "timings", {}))
            strict = bool(timings.get("strict_converged", getattr(sol, "converged", False)) and residual <= (2.5e-5 if tight else P.tol_eq))
            status, error, census, price = "ok", "", [], float(np.asarray(p_eq).reshape(-1)[0])
        except InfeasibleThetaError as exc:
            moments, loss, residual, timings, strict, status, error, census, price = {}, math.inf, math.inf, {}, False, "infeasible_theta", str(exc), list(exc.census), math.nan
        except Exception as exc:  # persist each failed proposal as a recoverable checkpoint
            moments, loss, residual, timings, strict, status, error, census, price = {}, math.inf, math.inf, {}, False, f"failed:{type(exc).__name__}", str(exc), [], math.nan
        record = {"case": eval_idx if tight_case is None else tight_case, "label": label, "arm": "E1", "status": status, "strict_converged": strict,
                  "rank_loss": loss, "market_residual": residual, "price": price, "theta": theta, "moments": moments,
                  "target_fit": target_fit(moments, targets, weights) if moments else [], "timings": timings,
                  "timing_diagnostics": {k: v for k, v in moments.items() if isinstance(v, (list, tuple, np.ndarray))},
                  "origin": origin, "unit_vector": np.asarray(unit), "feasibility_census": census, "error": error,
                  "elapsed_sec": time.perf_counter() - began}
        path = args.outdir / ("tight_cases.jsonl" if tight else "cases.jsonl")
        with path.open("a") as handle: handle.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
        if tight:
            (args.outdir / "latest_tight_case.json").write_text(json.dumps(jsonable(record), indent=2, sort_keys=True))
        if not tight:
            records.append(record); (args.outdir / "latest_completed_case.json").write_text(json.dumps(jsonable(record), indent=2, sort_keys=True))
            if score(record) < score(best): best = record; best_path.write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
            eval_idx += 1
        return record

    simplex = [(x0.copy(), evaluate(x0, "seed", {"phase": "seed"}))]
    for j in range(len(DOMAIN)):
        if not can_eval(): break
        step = float(args.initial_step) * (0.75 + 0.5 * rng.random()); trial = x0.copy()
        trial[j] = trial[j] + step if trial[j] + step <= 1.0 else max(0.0, trial[j] - step)
        simplex.append((trial, evaluate(trial, f"nm_init_{j:02d}", {"phase": "initial_simplex", "dim": j, "step": step})))
    alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, float(np.clip(args.shrink, 0.1, 0.9))
    iteration = 0
    while can_eval() and len(simplex) >= 2:
        simplex.sort(key=lambda item: score(item[1]))
        centroid = np.mean([unit for unit, _ in simplex[:-1]], axis=0)
        worst_unit, worst_record = simplex[-1]
        second_worst, best_loss = score(simplex[-2][1]), score(simplex[0][1])
        reflected = np.clip(centroid + alpha * (centroid - worst_unit), 0.0, 1.0)
        reflected_record = evaluate(reflected, f"nm_reflect_{iteration:04d}", {"phase": "reflect", "iteration": iteration})
        reflected_loss = score(reflected_record)
        if reflected_loss < best_loss and can_eval():
            expanded = np.clip(centroid + gamma * (reflected - centroid), 0.0, 1.0)
            expanded_record = evaluate(expanded, f"nm_expand_{iteration:04d}", {"phase": "expand", "iteration": iteration})
            simplex[-1] = (expanded, expanded_record) if score(expanded_record) < reflected_loss else (reflected, reflected_record)
        elif reflected_loss < second_worst:
            simplex[-1] = (reflected, reflected_record)
        else:
            if reflected_loss < score(worst_record):
                contracted, threshold, phase = np.clip(centroid + rho * (reflected - centroid), 0.0, 1.0), reflected_loss, "outside_contract"
            else:
                contracted, threshold, phase = np.clip(centroid - rho * (centroid - worst_unit), 0.0, 1.0), score(worst_record), "inside_contract"
            contract_record = evaluate(contracted, f"nm_contract_{iteration:04d}", {"phase": phase, "iteration": iteration}) if can_eval() else worst_record
            if score(contract_record) < threshold:
                simplex[-1] = (contracted, contract_record)
            else:
                best_unit, new_simplex = simplex[0][0].copy(), [simplex[0]]
                for slot, (unit, old_record) in enumerate(simplex[1:], start=1):
                    if not can_eval():
                        new_simplex.append((unit, old_record)); continue
                    shrunk = np.clip(best_unit + sigma * (unit - best_unit), 0.0, 1.0)
                    shrunk_record = evaluate(shrunk, f"nm_shrink_{iteration:04d}_{slot:02d}", {"phase": "shrink", "iteration": iteration, "slot": slot})
                    new_simplex.append((shrunk, shrunk_record))
                simplex = new_simplex
        iteration += 1
    tight_records = []
    if not args.smoke and best is not None:
        (args.outdir / "tight_cases.jsonl").write_text("")
        for repeat in range(2):
            tight_records.append(evaluate(np.asarray(best["unit_vector"]), f"tight_winner_repeat_{repeat + 1}", {"phase": "tight_winner_verification", "search_case": best["case"]}, tight=True, tight_case=repeat))
    strict_tight = [row for row in tight_records if row.get("strict_converged")]
    best_tight = min(strict_tight, key=lambda row: float(row["rank_loss"])) if strict_tight else None
    if best_tight is not None: (args.outdir / "best_tight.json").write_text(json.dumps(jsonable(best_tight), indent=2, sort_keys=True))
    repeat_check = None
    if len(tight_records) == 2 and all(row.get("moments") for row in tight_records):
        repeat_check = {"loss_abs_difference": abs(float(tight_records[0]["rank_loss"]) - float(tight_records[1]["rank_loss"])),
                        "max_abs_moment_difference": max(abs(float(tight_records[0]["moments"][n]) - float(tight_records[1]["moments"][n])) for n in targets),
                        "both_strict": all(row.get("strict_converged") for row in tight_records)}
    summary = {"status": "complete" if eval_idx >= args.max_evals else "time_budget_reached", "best": best, "best_search": best,
               "best_tight": best_tight, "tight_repeats": tight_records, "tight_repeat_check": repeat_check,
               "n_cases_completed": len(records), "n_strict": sum(bool(r.get("strict_converged")) for r in records),
               "elapsed_sec": time.perf_counter() - started, "iterations_completed": iteration, "metadata": metadata}
    (args.outdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
