"""Small Markov-income Howard speed audit for the test package."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from . import solver
from .calibration import base_overrides, extract_moments
from .local_panel import income_process_overrides


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare full-only and Markov Howard solves.")
    parser.add_argument("--J", type=int, default=16)
    parser.add_argument("--Nb", type=int, default=60)
    parser.add_argument("--income-states", type=int, default=5)
    parser.add_argument("--n-house", type=int, default=5)
    parser.add_argument("--max-iter-eq", type=int, default=3)
    parser.add_argument("--howard-freq", type=int, default=5)
    parser.add_argument("--howard-warmup", type=int, default=1)
    parser.add_argument("--polish-iter", type=int, default=3)
    parser.add_argument("--refine-mode", choices=("full", "eval"), default="full")
    parser.add_argument("--variant", choices=("both", "full", "howard"), default="both")
    parser.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_howard_test_speed_audit_20260622"))
    parser.add_argument("--label", type=str, default="")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    label = args.label or f"j{args.J}_nb{args.Nb}_nz{args.income_states}_h{args.n_house}_maxiter{args.max_iter_eq}"
    prefix = "comparison" if args.variant == "both" else args.variant
    outpath = args.outdir / f"{prefix}_polished_{label}.json"

    common = {
        **base_overrides(J=args.J, Nb=args.Nb, n_house=args.n_house, max_iter_eq=args.max_iter_eq),
        **income_process_overrides(args.income_states),
        "howard_freq": int(args.howard_freq),
        "markov_howard_warmup": int(args.howard_warmup),
        "markov_howard_refine_full_polish": True,
        "markov_howard_refine_full_polish_iter": int(args.polish_iter),
    }
    full_only = {
        **common,
        "force_full_bellman": True,
        "use_markov_howard": False,
        "use_markov_howard_refine": False,
    }
    markov_howard = {
        **common,
        "force_full_bellman": False,
        "use_markov_howard": True,
        "use_markov_howard_refine": str(args.refine_mode) == "eval",
    }

    result = {
        "status": "running",
        "package": "intergen_housing_fertility_howard_test",
        "command": " ".join(sys.argv),
        "configuration": {
            "J": int(args.J),
            "Nb": int(args.Nb),
            "income_states": int(args.income_states),
            "n_house": int(args.n_house),
            "H_own": jsonable(common["H_own"]),
            "max_iter_eq": int(args.max_iter_eq),
            "howard_freq": int(args.howard_freq),
            "howard_warmup": int(args.howard_warmup),
            "scalar_market_refine": bool(common["scalar_market_refine"]),
            "scalar_market_refine_method": str(common["scalar_market_refine_method"]),
            "scalar_market_refine_iter": int(common["scalar_market_refine_iter"]),
            "markov_howard_refine_full_polish": True,
            "markov_howard_refine_full_polish_iter": int(args.polish_iter),
            "markov_howard_refine_mode": str(args.refine_mode),
            "z_grid": jsonable(common["z_grid"]),
            "z_weights": jsonable(common["z_weights"]),
        },
    }
    if args.variant in {"both", "full"}:
        result["full_only"] = run_case("full_only", full_only, verbose=not args.quiet)
    if args.variant in {"both", "howard"}:
        result["markov_howard"] = run_case("markov_howard", markov_howard, verbose=not args.quiet)
    if args.variant == "both":
        result["deltas_howard_minus_full"] = summarize_deltas(result["markov_howard"], result["full_only"])
    result["status"] = "ok"

    outpath.write_text(json.dumps(jsonable(result), indent=2, sort_keys=True))
    report = {"status": "ok", "outpath": str(outpath)}
    if "deltas_howard_minus_full" in result:
        report["deltas_howard_minus_full"] = result["deltas_howard_minus_full"]
    print(json.dumps(jsonable(report), indent=2, sort_keys=True))


def run_case(label: str, overrides: dict[str, Any], *, verbose: bool) -> dict[str, Any]:
    fixed_price_calls: list[dict[str, Any]] = []
    original = solver.solve_markov_income_at_prices

    def traced_solve_markov_income_at_prices(p_eq, P, b_grid, *args, **kwargs):
        t0 = time.perf_counter()
        sol = original(p_eq, P, b_grid, *args, **kwargs)
        elapsed = time.perf_counter() - t0
        timings = getattr(sol, "timings", {})
        fixed_price_calls.append(
            {
                "call": len(fixed_price_calls) + 1,
                "price": jsonable(np.asarray(p_eq, dtype=float).reshape(-1)),
                "elapsed_sec": float(elapsed),
                "fast_stats": bool(kwargs.get("fast_stats", False)),
                "eval_mode_requested": bool(kwargs.get("eval_mode", False)),
                "bellman_mode": str(timings.get("bellman_mode", "")),
                "timings": jsonable(timings),
                "housing_demand": jsonable(getattr(sol, "housing_demand", np.array([]))),
                "housing_supply": jsonable(getattr(sol, "housing_supply", np.array([]))),
                "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
            }
        )
        return sol

    solver.solve_markov_income_at_prices = traced_solve_markov_income_at_prices
    t0 = time.perf_counter()
    try:
        sol, P, p_eq = solver.run_model_cp_dt(overrides, verbose=verbose)
    finally:
        solver.solve_markov_income_at_prices = original
    wall = time.perf_counter() - t0

    timings = getattr(sol, "timings", {})
    return {
        "label": label,
        "configuration": {
            "force_full_bellman": bool(overrides.get("force_full_bellman", False)),
            "use_markov_howard": bool(overrides.get("use_markov_howard", False)),
            "use_markov_howard_refine": bool(overrides.get("use_markov_howard_refine", False)),
        },
        "status": "ok",
        "wall_time_sec": float(wall),
        "p_eq": jsonable(p_eq),
        "market_residual": float(getattr(sol, "best_max_abs_rel_excess", np.nan)),
        "moments": jsonable(extract_moments(sol, P)),
        "final_solution_timings": jsonable(timings),
        "fixed_price_call_count": len(fixed_price_calls),
        "aggregate_fixed_price_timings": aggregate_call_timings(fixed_price_calls),
        "fixed_price_calls": fixed_price_calls,
    }


def aggregate_call_timings(calls: list[dict[str, Any]]) -> dict[str, Any]:
    timings = [c.get("timings", {}) for c in calls]
    return {
        "fixed_price_elapsed_sec": float(sum(float(c.get("elapsed_sec", 0.0)) for c in calls)),
        "bellman_full_sec": float(sum(float(t.get("bellman_full", 0.0)) for t in timings)),
        "bellman_eval_sec": float(sum(float(t.get("bellman_eval", 0.0)) for t in timings)),
        "distribution_sec": float(sum(float(t.get("distribution", 0.0)) for t in timings)),
        "n_full": int(sum(int(t.get("n_full", 0)) for t in timings)),
        "n_eval": int(sum(int(t.get("n_eval", 0)) for t in timings)),
        "n_dist": int(sum(int(t.get("n_dist", 0)) for t in timings)),
        "fast_stats_calls": int(sum(1 for c in calls if c.get("fast_stats"))),
        "full_stats_calls": int(sum(1 for c in calls if not c.get("fast_stats"))),
    }


def summarize_deltas(howard: dict[str, Any], full: dict[str, Any]) -> dict[str, Any]:
    out = {
        "wall_time_sec": float(howard["wall_time_sec"]) - float(full["wall_time_sec"]),
        "market_residual": float(howard["market_residual"]) - float(full["market_residual"]),
        "price_max_abs": max_abs_diff(howard.get("p_eq", []), full.get("p_eq", [])),
        "fixed_price_calls": int(howard["fixed_price_call_count"]) - int(full["fixed_price_call_count"]),
    }
    h_agg = howard.get("aggregate_fixed_price_timings", {})
    f_agg = full.get("aggregate_fixed_price_timings", {})
    for key in ("bellman_full_sec", "bellman_eval_sec", "distribution_sec", "n_full", "n_eval", "n_dist"):
        out[key] = float(h_agg.get(key, 0.0)) - float(f_agg.get(key, 0.0))
    h_mom = howard.get("moments", {})
    f_mom = full.get("moments", {})
    out["moments_max_abs"] = max(
        (abs(float(h_mom[k]) - float(f_mom[k])) for k in h_mom.keys() & f_mom.keys() if np.isfinite(float(h_mom[k])) and np.isfinite(float(f_mom[k]))),
        default=float("nan"),
    )
    return out


def max_abs_diff(a: Any, b: Any) -> float:
    aa = np.asarray(a, dtype=float).reshape(-1)
    bb = np.asarray(b, dtype=float).reshape(-1)
    if aa.shape != bb.shape:
        return float("nan")
    return float(np.max(np.abs(aa - bb))) if aa.size else 0.0


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


if __name__ == "__main__":
    main()
