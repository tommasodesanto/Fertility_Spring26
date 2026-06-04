"""Command line tools for the intergenerational housing/fertility model."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from .calibration import run_small_calibration
from .diagnostics import write_diagnostics
from .solver import run_model_cp_dt


def main() -> None:
    parser = argparse.ArgumentParser(description="Intergenerational housing/fertility model tools")
    sub = parser.add_subparsers(dest="cmd", required=True)

    smoke = sub.add_parser("smoke", help="Run a small fixed-price smoke solve")
    smoke.add_argument("--J", type=int, default=8)
    smoke.add_argument("--Nb", type=int, default=20)
    smoke.add_argument("--n-house", type=int, default=3)
    smoke.add_argument("--quiet", action="store_true")
    smoke.add_argument("--json", type=Path, default=None)

    solve = sub.add_parser("solve", help="Run the one-market housing price iteration")
    solve.add_argument("--max-iter-eq", type=int, default=20)
    solve.add_argument("--J", type=int, default=16)
    solve.add_argument("--Nb", type=int, default=60)
    solve.add_argument("--n-house", type=int, default=6)
    solve.add_argument("--quiet", action="store_true")
    solve.add_argument("--json", type=Path, default=None)

    diag = sub.add_parser("diagnostics", help="Solve and write a diagnostic packet")
    diag.add_argument("--fixed-prices", action="store_true")
    diag.add_argument("--max-iter-eq", type=int, default=20)
    diag.add_argument("--J", type=int, default=16)
    diag.add_argument("--Nb", type=int, default=60)
    diag.add_argument("--n-house", type=int, default=6)
    diag.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_smoke"))
    diag.add_argument("--quiet", action="store_true")
    diag.add_argument("--json", type=Path, default=None)

    cal = sub.add_parser("calibrate-small", help="Run a small diagnostic local calibration")
    cal.add_argument("--cases", type=int, default=24)
    cal.add_argument("--seed", type=int, default=1234)
    cal.add_argument("--J", type=int, default=12)
    cal.add_argument("--Nb", type=int, default=40)
    cal.add_argument("--n-house", type=int, default=4)
    cal.add_argument("--max-iter-eq", type=int, default=35)
    cal.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_small_calibration"))

    args = parser.parse_args()
    if args.cmd == "calibrate-small":
        summary = run_small_calibration(
            args.outdir,
            n_cases=int(args.cases),
            seed=int(args.seed),
            J=int(args.J),
            Nb=int(args.Nb),
            n_house=int(args.n_house),
            max_iter_eq=int(args.max_iter_eq),
        )
        print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))
        return

    if args.cmd == "smoke":
        overrides = smoke_overrides(args)
    elif args.cmd == "solve":
        overrides = one_market_overrides(run_size_overrides(args))
    elif args.cmd == "diagnostics":
        base: dict[str, Any] = run_size_overrides(args)
        if args.fixed_prices:
            base.update({"solve_mode": "pe", "p_fixed": np.array([1.0]), "w_fixed": np.array([1.0]), "entry_shares_fixed": np.array([1.0])})
        overrides = one_market_overrides(base)
    else:
        raise ValueError(args.cmd)

    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(overrides, verbose=not args.quiet)
    elapsed = time.perf_counter() - t0
    report = summarize_solution(sol, P, p_eq, elapsed)

    if args.cmd == "diagnostics":
        write_diagnostics(sol, P, args.outdir)
        report["outdir"] = str(args.outdir)

    text = json.dumps(_jsonable(report), indent=2, sort_keys=True)
    print(text)
    if getattr(args, "json", None) is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text)


def one_market_overrides(extra: dict[str, Any] | None = None) -> dict[str, Any]:
    overrides: dict[str, Any] = {
        "I": 1,
        "w_hat": np.array([1.0]),
        "E_loc": np.array([0.0]),
        "N_0": np.array([1.0]),
        "entry_shares": np.array([1.0]),
        "entry_by_loc": np.array([1.0 / 16.0]),
        "r_bar": np.array([0.16]),
        "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]),
        "stage_durations": np.array([1.0]),
        "use_pti_constraint": True,
        "pti_limit": 0.30,
    }
    if extra:
        overrides.update(extra)
    return overrides


def smoke_overrides(args: argparse.Namespace) -> dict[str, Any]:
    j = int(args.J)
    n_house = int(args.n_house)
    return one_market_overrides(
        {
            "J": j,
            "J_R": max(2, min(j - 2, 6)),
            "entry_by_loc": np.array([1.0 / j]),
            "A_f_start": 2,
            "A_f_end": min(4, j),
            "Nb": int(args.Nb),
            "n_house": n_house,
            "H_own": np.linspace(2.0, 8.0, n_house),
            "max_iter_eq": 2,
            "tol_eq": 1e-2,
            "force_full_bellman": True,
            "solve_mode": "pe",
            "p_fixed": np.array([1.0]),
            "w_fixed": np.array([1.0]),
            "entry_shares_fixed": np.array([1.0]),
        }
    )


def run_size_overrides(args: argparse.Namespace) -> dict[str, Any]:
    j = int(args.J)
    n_house = int(args.n_house)
    return {
        "J": j,
        "J_R": max(2, min(j - 2, 11)),
        "entry_by_loc": np.array([1.0 / j]),
        "A_f_start": 2,
        "A_f_end": min(6, j),
        "Nb": int(args.Nb),
        "n_house": n_house,
        "H_own": np.linspace(2.0, 11.0, n_house),
        "max_iter_eq": int(args.max_iter_eq),
    }


def summarize_solution(sol, P, p_eq, elapsed: float) -> dict[str, Any]:
    return {
        "elapsed_sec": elapsed,
        "J": P.J,
        "period_years": getattr(P, "period_years", P.da),
        "n_child_stages": P.n_child_stages,
        "markets": P.I,
        "income_process": str(getattr(P, "income_type_transition", "none")),
        "income_states": getattr(sol, "type_values", getattr(P, "z_grid", np.array([1.0]))),
        "income_state_weights": getattr(sol, "type_weights", getattr(P, "z_weights", np.array([1.0]))),
        "income_transition": getattr(sol, "income_transition", getattr(P, "Pi_z", np.eye(1))),
        "income_types": getattr(sol, "type_values", getattr(P, "z_grid", np.array([1.0]))),
        "income_type_weights": getattr(sol, "type_weights", getattr(P, "z_weights", np.array([1.0]))),
        "pti_constraint": bool(getattr(P, "use_pti_constraint", False)),
        "pti_limit": float(getattr(P, "pti_limit", np.nan)),
        "p_eq": p_eq,
        "owner_user_cost": sol.owner_user_cost,
        "housing_demand": sol.housing_demand,
        "housing_supply": sol.housing_supply,
        "aggregate_housing_demand": sol.aggregate_housing_demand,
        "aggregate_housing_supply": sol.aggregate_housing_supply,
        "aggregate_housing_excess": sol.aggregate_housing_excess,
        "best_max_abs_rel_excess": sol.best_max_abs_rel_excess,
        "own_rate": sol.own_rate,
        "young_owner_rate": sol.young_owner_rate,
        "old_owner_rate": sol.old_owner_rate,
        "mean_completed_fertility": sol.mean_completed_fertility,
        "own_rate_by_income_type": getattr(sol, "own_rate_by_income_type", np.array([])),
        "mean_fertility_by_income_type": getattr(sol, "mean_fertility_by_income_type", np.array([])),
        "housing_demand_by_income_type": getattr(sol, "housing_demand_by_income_type", np.array([])),
        "childless_rate": sol.childless_rate,
        "mean_age_first_birth": getattr(sol, "mean_age_first_birth", np.nan),
        "timings": getattr(sol, "timings", {}),
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


if __name__ == "__main__":
    main()
