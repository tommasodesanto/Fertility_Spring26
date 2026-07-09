"""Command line tools for the intergenerational housing/fertility model."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from .calibration import TARGET_SETS, run_informed_smoke, run_small_calibration
from .diagnostics import write_diagnostics
from .local_panel import run_global_de_panel, run_local_panel, run_local_polish
from .production_profile import (
    PRODUCTION_INCOME_STATES,
    PRODUCTION_J,
    PRODUCTION_PROFILE_NAME,
    PRODUCTION_SEARCH_NB,
    PRODUCTION_TARGET_SET,
)
from .solver import run_model_cp_dt


PERIOD_YEARS = 4.0
AGE_START = 18.0
FERTILITY_START_AGE = 18.0
FERTILITY_END_AGE = 42.0
RETIREMENT_AGE = 66.0
CHILD_MATURITY_AGE = 18.0


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
    solve.add_argument("--J", type=int, default=17)
    solve.add_argument("--Nb", type=int, default=60)
    solve.add_argument("--n-house", type=int, default=6)
    solve.add_argument("--quiet", action="store_true")
    solve.add_argument("--json", type=Path, default=None)

    diag = sub.add_parser("diagnostics", help="Solve and write a diagnostic packet")
    diag.add_argument("--fixed-prices", action="store_true")
    diag.add_argument("--max-iter-eq", type=int, default=20)
    diag.add_argument("--J", type=int, default=17)
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
    cal.add_argument("--target-set", choices=sorted(TARGET_SETS), default="old_nonlocation")
    cal.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_small_calibration"))

    informed = sub.add_parser("informed-smoke", help="Run a deterministic parameter-led smoke panel")
    informed.add_argument("--J", type=int, default=17)
    informed.add_argument("--Nb", type=int, default=40)
    informed.add_argument("--n-house", type=int, default=6)
    informed.add_argument("--max-iter-eq", type=int, default=25)
    informed.add_argument("--labels", type=str, default="")
    informed.add_argument("--case-limit", type=int, default=None)
    informed.add_argument("--fixed-price", type=float, default=None)
    informed.add_argument("--quiet", action="store_true")
    informed.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_informed_smoke"))

    panel = sub.add_parser("local-panel", help="Run a bounded multicore local diagnostic panel")
    panel.add_argument("--cases", type=int, default=144)
    panel.add_argument("--seed", type=int, default=20260608)
    panel.add_argument("--J", type=int, default=17)
    panel.add_argument("--Nb", type=int, default=60)
    panel.add_argument("--income-states", type=int, default=5)
    panel.add_argument("--n-house", type=int, default=6)
    panel.add_argument("--max-iter-eq", type=int, default=25)
    panel.add_argument("--workers", type=int, default=6)
    panel.add_argument("--minutes", type=float, default=30.0)
    panel.add_argument("--diagnostic-best", type=int, default=3)
    panel.add_argument("--target-set", choices=sorted(TARGET_SETS), default="candidate_no_timing_v0")
    panel.add_argument("--random-only", action="store_true", help="Skip deterministic anchor cases in local-panel draws")
    panel.add_argument("--seed-theta-json", type=Path, default=None, help="JSON file with a top-level theta object to seed the panel.")
    panel.add_argument("--quiet", action="store_true")
    panel.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_local_panel"))

    global_de = sub.add_parser("global-de-panel", help="Run a bounded differential-evolution global panel")
    global_de.add_argument("--max-evals", type=int, default=240)
    global_de.add_argument("--seed", type=int, default=20260609)
    global_de.add_argument("--J", type=int, default=17)
    global_de.add_argument("--Nb", type=int, default=60)
    global_de.add_argument("--income-states", type=int, default=5)
    global_de.add_argument("--n-house", type=int, default=6)
    global_de.add_argument("--max-iter-eq", type=int, default=25)
    global_de.add_argument("--minutes", type=float, default=115.0)
    global_de.add_argument("--pop-size", type=int, default=20)
    global_de.add_argument("--mutation", type=float, default=0.85)
    global_de.add_argument("--crossover", type=float, default=0.70)
    global_de.add_argument("--target-set", choices=sorted(TARGET_SETS), default="candidate_no_timing_v0")
    global_de.add_argument("--seed-theta-json", type=Path, default=None, help="JSON file with a top-level theta object to seed the DE population.")
    global_de.add_argument("--quiet", action="store_true")
    global_de.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_global_de"))

    polish = sub.add_parser("local-polish", help="Run a true derivative-free local polish around a seed theta")
    polish.add_argument("--max-evals", type=int, default=240)
    polish.add_argument("--seed", type=int, default=20260628)
    polish.add_argument("--J", type=int, default=PRODUCTION_J)
    polish.add_argument("--Nb", type=int, default=PRODUCTION_SEARCH_NB)
    polish.add_argument("--income-states", type=int, default=PRODUCTION_INCOME_STATES)
    polish.add_argument("--n-house", type=int, default=5)
    polish.add_argument("--max-iter-eq", type=int, default=25)
    polish.add_argument("--minutes", type=float, default=115.0)
    polish.add_argument("--target-set", choices=sorted(TARGET_SETS), default=PRODUCTION_TARGET_SET)
    polish.add_argument("--seed-theta-json", type=Path, required=True, help="JSON file with a top-level theta object to polish.")
    polish.add_argument("--method", choices=["nelder-mead", "pattern"], default="nelder-mead")
    polish.add_argument("--initial-step", type=float, default=0.06, help="Initial step in normalized bounded parameter space.")
    polish.add_argument("--min-step", type=float, default=0.003, help="Minimum step in normalized bounded parameter space.")
    polish.add_argument("--shrink", type=float, default=0.5)
    polish.add_argument(
        "--profile",
        choices=[PRODUCTION_PROFILE_NAME, "none"],
        default=PRODUCTION_PROFILE_NAME,
        help="Named source-controlled model/search profile; use 'none' only for legacy diagnostics.",
    )
    polish.add_argument("--fixed-beta-annual", type=float, default=None)
    polish.add_argument("--fixed-theta-n", type=float, default=None)
    polish.add_argument("--fixed-chi", type=float, default=None)
    polish.add_argument("--quiet", action="store_true")
    polish.add_argument("--outdir", type=Path, default=Path("../../output/model/intergen_housing_fertility_local_polish"))

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
            target_set=str(args.target_set),
        )
        print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))
        return
    if args.cmd == "informed-smoke":
        summary = run_informed_smoke(
            args.outdir,
            J=int(args.J),
            Nb=int(args.Nb),
            n_house=int(args.n_house),
            max_iter_eq=int(args.max_iter_eq),
            labels=[x.strip() for x in str(args.labels).split(",") if x.strip()],
            case_limit=args.case_limit,
            fixed_price=args.fixed_price,
            progress=not bool(args.quiet),
        )
        print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))
        return
    if args.cmd == "local-panel":
        summary = run_local_panel(
            args.outdir,
            n_cases=int(args.cases),
            seed=int(args.seed),
            J=int(args.J),
            Nb=int(args.Nb),
            income_states=int(args.income_states),
            n_house=int(args.n_house),
            max_iter_eq=int(args.max_iter_eq),
            workers=int(args.workers),
            minutes=float(args.minutes),
            diagnostic_best=int(args.diagnostic_best),
            target_set=str(args.target_set),
            include_anchors=not bool(args.random_only),
            seed_theta=load_seed_theta(args.seed_theta_json),
            progress=not bool(args.quiet),
        )
        print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))
        return
    if args.cmd == "global-de-panel":
        summary = run_global_de_panel(
            args.outdir,
            max_evals=int(args.max_evals),
            seed=int(args.seed),
            J=int(args.J),
            Nb=int(args.Nb),
            income_states=int(args.income_states),
            n_house=int(args.n_house),
            max_iter_eq=int(args.max_iter_eq),
            minutes=float(args.minutes),
            pop_size=int(args.pop_size),
            mutation=float(args.mutation),
            crossover=float(args.crossover),
            target_set=str(args.target_set),
            seed_theta=load_seed_theta(args.seed_theta_json),
            progress=not bool(args.quiet),
        )
        print(json.dumps(_jsonable(summary), indent=2, sort_keys=True))
        return
    if args.cmd == "local-polish":
        fixed_theta: dict[str, float] = {}
        if args.fixed_beta_annual is not None:
            fixed_theta["beta"] = float(args.fixed_beta_annual) ** PERIOD_YEARS
        if args.fixed_theta_n is not None:
            fixed_theta["theta_n"] = float(args.fixed_theta_n)
        if args.fixed_chi is not None:
            fixed_theta["chi"] = float(args.fixed_chi)
        summary = run_local_polish(
            args.outdir,
            max_evals=int(args.max_evals),
            seed=int(args.seed),
            J=int(args.J),
            Nb=int(args.Nb),
            income_states=int(args.income_states),
            n_house=int(args.n_house),
            max_iter_eq=int(args.max_iter_eq),
            minutes=float(args.minutes),
            target_set=str(args.target_set),
            seed_theta=load_seed_theta(args.seed_theta_json),
            method=str(args.method),
            initial_step=float(args.initial_step),
            min_step=float(args.min_step),
            shrink=float(args.shrink),
            profile_name=None if str(args.profile) == "none" else str(args.profile),
            fixed_theta=fixed_theta,
            progress=not bool(args.quiet),
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
        "use_pti_constraint": False,
        "pti_limit": 0.30,
    }
    if extra:
        overrides.update(extra)
    return overrides


def load_seed_theta(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    with Path(path).open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    theta = payload.get("theta", payload)
    if not isinstance(theta, dict):
        raise ValueError(f"{path} does not contain a theta object")
    return dict(theta)


def smoke_overrides(args: argparse.Namespace) -> dict[str, Any]:
    j = int(args.J)
    n_house = int(args.n_house)
    lifecycle = lifecycle_overrides(j)
    return one_market_overrides(
        {
            **lifecycle,
            "J": j,
            "entry_by_loc": np.array([1.0 / j]),
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
    lifecycle = lifecycle_overrides(j)
    return {
        **lifecycle,
        "J": j,
        "entry_by_loc": np.array([1.0 / j]),
        "Nb": int(args.Nb),
        "n_house": n_house,
        "H_own": np.linspace(2.0, 11.0, n_house),
        "max_iter_eq": int(args.max_iter_eq),
    }


def lifecycle_overrides(J: int) -> dict[str, Any]:
    j = int(J)
    fertile_start = age_to_period_number(FERTILITY_START_AGE)
    fertile_end = age_to_period_number(FERTILITY_END_AGE)
    retirement_idx = int(round((RETIREMENT_AGE - AGE_START) / PERIOD_YEARS))
    return {
        "period_years": PERIOD_YEARS,
        "da": PERIOD_YEARS,
        "age_start": AGE_START,
        "A_m": CHILD_MATURITY_AGE,
        "stage_durations": np.array([CHILD_MATURITY_AGE / PERIOD_YEARS]),
        "J_R": max(2, min(j - 1, retirement_idx)),
        "A_f_start": max(1, min(j, fertile_start)),
        "A_f_end": max(1, min(j, fertile_end)),
    }


def age_to_period_number(age: float) -> int:
    return int(round((float(age) - AGE_START) / PERIOD_YEARS)) + 1


def summarize_solution(sol, P, p_eq, elapsed: float) -> dict[str, Any]:
    return {
        "elapsed_sec": elapsed,
        "J": P.J,
        "period_years": getattr(P, "period_years", P.da),
        "age_start": getattr(P, "age_start", np.nan),
        "retirement_age": getattr(P, "age_start", np.nan) + getattr(P, "J_R", np.nan) * getattr(P, "da", np.nan),
        "fertility_start_age": getattr(P, "age_start", np.nan) + (getattr(P, "A_f_start", np.nan) - 1) * getattr(P, "da", np.nan),
        "fertility_end_age": getattr(P, "age_start", np.nan) + (getattr(P, "A_f_end", np.nan) - 1) * getattr(P, "da", np.nan),
        "child_maturity_age": getattr(P, "A_m", np.nan),
        "stage_durations": getattr(P, "stage_durations", np.array([])),
        "expected_child_years": np.asarray(getattr(P, "stage_durations", np.array([np.nan])), dtype=float) * getattr(P, "period_years", getattr(P, "da", np.nan)),
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
        "aggregate_own_rate": sol.own_rate,
        "own_rate_3055": getattr(sol, "own_rate_3055", np.nan),
        "target_own_rate": getattr(sol, "own_rate_3055", np.nan),
        "young_owner_rate": sol.young_owner_rate,
        "old_owner_rate": sol.old_owner_rate,
        "tfr": 2.0 * sol.mean_completed_fertility,
        "mean_completed_fertility": sol.mean_completed_fertility,
        "tfr_by_income_type": 2.0 * getattr(sol, "mean_fertility_by_income_type", np.array([])),
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
