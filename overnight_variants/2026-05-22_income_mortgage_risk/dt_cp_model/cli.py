"""Command line tools for the Python DT model port."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .objective import smm_objective_dt
from .parameters import build_calibration_setup
from .evaluate import solve_theta
from .solver import accounting_population_scale, run_model_cp_dt
from .utils import make_grid


def main() -> None:
    parser = argparse.ArgumentParser(description="Python DT center-periphery model tools")
    sub = parser.add_subparsers(dest="cmd", required=True)

    smoke = sub.add_parser("smoke", help="Run a tiny fixed-price smoke solve")
    smoke.add_argument("--J", type=int, default=8)
    smoke.add_argument("--Nb", type=int, default=20)
    smoke.add_argument("--n-house", type=int, default=2)
    smoke.add_argument("--max-iter-eq", type=int, default=2)
    smoke.add_argument("--quiet", action="store_true")

    bench = sub.add_parser("benchmark", help="Run a GE benchmark solve")
    bench.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    bench.add_argument("--max-iter-eq", type=int, default=6)
    bench.add_argument("--force-full", action="store_true")
    bench.add_argument("--quiet", action="store_true")
    bench.add_argument("--json", type=Path, default=None)

    obj = sub.add_parser("objective", help="Evaluate the Stage A SMM objective at x0")
    obj.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    obj.add_argument("--max-iter-eq", type=int, default=4)
    obj.add_argument("--inv-iters", type=int, default=1)
    obj.add_argument("--quiet", action="store_true")

    theta = sub.add_parser("solve-theta", help="Solve one GE model from a theta vector, with no SMM loss or inversion")
    theta.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    theta.add_argument("--theta", default="x0", help="Comma-separated theta vector, or 'x0'")
    theta.add_argument("--max-iter-eq", type=int, default=None)
    theta.add_argument("--quiet", action="store_true")
    theta.add_argument("--json", type=Path, default=None)

    scale = sub.add_parser("accounting-scale", help="Compute fast implied population scale from a normalized solve")
    scale.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    scale.add_argument("--max-iter-eq", type=int, default=1)
    scale.add_argument("--force-full", action="store_true")
    scale.add_argument("--outside-entry-flow", type=float, default=None)
    scale.add_argument("--local-birth-weight", type=float, default=1.0)
    scale.add_argument("--outside-value", type=float, default=None)
    scale.add_argument("--kappa-entry", type=float, default=2.0)
    scale.add_argument("--quiet", action="store_true")
    scale.add_argument("--json", type=Path, default=None)

    scaled = sub.add_parser("scaled-equilibrium", help="Solve GE with accounting-scale housing demand")
    scaled.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    scaled.add_argument("--max-iter-eq", type=int, default=60)
    scaled.add_argument("--baseline-max-iter-eq", type=int, default=120)
    scaled.add_argument("--outside-entry-flow", type=float, default=None)
    scaled.add_argument("--local-birth-weight", type=float, default=1.0)
    scaled.add_argument("--outside-value", type=float, default=None)
    scaled.add_argument("--kappa-entry", type=float, default=2.0)
    scaled.add_argument("--force-full", action="store_true")
    scaled.add_argument("--quiet", action="store_true")
    scaled.add_argument("--json", type=Path, default=None)

    args = parser.parse_args()
    if args.cmd == "smoke":
        run_smoke(args)
    elif args.cmd == "benchmark":
        run_benchmark(args)
    elif args.cmd == "objective":
        run_objective(args)
    elif args.cmd == "solve-theta":
        run_solve_theta(args)
    elif args.cmd == "accounting-scale":
        run_accounting_scale(args)
    elif args.cmd == "scaled-equilibrium":
        run_scaled_equilibrium(args)


def run_smoke(args) -> None:
    override = {
        "J": args.J,
        "J_R": max(2, args.J - 2),
        "A_f_end": min(4, args.J),
        "Nb": args.Nb,
        "n_house": args.n_house,
        "H_own": np.linspace(4.0, 8.0, args.n_house),
        "max_iter_eq": args.max_iter_eq,
        "tol_eq": 1e-2,
        "force_full_bellman": True,
        "solve_mode": "pe",
        "p_fixed": np.array([1.0, 1.2]),
        "w_fixed": np.array([1.0, 1.0]),
        "entry_shares_fixed": np.array([0.55, 0.45]),
    }
    t0 = time.perf_counter()
    sol, P, p_eq = run_model_cp_dt(override, verbose=not args.quiet)
    elapsed = time.perf_counter() - t0
    report = summarize_solution(sol, p_eq, elapsed)
    print(json.dumps(report, indent=2, sort_keys=True))


def run_benchmark(args) -> None:
    setup = build_calibration_setup(args.setup)
    P = setup.P_base
    P.max_iter_eq = args.max_iter_eq
    if args.force_full:
        P.force_full_bellman = True
    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=not args.quiet)
    elapsed = time.perf_counter() - t0
    report = summarize_solution(sol, p_eq, elapsed)
    report["setup"] = args.setup
    report["max_iter_eq"] = args.max_iter_eq
    report["timings"] = getattr(sol, "timings", {})
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, sort_keys=True))


def run_objective(args) -> None:
    setup = build_calibration_setup(args.setup)
    setup.P_base.max_iter_eq = args.max_iter_eq
    setup.inversion_targets["max_iter"] = args.inv_iters
    t0 = time.perf_counter()
    loss, moments, _, converged, _, _ = smm_objective_dt(
        setup.x0, setup.targets, setup.weights, setup.P_base, setup.inversion_targets, verbose=not args.quiet
    )
    elapsed = time.perf_counter() - t0
    report = {
        "loss": loss,
        "converged": converged,
        "elapsed_sec": elapsed,
        "moments": {k: _jsonable(v) for k, v in vars(moments).items()},
    }
    print(json.dumps(report, indent=2, sort_keys=True))


def run_solve_theta(args) -> None:
    setup = build_calibration_setup(args.setup)
    if args.theta.strip().lower() == "x0":
        theta = setup.x0
    else:
        theta = np.array([float(x) for x in args.theta.split(",")])
    t0 = time.perf_counter()
    sol, _, p_eq = solve_theta(theta, setup_mode=args.setup, max_iter_eq=args.max_iter_eq, verbose=not args.quiet)
    elapsed = time.perf_counter() - t0
    report = summarize_solution(sol, p_eq, elapsed)
    report["setup"] = args.setup
    report["max_iter_eq"] = args.max_iter_eq
    report["timings"] = getattr(sol, "timings", {})
    report["theta"] = [float(x) for x in theta]
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, sort_keys=True))


def run_accounting_scale(args) -> None:
    setup = build_calibration_setup(args.setup)
    P = setup.P_base
    P.max_iter_eq = args.max_iter_eq
    if args.force_full:
        P.force_full_bellman = True
    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=not args.quiet)
    b_grid = make_grid(P_out)
    scale = accounting_population_scale(
        sol,
        P_out,
        b_grid,
        outside_entry_flow=args.outside_entry_flow,
        local_birth_entry_weight=args.local_birth_weight,
        outside_value=args.outside_value,
        kappa_entry=args.kappa_entry,
        calibrate_outside_value=args.outside_value is None,
    )
    elapsed = time.perf_counter() - t0
    report = summarize_solution(sol, p_eq, elapsed)
    report["setup"] = args.setup
    report["max_iter_eq"] = args.max_iter_eq
    report["accounting_scale"] = _jsonable(vars(scale))
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, sort_keys=True))


def run_scaled_equilibrium(args) -> None:
    setup0 = build_calibration_setup(args.setup)
    P0 = setup0.P_base
    P0.max_iter_eq = args.baseline_max_iter_eq
    sol0, P0_out, p0 = run_model_cp_dt(P0, verbose=False)
    b_grid0 = make_grid(P0_out)
    scale0 = accounting_population_scale(
        sol0,
        P0_out,
        b_grid0,
        outside_entry_flow=args.outside_entry_flow,
        local_birth_entry_weight=args.local_birth_weight,
        outside_value=args.outside_value,
        kappa_entry=args.kappa_entry,
        calibrate_outside_value=args.outside_value is None,
    )

    setup = build_calibration_setup(args.setup)
    P = setup.P_base
    P.max_iter_eq = args.max_iter_eq
    P.population_closure = "accounting_scale_prices"
    P.outside_value = scale0.outside_value
    P.outside_value_is_calibrated = True
    P.kappa_entry = args.kappa_entry
    P.local_birth_entry_weight = args.local_birth_weight
    if args.force_full:
        P.force_full_bellman = True
    if args.outside_entry_flow is not None:
        P.outside_entry_flow = args.outside_entry_flow
    P.collect_ge_trace = True

    t0 = time.perf_counter()
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=not args.quiet)
    elapsed = time.perf_counter() - t0
    report = summarize_solution(sol, p_eq, elapsed)
    report["setup"] = args.setup
    report["max_iter_eq"] = args.max_iter_eq
    report["population_closure"] = P_out.population_closure
    report["timings"] = getattr(sol, "timings", {})
    report["outside_value_calibration"] = {
        "outside_value": float(scale0.outside_value),
        "baseline_p_eq": [float(x) for x in p0],
        "baseline_tfr": float(2 * sol0.mean_parity),
        "baseline_scale": _jsonable(vars(scale0)),
    }
    if hasattr(sol, "accounting_scale"):
        report["accounting_scale"] = _jsonable(vars(sol.accounting_scale))
    if hasattr(sol, "ge_trace"):
        report["first_trace"] = sol.ge_trace[:3]
        report["last_trace"] = sol.ge_trace[-8:]
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, sort_keys=True))


def summarize_solution(sol: SimpleNamespace, p_eq: np.ndarray, elapsed: float) -> dict:
    return {
        "elapsed_sec": elapsed,
        "p_eq": [float(x) for x in p_eq],
        "tfr": float(2 * sol.mean_parity),
        "own_rate": float(sol.own_rate),
        "pop_share": [float(x) for x in sol.pop_share],
        "mean_age_first_birth": float(getattr(sol, "mean_age_first_birth", np.nan)),
        "migration_rate_2245": float(getattr(sol, "migration_rate_2245", np.nan)),
        "prime_childless_renter_median_rooms": float(getattr(sol, "prime_childless_renter_median_rooms", np.nan)),
        "prime_childless_owner_median_rooms": float(getattr(sol, "prime_childless_owner_median_rooms", np.nan)),
        "housing_increment_0to1": float(getattr(sol, "housing_increment_0to1_eventstudy_t3", np.nan)),
        "housing_increment_1to2": float(getattr(sol, "housing_increment_1to2_proxy_t3", np.nan)),
        "total_mass": float(getattr(sol, "total_mass", np.nan)),
    }


def _jsonable(v):
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    return v


if __name__ == "__main__":
    main()
