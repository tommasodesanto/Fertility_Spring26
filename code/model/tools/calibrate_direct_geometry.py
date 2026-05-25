#!/usr/bin/env python3
"""Run a direct no-inversion outside-option calibration job.

Each invocation is one sequential optimizer worker. On the cluster we launch an
array of independent workers with different seeds and collect the best JSON
records afterwards.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dt_cp_model.direct_calibration import (  # noqa: E402
    DirectCalibrationSetup,
    build_direct_calibration_setup,
    direct_theta_record,
    evaluate_direct_theta,
)


def main() -> None:
    args = parse_args()
    setup = build_direct_calibration_setup(
        args.setup,
        geo_weight=args.geo_weight,
        population_closure=args.population_closure,
        scale_target=args.scale_target,
        scale_weight=args.scale_weight,
        outside_value_x0=args.outside_value_x0,
        outside_flow_x0=args.outside_flow_x0,
        renewal_retention=args.renewal_retention,
        target_city_entry_prob=args.target_city_entry_prob,
        kappa_entry=args.kappa_entry,
        housing_product_market=args.housing_product_market,
    )
    apply_run_overrides(setup, args)
    lb, ub = apply_bound_profile(setup, args.bound_profile)
    seed = int(args.seed_base + 1009 * args.job_id)
    rng = np.random.default_rng(seed)

    job_dir = args.results_dir / f"job_{args.job_id:03d}"
    job_dir.mkdir(parents=True, exist_ok=True)
    eval_path = job_dir / "evaluations.jsonl"
    best_path = job_dir / "best.json"
    status_path = job_dir / "status.json"
    config_path = job_dir / "config.json"

    config = {
        "job_id": args.job_id,
        "run_tag": args.run_tag,
        "seed": seed,
        "setup": args.setup,
        "bound_profile": args.bound_profile,
        "max_iter_eq": args.max_iter_eq,
        "max_evals": args.max_evals,
        "budget_sec": args.budget_sec,
        "geo_weight": args.geo_weight,
        "population_closure": args.population_closure,
        "scale_target": args.scale_target,
        "scale_weight": args.scale_weight,
        "outside_value_x0": args.outside_value_x0,
        "outside_flow_x0": args.outside_flow_x0,
        "renewal_retention": args.renewal_retention,
        "target_city_entry_prob": args.target_city_entry_prob,
        "kappa_entry": args.kappa_entry,
        "housing_product_market": args.housing_product_market,
        "hR_max": float(setup.P_base.hR_max),
        "H_own": [float(x) for x in np.asarray(setup.P_base.H_own).reshape(-1)],
        "eq_penalty_weight": args.eq_penalty_weight,
        "max_tfr": args.max_tfr,
        "theta_names": setup.names,
        "lb": lb.tolist(),
        "ub": ub.tolist(),
        "x0": setup.x0.tolist(),
    }
    write_json(config_path, config)

    best_record: dict[str, Any] | None = None
    if args.resume and best_path.exists():
        best_record = json.loads(best_path.read_text())
        print(f"[resume] loaded best loss {best_record['loss']:.6g} from {best_path}", flush=True)
    elif eval_path.exists():
        eval_path.unlink()

    best_theta = np.array(best_record["theta"], dtype=float) if best_record else None
    best_loss = float(best_record["loss"]) if best_record else math.inf
    sigma = args.initial_scale * (ub - lb)
    seed_bank = build_seed_bank(setup, lb, ub, args.job_id)
    t_start = time.perf_counter()
    n_eval = count_jsonl(eval_path) if args.resume and eval_path.exists() else 0
    no_improve = 0

    print_header(args, job_dir, seed, lb, ub, setup)

    while n_eval < args.max_evals:
        elapsed = time.perf_counter() - t_start
        if elapsed >= args.budget_sec:
            break

        theta = propose_theta(n_eval, seed_bank, best_theta, lb, ub, sigma, rng, args.global_prob)
        result = evaluate_direct_theta(
            theta,
            setup,
            max_iter_eq=args.max_iter_eq,
            force_full=args.force_full,
            eq_penalty_weight=args.eq_penalty_weight,
            max_tfr=args.max_tfr,
            verbose=args.verbose_solver,
        )
        n_eval += 1

        improved = result.loss < best_loss
        if improved:
            best_loss = result.loss
            best_theta = result.theta.copy()
            no_improve = 0
            sigma = np.minimum(0.50 * (ub - lb), args.grow * sigma)
        else:
            no_improve += 1
            if no_improve > 0 and no_improve % args.stall_window == 0:
                sigma = np.maximum(args.min_scale * (ub - lb), args.shrink * sigma)

        record = make_record(result, setup, args, n_eval, elapsed=time.perf_counter() - t_start, improved=improved)
        append_jsonl(eval_path, record)
        if improved or best_record is None:
            best_record = record
            write_json(best_path, record)

        status = {
            "job_id": args.job_id,
            "run_tag": args.run_tag,
            "n_eval": n_eval,
            "elapsed_sec": time.perf_counter() - t_start,
            "best_loss": best_loss,
            "best_theta": best_theta.tolist() if best_theta is not None else None,
            "last_loss": result.loss,
            "last_elapsed_sec": result.elapsed_sec,
            "last_error": result.error,
            "sigma": sigma.tolist(),
        }
        write_json(status_path, status)

        if improved or n_eval == 1 or n_eval % args.log_every == 0:
            tag = "BEST" if improved else "eval"
            print(
                f"[{tag} {n_eval:04d}] loss={result.loss:.6g} solve={result.elapsed_sec:.1f}s "
                f"TFR={record['moments'].get('tfr', math.nan):.3f} "
                f"own={record['moments'].get('own_rate', math.nan):.3f} "
                f"N={record['moments'].get('implied_total_population', math.nan):.3f} "
                f"popC={record['moments'].get('inv_pop_share_C', math.nan):.3f} "
                f"rr={record['moments'].get('inv_rent_ratio_C_over_P', math.nan):.3f}",
                flush=True,
            )

    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    status["finished"] = True
    status["finish_reason"] = "max_evals" if n_eval >= args.max_evals else "time_budget"
    status["elapsed_sec"] = time.perf_counter() - t_start
    write_json(status_path, status)
    print(f"finished job {args.job_id}: evals={n_eval}, best_loss={best_loss:.6g}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Direct no-inversion outside-option calibration worker.")
    parser.add_argument("--job-id", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", "1")))
    parser.add_argument("--run-tag", default=os.environ.get("DT_DIRECT_RUN_TAG", time.strftime("direct_geo_%Y%m%d_%H%M%S")))
    parser.add_argument("--setup", choices=["fast", "benchmark"], default=os.environ.get("DT_DIRECT_SETUP", "benchmark"))
    parser.add_argument("--bound-profile", choices=["scout", "wide", "global"], default=os.environ.get("DT_DIRECT_BOUNDS", "scout"))
    parser.add_argument("--results-dir", type=Path, default=Path(os.environ.get("DT_DIRECT_RESULTS_DIR", "results_python_direct_geometry")))
    parser.add_argument("--max-iter-eq", type=int, default=int(os.environ.get("DT_DIRECT_MAX_ITER_EQ", "0")))
    parser.add_argument("--max-evals", type=int, default=int(os.environ.get("DT_DIRECT_MAX_EVALS", "1000000")))
    parser.add_argument("--budget-sec", type=float, default=float(os.environ.get("DT_DIRECT_BUDGET_SEC", "41400")))
    parser.add_argument("--seed-base", type=int, default=int(os.environ.get("DT_DIRECT_SEED_BASE", "91000")))
    parser.add_argument("--geo-weight", type=float, default=float(os.environ.get("DT_DIRECT_GEO_WEIGHT", "100")))
    parser.add_argument(
        "--population-closure",
        choices=[
            "renewal_valve_calibrated",
            "renewal_valve",
            "accounting_scale_prices",
            "outside_option_benchmark_normalized",
            "outside_option",
            "normalized",
        ],
        default=os.environ.get("DT_DIRECT_POPULATION_CLOSURE", "outside_option_benchmark_normalized"),
    )
    parser.add_argument("--scale-target", type=float, default=float(os.environ.get("DT_DIRECT_SCALE_TARGET", "1.0")))
    parser.add_argument("--scale-weight", type=float, default=float(os.environ.get("DT_DIRECT_SCALE_WEIGHT", "100")))
    parser.add_argument(
        "--outside-value-x0",
        type=float,
        default=float(os.environ.get("DT_DIRECT_OUTSIDE_VALUE_X0", "-41.95")),
    )
    parser.add_argument(
        "--outside-flow-x0",
        type=float,
        default=float(os.environ.get("DT_DIRECT_OUTSIDE_FLOW_X0", "0.003")),
    )
    parser.add_argument(
        "--renewal-retention",
        type=float,
        default=float(os.environ.get("DT_DIRECT_RENEWAL_RETENTION", "1.0")),
    )
    parser.add_argument(
        "--target-city-entry-prob",
        type=float,
        default=float(os.environ.get("DT_DIRECT_TARGET_CITY_ENTRY_PROB", "0.89")),
    )
    parser.add_argument(
        "--kappa-entry",
        type=float,
        default=(
            float(os.environ["DT_DIRECT_KAPPA_ENTRY"])
            if "DT_DIRECT_KAPPA_ENTRY" in os.environ
            else None
        ),
    )
    parser.add_argument(
        "--housing-product-market",
        action="store_true",
        default=os.environ.get("DT_DIRECT_HOUSING_PRODUCT_MARKET", "0").lower() in ("1", "true", "yes"),
    )
    parser.add_argument(
        "--h-own-min",
        type=float,
        default=(
            float(os.environ["DT_DIRECT_H_OWN_MIN"])
            if "DT_DIRECT_H_OWN_MIN" in os.environ
            else None
        ),
    )
    parser.add_argument("--hR-max", default=os.environ.get("DT_DIRECT_HR_MAX"))
    parser.add_argument("--eq-penalty-weight", type=float, default=float(os.environ.get("DT_DIRECT_EQ_PENALTY_WEIGHT", "0")))
    parser.add_argument(
        "--max-tfr",
        type=float,
        default=(
            float(os.environ["DT_DIRECT_MAX_TFR"])
            if "DT_DIRECT_MAX_TFR" in os.environ
            else None
        ),
        help="Optional hard cap: candidates with TFR >= this value receive loss 1e6.",
    )
    parser.add_argument("--initial-scale", type=float, default=float(os.environ.get("DT_DIRECT_INITIAL_SCALE", "0.18")))
    parser.add_argument("--min-scale", type=float, default=float(os.environ.get("DT_DIRECT_MIN_SCALE", "0.015")))
    parser.add_argument("--global-prob", type=float, default=float(os.environ.get("DT_DIRECT_GLOBAL_PROB", "0.12")))
    parser.add_argument("--stall-window", type=int, default=int(os.environ.get("DT_DIRECT_STALL_WINDOW", "12")))
    parser.add_argument("--shrink", type=float, default=float(os.environ.get("DT_DIRECT_SHRINK", "0.72")))
    parser.add_argument("--grow", type=float, default=float(os.environ.get("DT_DIRECT_GROW", "1.04")))
    parser.add_argument("--log-every", type=int, default=int(os.environ.get("DT_DIRECT_LOG_EVERY", "5")))
    parser.add_argument("--force-full", action="store_true")
    parser.add_argument("--verbose-solver", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.max_iter_eq <= 0:
        args.max_iter_eq = None
    return args


def apply_run_overrides(setup: DirectCalibrationSetup, args: argparse.Namespace) -> None:
    if args.h_own_min is not None:
        current = np.asarray(setup.P_base.H_own, dtype=float).reshape(-1)
        h_min = float(args.h_own_min)
        h_max = float(np.max(current))
        if not np.isfinite(h_min) or h_min < 0 or h_min >= h_max:
            raise ValueError(f"Invalid owner-grid minimum override: {args.h_own_min}")
        setup.P_base.H_own = np.linspace(h_min, h_max, len(current))
        setup.P_base.n_house = len(current)
        setup.P_base.h_own_min = h_min
        setup.P_base.h_own_max = h_max
    if args.hR_max is None or str(args.hR_max).strip() == "":
        return
    raw = str(args.hR_max).strip().lower()
    if raw in ("owner_max", "max_owner", "h_own_max", "common", "due_common"):
        hR_max = float(np.max(np.asarray(setup.P_base.H_own, dtype=float)))
    else:
        hR_max = float(raw)
    if not np.isfinite(hR_max) or hR_max <= 0:
        raise ValueError(f"Invalid hR_max override: {args.hR_max}")
    setup.P_base.hR_max = hR_max


def apply_bound_profile(setup: DirectCalibrationSetup, profile: str) -> tuple[np.ndarray, np.ndarray]:
    lb = setup.lb.copy()
    ub = setup.ub.copy()
    if profile == "global":
        return lb, ub
    if profile == "wide":
        profile_bounds = {
            "beta": (0.920, 0.980),
            "b_entry_fixed": (0.00, 0.70),
            "psi_child": (0.030, 0.120),
            "h_bar_jump": (0.35, 2.50),
            "h_bar_n": (0.20, 1.00),
            "c_bar_n": (0.06, 0.23),
            "kappa_fert": (1.30, 7.00),
            "chi": (1.02, 1.10),
            "kappa_loc": (0.55, 3.50),
            "mu_move": (0.00, 8.00),
            "theta0": (0.10, 1.40),
            "theta_n": (0.00, 0.70),
            "h_bar_0": (2.00, 4.50),
            "E_C": (-1.00, 1.80),
            "r_bar_C": (0.020, 0.240),
            "outside_value": (-70.00, -25.00),
            "outside_entry_flow": (0.0005, 0.0300),
        }
    elif profile == "scout":
        profile_bounds = {
            "beta": (0.925, 0.970),
            "b_entry_fixed": (0.00, 0.50),
            "psi_child": (0.040, 0.120),
            "h_bar_jump": (0.50, 2.50),
            "h_bar_n": (0.30, 1.00),
            "c_bar_n": (0.07, 0.20),
            "kappa_fert": (1.50, 6.00),
            "chi": (1.03, 1.10),
            "kappa_loc": (0.70, 3.00),
            "mu_move": (0.00, 5.00),
            "theta0": (0.20, 1.20),
            "theta_n": (0.00, 0.60),
            "h_bar_0": (2.50, 4.50),
            "E_C": (-0.50, 1.20),
            "r_bar_C": (0.040, 0.180),
            "outside_value": (-58.00, -32.00),
            "outside_entry_flow": (0.0010, 0.0200),
        }
    else:
        raise ValueError(f"unknown bound profile: {profile}")
    name_to_idx = {name: idx for idx, name in enumerate(setup.names)}
    for name, (lo, hi) in profile_bounds.items():
        if name not in name_to_idx:
            continue
        idx = name_to_idx[name]
        lb[idx] = max(lb[idx], float(lo))
        ub[idx] = min(ub[idx], float(hi))
    return lb, ub


def build_seed_bank(setup: DirectCalibrationSetup, lb: np.ndarray, ub: np.ndarray, job_id: int) -> list[np.ndarray]:
    base = setup.x0.copy()
    names = {name: idx for idx, name in enumerate(setup.names)}

    def seed(**updates: float) -> np.ndarray:
        theta = base.copy()
        for key, value in updates.items():
            if key in names:
                theta[names[key]] = float(value)
        return np.clip(theta, lb, ub)

    seeds = [
        seed(),
        seed(beta=0.940, psi_child=0.120, h_bar_jump=2.30, h_bar_n=1.00, c_bar_n=0.10, chi=1.09, kappa_loc=1.50, mu_move=0.04, outside_value=-42.0, outside_entry_flow=0.003),
        seed(beta=0.940, psi_child=0.115, h_bar_jump=2.30, h_bar_n=1.00, c_bar_n=0.10, chi=1.08, kappa_loc=1.50, mu_move=0.04, E_C=0.45, r_bar_C=0.08, outside_value=-45.0, outside_entry_flow=0.004),
        seed(beta=0.955, psi_child=0.085, h_bar_jump=1.70, h_bar_n=0.70, c_bar_n=0.10, kappa_fert=3.50, kappa_loc=2.00, mu_move=2.00, outside_value=-42.0, outside_entry_flow=0.006),
        seed(beta=0.948, psi_child=0.070, h_bar_jump=1.20, h_bar_n=0.60, c_bar_n=0.12, kappa_fert=4.00, kappa_loc=1.20, mu_move=1.00, E_C=0.15, r_bar_C=0.11, outside_value=-39.0, outside_entry_flow=0.008),
        seed(beta=0.965, psi_child=0.050, h_bar_jump=0.90, h_bar_n=0.45, c_bar_n=0.08, kappa_fert=2.50, chi=1.05, kappa_loc=2.50, mu_move=3.00, E_C=0.75, r_bar_C=0.07, outside_value=-46.0, outside_entry_flow=0.010),
    ]
    shift = (job_id - 1) % len(seeds)
    return seeds[shift:] + seeds[:shift]


def propose_theta(
    n_eval: int,
    seed_bank: list[np.ndarray],
    best_theta: np.ndarray | None,
    lb: np.ndarray,
    ub: np.ndarray,
    sigma: np.ndarray,
    rng: np.random.Generator,
    global_prob: float,
) -> np.ndarray:
    if n_eval < len(seed_bank):
        return seed_bank[n_eval].copy()
    if best_theta is None or rng.random() < global_prob:
        return rng.uniform(lb, ub)
    theta = best_theta + rng.normal(0.0, sigma)
    if n_eval % (len(best_theta) + 1) == 0:
        theta = best_theta.copy()
        idx = int(rng.integers(0, len(theta)))
        theta[idx] += rng.normal(0.0, 2.0 * sigma[idx])
    return reflect_bounds(theta, lb, ub)


def reflect_bounds(theta: np.ndarray, lb: np.ndarray, ub: np.ndarray) -> np.ndarray:
    out = np.asarray(theta, dtype=float).copy()
    width = ub - lb
    for k in range(len(out)):
        if width[k] <= 0:
            out[k] = lb[k]
            continue
        while out[k] < lb[k] or out[k] > ub[k]:
            if out[k] < lb[k]:
                out[k] = lb[k] + (lb[k] - out[k])
            if out[k] > ub[k]:
                out[k] = ub[k] - (out[k] - ub[k])
        out[k] = min(max(out[k], lb[k]), ub[k])
    return out


def make_record(result, setup: DirectCalibrationSetup, args: argparse.Namespace, eval_id: int, elapsed: float, improved: bool) -> dict[str, Any]:
    return {
        "run_tag": args.run_tag,
        "job_id": args.job_id,
        "eval_id": eval_id,
        "elapsed_sec_total": elapsed,
        "loss": float(result.loss),
        "improved": bool(improved),
        "solve_ok": bool(result.solve_ok),
        "converged": bool(result.converged),
        "error": result.error,
        "solve_elapsed_sec": float(result.elapsed_sec),
        "theta": [float(x) for x in result.theta],
        "parameters": direct_theta_record(result.theta, setup.names),
        "moments": result.moments,
        "p_eq": result.p_eq,
        "timings": result.timings,
    }


def print_header(
    args: argparse.Namespace,
    job_dir: Path,
    seed: int,
    lb: np.ndarray,
    ub: np.ndarray,
    setup: DirectCalibrationSetup,
) -> None:
    print("=" * 72)
    print("Direct no-inversion outside-option Python calibration")
    print(f"job_id={args.job_id} run_tag={args.run_tag} seed={seed}")
    print(
        f"setup={args.setup} bounds={args.bound_profile} "
        f"closure={args.population_closure} max_iter_eq={args.max_iter_eq}"
    )
    if args.population_closure in ("renewal_valve_calibrated", "outside_option_benchmark_normalized"):
        print(f"scale_target={args.scale_target} imposed mechanically; scale_weight inactive")
    else:
        print(f"scale_target={args.scale_target} scale_weight={args.scale_weight}")
    print(f"target_city_entry_prob={args.target_city_entry_prob} kappa_entry={args.kappa_entry}")
    print(f"housing_product_market={args.housing_product_market}")
    print(f"hR_max={setup.P_base.hR_max} H_own={np.asarray(setup.P_base.H_own).tolist()}")
    print(f"max_tfr={args.max_tfr}")
    print(f"renewal_retention={args.renewal_retention}")
    print(f"budget_sec={args.budget_sec:.0f} max_evals={args.max_evals} job_dir={job_dir}")
    print(f"dim={len(lb)}")
    print("=" * 72, flush=True)


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True))


def append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(obj, sort_keys=True) + "\n")


def count_jsonl(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


if __name__ == "__main__":
    main()
