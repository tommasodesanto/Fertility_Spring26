#!/usr/bin/env python3
"""Parallel global search for the V5 Nz=5 outside-closure prototype."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from dt_cp_model.direct_calibration import compute_smm_loss
from dt_cp_model.objective import extract_moments
from dt_cp_model.solver import make_grid, solve_equilibrium_hank_z

from run_income_mortgage_risk_v2_scenarios import BENCHMARK_P, MOMENT_ORDER, namespace_float_dict
from run_income_mortgage_risk_v5_hank_z_outside_closure import prepare_parameters


PARAM_SPECS = [
    ("beta", 0.900, 0.990),
    ("b_entry_fixed", 0.000, 1.000),
    ("psi_child", 0.020, 0.120),
    ("h_bar_jump", 0.200, 2.500),
    ("h_bar_n", 0.100, 1.000),
    ("c_bar_n", 0.050, 0.250),
    ("kappa_fert", 1.250, 8.000),
    ("chi", 1.010, 1.100),
    ("kappa_loc", 0.500, 4.000),
    ("mu_move", 0.000, 10.000),
    ("theta0", 0.100, 1.500),
    ("theta_n", 0.000, 0.750),
    ("h_bar_0", 1.500, 4.500),
    ("E_C", -1.500, 2.500),
    ("r_bar_C", 0.010, 0.300),
    ("alpha_cons", 0.550, 0.900),
    ("phi", 0.600, 0.950),
    ("hR_max", 5.000, 9.000),
    ("h_own_max", 7.000, 12.000),
]

PARAM_NAMES = [name for name, _, _ in PARAM_SPECS]
LB = np.array([lo for _, lo, _ in PARAM_SPECS], dtype=float)
UB = np.array([hi for _, _, hi in PARAM_SPECS], dtype=float)


def configure_closure(P: SimpleNamespace, config: dict[str, Any]) -> SimpleNamespace:
    P.population_closure = "outside_option_benchmark_normalized"
    P.kappa_entry = float(config["kappa_entry"])
    P.local_birth_entry_weight = 1.0
    P.target_city_entry_prob = float(config["target_city_entry_prob"])
    P.calibrate_outside_value_to_entry_prob = True
    P.outside_benchmark_target_total_population = float(getattr(P, "N_target", 1.0))
    P.allow_uncalibrated_outside_value = False
    P.collect_ge_trace = True
    return P


def vector_from_params(P: SimpleNamespace) -> np.ndarray:
    values = []
    for name in PARAM_NAMES:
        if name == "E_C":
            values.append(float(P.E_loc[1]))
        elif name == "r_bar_C":
            values.append(float(P.r_bar[1]))
        elif name == "phi":
            values.append(float(np.asarray(P.phi).reshape(-1)[0]))
        else:
            values.append(float(getattr(P, name)))
    return np.clip(np.array(values, dtype=float), LB, UB)


def params_from_vector(theta: np.ndarray) -> dict[str, float]:
    return {name: float(value) for name, value in zip(PARAM_NAMES, np.asarray(theta, dtype=float))}


def apply_vector(P: SimpleNamespace, theta: np.ndarray) -> SimpleNamespace:
    params = params_from_vector(theta)
    for name, value in params.items():
        if name == "E_C":
            P.E_loc = np.asarray(P.E_loc, dtype=float).copy()
            P.E_loc[1] = value
        elif name == "r_bar_C":
            P.r_bar = np.asarray(P.r_bar, dtype=float).copy()
            P.r_bar[1] = value
        elif name == "phi":
            P.phi = value * np.ones(P.n_parity)
        elif name == "h_own_max":
            P.h_own_max = value
            P.H_own = np.linspace(float(P.h_own_min), float(P.h_own_max), int(P.n_house))
        else:
            setattr(P, name, value)
    P.eps_fert = float(P.kappa_fert)
    P.eps_loc = float(P.kappa_loc)
    P.R_gross = 1.0 + float(P.q)
    P.rho = 1.0 / float(P.beta) - 1.0
    P.rho_hat = P.rho
    P.user_cost_rate = float(P.q) + float(P.delta) + float(P.tau_H)
    return P


def evaluate_candidate(eval_id: int, theta_raw: list[float], config: dict[str, Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    theta = np.asarray(theta_raw, dtype=float)
    record: dict[str, Any] = {
        "eval_id": int(eval_id),
        "worker_pid": int(os.getpid()),
        "theta": [float(x) for x in theta],
        "parameters": params_from_vector(theta),
        "ok": False,
        "error": None,
    }
    try:
        P, targets, weights = prepare_parameters(
            int(config["nb"]),
            int(config["max_iter_eq"]),
            float(config["tol_eq"]),
        )
        P = configure_closure(P, config)
        P = apply_vector(P, theta)
        b_grid = make_grid(P)
        sol, P_out, p_eq = solve_equilibrium_hank_z(
            BENCHMARK_P,
            P,
            b_grid,
            nz=int(config["nz"]),
            rho_z=float(config["rho_z"]),
            sigma_z=float(config["sigma_z"]),
            verbose=False,
        )
        moments_ns = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        moments_ns.inv_pop_share_C = float(sol.pop_share[1])
        moments_ns.inv_rent_ratio_C_over_P = float(
            (P_out.user_cost_rate * p_eq[1]) / max(P_out.user_cost_rate * p_eq[0], 1e-12)
        )
        moments = namespace_float_dict(moments_ns)
        moment_loss = compute_smm_loss(moments, targets, weights)
        timings = dict(getattr(sol, "timings", {}))
        scale = getattr(sol, "accounting_scale", SimpleNamespace())
        final_err = float(timings.get("final_eq_error", np.nan))
        accepted = bool(timings.get("accepted", False))
        finite_scale = bool(getattr(scale, "finite_stationary_scale", False))
        scale_factor = float(getattr(scale, "scale_factor", np.nan))
        outside_flow = float(getattr(scale, "outside_entry_flow", np.nan))
        penalty = compute_penalty(
            accepted=accepted,
            final_eq_error=final_err,
            tol_eq=float(config["tol_eq"]),
            finite_scale=finite_scale,
            scale_factor=scale_factor,
            outside_flow=outside_flow,
            eq_penalty=float(config["eq_penalty"]),
            scale_penalty=float(config["scale_penalty"]),
        )
        objective_loss = float(moment_loss + penalty)
        record.update(
            {
                "ok": True,
                "solve_ok": True,
                "accepted": accepted,
                "strict_converged": bool(timings.get("strict_converged", False)),
                "convergence_reason": timings.get("convergence_reason"),
                "iterations_completed": int(timings.get("iterations_completed", 0)),
                "best_eq_error": float(timings.get("best_eq_error", np.nan)),
                "final_eq_error": final_err,
                "moment_loss": float(moment_loss),
                "penalty": float(penalty),
                "objective_loss": objective_loss,
                "moments": moments,
                "targets": targets,
                "weights": weights,
                "p_eq": [float(x) for x in p_eq],
                "timings": timings,
                "scale": {
                    "finite_stationary_scale": finite_scale,
                    "scale_factor": scale_factor,
                    "city_entry_prob_total": float(getattr(scale, "city_entry_prob_total", np.nan)),
                    "outside_entry_prob": float(getattr(scale, "outside_entry_prob", np.nan)),
                    "outside_entry_flow": outside_flow,
                    "outside_value": float(getattr(scale, "outside_value", np.nan)),
                    "denominator": float(getattr(scale, "denominator", np.nan)),
                },
            }
        )
    except Exception as exc:
        record.update(
            {
                "solve_ok": False,
                "accepted": False,
                "moment_loss": math.inf,
                "penalty": math.inf,
                "objective_loss": math.inf,
                "error": repr(exc),
            }
        )
    record["elapsed_sec"] = float(time.perf_counter() - t0)
    return to_builtin(record)


def compute_penalty(
    *,
    accepted: bool,
    final_eq_error: float,
    tol_eq: float,
    finite_scale: bool,
    scale_factor: float,
    outside_flow: float,
    eq_penalty: float,
    scale_penalty: float,
) -> float:
    penalty = 0.0
    if (not accepted) or (not np.isfinite(final_eq_error)):
        penalty += eq_penalty
    elif final_eq_error > tol_eq:
        penalty += eq_penalty * (final_eq_error / max(tol_eq, 1e-12) - 1.0) ** 2
    if (not finite_scale) or (not np.isfinite(scale_factor)) or (not np.isfinite(outside_flow)) or outside_flow <= 0.0:
        penalty += scale_penalty
        if np.isfinite(outside_flow) and outside_flow <= 0.0:
            penalty += scale_penalty * min(1.0, abs(outside_flow) / 0.01)
    return float(penalty)


def seed_bank(config: dict[str, Any]) -> list[np.ndarray]:
    P, _, _ = prepare_parameters(int(config["nb"]), int(config["max_iter_eq"]), float(config["tol_eq"]))
    P = configure_closure(P, config)
    base = vector_from_params(P)
    idx = {name: k for k, name in enumerate(PARAM_NAMES)}

    def seed(**updates: float) -> np.ndarray:
        theta = base.copy()
        for name, value in updates.items():
            theta[idx[name]] = float(value)
        return np.clip(theta, LB, UB)

    return [
        seed(),
        seed(alpha_cons=0.80, kappa_fert=4.00, phi=0.90),
        seed(alpha_cons=0.80),
        seed(kappa_fert=4.00, phi=0.90),
        seed(kappa_fert=5.00),
        seed(phi=0.90),
        seed(beta=0.965),
        seed(beta=0.972, phi=0.90, kappa_fert=4.00),
        seed(E_C=0.25, r_bar_C=0.055),
        seed(psi_child=0.12),
        seed(c_bar_n=0.075, h_bar_jump=0.20, h_bar_n=0.55),
        seed(theta0=0.30, theta_n=0.30),
    ]


def propose_candidate(
    eval_id: int,
    seeds: list[np.ndarray],
    best_theta: np.ndarray | None,
    sigma: np.ndarray,
    rng: np.random.Generator,
    global_prob: float,
) -> np.ndarray:
    if eval_id <= len(seeds):
        return seeds[eval_id - 1].copy()
    if best_theta is None or rng.random() < global_prob:
        return rng.uniform(LB, UB)
    theta = best_theta + rng.normal(0.0, sigma)
    if rng.random() < 0.20:
        k = int(rng.integers(0, len(theta)))
        theta[k] = rng.uniform(LB[k], UB[k])
    return reflect_bounds(theta)


def reflect_bounds(theta: np.ndarray) -> np.ndarray:
    out = np.asarray(theta, dtype=float).copy()
    width = UB - LB
    for k in range(len(out)):
        if width[k] <= 0:
            out[k] = LB[k]
            continue
        while out[k] < LB[k] or out[k] > UB[k]:
            if out[k] < LB[k]:
                out[k] = LB[k] + (LB[k] - out[k])
            if out[k] > UB[k]:
                out[k] = UB[k] - (out[k] - UB[k])
        out[k] = min(max(out[k], LB[k]), UB[k])
    return out


def main() -> None:
    args = parse_args()
    run_dir = args.results_dir / args.run_tag
    run_dir.mkdir(parents=True, exist_ok=True)
    config = make_config(args, run_dir)
    write_json(run_dir / "config.json", config)
    write_json(run_dir / "bounds.json", {"names": PARAM_NAMES, "lb": LB.tolist(), "ub": UB.tolist()})

    eval_jsonl = run_dir / "evaluations.jsonl"
    eval_csv = run_dir / "evaluations.csv"
    best_path = run_dir / "best.json"
    latest_path = run_dir / "latest.json"
    status_path = run_dir / "status.json"
    heartbeat_path = run_dir / "heartbeat.txt"
    summary_path = run_dir / "best_summary.md"

    if not args.resume:
        for path in [eval_jsonl, eval_csv, best_path, latest_path, status_path, heartbeat_path, summary_path]:
            if path.exists():
                path.unlink()

    rng = np.random.default_rng(int(args.seed))
    seeds = seed_bank(config)
    best_record = load_json(best_path) if args.resume and best_path.exists() else None
    best_theta = np.asarray(best_record["theta"], dtype=float) if best_record else None
    best_loss = float(best_record["objective_loss"]) if best_record else math.inf
    sigma = float(args.local_scale) * (UB - LB)
    submitted = count_jsonl(eval_jsonl) if args.resume and eval_jsonl.exists() else 0
    completed = submitted
    active: dict[Any, int] = {}
    t0 = time.perf_counter()
    last_heartbeat = 0.0

    print_header(args, run_dir)
    write_status(
        status_path,
        heartbeat_path,
        run_tag=args.run_tag,
        submitted=submitted,
        completed=completed,
        active=0,
        elapsed=time.perf_counter() - t0,
        best_loss=best_loss,
        best_record=best_record,
        finished=False,
        finish_reason="running",
    )

    with ProcessPoolExecutor(max_workers=int(args.workers)) as pool:
        while True:
            elapsed = time.perf_counter() - t0
            can_submit = elapsed < float(args.budget_sec) and submitted < int(args.max_evals)
            while can_submit and len(active) < int(args.workers):
                submitted += 1
                theta = propose_candidate(
                    submitted,
                    seeds,
                    best_theta,
                    sigma,
                    rng,
                    float(args.global_prob),
                )
                future = pool.submit(evaluate_candidate, submitted, theta.tolist(), config)
                active[future] = submitted
                elapsed = time.perf_counter() - t0
                can_submit = elapsed < float(args.budget_sec) and submitted < int(args.max_evals)

            if not active:
                break

            done, _ = wait(active.keys(), timeout=float(args.poll_sec), return_when=FIRST_COMPLETED)
            if not done:
                now = time.perf_counter()
                if now - last_heartbeat >= float(args.heartbeat_sec):
                    write_status(
                        status_path,
                        heartbeat_path,
                        run_tag=args.run_tag,
                        submitted=submitted,
                        completed=completed,
                        active=len(active),
                        elapsed=now - t0,
                        best_loss=best_loss,
                        best_record=best_record,
                        finished=False,
                        finish_reason="running",
                    )
                    last_heartbeat = now
                continue

            for future in done:
                active.pop(future)
                result = future.result()
                completed += 1
                result["run_tag"] = args.run_tag
                result["completed_id"] = completed
                result["total_elapsed_sec"] = float(time.perf_counter() - t0)
                objective_loss = float(result.get("objective_loss", math.inf))
                improved = objective_loss < best_loss
                result["improved"] = bool(improved)
                append_jsonl(eval_jsonl, result)
                append_csv(eval_csv, result)
                write_json(latest_path, result)
                if improved:
                    best_record = result
                    best_loss = objective_loss
                    best_theta = np.asarray(result["theta"], dtype=float)
                    sigma = np.minimum(0.35 * (UB - LB), float(args.grow) * sigma)
                    write_json(best_path, result)
                    write_best_summary(summary_path, result)
                elif completed % int(args.stall_window) == 0:
                    sigma = np.maximum(float(args.min_local_scale) * (UB - LB), float(args.shrink) * sigma)

                if improved or completed == 1 or completed % int(args.log_every) == 0:
                    print_progress(completed, submitted, result, best_loss, improved)

                write_status(
                    status_path,
                    heartbeat_path,
                    run_tag=args.run_tag,
                    submitted=submitted,
                    completed=completed,
                    active=len(active),
                    elapsed=time.perf_counter() - t0,
                    best_loss=best_loss,
                    best_record=best_record,
                    finished=False,
                    finish_reason="running",
                )

            if time.perf_counter() - t0 >= float(args.budget_sec) and not active:
                break
            if submitted >= int(args.max_evals) and not active:
                break

    finish_reason = "max_evals" if submitted >= int(args.max_evals) else "time_budget"
    write_status(
        status_path,
        heartbeat_path,
        run_tag=args.run_tag,
        submitted=submitted,
        completed=completed,
        active=0,
        elapsed=time.perf_counter() - t0,
        best_loss=best_loss,
        best_record=best_record,
        finished=True,
        finish_reason=finish_reason,
    )
    print(f"finished {args.run_tag}: completed={completed}, best_objective_loss={best_loss:.6g}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tag", default=time.strftime("v5_nz5_global_%Y%m%d_%H%M%S"))
    parser.add_argument("--results-dir", type=Path, default=Path("global_search_v5_nz5"))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--budget-sec", type=float, default=8 * 60 * 60)
    parser.add_argument("--max-evals", type=int, default=1000000)
    parser.add_argument("--seed", type=int, default=20260522)
    parser.add_argument("--nb", type=int, default=30)
    parser.add_argument("--nz", type=int, default=5)
    parser.add_argument("--rho-z", type=float, default=0.95)
    parser.add_argument("--sigma-z", type=float, default=0.35)
    parser.add_argument("--target-city-entry-prob", type=float, default=0.9)
    parser.add_argument("--kappa-entry", type=float, default=1_000_000.0)
    parser.add_argument("--max-iter-eq", type=int, default=35)
    parser.add_argument("--tol-eq", type=float, default=5e-4)
    parser.add_argument("--global-prob", type=float, default=0.75)
    parser.add_argument("--local-scale", type=float, default=0.18)
    parser.add_argument("--min-local-scale", type=float, default=0.02)
    parser.add_argument("--grow", type=float, default=1.04)
    parser.add_argument("--shrink", type=float, default=0.72)
    parser.add_argument("--stall-window", type=int, default=24)
    parser.add_argument("--eq-penalty", type=float, default=500.0)
    parser.add_argument("--scale-penalty", type=float, default=500.0)
    parser.add_argument("--poll-sec", type=float, default=5.0)
    parser.add_argument("--heartbeat-sec", type=float, default=60.0)
    parser.add_argument("--log-every", type=int, default=10)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def make_config(args: argparse.Namespace, run_dir: Path) -> dict[str, Any]:
    return {
        "run_tag": args.run_tag,
        "run_dir": str(run_dir),
        "nb": int(args.nb),
        "nz": int(args.nz),
        "rho_z": float(args.rho_z),
        "sigma_z": float(args.sigma_z),
        "target_city_entry_prob": float(args.target_city_entry_prob),
        "kappa_entry": float(args.kappa_entry),
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "eq_penalty": float(args.eq_penalty),
        "scale_penalty": float(args.scale_penalty),
        "param_names": PARAM_NAMES,
        "lb": LB.tolist(),
        "ub": UB.tolist(),
        "moment_order": MOMENT_ORDER,
        "objective": "full 19-moment weighted SMM loss plus GE/scale penalties",
    }


def append_csv(path: Path, result: dict[str, Any]) -> None:
    exists = path.exists()
    row = flat_result(result)
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def flat_result(result: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "eval_id": result.get("eval_id"),
        "completed_id": result.get("completed_id"),
        "worker_pid": result.get("worker_pid"),
        "ok": result.get("ok"),
        "solve_ok": result.get("solve_ok"),
        "accepted": result.get("accepted"),
        "strict_converged": result.get("strict_converged"),
        "convergence_reason": result.get("convergence_reason"),
        "iterations_completed": result.get("iterations_completed"),
        "elapsed_sec": result.get("elapsed_sec"),
        "total_elapsed_sec": result.get("total_elapsed_sec"),
        "objective_loss": result.get("objective_loss"),
        "moment_loss": result.get("moment_loss"),
        "penalty": result.get("penalty"),
        "improved": result.get("improved"),
        "final_eq_error": result.get("final_eq_error"),
        "p_P": result.get("p_eq", [np.nan, np.nan])[0] if isinstance(result.get("p_eq"), list) else np.nan,
        "p_C": result.get("p_eq", [np.nan, np.nan])[1] if isinstance(result.get("p_eq"), list) else np.nan,
        "scale_factor": result.get("scale", {}).get("scale_factor") if isinstance(result.get("scale"), dict) else np.nan,
        "outside_entry_flow": result.get("scale", {}).get("outside_entry_flow") if isinstance(result.get("scale"), dict) else np.nan,
        "outside_value": result.get("scale", {}).get("outside_value") if isinstance(result.get("scale"), dict) else np.nan,
        "error": result.get("error"),
    }
    params = result.get("parameters", {})
    if isinstance(params, dict):
        for name in PARAM_NAMES:
            row[f"param_{name}"] = params.get(name, np.nan)
    moments = result.get("moments", {})
    if isinstance(moments, dict):
        for name in MOMENT_ORDER:
            row[name] = moments.get(name, np.nan)
    return row


def write_best_summary(path: Path, result: dict[str, Any]) -> None:
    moments = result.get("moments", {})
    targets = result.get("targets", {})
    params = result.get("parameters", {})
    lines = [
        "# V5 Nz=5 Global Search Best",
        "",
        f"- eval id: `{result.get('eval_id')}`",
        f"- objective loss: `{float(result.get('objective_loss', math.inf)):.6g}`",
        f"- moment loss: `{float(result.get('moment_loss', math.inf)):.6g}`",
        f"- penalty: `{float(result.get('penalty', math.inf)):.6g}`",
        f"- accepted: `{result.get('accepted')}`",
        f"- final GE error: `{float(result.get('final_eq_error', math.nan)):.6g}`",
        f"- solve seconds: `{float(result.get('elapsed_sec', math.nan)):.2f}`",
        f"- prices: `{result.get('p_eq')}`",
        "",
        "## Parameters",
        "",
        "| Parameter | Value |",
        "|---|---:|",
    ]
    for name in PARAM_NAMES:
        lines.append(f"| `{name}` | {float(params.get(name, math.nan)):.6g} |")
    lines.extend(["", "## Moments", "", "| Moment | Target | Model | Gap |", "|---|---:|---:|---:|"])
    for name in MOMENT_ORDER:
        target = float(targets.get(name, math.nan))
        model = float(moments.get(name, math.nan))
        lines.append(f"| `{name}` | {target:.6g} | {model:.6g} | {model - target:.6g} |")
    path.write_text("\n".join(lines) + "\n")


def write_status(
    status_path: Path,
    heartbeat_path: Path,
    *,
    run_tag: str,
    submitted: int,
    completed: int,
    active: int,
    elapsed: float,
    best_loss: float,
    best_record: dict[str, Any] | None,
    finished: bool,
    finish_reason: str,
) -> None:
    status = {
        "run_tag": run_tag,
        "submitted": int(submitted),
        "completed": int(completed),
        "active": int(active),
        "elapsed_sec": float(elapsed),
        "best_objective_loss": float(best_loss),
        "best_eval_id": best_record.get("eval_id") if best_record else None,
        "best_moment_loss": best_record.get("moment_loss") if best_record else None,
        "finished": bool(finished),
        "finish_reason": finish_reason,
        "updated_at_epoch": time.time(),
    }
    write_json(status_path, status)
    heartbeat_path.write_text(
        "run_tag={run_tag}\nsubmitted={submitted}\ncompleted={completed}\nactive={active}\n"
        "elapsed_sec={elapsed:.2f}\nbest_objective_loss={best_loss:.6g}\nfinished={finished}\n"
        "finish_reason={finish_reason}\nupdated_at_epoch={updated:.3f}\n".format(
            run_tag=run_tag,
            submitted=submitted,
            completed=completed,
            active=active,
            elapsed=elapsed,
            best_loss=best_loss,
            finished=finished,
            finish_reason=finish_reason,
            updated=time.time(),
        )
    )


def print_header(args: argparse.Namespace, run_dir: Path) -> None:
    expected_low = int(max(1, args.workers) * args.budget_sec / 300)
    expected_high = int(max(1, args.workers) * args.budget_sec / 100)
    print("=" * 72, flush=True)
    print("V5 Nz=5 parallel global search", flush=True)
    print(f"run_tag={args.run_tag}", flush=True)
    print(f"run_dir={run_dir}", flush=True)
    print(f"workers={args.workers} budget_sec={args.budget_sec:.0f} max_evals={args.max_evals}", flush=True)
    print(f"expected completed evals roughly {expected_low}-{expected_high}", flush=True)
    print(f"objective: all {len(MOMENT_ORDER)} moments plus GE/scale penalties", flush=True)
    print("=" * 72, flush=True)


def print_progress(completed: int, submitted: int, result: dict[str, Any], best_loss: float, improved: bool) -> None:
    moments = result.get("moments", {}) if isinstance(result.get("moments"), dict) else {}
    tag = "BEST" if improved else "eval"
    print(
        f"[{tag} completed={completed:04d} submitted={submitted:04d}] "
        f"obj={float(result.get('objective_loss', math.inf)):.6g} "
        f"mom={float(result.get('moment_loss', math.inf)):.6g} "
        f"pen={float(result.get('penalty', math.inf)):.3g} "
        f"best={best_loss:.6g} "
        f"solve={float(result.get('elapsed_sec', math.nan)):.1f}s "
        f"TFR={float(moments.get('tfr', math.nan)):.3f} "
        f"childless={float(moments.get('childless_rate', math.nan)):.3f} "
        f"own={float(moments.get('own_rate', math.nan)):.3f} "
        f"grad={float(moments.get('own_gradient', math.nan)):.3f} "
        f"wealth={float(moments.get('young_liquid_wealth_to_income', math.nan)):.3f}",
        flush=True,
    )


def to_builtin(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): to_builtin(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_builtin(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [to_builtin(v) for v in obj.tolist()]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_builtin(obj), indent=2, sort_keys=True))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def append_jsonl(path: Path, obj: dict[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(to_builtin(obj), sort_keys=True) + "\n")


def count_jsonl(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


if __name__ == "__main__":
    main()
