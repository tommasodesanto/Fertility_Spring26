"""Simple finite-difference gradient descent on the SMM objective.

Skips the geography inversion: uses fixed bridge-benchmark E_C, r_bar_C.
Targets the live calibration table directly.

Usage:
    PYTHONPATH=. .venv/bin/python tools/gradient_descent.py [--max-iters N]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.theta import apply_theta


BRIDGE_THETA = np.array([
    0.935, 0.1226, 0.0692, 0.9000, 0.9084, 0.0871, 2.3, 1.0464, 1.5, 5.0, 0.53, 0.25, 2.59,
])

BRIDGE_OVERRIDES = {
    "E_loc": np.array([0.0, 0.1318]),
    "r_bar": np.array([0.04, 0.0762]),
    "H_own": np.array([5.3725, 5.84, 6.3075, 6.875, 7.625, 8.68]),
    "hR_max": 5.1,
}

# Subset of theta dims to search; the rest stay at BRIDGE_THETA values.
# Names (full theta order): beta, b_entry_fixed, psi_child, h_bar_jump, h_bar_n,
#                           c_bar_n, kappa_fert, chi, kappa_loc, mu_move,
#                           theta0, theta_n, h_bar_0
SEARCH_NAMES = ["beta", "psi_child", "h_bar_jump", "h_bar_n", "c_bar_n", "kappa_fert", "chi"]
LB = np.array([0.90, 0.00, 0.30, 0.30, 0.02, 1.0, 0.95])
UB = np.array([0.97, 0.20, 2.00, 1.50, 0.30, 5.0, 1.20])


def evaluator(P_base, names, targets, weights, max_iter_eq=120):
    """Build a closure that takes a theta dict and returns (loss, moments)."""
    name_to_idx = {n: i for i, n in enumerate(names)}
    search_idx = np.array([name_to_idx[n] for n in SEARCH_NAMES])

    def f(theta_search: np.ndarray) -> tuple[float, SimpleNamespace, np.ndarray]:
        theta_full = BRIDGE_THETA.copy()
        theta_full[search_idx] = theta_search
        P = apply_theta(P_base, theta_full, names)
        for k, v in BRIDGE_OVERRIDES.items():
            setattr(P, k, v)
        P.max_iter_eq = max_iter_eq
        t0 = time.perf_counter()
        try:
            sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
        except Exception:
            return 1e6, SimpleNamespace(), np.array([np.nan, np.nan])
        eval_time = time.perf_counter() - t0
        moments = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)
        loss = 0.0
        for tname, tval in targets.items():
            mv = float(getattr(moments, tname, np.nan))
            if not np.isfinite(mv):
                return 1e6, moments, p_eq
            denom = abs(tval) if abs(tval) > 1e-6 else 1.0
            dev = (mv - tval) / denom
            loss += weights[tname] * dev * dev
        moments._eval_time = eval_time
        return float(loss), moments, p_eq

    return f


def gradient_descent(f_eval, theta0, lb, ub, max_iters=15, lr0=0.05, fd_rel=0.04, verbose=True):
    """Forward-FD gradient descent with backtracking line search and box constraints."""
    n = theta0.size
    theta = theta0.copy()
    f, m0, _ = f_eval(theta)
    eval_t = getattr(m0, "_eval_time", float("nan"))
    n_evals = 1
    lr = lr0
    history = [{"iter": 0, "n_evals": n_evals, "loss": f, "theta": theta.copy()}]
    if verbose:
        print(f"iter 0  evals=  1  loss={f:8.4f}  t/eval={eval_t:5.1f}s  theta={np.array2string(theta, precision=4)}", flush=True)

    for it in range(1, max_iters + 1):
        # Forward finite-difference gradient (n evals)
        h = fd_rel * np.maximum(np.abs(theta), 1.0) * (ub - lb) / 2.0
        h = np.clip(h, 1e-4, None)
        grad = np.zeros(n)
        for i in range(n):
            theta_p = theta.copy()
            theta_p[i] = min(theta[i] + h[i], ub[i])
            actual_h = theta_p[i] - theta[i]
            if actual_h <= 0:
                theta_p[i] = max(theta[i] - h[i], lb[i])
                actual_h = theta_p[i] - theta[i]
            f_p, m_p, _ = f_eval(theta_p)
            n_evals += 1
            grad[i] = (f_p - f) / actual_h if abs(actual_h) > 1e-12 else 0.0
            if verbose:
                eval_t = getattr(m_p, "_eval_time", float("nan"))
                print(f"  fd[{i}]={grad[i]:+9.2f}  t/eval={eval_t:5.1f}s", flush=True)
        gnorm = np.linalg.norm(grad)
        if gnorm < 1e-6:
            if verbose:
                print(f"iter {it}: gradient ~ 0, stop.")
            break

        scale = (ub - lb) / 2.0
        d_raw = -grad
        d_scaled = d_raw / gnorm
        d = scale * d_scaled / max(np.abs(scale * d_scaled).max(), 1e-12)

        accepted = False
        trial_lr = lr
        for _ls in range(5):
            theta_new = np.clip(theta + trial_lr * d, lb, ub)
            f_new, m_new, _ = f_eval(theta_new)
            n_evals += 1
            if verbose:
                eval_t = getattr(m_new, "_eval_time", float("nan"))
                print(f"  ls lr={trial_lr:.4f}  loss={f_new:8.4f}  t/eval={eval_t:5.1f}s", flush=True)
            if f_new < f - 1e-6:
                theta = theta_new
                f = f_new
                lr = min(trial_lr * 1.4, lr0 * 4.0)
                accepted = True
                break
            trial_lr *= 0.5

        history.append({"iter": it, "n_evals": n_evals, "loss": f, "theta": theta.copy(),
                        "grad_norm": float(gnorm), "lr": float(trial_lr if accepted else 0.0)})
        if verbose:
            tag = "*" if accepted else " "
            print(f"iter {it:2d}{tag} evals={n_evals:3d}  loss={f:8.4f}  lr={trial_lr:.4f}  "
                  f"theta={np.array2string(theta, precision=4)}", flush=True)
        if not accepted:
            if verbose:
                print(f"iter {it}: line search failed, stopping.")
            break

    return theta, f, history, n_evals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-iters", type=int, default=12)
    parser.add_argument("--lr0", type=float, default=0.05)
    parser.add_argument("--max-iter-eq", type=int, default=120)
    parser.add_argument("--setup", default="benchmark")
    parser.add_argument("--out", type=Path, default=Path("benchmarks/grad_descent_bridge.json"))
    args = parser.parse_args()

    setup = build_calibration_setup(args.setup)
    f_eval = evaluator(setup.P_base, setup.names, setup.targets, setup.weights, max_iter_eq=args.max_iter_eq)

    name_to_idx = {n: i for i, n in enumerate(setup.names)}
    search_idx = np.array([name_to_idx[n] for n in SEARCH_NAMES])
    theta0 = BRIDGE_THETA[search_idx].copy()
    lb = np.maximum(LB, theta0 - 0.5 * (UB - LB))  # local box for safety
    ub = np.minimum(UB, theta0 + 0.5 * (UB - LB))

    print(f"Search dims (subset): {SEARCH_NAMES}")
    print(f"Start theta:          {theta0}")
    print(f"Lower bounds:         {lb}")
    print(f"Upper bounds:         {ub}")
    print()

    t0 = time.perf_counter()
    theta_best, f_best, history, n_evals = gradient_descent(
        f_eval, theta0, lb, ub, max_iters=args.max_iters, lr0=args.lr0
    )
    elapsed = time.perf_counter() - t0

    print()
    print(f"Done. Total evals = {n_evals}, wall time = {elapsed:.1f}s ({elapsed/n_evals:.2f}s/eval)")
    print(f"Best loss: {f_best:.4f}  (started at {history[0]['loss']:.4f})")
    print(f"Best theta: {dict(zip(SEARCH_NAMES, theta_best.round(4).tolist()))}")

    # Reconstruct full theta and dump moments at best
    full_best = BRIDGE_THETA.copy()
    full_best[search_idx] = theta_best
    f_chk, moments, p_eq = f_eval(theta_best)
    print()
    print(f"{'moment':40s} {'target':>10s} {'best':>10s} {'start':>10s}")
    print("-" * 75)
    f_start, m_start, _ = f_eval(theta0)  # one extra eval to print start moments
    for tname, tval in setup.targets.items():
        b = float(getattr(moments, tname, np.nan))
        s = float(getattr(m_start, tname, np.nan))
        print(f"{tname:40s} {tval:10.4f} {b:10.4f} {s:10.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "search_names": SEARCH_NAMES,
        "best_theta_search": theta_best.tolist(),
        "best_theta_full": full_best.tolist(),
        "best_loss": float(f_best),
        "start_loss": float(history[0]["loss"]),
        "n_evals": int(n_evals),
        "elapsed_sec": float(elapsed),
        "history": [
            {"iter": h["iter"], "n_evals": h["n_evals"], "loss": float(h["loss"]),
             "theta": h["theta"].tolist()}
            for h in history
        ],
    }
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
