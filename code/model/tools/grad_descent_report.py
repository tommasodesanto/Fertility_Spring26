"""Read the gradient_descent.json output and print a clean comparison."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from dt_cp_model.objective import extract_moments
from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.theta import apply_theta


BRIDGE_OVERRIDES = {
    "E_loc": np.array([0.0, 0.1318]),
    "r_bar": np.array([0.04, 0.0762]),
    "H_own": np.array([5.3725, 5.84, 6.3075, 6.875, 7.625, 8.68]),
    "hR_max": 5.1,
}


def solve_and_get_moments(theta_full: np.ndarray, P_base, names, max_iter_eq: int = 120):
    P = apply_theta(P_base, theta_full, names)
    for k, v in BRIDGE_OVERRIDES.items():
        setattr(P, k, v)
    P.max_iter_eq = max_iter_eq
    sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
    return extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True), p_eq, sol


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gd-json", type=Path, default=Path("benchmarks/grad_descent_bench.json"))
    parser.add_argument("--setup", default="benchmark")
    args = parser.parse_args()

    payload = json.loads(args.gd_json.read_text())
    theta_start_full = np.array(payload["history"][0]["theta"])  # this is the search-subset
    theta_best_full = np.array(payload["best_theta_full"])
    theta_start_full_full = np.array(payload["best_theta_full"])  # placeholder
    # Reconstruct full theta at start using the search idx logic
    setup = build_calibration_setup(args.setup)
    name_to_idx = {n: i for i, n in enumerate(setup.names)}
    search_idx = np.array([name_to_idx[n] for n in payload["search_names"]])
    BRIDGE_THETA = np.array([0.935, 0.1226, 0.0692, 0.9, 0.9084, 0.0871, 2.3, 1.0464, 1.5, 5.0, 0.53, 0.25, 2.59])
    theta_start_full_arr = BRIDGE_THETA.copy()
    theta_start_full_arr[search_idx] = np.array(payload["history"][0]["theta"])

    print("Loss trajectory:")
    print(f"{'iter':>4s} {'evals':>6s} {'loss':>10s}")
    for h in payload["history"]:
        print(f"{h['iter']:4d} {h['n_evals']:6d} {h['loss']:10.4f}")
    print()
    print(f"Total evals: {payload['n_evals']}")
    print(f"Wall time:   {payload['elapsed_sec']:.1f}s ({payload['elapsed_sec']/payload['n_evals']:.1f}s/eval)")
    print()

    print("Theta change (search subset only):")
    print(f"{'param':20s} {'start':>10s} {'best':>10s} {'change':>10s}")
    for i, name in enumerate(payload["search_names"]):
        s = float(payload["history"][0]["theta"][i])
        b = float(payload["best_theta_search"][i])
        print(f"{name:20s} {s:10.4f} {b:10.4f} {b-s:+10.4f}")
    print()

    print("Solving model at start and best to compare moments...")
    m_start, _, _ = solve_and_get_moments(theta_start_full_arr, setup.P_base, setup.names)
    m_best, _, _ = solve_and_get_moments(theta_best_full, setup.P_base, setup.names)

    print()
    print(f"{'moment':40s} {'target':>10s} {'start':>10s} {'best':>10s} {'wt*err^2 best':>15s}")
    print("-" * 90)
    total_start = 0.0
    total_best = 0.0
    for tname, tval in setup.targets.items():
        s = float(getattr(m_start, tname, np.nan))
        b = float(getattr(m_best, tname, np.nan))
        denom = abs(tval) if abs(tval) > 1e-6 else 1.0
        contrib_b = setup.weights[tname] * ((b - tval) / denom) ** 2
        contrib_s = setup.weights[tname] * ((s - tval) / denom) ** 2
        total_start += contrib_s
        total_best += contrib_b
        print(f"{tname:40s} {tval:10.4f} {s:10.4f} {b:10.4f} {contrib_b:15.4f}")
    print("-" * 90)
    print(f"{'TOTAL SMM LOSS':40s} {'':>10s} {total_start:10.4f} {total_best:10.4f} {'':>15s}")


if __name__ == "__main__":
    main()
