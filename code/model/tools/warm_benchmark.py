"""Run two Python model benchmarks in one process to separate JIT warmup cost."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.evaluate import apply_theta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--setup", choices=["fast", "benchmark"], default="fast")
    parser.add_argument("--theta", default="", help="Use 'x0' to apply calibration x0 before solving.")
    parser.add_argument("--max-iter-eq", type=int, default=1)
    parser.add_argument("--force-full", action="store_true")
    args = parser.parse_args()

    reports = []
    for run in range(2):
        setup = build_calibration_setup(args.setup)
        P = setup.P_base
        if args.theta.strip().lower() == "x0":
            P = apply_theta(P, setup.x0, setup.names)
        P.max_iter_eq = args.max_iter_eq
        if args.force_full:
            P.force_full_bellman = True
        t0 = time.perf_counter()
        sol, _, p_eq = run_model_cp_dt(P, verbose=False)
        reports.append(
            {
                "run": run + 1,
                "elapsed_sec": time.perf_counter() - t0,
                "timings": getattr(sol, "timings", {}),
                "tfr": float(2 * sol.mean_parity),
                "own_rate": float(sol.own_rate),
                "pop_share": [float(x) for x in sol.pop_share],
                "p_eq": [float(x) for x in p_eq],
                "total_mass": float(sol.total_mass),
            }
        )
    print(json.dumps(reports, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
