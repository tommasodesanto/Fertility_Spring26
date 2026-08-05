#!/usr/bin/env python3
"""Run the local E5F strict exact-repeat smoke (no cluster submission)."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SEED = ROOT / "output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-record", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.environ.update(
        E3_L4="1",
        E5="1",
        E5_MATURATION_REPAIR="1",
        E5F="1",
        E2_SEED_RECORD=str(args.seed_record),
        NUMBA_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
    )

    from intergen_eqscale_seq_optimized import run_e1_chain as chain
    from intergen_eqscale_seq_optimized.e5_profile import e5_target_system

    chain.load_runtime()
    target_system = e5_target_system()
    targets, weights = target_system.targets_dict(), target_system.weights_dict()
    runtime_args = SimpleNamespace(J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5)
    overrides = chain.common_overrides(runtime_args)
    seed = chain.build_seed_theta(args.seed_record)
    rows: list[dict[str, object]] = []

    for floor in (0.10, 0.60, 1.80):
        theta = {**seed, "hbar_child_rooms": floor}
        repeats: list[dict[str, object]] = []
        for repeat in (1, 2):
            started = time.perf_counter()
            solution, parameters, price = chain.run_model_cp_dt(
                {**overrides, **theta, "max_iter_eq": 40, "tol_eq": 2.5e-5},
                verbose=False,
            )
            moments = chain.extract_moments(solution, parameters)
            loss = chain.objective_loss(moments, targets, weights, "canonical")
            residual = float(getattr(solution, "best_max_abs_rel_excess", np.inf))
            strict = bool(
                getattr(solution, "timings", {}).get(
                    "strict_converged", getattr(solution, "converged", False)
                )
                and residual <= 2.5e-5
            )
            record = {
                "repeat": repeat,
                "elapsed_sec": time.perf_counter() - started,
                "strict_converged": strict,
                "market_residual": residual,
                "loss": float(loss),
                "price": np.asarray(price, dtype=float).tolist(),
                "moments": {name: float(moments[name]) for name in targets},
            }
            repeats.append(record)
            print(
                f"floor={floor:.2f} repeat={repeat} "
                f"strict={strict} loss={loss:.12g} residual={residual:.3e} "
                f"seconds={record['elapsed_sec']:.2f}",
                flush=True,
            )
        moment_gap = max(
            abs(float(repeats[0]["moments"][name]) - float(repeats[1]["moments"][name]))
            for name in targets
        )
        loss_gap = abs(float(repeats[0]["loss"]) - float(repeats[1]["loss"]))
        if loss_gap != 0.0 or moment_gap != 0.0:
            raise RuntimeError(
                f"floor={floor:g} failed exact repeat: loss gap={loss_gap}, moment gap={moment_gap}"
            )
        if not all(bool(row["strict_converged"]) for row in repeats):
            raise RuntimeError(f"floor={floor:g} failed strict convergence")
        rows.append(
            {
                "hbar_child_rooms": floor,
                "loss_abs_difference": loss_gap,
                "max_abs_moment_difference": moment_gap,
                "repeats": repeats,
            }
        )

    payload = {
        "status": "passed",
        "arm": "E5F",
        "J": 17,
        "Nb": 120,
        "max_iter_eq": 40,
        "tol_eq": 2.5e-5,
        "draws": rows,
    }
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
