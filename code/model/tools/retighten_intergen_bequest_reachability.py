#!/usr/bin/env python3
"""Re-solve one reachability-cell winner with a stricter equilibrium budget."""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from intergen_housing_fertility.calibration import extract_moments
from intergen_housing_fertility.local_panel import jsonable
from intergen_housing_fertility.solver import run_model_cp_dt
from tools.profile_intergen_bequest_reachability import DISPERSION, FAMILY_GAP, MEDIAN
from tools.run_intergen_bequest_exit_chain import arm_contract, common_overrides, target_system


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-dir", type=Path, required=True)
    parser.add_argument("--max-iter-eq", type=int, default=100)
    parser.add_argument("--tol-eq", type=float, default=2.5e-5)
    parser.add_argument("--scalar-market-refine-iter", type=int, default=64)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary_path = args.cell_dir / "summary.json"
    summary: dict[str, Any] = json.loads(summary_path.read_text())
    meta = summary["metadata"]
    source = summary.get("best_search")
    if not isinstance(source, dict) or not isinstance(source.get("theta"), dict):
        raise RuntimeError("cell summary has no search winner")
    theta = {str(k): float(v) for k, v in source["theta"].items() if isinstance(v, (int, float))}
    theta1 = float(meta["theta1"])
    contract = SimpleNamespace(
        arm="M3",
        ltv_terminal=0.0,
        theta1=theta1,
        fixed_theta0=None,
        J=int(meta["J"]),
        Nb=int(meta["Nb"]),
        max_iter_eq=int(args.max_iter_eq),
        tol_eq=float(args.tol_eq),
    )
    _, _, mechanism = arm_contract(contract)
    overrides = common_overrides(contract, mechanism)
    overrides.update(
        max_iter_eq=int(args.max_iter_eq),
        tol_eq=float(args.tol_eq),
        scalar_market_refine_iter=int(args.scalar_market_refine_iter),
        scalar_market_refine_max_expand=12,
    )
    targets, weights = target_system("candidate_replacement_bequest_internal_v1")
    names = (MEDIAN, DISPERSION, FAMILY_GAP)

    records: list[dict[str, Any]] = []
    for repeat in range(2):
        t0 = time.perf_counter()
        sol, P, p_eq = run_model_cp_dt({**overrides, **theta}, verbose=False)
        moments = extract_moments(sol, P)
        residual = float(getattr(sol, "best_max_abs_rel_excess", math.inf))
        timings = dict(getattr(sol, "timings", {}))
        strict = bool(
            timings.get("strict_converged", getattr(sol, "converged", False))
            and math.isfinite(residual)
            and residual <= float(args.tol_eq)
        )
        gaps = {name: float(moments[name]) - float(targets[name]) for name in names}
        loss = float(weights[MEDIAN]) * gaps[MEDIAN] ** 2 + float(weights[FAMILY_GAP]) * gaps[FAMILY_GAP] ** 2
        record = {
            "case": repeat,
            "label": f"retight_repeat_{repeat+1}",
            "status": "ok",
            "strict_converged": strict,
            "profile_loss": loss,
            "market_residual": residual,
            "price": float(np.asarray(p_eq).reshape(-1)[0]),
            "theta": theta,
            "moments": {name: moments[name] for name in names},
            "gaps": gaps,
            "unit_vector": source.get("unit_vector"),
            "origin": {"phase": "retighten", "source_case": source.get("case")},
            "timings": timings,
            "error": "",
            "evaluator": {"max_iter_eq": int(args.max_iter_eq), "tol_eq": float(args.tol_eq)},
            "elapsed_sec": time.perf_counter() - t0,
        }
        records.append(record)
        (args.cell_dir / "retight_latest.json").write_text(json.dumps(jsonable(record), indent=2, sort_keys=True))
        print(
            f"cell={meta['cell']} repeat={repeat+1}/2 strict={strict} "
            f"loss={loss:.6g} resid={residual:.2e}",
            flush=True,
        )
        if not strict:
            raise RuntimeError(f"cell {meta['cell']} did not converge at max_iter_eq={args.max_iter_eq}")

    repeat_check = {
        "both_strict": True,
        "loss_abs_difference": abs(records[0]["profile_loss"] - records[1]["profile_loss"]),
        "max_abs_moment_difference": max(
            abs(float(records[0]["moments"][name]) - float(records[1]["moments"][name]))
            for name in names
        ),
    }
    if repeat_check["loss_abs_difference"] > 1e-10 or repeat_check["max_abs_moment_difference"] > 1e-10:
        raise RuntimeError("retightened repeats are not deterministic")
    best = min(records, key=lambda row: row["profile_loss"])
    (args.cell_dir / "tight_cases.jsonl").write_text(
        "".join(json.dumps(jsonable(row), sort_keys=True) + "\n" for row in records)
    )
    (args.cell_dir / "best_tight.json").write_text(json.dumps(jsonable(best), indent=2, sort_keys=True))
    summary["best_tight"] = best
    summary["tight_repeat_check"] = repeat_check
    summary["retightening"] = {
        "max_iter_eq": int(args.max_iter_eq),
        "tol_eq": float(args.tol_eq),
        "scalar_market_refine_iter": int(args.scalar_market_refine_iter),
        "completed": True,
    }
    summary_path.write_text(json.dumps(jsonable(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
