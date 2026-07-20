#!/usr/bin/env python3
"""Run the two-cell sequential-fertility fork smoke.

Run from the repository root with:
``PYTHONPATH=code/model code/model/.venv/bin/python3 -m intergen_eqscale_seq_optimized.run_seq_smoke``.
This driver is intentionally not run by the implementation verification.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from intergen_eqscale_seq_optimized.calibration import diagnostic_loss, extract_moments
from intergen_eqscale_seq_optimized.local_panel import jsonable
from intergen_eqscale_seq_optimized.solver import InfeasibleThetaError, run_model_cp_dt
from tools.run_intergen_bequest_exit_chain import (
    arm_contract,
    common_overrides,
    load_theta,
    target_system,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = ROOT / "output/model/eqscale_seq_fork_smoke_20260718"
M5_RESULTS = ROOT / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json"
GATE0_RESULT = ROOT / "output/model/intergen_m5_draft_refresh_20260717/gate0/gate0_result.json"
TARGET_SET = "candidate_replacement_income_disciplined_v1"

CELLS: tuple[dict[str, Any], ...] = (
    {
        "index": 0,
        "name": "nested_baseline",
        "fecundity_omega1": 0.0,
        "fecundity_omega2": 0.0,
        "fecundity_terminal_age": 45.0,
        "label": "nested baseline: default fecundity schedule",
    },
    {
        "index": 1,
        "name": "clock_smoke",
        "fecundity_omega1": 0.02,
        "fecundity_omega2": 0.134,
        "fecundity_terminal_age": 45.0,
        "label": "ILLUSTRATIVE fecundity schedule; not calibrated",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", default="0,1", help="Comma-separated cell indices.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def target_fit(moments: dict[str, Any], targets: dict[str, float], weights: dict[str, float]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for name, target in targets.items():
        model = float(moments.get(name, math.nan))
        gap = model - float(target)
        weight = float(weights[name])
        rows.append(
            {
                "moment": name,
                "target": float(target),
                "model": model,
                "gap": gap,
                "weight": weight,
                "loss_contribution": weight * gap * gap,
            }
        )
    return rows


def compare_gate0(loss: float | None, market_residual: float | None) -> dict[str, Any]:
    expected = json.loads(GATE0_RESULT.read_text())
    expected_loss = float(expected["loss"])
    expected_residual = float(expected["market_residual"])
    status = "MATCH" if loss == expected_loss and market_residual == expected_residual else "MISMATCH"
    print(
        f"gate0 {status}: loss current={loss!r} expected={expected_loss!r}; "
        f"market_residual current={market_residual!r} expected={expected_residual!r}",
        flush=True,
    )
    return {
        "status": status,
        "loss_current": loss,
        "loss_expected": expected_loss,
        "market_residual_current": market_residual,
        "market_residual_expected": expected_residual,
    }


def main() -> None:
    cli = parse_args()
    selected = [int(value.strip()) for value in cli.cells.split(",") if value.strip()]
    by_index = {int(cell["index"]): cell for cell in CELLS}
    if not selected or set(selected) - set(by_index):
        raise ValueError("--cells must select a nonempty subset of 0,1")

    args = argparse.Namespace(
        arm="M5", J=17, Nb=120, max_iter_eq=40, tol_eq=2.5e-5,
        ltv_terminal=0.4, theta1=0.25, seed_theta0=0.30, seed_theta_n=0.75,
        seed_kappa=0.0, fixed_theta0=None,
    )
    active, fixed, mechanism = arm_contract(args)
    overrides = {**common_overrides(args, mechanism), "max_iter_eq": 40, "tol_eq": 2.5e-5}
    theta = load_theta(M5_RESULTS, seed_arm="M5")
    targets, weights = target_system(TARGET_SET)
    if len(targets) != 15:
        raise ValueError(f"sequential-fertility smoke requires 15 targets, found {len(targets)}")

    cli.outdir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    with (cli.outdir / "results.jsonl").open("a") as stream:
        for index in selected:
            cell = by_index[index]
            started = time.perf_counter()
            try:
                sol, P, _ = run_model_cp_dt({**overrides, **theta, **{k: cell[k] for k in cell if k.startswith("fecundity_")}}, verbose=False)
                moments = extract_moments(sol, P)
                loss = float(diagnostic_loss(moments, targets=targets, weights=weights))
                record: dict[str, Any] = {
                    "cell": cell,
                    "status": "ok",
                    "elapsed_seconds": float(time.perf_counter() - started),
                    "market_residual": float(sol.best_max_abs_rel_excess),
                    "target_fit": target_fit(moments, targets, weights),
                    "diagnostic_loss": loss,
                    "moments": moments,
                    "timing_moments": {
                        name: getattr(sol, name)
                        for name in (
                            "attempt_hazard_by_age", "first_birth_hazard_by_age",
                            "first_birth_age_distribution", "share_first_births_age30plus",
                            "childless_chosen_45", "childless_clock_45",
                        )
                    },
                    "tfr": float(moments["tfr"]),
                    "childless_rate": float(moments["childless_rate"]),
                    "mean_age_first_birth": float(moments["mean_age_first_birth"]),
                }
            except InfeasibleThetaError as exc:
                record = {
                    "cell": cell, "status": "infeasible", "stage": exc.stage,
                    "dead_mass": exc.dead_mass, "census_head": exc.census[:5],
                    "elapsed_seconds": float(time.perf_counter() - started),
                    "market_residual": None, "target_fit": [], "diagnostic_loss": None,
                    "moments": {}, "timing_moments": {},
                }
            if index == 0:
                record["gate0_comparison"] = compare_gate0(
                    record.get("diagnostic_loss"), record.get("market_residual")
                )
            records.append(record)
            stream.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
            stream.flush()
            print(f"cell {index} {cell['name']}: status={record['status']} elapsed={record['elapsed_seconds']:.3f}s", flush=True)

    final = {
        "selected_cells": selected, "target_set": TARGET_SET,
        "active_parameters": active, "fixed_parameters": fixed,
        "mechanism": mechanism, "results": records,
    }
    (cli.outdir / "results.json").write_text(json.dumps(jsonable(final), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
