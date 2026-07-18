#!/usr/bin/env python3
"""Run the fixed-M5 income-risk, transfer-floor, and subsistence frontier."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from intergen_housing_fertility.calibration import (
    diagnostic_loss,
    extract_moments,
)
from intergen_housing_fertility.local_panel import (
    income_process_overrides,
    jsonable,
)
from intergen_housing_fertility.production_profile import (
    PRODUCTION_J,
    PRODUCTION_SEARCH_NB,
)
from intergen_housing_fertility.solver import (
    InfeasibleThetaError,
    run_model_cp_dt,
)
from tools.run_intergen_bequest_exit_chain import (
    arm_contract,
    common_overrides,
    load_theta,
    target_system,
)
from tools.run_intergen_transfer_floor_probe import (
    receipt_statistics,
    target_fit,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = ROOT / "output/model/income_feasibility_frontier_20260718"
M5_RESULTS = ROOT / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json"
TARGET_SET = "candidate_replacement_income_disciplined_v1"
ANNUAL_RHO = 0.9601845894041878
PERIOD_YEARS = 4.0

SIGMAS: list[float] = [0.0645, 0.09, 0.12, 0.15, 0.18, 0.20, 0.22]
VARIANTS: tuple[dict[str, Any], ...] = (
    {
        "variant": "A",
        "G0_period": 0.0,
        "Gn_period": 0.0,
        "c_bar_0_period": None,
    },
    {
        "variant": "B",
        "G0_period": 0.0,
        "Gn_period": 0.0,
        "c_bar_0_period": 0.40,
    },
    {
        "variant": "C",
        "G0_period": 0.52,
        "Gn_period": 0.40,
        "c_bar_0_period": None,
    },
    {
        "variant": "D",
        "G0_period": 0.52,
        "Gn_period": 0.40,
        "c_bar_0_period": 0.40,
    },
)
PROBE_CROSS_CHECKS: dict[tuple[str, float], dict[str, Any]] = {
    ("A", 0.0645): {
        "probe_cell": 0,
        "expectation": "gate0 bitwise",
    },
    ("C", 0.0645): {
        "probe_cell": 1,
        "expectation": "same solution, zero receipt",
    },
    ("D", 0.20): {
        "probe_cell": 4,
        "expectation": "same solution",
    },
}


def sigma_label(sigma_annual: float) -> str:
    return "0.0645" if sigma_annual == 0.0645 else f"{sigma_annual:.2f}"


def build_cells() -> tuple[dict[str, Any], ...]:
    cells: list[dict[str, Any]] = []
    for variant in VARIANTS:
        for sigma_annual in SIGMAS:
            variant_name = str(variant["variant"])
            cells.append(
                {
                    "index": len(cells),
                    "name": f"{variant_name}_s{sigma_label(sigma_annual)}",
                    "variant": variant_name,
                    "sigma_annual": sigma_annual,
                    "G0_period": float(variant["G0_period"]),
                    "Gn_period": float(variant["Gn_period"]),
                    "c_bar_0_period": variant["c_bar_0_period"],
                    "probe_cross_check": PROBE_CROSS_CHECKS.get(
                        (variant_name, sigma_annual)
                    ),
                }
            )
    return tuple(cells)


CELLS = build_cells()
CELL_BY_INDEX = {int(cell["index"]): cell for cell in CELLS}
CELL_INDEX = {
    (str(cell["variant"]), float(cell["sigma_annual"])): int(cell["index"])
    for cell in CELLS
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cells",
        default=",".join(str(cell["index"]) for cell in CELLS),
        help="Comma-separated cell indices (default: all).",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--list",
        action="store_true",
        help="List the 28-cell frontier and exit without solving.",
    )
    return parser.parse_args()


def parse_cell_indices(raw: str) -> list[int]:
    try:
        indices = [int(piece.strip()) for piece in str(raw).split(",") if piece.strip()]
    except ValueError as exc:
        raise ValueError("--cells must be a comma-separated list of integers") from exc
    if not indices:
        raise ValueError("--cells selects no frontier cells")
    if len(indices) != len(set(indices)):
        raise ValueError("--cells may not contain duplicate indices")
    invalid = sorted(set(indices) - set(CELL_BY_INDEX))
    if invalid:
        raise ValueError(f"unknown frontier cell indices: {invalid}")
    return indices


def list_cells() -> None:
    for cell in CELLS:
        c_bar = (
            "theta"
            if cell["c_bar_0_period"] is None
            else f"{float(cell['c_bar_0_period']):.2f}"
        )
        income = "base" if float(cell["sigma_annual"]) == SIGMAS[0] else "reapplied"
        print(
            f"{int(cell['index']):2d} -> (variant={cell['variant']}, "
            f"sigma={sigma_label(float(cell['sigma_annual']))}, "
            f"G0={float(cell['G0_period']):.2f}, "
            f"Gn={float(cell['Gn_period']):.2f}, "
            f"c_bar_0={c_bar}, income={income})"
        )


def cell_parameters(spec: dict[str, Any], theta: dict[str, float]) -> dict[str, Any]:
    c_bar_period = (
        float(theta["c_bar_0"])
        if spec["c_bar_0_period"] is None
        else float(spec["c_bar_0_period"])
    )
    return {
        "index": int(spec["index"]),
        "name": str(spec["name"]),
        "variant": str(spec["variant"]),
        "sigma_annual": float(spec["sigma_annual"]),
        "rho_annual": ANNUAL_RHO,
        "c_bar_0_annual": c_bar_period / PERIOD_YEARS,
        "c_bar_0_period": c_bar_period,
        "G0_annual": float(spec["G0_period"]) / PERIOD_YEARS,
        "G0_period": float(spec["G0_period"]),
        "Gn_annual": float(spec["Gn_period"]) / PERIOD_YEARS,
        "Gn_period": float(spec["Gn_period"]),
        "income_overrides_reapplied": float(spec["sigma_annual"]) != SIGMAS[0],
        "probe_cross_check": spec["probe_cross_check"],
    }


def frontier_entry(record: dict[str, Any]) -> str:
    if record["status"] == "ok":
        return f"ok loss={float(record['diagnostic_loss']):.2f}"
    stage = str(record["stage"]).replace("|", "\\|")
    return f"DEAD@{stage}"


def write_frontier_table(
    path: Path,
    records_by_index: dict[int, dict[str, Any]],
) -> None:
    lines = [
        "# Income feasibility frontier",
        "",
        "| sigma_annual | A | B | C | D |",
        "| ---: | :--- | :--- | :--- | :--- |",
    ]
    for sigma_annual in SIGMAS:
        entries = []
        for variant in ("A", "B", "C", "D"):
            index = CELL_INDEX[(variant, sigma_annual)]
            record = records_by_index.get(index)
            entries.append("—" if record is None else frontier_entry(record))
        lines.append(
            f"| {sigma_label(sigma_annual)} | " + " | ".join(entries) + " |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    cli = parse_args()
    if cli.list:
        list_cells()
        return

    selected = parse_cell_indices(cli.cells)
    args = argparse.Namespace(
        arm="M5",
        J=PRODUCTION_J,
        Nb=PRODUCTION_SEARCH_NB,
        max_iter_eq=40,
        tol_eq=2.5e-5,
        ltv_terminal=0.4,
        theta1=0.25,
        seed_theta0=0.30,
        seed_theta_n=0.75,
        seed_kappa=0.0,
        fixed_theta0=None,
    )
    active, fixed, mechanism = arm_contract(args)
    base = common_overrides(args, mechanism)
    tight = {**base, "max_iter_eq": 40, "tol_eq": 2.5e-5}
    theta = load_theta(M5_RESULTS, seed_arm="M5")
    targets, weights = target_system(TARGET_SET)
    if len(targets) != 15:
        raise ValueError(f"income frontier requires 15 targets, found {len(targets)}")

    cli.outdir.mkdir(parents=True, exist_ok=True)
    results_path = cli.outdir / "results.jsonl"
    frontier_path = cli.outdir / "frontier_table.md"
    records: list[dict[str, Any]] = []
    records_by_index: dict[int, dict[str, Any]] = {}

    with results_path.open("a") as stream:
        for index in selected:
            spec = CELL_BY_INDEX[index]
            params = cell_parameters(spec, theta)
            overrides = {**tight}
            if float(spec["sigma_annual"]) != SIGMAS[0]:
                overrides.update(
                    income_process_overrides(
                        5,
                        "rouwenhorst",
                        float(spec["sigma_annual"]),
                        ANNUAL_RHO,
                    )
                )
            overrides["transfer_floor_G0"] = float(spec["G0_period"])
            overrides["transfer_floor_Gn"] = float(spec["Gn_period"])
            theta_cell = {**theta}
            if spec["c_bar_0_period"] is not None:
                theta_cell["c_bar_0"] = float(spec["c_bar_0_period"])

            started = time.perf_counter()
            try:
                sol, P, _ = run_model_cp_dt({**overrides, **theta_cell}, verbose=False)
                moments = {
                    str(name): float(value)
                    for name, value in extract_moments(sol, P).items()
                }
                loss = float(diagnostic_loss(moments, targets=targets, weights=weights))
                record: dict[str, Any] = {
                    "cell": params,
                    "status": "ok",
                    "elapsed_seconds": float(time.perf_counter() - started),
                    "market_residual": float(sol.best_max_abs_rel_excess),
                    "target_fit": target_fit(moments, targets, weights),
                    "diagnostic_loss": loss,
                    "moments": moments,
                    "entry_censored_share": float(
                        getattr(sol, "entry_censored_share", -1.0)
                    ),
                    "receipt_statistics": receipt_statistics(sol, P),
                }
            except InfeasibleThetaError as exc:
                loss = None
                record = {
                    "cell": params,
                    "status": "infeasible",
                    "stage": exc.stage,
                    "dead_mass": exc.dead_mass,
                    "census_head": exc.census[:5],
                    "elapsed_seconds": float(time.perf_counter() - started),
                    "market_residual": None,
                    "target_fit": [],
                    "diagnostic_loss": None,
                    "moments": {},
                    "entry_censored_share": -1.0,
                    "receipt_statistics": {},
                }

            records.append(record)
            records_by_index[index] = record
            stream.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
            stream.flush()
            write_frontier_table(frontier_path, records_by_index)
            loss_text = "n/a" if loss is None else repr(loss)
            print(
                f"cell {index} {spec['name']}: status={record['status']} "
                f"loss={loss_text} elapsed={record['elapsed_seconds']:.3f}s",
                flush=True,
            )

    final = {
        "selected_cells": selected,
        "target_set": TARGET_SET,
        "active_parameters": active,
        "fixed_parameters": fixed,
        "mechanism": mechanism,
        "results": records,
    }
    (cli.outdir / "results.json").write_text(
        json.dumps(jsonable(final), indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
