#!/usr/bin/env python3
"""Run the fixed-theta eight-cell means-tested transfer-floor probe."""

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

import numpy as np

from intergen_housing_fertility.calibration import (
    diagnostic_loss,
    extract_moments,
    get_target_set,
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
    income_at_state,
    income_transition_values,
    run_model_cp_dt,
)
from tools.run_intergen_bequest_exit_chain import (
    MATCHED_ANNUAL_INNOVATION_SD,
    arm_contract,
    common_overrides,
    load_theta,
    survival_schedule,
    target_system,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = ROOT / "output/model/transfer_floor_probe_20260718"
M5_RESULTS = ROOT / "output/model/intergen_income_disciplined_recalibration_20260716/report/results.json"
GATE0_RESULT = ROOT / "output/model/intergen_m5_draft_refresh_20260717/gate0/gate0_result.json"
TARGET_SET = "candidate_replacement_income_disciplined_v1"
ANNUAL_RHO = 0.9601845894041878
PERIOD_YEARS = 4.0


CELLS: tuple[dict[str, Any], ...] = (
    {
        "index": 0,
        "name": "base_repro",
        "sigma_annual": None,
        "c_bar_0_period": None,
        "G0_period": 0.0,
        "Gn_period": 0.0,
        "expectation": "must reproduce gate0",
    },
    {
        "index": 1,
        "name": "floor_only",
        "sigma_annual": None,
        "c_bar_0_period": None,
        "G0_period": 0.52,
        "Gn_period": 0.40,
        "expectation": "receipt confined to bottom-income states and poor retirees; moments near baseline",
    },
    {
        "index": 2,
        "name": "lowc_nofloor",
        "sigma_annual": 0.20,
        "c_bar_0_period": 0.40,
        "G0_period": 0.0,
        "Gn_period": 0.0,
        "expectation": "does low c_bar_0 alone restore feasibility?",
    },
    {
        "index": 3,
        "name": "s12_lowc",
        "sigma_annual": 0.12,
        "c_bar_0_period": 0.40,
        "G0_period": 0.52,
        "Gn_period": 0.40,
        "expectation": "main grid",
    },
    {
        "index": 4,
        "name": "s20_lowc",
        "sigma_annual": 0.20,
        "c_bar_0_period": 0.40,
        "G0_period": 0.52,
        "Gn_period": 0.40,
        "expectation": "main grid",
    },
    {
        "index": 5,
        "name": "s12_midc",
        "sigma_annual": 0.12,
        "c_bar_0_period": 0.56,
        "G0_period": 0.72,
        "Gn_period": 0.40,
        "expectation": "main grid",
    },
    {
        "index": 6,
        "name": "s20_midc",
        "sigma_annual": 0.20,
        "c_bar_0_period": 0.56,
        "G0_period": 0.72,
        "Gn_period": 0.40,
        "expectation": "main grid",
    },
    {
        "index": 7,
        "name": "gate_check_smallG",
        "sigma_annual": 0.20,
        "c_bar_0_period": None,
        "G0_period": 0.52,
        "Gn_period": 0.40,
        "expectation": "must raise InfeasibleThetaError",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cells",
        default=",".join(str(cell["index"]) for cell in CELLS),
        help="Comma-separated cell indices (default: all).",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def parse_cell_indices(raw: str) -> list[int]:
    try:
        indices = [int(piece.strip()) for piece in str(raw).split(",") if piece.strip()]
    except ValueError as exc:
        raise ValueError("--cells must be a comma-separated list of integers") from exc
    if not indices:
        raise ValueError("--cells selects no probe cells")
    if len(indices) != len(set(indices)):
        raise ValueError("--cells may not contain duplicate indices")
    valid = {int(cell["index"]) for cell in CELLS}
    invalid = sorted(set(indices) - valid)
    if invalid:
        raise ValueError(f"unknown probe cell indices: {invalid}")
    return indices


def target_fit(
    moments: dict[str, float],
    targets: dict[str, float],
    weights: dict[str, float],
) -> list[dict[str, float | str]]:
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


def receipt_statistics(sol: Any, P: Any) -> dict[str, float]:
    distribution = np.asarray(sol.g, dtype=float)
    if distribution.ndim != 7:
        raise ValueError(f"transfer-floor probe requires rank-7 Markov distribution, got {distribution.shape}")
    b_grid = np.asarray(sol.b_grid, dtype=float)
    Nb, _, I, J, Nz, npar, ncs = distribution.shape
    if b_grid.shape != (Nb,):
        raise ValueError("solution asset grid does not match the distribution")

    z_grid, _, _ = income_transition_values(P)
    if len(z_grid) != Nz:
        raise ValueError("income-state grid does not match the distribution")
    income = np.empty((I, J, Nz), dtype=float)
    for i in range(I):
        for j in range(J):
            for zz, z_value in enumerate(z_grid):
                income[i, j, zz] = income_at_state(P, i, j, float(z_value))

    parity = np.arange(npar, dtype=float)[:, None]
    child_state = np.arange(ncs)[None, :]
    has_children = (child_state >= 1) & (child_state < int(P.n_child_stages) + 1)
    g0 = float(getattr(P, "transfer_floor_G0", 0.0))
    gn = float(getattr(P, "transfer_floor_Gn", 0.0))
    g_bar = g0 + np.where(has_children, gn * parity, 0.0)

    test_resources = (
        float(P.R_gross) * np.maximum(b_grid, 0.0)[:, None, None, None, None, None, None]
        + income[None, None, :, :, :, None, None]
    )
    guarantee = g_bar[None, None, None, None, None, :, :]
    transfer = np.clip(guarantee - test_resources, 0.0, guarantee)
    recipient_mass = distribution * (transfer > 0.0)
    outlays = distribution * transfer
    income_mass = distribution * income[None, None, :, :, :, None, None]

    total_mass = float(np.sum(distribution))
    gross_income_total = float(np.sum(income_mass))
    denom = max(total_mass, 1e-300)
    ages = float(P.age_start) + np.arange(J, dtype=float) * float(P.da)

    def mass_share(value: np.ndarray) -> float:
        return float(np.sum(value) / denom)

    return {
        "receipt_mass_share": mass_share(recipient_mass),
        "outlays_total": float(np.sum(outlays)),
        "outlays_over_income": float(np.sum(outlays) / max(gross_income_total, 1e-300)),
        "receipt_mass_share_renter": mass_share(recipient_mass[:, 0]),
        "receipt_mass_share_owner": mass_share(recipient_mass[:, 1:]),
        "receipt_mass_share_b_negative": mass_share(recipient_mass[b_grid < 0.0]),
        "receipt_mass_share_b_nonnegative": mass_share(recipient_mass[b_grid >= 0.0]),
        "receipt_mass_share_age_18_33": mass_share(recipient_mass[:, :, :, (ages >= 18.0) & (ages <= 33.0)]),
        "receipt_mass_share_age_34_65": mass_share(recipient_mass[:, :, :, (ages >= 34.0) & (ages <= 65.0)]),
        "receipt_mass_share_age_66_plus": mass_share(recipient_mass[:, :, :, ages >= 66.0]),
    }


def cell_parameters(spec: dict[str, Any], theta: dict[str, float]) -> dict[str, Any]:
    c_bar_period = (
        float(theta["c_bar_0"])
        if spec["c_bar_0_period"] is None
        else float(spec["c_bar_0_period"])
    )
    sigma_annual = (
        float(MATCHED_ANNUAL_INNOVATION_SD)
        if spec["sigma_annual"] is None
        else float(spec["sigma_annual"])
    )
    return {
        "index": int(spec["index"]),
        "name": str(spec["name"]),
        "expectation": str(spec["expectation"]),
        "sigma_annual": sigma_annual,
        "rho_annual": ANNUAL_RHO,
        "c_bar_0_annual": c_bar_period / PERIOD_YEARS,
        "c_bar_0_period": c_bar_period,
        "G0_annual": float(spec["G0_period"]) / PERIOD_YEARS,
        "G0_period": float(spec["G0_period"]),
        "Gn_annual": float(spec["Gn_period"]) / PERIOD_YEARS,
        "Gn_period": float(spec["Gn_period"]),
        "income_overrides_reapplied": spec["sigma_annual"] is not None,
    }


def compare_gate0(loss: float | None, market_residual: float | None) -> dict[str, Any]:
    expected = json.loads(GATE0_RESULT.read_text())
    expected_loss = float(expected["loss"])
    expected_residual = float(expected["market_residual"])
    matches = loss == expected_loss and market_residual == expected_residual
    label = "MATCH" if matches else "MISMATCH"
    print(
        f"gate0 {label}: loss current={loss!r} expected={expected_loss!r}; "
        f"market_residual current={market_residual!r} expected={expected_residual!r}",
        flush=True,
    )
    return {
        "status": label,
        "loss_current": loss,
        "loss_expected": expected_loss,
        "market_residual_current": market_residual,
        "market_residual_expected": expected_residual,
    }


def main() -> None:
    cli = parse_args()
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
        raise ValueError(f"transfer-floor probe requires 15 targets, found {len(targets)}")

    cli.outdir.mkdir(parents=True, exist_ok=True)
    results_path = cli.outdir / "results.jsonl"
    records: list[dict[str, Any]] = []
    by_index = {int(cell["index"]): cell for cell in CELLS}

    with results_path.open("a") as stream:
        for index in selected:
            spec = by_index[index]
            params = cell_parameters(spec, theta)
            overrides = {**tight}
            if spec["sigma_annual"] is not None:
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
                moments = {str(name): float(value) for name, value in extract_moments(sol, P).items()}
                loss = float(diagnostic_loss(moments, targets=targets, weights=weights))
                market_residual = float(sol.best_max_abs_rel_excess)
                record: dict[str, Any] = {
                    "cell": params,
                    "status": "ok",
                    "elapsed_seconds": float(time.perf_counter() - started),
                    "market_residual": market_residual,
                    "target_fit": target_fit(moments, targets, weights),
                    "diagnostic_loss": loss,
                    "moments": moments,
                    "entry_censored_share": float(getattr(sol, "entry_censored_share", -1.0)),
                    "receipt_statistics": receipt_statistics(sol, P),
                }
            except InfeasibleThetaError as exc:
                loss = None
                market_residual = None
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

            if index == 0:
                record["gate0_comparison"] = compare_gate0(loss, market_residual)
            records.append(record)
            stream.write(json.dumps(jsonable(record), sort_keys=True) + "\n")
            stream.flush()
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
