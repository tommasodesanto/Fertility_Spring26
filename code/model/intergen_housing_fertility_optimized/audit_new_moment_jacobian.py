#!/usr/bin/env python3
"""Local identification audit for the provisional one-shot 14-by-14 system.

This is not a calibration.  It perturbs each free parameter in the transformed
search coordinate, recomputes all fourteen empirical auxiliary moments, and
writes a checkpointed normalized Jacobian.  A search must remain disabled
unless the complete audit is finite and full rank at the declared tolerance.
"""

from __future__ import annotations

import argparse
import csv
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

from .calibration import extract_moments
from .calibration_search import DOMAIN, theta_from_unit, unit_from_theta
from .new_moment_profile import (
    NEW_MOMENT_SEED,
    new_moment_overrides,
    new_moment_target_system,
)
from .solver import run_model_cp_dt


# Rows and columns are deliberately paired in the economic mapping proposed in
# the calibration ledger.  Full-system rank, rather than diagonal dominance,
# is the identification gate because several blocks are necessarily joint.
MOMENT_PARAMETER_MAP: tuple[tuple[str, str], ...] = (
    ("aggregate_wealth_to_annual_after_tax_earnings", "beta_annual"),
    ("annual_bequest_flow_to_aggregate_wealth", "theta0"),
    ("old_total_estate_wealth_to_annual_income_p90_p50_7684", "theta1"),
    ("childless_renter_rent_expenditure_slope", "alpha_cons"),
    ("childless_renter_intercept_at_mean_price_model_units", "c_bar_0"),
    ("bottom_quintile_childless_renter_mean_rooms", "h_bar_0"),
    ("one_shot_parity_consumption_coefficient", "c_bar_n"),
    ("housing_increment_0to1", "h_bar_jump"),
    ("prime30_55_parent_3plus_minus_1to2_mean_rooms", "h_bar_n"),
    ("tfr", "psi_child"),
    ("childless_rate", "kappa_fert"),
    ("own_rate", "chi"),
    ("model_feasible_four_year_tenure_brier", "tenure_choice_kappa"),
    ("aggregate_mean_occupied_rooms_18_85", "H0"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--unit-step", type=float, default=0.025)
    parser.add_argument("--rank-relative-tol", type=float, default=1e-6)
    parser.add_argument(
        "--max-parameters",
        type=int,
        default=len(MOMENT_PARAMETER_MAP),
        help="Audit only the first N mapped parameters; intended for loop smoke tests.",
    )
    parser.add_argument("--tight", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    system = new_moment_target_system()
    moment_names = [name for name, _parameter in MOMENT_PARAMETER_MAP]
    if tuple(moment_names) != system.moment_names:
        raise RuntimeError("Jacobian row order differs from the provisional target contract")
    domain_names = [name for name, _lower, _upper, _kind in DOMAIN]
    parameter_names = [parameter for _moment, parameter in MOMENT_PARAMETER_MAP]
    if set(parameter_names) != set(domain_names):
        raise RuntimeError("Jacobian columns differ from the live 14-parameter domain")
    selected_parameters = parameter_names[: max(0, min(int(args.max_parameters), len(parameter_names)))]
    if not selected_parameters:
        raise ValueError("--max-parameters must select at least one parameter")

    config = {
        "status": "identification_diagnostic_not_calibration",
        "tight": bool(args.tight),
        "unit_step": float(args.unit_step),
        "rank_relative_tolerance": float(args.rank_relative_tol),
        "moments": moment_names,
        "parameters": selected_parameters,
        "complete_14_by_14": len(selected_parameters) == len(parameter_names),
        "mapping": [
            {"moment": moment, "assigned_parameter": parameter}
            for moment, parameter in MOMENT_PARAMETER_MAP
        ],
    }
    write_json(outdir / "config.json", config)
    cases_path = outdir / "cases.jsonl"
    cases_path.write_text("")

    baseline_unit = unit_from_theta(NEW_MOMENT_SEED)
    domain_index = {name: index for index, name in enumerate(domain_names)}
    start = time.perf_counter()
    baseline = solve_case("baseline", baseline_unit, args.tight, system)
    append_jsonl(cases_path, baseline)
    write_json(outdir / "baseline.json", baseline)
    print_case(baseline, 1, 1 + 2 * len(selected_parameters))

    sides: dict[str, dict[str, dict[str, Any]]] = {}
    completed = 1
    for parameter in selected_parameters:
        index = domain_index[parameter]
        lower_unit = max(0.0, float(baseline_unit[index]) - float(args.unit_step))
        upper_unit = min(1.0, float(baseline_unit[index]) + float(args.unit_step))
        if not upper_unit > lower_unit:
            raise RuntimeError(f"no transformed perturbation available for {parameter}")
        sides[parameter] = {}
        for side, value in (("minus", lower_unit), ("plus", upper_unit)):
            unit = baseline_unit.copy()
            unit[index] = value
            record = solve_case(f"{parameter}_{side}", unit, args.tight, system)
            record["perturbation"] = {
                "parameter": parameter,
                "side": side,
                "base_unit": float(baseline_unit[index]),
                "case_unit": float(value),
                "unit_step": float(value - baseline_unit[index]),
            }
            sides[parameter][side] = record
            append_jsonl(cases_path, record)
            completed += 1
            write_json(
                outdir / "latest_completed_case.json",
                {"completed": completed, "total": 1 + 2 * len(selected_parameters), **record},
            )
            print_case(record, completed, 1 + 2 * len(selected_parameters))

    matrix, rows = build_jacobian(
        baseline,
        sides,
        selected_parameters,
        moment_names,
        system.targets_dict(),
    )
    write_csv(outdir / "jacobian_long.csv", rows)
    singular = np.linalg.svd(matrix, compute_uv=False)
    largest = float(singular[0]) if singular.size else 0.0
    relative = singular / largest if largest > 0.0 else np.zeros_like(singular)
    rank = int(np.sum(relative > float(args.rank_relative_tol)))
    condition = float(largest / singular[-1]) if singular.size and singular[-1] > 0.0 else math.inf
    singular_rows = [
        {
            "index": index + 1,
            "singular_value": float(value),
            "relative_to_largest": float(rel),
            "rank": rank,
            "condition_number": condition,
        }
        for index, (value, rel) in enumerate(zip(singular, relative))
    ]
    write_csv(outdir / "singular_values.csv", singular_rows)
    mapping_rows = mapping_diagnostics(matrix, moment_names, selected_parameters)
    write_csv(outdir / "mapping_diagnostics.csv", mapping_rows)

    all_records = [baseline] + [sides[p][s] for p in selected_parameters for s in ("minus", "plus")]
    all_solved = all(record["status"] == "ok" for record in all_records)
    all_finite = bool(np.all(np.isfinite(matrix)))
    full = len(selected_parameters) == len(parameter_names)
    column_norms = np.linalg.norm(matrix, axis=0)
    no_zero_columns = bool(np.all(column_norms > 1e-5))
    full_rank = bool(full and rank == len(parameter_names))
    gate_pass = bool(all_solved and all_finite and no_zero_columns and full_rank)
    summary = {
        "status": "pass" if gate_pass else "fail_or_incomplete",
        "gate_pass": gate_pass,
        "complete_14_by_14": full,
        "all_cases_solved": all_solved,
        "all_jacobian_entries_finite": all_finite,
        "no_near_zero_columns": no_zero_columns,
        "rank": rank,
        "required_rank": len(parameter_names),
        "condition_number": condition,
        "singular_values": [float(value) for value in singular],
        "relative_singular_values": [float(value) for value in relative],
        "column_norms": {
            parameter: float(column_norms[index])
            for index, parameter in enumerate(selected_parameters)
        },
        "baseline_loss": float(baseline["loss"]),
        "baseline_market_residual": float(baseline["market_residual"]),
        "elapsed_seconds": float(time.perf_counter() - start),
    }
    write_json(outdir / "summary.json", summary)
    print(
        f"audit {summary['status']}: rank={rank}/{len(parameter_names)}, "
        f"condition={condition:.3g}, elapsed={summary['elapsed_seconds'] / 60.0:.1f} min",
        flush=True,
    )


def solve_case(
    label: str,
    unit: np.ndarray,
    tight: bool,
    system,
) -> dict[str, Any]:
    theta = theta_from_unit(unit)
    start = time.perf_counter()
    try:
        sol, parameters, price = run_model_cp_dt(
            {**new_moment_overrides(tight=tight), **theta},
            verbose=False,
        )
        moments_all = extract_moments(sol, parameters)
        moments = {
            name: float(moments_all.get(name, np.nan))
            for name in system.moment_names
        }
        loss = float(system.loss(moments))
        status = "ok" if all(np.isfinite(value) for value in moments.values()) else "nonfinite_moment"
        residual = float(getattr(sol, "best_max_abs_rel_excess", np.nan))
        price_value = float(np.asarray(price, dtype=float).reshape(-1)[0])
    except Exception as error:  # Keep failed perturbations in the audit packet.
        moments = {name: math.nan for name in system.moment_names}
        loss = math.inf
        status = f"failed: {type(error).__name__}: {error}"
        residual = math.inf
        price_value = math.nan
    return {
        "case": label,
        "status": status,
        "theta": theta,
        "unit": [float(value) for value in unit],
        "moments": moments,
        "loss": loss,
        "market_residual": residual,
        "price": price_value,
        "elapsed_seconds": float(time.perf_counter() - start),
    }


def build_jacobian(
    baseline: dict[str, Any],
    sides: dict[str, dict[str, dict[str, Any]]],
    parameters: list[str],
    moments: list[str],
    targets: dict[str, float],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    matrix = np.full((len(moments), len(parameters)), np.nan)
    rows: list[dict[str, Any]] = []
    for column, parameter in enumerate(parameters):
        minus = sides[parameter]["minus"]
        plus = sides[parameter]["plus"]
        parameter_index = [name for name, _lower, _upper, _kind in DOMAIN].index(parameter)
        denominator = float(plus["unit"][parameter_index]) - float(minus["unit"][parameter_index])
        for row, moment in enumerate(moments):
            minus_value = float(minus["moments"][moment])
            plus_value = float(plus["moments"][moment])
            target = float(targets[moment])
            derivative = (plus_value - minus_value) / denominator / abs(target)
            matrix[row, column] = derivative
            rows.append(
                {
                    "moment": moment,
                    "parameter": parameter,
                    "target": target,
                    "baseline_moment": float(baseline["moments"][moment]),
                    "minus_moment": minus_value,
                    "plus_moment": plus_value,
                    "derivative_of_moment_over_target_wrt_unit": derivative,
                }
            )
    return matrix, rows


def mapping_diagnostics(
    matrix: np.ndarray,
    moments: list[str],
    parameters: list[str],
) -> list[dict[str, Any]]:
    assigned = dict(MOMENT_PARAMETER_MAP)
    rows: list[dict[str, Any]] = []
    for row, moment in enumerate(moments):
        effects = np.abs(matrix[row, :])
        order = np.argsort(-effects)
        top_index = int(order[0])
        assigned_parameter = assigned[moment]
        if assigned_parameter in parameters:
            assigned_index = parameters.index(assigned_parameter)
            assigned_effect = float(effects[assigned_index])
            assigned_rank = int(np.where(order == assigned_index)[0][0]) + 1
        else:
            assigned_effect = math.nan
            assigned_rank = 0
        rows.append(
            {
                "moment": moment,
                "assigned_parameter": assigned_parameter,
                "assigned_abs_effect": assigned_effect,
                "assigned_rank": assigned_rank,
                "top_parameter": parameters[top_index],
                "top_abs_effect": float(effects[top_index]),
                "assigned_share_of_row_l1": assigned_effect / max(float(np.sum(effects)), 1e-300)
                if np.isfinite(assigned_effect)
                else math.nan,
            }
        )
    return rows


def print_case(record: dict[str, Any], completed: int, total: int) -> None:
    print(
        f"case {completed}/{total} {record['case']}: status={record['status']}, "
        f"loss={record['loss']:.4g}, residual={record['market_residual']:.2e}, "
        f"elapsed={record['elapsed_seconds']:.1f}s",
        flush=True,
    )


def append_jsonl(path: Path, payload: Any) -> None:
    with path.open("a") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
