#!/usr/bin/env python3
"""Run the funded extension of the established Phase-9b entry/scale exercise.

Phase 9b recovers the baseline outside-entry objects, holds them fixed across
policies, and clears the housing market after scaling demand by the endogenous
population mass. This driver preserves that protocol and jointly solves the
house price and lump-sum rebate required for fiscal balance.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from intergen_housing_fertility_optimized.calibration import extract_moments
from intergen_housing_fertility_optimized.solver import (
    entry_wealth_grid_weights,
    income_transition_values,
    run_model_cp_dt,
)
from run_intergen_funded_property_tax_test import (
    case_overrides,
    jsonable,
    load_profile,
    write_csv,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = (
    ROOT
    / "output/model/intergen_new_moment_unrestricted_overnight_20260723_w2/report/results.json"
)
DEFAULT_FIXED_DIR = ROOT / "output/model/intergen_current_m_figures_policy_20260724/funded_policy"
DEFAULT_OUTDIR = ROOT / "output/model/intergen_current_m_figures_policy_20260724/funded_policy_entry"
TARGET_QBAR = 0.5
ENTRY_TASTE_SCALE = 2.0
INFEASIBLE_VALUE_CUTOFF = -1.0e9


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--fixed-dir", type=Path, default=DEFAULT_FIXED_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def safe_logistic(value: float) -> float:
    clipped = float(np.clip(value, -60.0, 60.0))
    return float(1.0 / (1.0 + math.exp(-clipped)))


def node_entry_probabilities(
    value_function: np.ndarray,
    wealth_grid: np.ndarray,
    parameters: Any,
    *,
    outside_value: float,
    taste_scale: float,
) -> SimpleNamespace:
    values = np.asarray(value_function, dtype=float)
    grid = np.asarray(wealth_grid, dtype=float).reshape(-1)
    rows: list[dict[str, Any]] = []
    qbar = 0.0
    feasible_weight = 0.0
    if values.ndim == 7:
        z_grid, z_weights, _ = income_transition_values(parameters)
        for z_index, (z_value, z_weight) in enumerate(zip(z_grid, z_weights)):
            indices, node_weights = entry_wealth_grid_weights(
                grid,
                parameters,
                i=0,
                j=0,
                z_value=float(z_value),
            )
            for node_index, node_weight in zip(indices, node_weights):
                value = float(values[int(node_index), 0, 0, 0, z_index, 0, 0])
                feasible = bool(np.isfinite(value) and value > INFEASIBLE_VALUE_CUTOFF)
                probability = (
                    safe_logistic(float(taste_scale) * (value - float(outside_value)))
                    if feasible
                    else 0.0
                )
                combined_weight = float(z_weight) * float(node_weight)
                qbar += combined_weight * probability
                feasible_weight += combined_weight * float(feasible)
                rows.append(
                    {
                        "z_index": int(z_index),
                        "z": float(z_value),
                        "z_weight": float(z_weight),
                        "wealth_index": int(node_index),
                        "wealth": float(grid[int(node_index)]),
                        "node_weight": float(node_weight),
                        "combined_weight": combined_weight,
                        "value": value,
                        "feasible": feasible,
                        "entry_probability": probability,
                    }
                )
    elif values.ndim == 6:
        indices, node_weights = entry_wealth_grid_weights(
            grid,
            parameters,
            i=0,
            j=0,
            z_value=1.0,
        )
        for node_index, node_weight in zip(indices, node_weights):
            value = float(values[int(node_index), 0, 0, 0, 0, 0])
            feasible = bool(np.isfinite(value) and value > INFEASIBLE_VALUE_CUTOFF)
            probability = (
                safe_logistic(float(taste_scale) * (value - float(outside_value)))
                if feasible
                else 0.0
            )
            qbar += float(node_weight) * probability
            feasible_weight += float(node_weight) * float(feasible)
            rows.append(
                {
                    "z_index": 0,
                    "z": 1.0,
                    "z_weight": 1.0,
                    "wealth_index": int(node_index),
                    "wealth": float(grid[int(node_index)]),
                    "node_weight": float(node_weight),
                    "combined_weight": float(node_weight),
                    "value": value,
                    "feasible": feasible,
                    "entry_probability": probability,
                }
            )
    else:
        raise ValueError(f"Unsupported value-function dimension: {values.ndim}")
    return SimpleNamespace(
        qbar=float(qbar),
        feasible_weight=float(feasible_weight),
        rows=rows,
    )


def calibrate_outside_objects(
    solution: Any,
    parameters: Any,
    *,
    target_qbar: float,
    taste_scale: float,
) -> tuple[dict[str, float], list[dict[str, Any]]]:
    wealth_grid = np.asarray(solution.b_grid, dtype=float)
    probe = node_entry_probabilities(
        solution.V,
        wealth_grid,
        parameters,
        outside_value=0.0,
        taste_scale=taste_scale,
    )
    feasible_values = [float(row["value"]) for row in probe.rows if bool(row["feasible"])]
    if not feasible_values or probe.feasible_weight + 1e-12 < target_qbar:
        raise RuntimeError(
            f"Entry target {target_qbar:.3f} exceeds feasible node mass {probe.feasible_weight:.3f}."
        )
    lower = min(feasible_values) - 20.0
    upper = max(feasible_values) + 20.0
    midpoint = 0.5 * (lower + upper)
    nodes = probe
    for _ in range(240):
        midpoint = 0.5 * (lower + upper)
        nodes = node_entry_probabilities(
            solution.V,
            wealth_grid,
            parameters,
            outside_value=midpoint,
            taste_scale=taste_scale,
        )
        gap = float(nodes.qbar) - float(target_qbar)
        if abs(gap) <= 1e-12 or abs(upper - lower) <= 1e-12:
            break
        if gap > 0.0:
            lower = midpoint
        else:
            upper = midpoint
    entry_flow = max(float(solution.entry_rate), 1e-14)
    mature_flow = float(solution.entrants_mature_total)
    local_weight = float(getattr(parameters, "local_birth_entry_weight", 1.0))
    outside_flow = entry_flow / float(nodes.qbar) - local_weight * mature_flow
    if outside_flow < 0.0:
        raise RuntimeError(f"Recovered outside entrant flow is negative: {outside_flow:.6g}")
    recovered = {
        "outside_value": float(midpoint),
        "outside_entry_flow": float(outside_flow),
        "target_qbar": float(target_qbar),
        "recovered_qbar": float(nodes.qbar),
        "entry_taste_scale": float(taste_scale),
        "local_birth_entry_weight": local_weight,
        "baseline_entry_flow": entry_flow,
        "baseline_mature_cityborn_flow": mature_flow,
        "feasible_entry_weight": float(nodes.feasible_weight),
    }
    return recovered, nodes.rows


def population_scale(
    solution: Any,
    parameters: Any,
    outside: dict[str, float],
    *,
    value_function: np.ndarray | None = None,
) -> tuple[SimpleNamespace, list[dict[str, Any]]]:
    values = (
        np.asarray(value_function, dtype=float)
        if value_function is not None
        else np.asarray(solution.V, dtype=float)
    )
    nodes = node_entry_probabilities(
        values,
        np.asarray(solution.b_grid, dtype=float),
        parameters,
        outside_value=float(outside["outside_value"]),
        taste_scale=float(outside["entry_taste_scale"]),
    )
    qbar = float(nodes.qbar)
    entry_flow = max(float(solution.entry_rate), 1e-14)
    mature_flow = float(solution.entrants_mature_total)
    local_weight = float(outside["local_birth_entry_weight"])
    denominator = 1.0 - qbar * local_weight * mature_flow / entry_flow
    finite = bool(qbar > 0.0 and denominator > 1e-12)
    if finite:
        implied_entry = qbar * float(outside["outside_entry_flow"]) / denominator
        scale = implied_entry / entry_flow
        potential_mass = (
            float(outside["outside_entry_flow"]) + local_weight * scale * mature_flow
        )
        entry_residual = implied_entry - qbar * potential_mass
    else:
        implied_entry = math.inf
        scale = math.inf
        potential_mass = math.inf
        entry_residual = math.nan
    return (
        SimpleNamespace(
            finite=finite,
            qbar=qbar,
            denominator=float(denominator),
            scale_factor=float(scale),
            implied_entry_total=float(implied_entry),
            implied_potential_entrant_mass=float(potential_mass),
            entry_residual=float(entry_residual),
            implied_total_population=float(solution.total_mass) * float(scale),
        ),
        nodes.rows,
    )


def load_fixed_results(
    fixed_dir: Path,
    *,
    theta: dict[str, float],
    smoke: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    metadata = json.loads((fixed_dir / "metadata.json").read_text(encoding="utf-8"))
    if bool(metadata.get("smoke")) != bool(smoke):
        raise ValueError("Fixed-result grid does not match --smoke.")
    fixed_theta = {str(key): float(value) for key, value in metadata["theta"].items()}
    if fixed_theta != theta:
        raise ValueError("Fixed-result theta does not match the selected current-M theta.")
    rows: list[dict[str, Any]] = []
    with (fixed_dir / "summary.csv").open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row: dict[str, Any] = {"case": raw["case"]}
            for key, value in raw.items():
                if key == "case":
                    continue
                if key in {"purchase_grant", "strict_converged"}:
                    row[key] = value == "True"
                else:
                    row[key] = float(value)
            rows.append(row)
    expected = [
        "rebated_tax1_baseline",
        "rebated_tax2",
        "rebated_tax2_grant0p4_Hge6",
    ]
    if [row["case"] for row in rows] != expected:
        raise ValueError(f"Unexpected fixed-result cases: {[row['case'] for row in rows]}")
    return rows, metadata


def solve_partial_case(
    base_overrides: dict[str, Any],
    *,
    annual_tax: float,
    grant: bool,
    transfer: float,
    price: float,
    smoke: bool,
) -> tuple[Any, Any]:
    overrides = case_overrides(
        base_overrides,
        annual_tax,
        grant,
        transfer,
        smoke,
    )
    overrides.update(
        solve_mode="pe",
        p_fixed=np.array([float(price)], dtype=float),
    )
    solution, parameters, _ = run_model_cp_dt(overrides, verbose=False)
    return solution, parameters


def solve_joint_case(
    base_overrides: dict[str, Any],
    outside: dict[str, float],
    *,
    label: str,
    annual_tax: float,
    grant: bool,
    initial_price: float,
    initial_transfer: float,
    smoke: bool,
    population_scale_function: Any | None = None,
) -> SimpleNamespace:
    cache: dict[tuple[float, float], SimpleNamespace] = {}
    scale_function = population_scale if population_scale_function is None else population_scale_function

    def evaluate(point: np.ndarray) -> SimpleNamespace:
        log_price = float(point[0])
        transfer = float(point[1])
        price = float(math.exp(log_price))
        key = (float(round(log_price, 10)), float(round(transfer, 10)))
        if key not in cache:
            solution, parameters = solve_partial_case(
                base_overrides,
                annual_tax=annual_tax,
                grant=grant,
                transfer=transfer,
                price=price,
                smoke=smoke,
            )
            scale, _ = scale_function(solution, parameters, outside)
            if not scale.finite:
                raise RuntimeError(f"{label}: entry/population closure is not finite.")
            demand = (
                float(scale.scale_factor)
                * float(np.asarray(solution.housing_demand, dtype=float).reshape(-1)[0])
            )
            supply = float(parameters.H0[0]) * (
                float(parameters.user_cost_rate) * price / float(parameters.r_bar[0])
            ) ** float(parameters.xi_supply[0])
            market_residual = (demand - supply) / max(abs(supply), 1e-14)
            raw_fiscal_residual = float(solution.property_tax_budget_residual)
            fiscal_scale = max(
                abs(float(solution.property_tax_revenue)),
                abs(float(solution.property_tax_transfer_outlays))
                + abs(float(solution.birth_entry_grant_outlays)),
                0.1,
            )
            fiscal_residual = raw_fiscal_residual / fiscal_scale
            cache[key] = SimpleNamespace(
                point=np.array([log_price, transfer], dtype=float),
                price=price,
                transfer=transfer,
                solution=solution,
                parameters=parameters,
                scale=scale,
                relative_excess=float(market_residual),
                scaled_fiscal_residual=(
                    float(scale.scale_factor) * raw_fiscal_residual
                ),
                residual=np.array(
                    [float(market_residual), float(fiscal_residual)],
                    dtype=float,
                ),
            )
            print(
                f"  {label}: T={transfer:.6f}, p={price:.6f}, "
                f"S={scale.scale_factor:.5f}, "
                f"market={market_residual:.2e}, fiscal={raw_fiscal_residual:.2e}",
                flush=True,
            )
        return cache[key]

    tolerance = 2.0e-4 if smoke else 2.5e-5
    point = np.array([math.log(float(initial_price)), float(initial_transfer)], dtype=float)
    result = evaluate(point)
    result.joint_iterations = 0
    if float(np.max(np.abs(result.residual))) > tolerance:
        steps = np.array([0.01, 0.01], dtype=float)
        jacobian = np.column_stack(
            [
                (evaluate(point + np.array([steps[0], 0.0])).residual - result.residual)
                / steps[0],
                (evaluate(point + np.array([0.0, steps[1]])).residual - result.residual)
                / steps[1],
            ]
        )
        for iteration in range(14):
            if float(np.max(np.abs(result.residual))) <= tolerance:
                break
            try:
                direction = np.linalg.solve(jacobian, -result.residual)
            except np.linalg.LinAlgError:
                direction = np.linalg.lstsq(
                    jacobian,
                    -result.residual,
                    rcond=None,
                )[0]
            direction = np.clip(direction, [-0.25, -0.15], [0.25, 0.15])
            current_norm = float(np.linalg.norm(result.residual))
            accepted: SimpleNamespace | None = None
            accepted_point: np.ndarray | None = None
            for damping in (1.0, 0.5, 0.25, 0.125):
                trial_point = point + damping * direction
                trial_point[0] = float(
                    np.clip(
                        trial_point[0],
                        math.log(float(base_overrides.get("p_min", 0.01))),
                        math.log(float(base_overrides.get("p_max", 30.0))),
                    )
                )
                trial_point[1] = float(np.clip(trial_point[1], 0.0, 8.0))
                trial = evaluate(trial_point)
                if float(np.linalg.norm(trial.residual)) < current_norm:
                    accepted = trial
                    accepted_point = trial_point
                    break
            if accepted is None or accepted_point is None:
                steps *= 0.5
                jacobian = np.column_stack(
                    [
                        (
                            evaluate(point + np.array([steps[0], 0.0])).residual
                            - result.residual
                        )
                        / steps[0],
                        (
                            evaluate(point + np.array([0.0, steps[1]])).residual
                            - result.residual
                        )
                        / steps[1],
                    ]
                )
                continue
            displacement = accepted_point - point
            residual_change = accepted.residual - result.residual
            denominator = float(displacement @ displacement)
            if denominator > 1e-14:
                jacobian += np.outer(
                    residual_change - jacobian @ displacement,
                    displacement,
                ) / denominator
            point = accepted_point
            result = accepted
            result.joint_iterations = iteration + 1
    if float(np.max(np.abs(result.residual))) > tolerance:
        raise RuntimeError(
            f"{label}: joint root did not converge after {len(cache)} model evaluations; "
            f"residuals={result.residual.tolist()}."
        )
    result.joint_model_evaluations = len(cache)
    result.price_bracket_expansions = 0
    return result


def make_plot(outdir: Path, rows: list[dict[str, Any]]) -> None:
    policies = rows[1:]
    labels = ["Rebated 2% tax", "Rebated 2% + grant"]
    y = np.arange(len(policies))[::-1]
    fixed_tfr = np.asarray([row["fixed_tfr_change_percent"] for row in policies])
    entry_tfr = np.asarray([row["entry_tfr_change_percent"] for row in policies])
    fixed_price = np.asarray([row["fixed_price_change_percent"] for row in policies])
    entry_price = np.asarray([row["entry_price_change_percent"] for row in policies])
    total_births = np.asarray([row["total_births_change_percent"] for row in policies])

    plt.rcdefaults()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
        }
    )
    navy = "#17365D"
    orange = "#C45A3B"
    green = "#2B6F68"
    gray = "#A8ADB3"
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(13.0, 5.2),
        gridspec_kw={"width_ratios": [1.22, 1.12, 0.90]},
    )

    for axis, fixed, entry, title, xlabel in (
        (
            axes[0],
            fixed_tfr,
            entry_tfr,
            "A. Completed fertility",
            "Percent change in TFR",
        ),
        (
            axes[1],
            fixed_price,
            entry_price,
            "B. Equilibrium house price",
            "Percent change",
        ),
    ):
        axis.axvline(0.0, color="#4B4B4B", linewidth=0.7)
        for yy, fixed_value, entry_value in zip(y, fixed, entry):
            axis.hlines(yy, fixed_value, entry_value, color=gray, linewidth=1.3)
        axis.scatter(fixed, y, color=navy, s=45, label="Population fixed", zorder=3)
        axis.scatter(
            entry,
            y,
            color=orange,
            marker="s",
            s=45,
            label="Entry adjusts",
            zorder=3,
        )
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel(xlabel)
    axes[0].set_yticks(y, labels)
    axes[1].set_yticks(y, ["", ""])

    axes[2].axvline(0.0, color="#4B4B4B", linewidth=0.7)
    axes[2].scatter(total_births, y, color=green, s=48, zorder=3)
    for yy, value in zip(y, total_births):
        axes[2].annotate(
            f"{value:+.2f}%",
            (value, yy),
            xytext=(6 if value >= 0.0 else -6, 0),
            textcoords="offset points",
            ha="left" if value >= 0.0 else "right",
            va="center",
        )
    axes[2].set_yticks(y, ["", ""])
    axes[2].set_title("C. Total births", loc="left", fontweight="bold")
    axes[2].set_xlabel("Percent change")

    handles, legend_labels = axes[0].get_legend_handles_labels()
    for axis in axes:
        axis.grid(axis="x", color="#D9D9D9", linewidth=0.5)
        axis.set_axisbelow(True)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
    note = (
        "Notes: Fixed current-M parameters. All policies are fiscally closed relative to "
        "the rebated 1% baseline. The entry-adjusting equilibrium holds the calibrated "
        "outside option fixed and lets population scale and the house price move jointly."
    )
    fig.subplots_adjust(left=0.14, right=0.985, top=0.82, bottom=0.23, wspace=0.25)
    fig.legend(
        handles,
        legend_labels,
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.98),
        ncol=2,
    )
    fig.text(0.14, 0.07, note, ha="left", va="bottom", fontsize=8.3, color="#333333")
    fig.savefig(
        outdir / "policy_general_equilibrium_feedback_nb120.png",
        dpi=240,
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    source = args.source.resolve()
    fixed_dir = args.fixed_dir.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    profile_label, theta, base_overrides, target_system, _ = load_profile(
        "new-moment",
        source,
        args.smoke,
    )
    fixed_rows, fixed_metadata = load_fixed_results(
        fixed_dir,
        theta=theta,
        smoke=args.smoke,
    )
    cases = (
        ("rebated_tax1_baseline", 0.01, False),
        ("rebated_tax2", 0.02, False),
        ("rebated_tax2_grant0p4_Hge6", 0.02, True),
    )
    metadata = {
        "status": "running",
        "profile": profile_label,
        "source": str(source),
        "fixed_results": str(fixed_dir),
        "theta": theta,
        "smoke": bool(args.smoke),
        "entry_protocol": {
            "source": "funded extension of phase9b_entry_margin_scale_loop.py",
            "aggregation": "node-level logit over entrant income and wealth nodes",
            "target_baseline_city_entry_probability": TARGET_QBAR,
            "entry_taste_scale": ENTRY_TASTE_SCALE,
            "outside_objects": "recovered at funded baseline and fixed across policies",
        },
        "fiscal_rule": (
            "scaled property-tax revenue = scaled equal lump-sum transfers "
            "+ scaled targeted purchase-grant outlays"
        ),
    }
    (outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    fixed_baseline = fixed_rows[0]
    baseline_solution, baseline_parameters = solve_partial_case(
        base_overrides,
        annual_tax=0.01,
        grant=False,
        transfer=float(fixed_baseline["lump_sum_transfer_period_units"]),
        price=float(fixed_baseline["price"]),
        smoke=args.smoke,
    )
    outside, entry_rows = calibrate_outside_objects(
        baseline_solution,
        baseline_parameters,
        target_qbar=TARGET_QBAR,
        taste_scale=ENTRY_TASTE_SCALE,
    )
    identity_scale, _ = population_scale(
        baseline_solution,
        baseline_parameters,
        outside,
    )
    outside["identity_scale_factor"] = float(identity_scale.scale_factor)
    outside["identity_entry_residual"] = float(identity_scale.entry_residual)
    (outdir / "benchmark_outside_objects.json").write_text(
        json.dumps(jsonable(outside), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv(outdir / "benchmark_entry_nodes.csv", entry_rows)

    results: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    for index, ((label, annual_tax, grant), fixed) in enumerate(
        zip(cases, fixed_rows),
        start=1,
    ):
        case_started = time.perf_counter()
        print(f"case {index}/{len(cases)}: {label}", flush=True)
        joint = solve_joint_case(
            base_overrides,
            outside,
            label=label,
            annual_tax=annual_tax,
            grant=grant,
            initial_price=float(fixed["price"]),
            initial_transfer=float(fixed["lump_sum_transfer_period_units"]),
            smoke=args.smoke,
        )
        solution = joint.solution
        parameters = joint.parameters
        scale = joint.scale
        moments = extract_moments(solution, parameters)
        result = {
            "case": label,
            "annual_property_tax_rate": annual_tax,
            "purchase_grant": grant,
            "fixed_tfr": float(fixed["tfr"]),
            "entry_tfr": float(moments["tfr"]),
            "fixed_price": float(fixed["price"]),
            "entry_price": float(joint.price),
            "population_scale": float(scale.scale_factor),
            "entry_probability": float(scale.qbar),
            "normalized_total_births": float(solution.total_births_kfe),
            "scaled_total_births": (
                float(scale.scale_factor) * float(solution.total_births_kfe)
            ),
            "lump_sum_transfer_period_units": float(joint.transfer),
            "lump_sum_transfer_annual_income_units": (
                float(joint.transfer) / float(parameters.period_years)
            ),
            "scaled_property_tax_revenue": (
                float(scale.scale_factor) * float(solution.property_tax_revenue)
            ),
            "scaled_grant_outlays": (
                float(scale.scale_factor) * float(solution.birth_entry_grant_outlays)
            ),
            "scaled_transfer_outlays": (
                float(scale.scale_factor) * float(solution.property_tax_transfer_outlays)
            ),
            "scaled_fiscal_residual": float(joint.scaled_fiscal_residual),
            "scaled_market_residual": float(joint.relative_excess),
            "joint_model_evaluations": int(joint.joint_model_evaluations),
            "joint_iterations": int(joint.joint_iterations),
            "elapsed_seconds": time.perf_counter() - case_started,
        }
        results.append(result)
        for name, target, weight in zip(
            target_system.moment_names,
            target_system.target_values,
            target_system.weights,
        ):
            model = float(moments[name])
            gap = model - float(target)
            target_rows.append(
                {
                    "case": label,
                    "moment": name,
                    "target": float(target),
                    "model": model,
                    "gap": gap,
                    "weight": float(weight),
                    "loss_contribution": float(weight) * gap * gap,
                }
            )
        write_csv(outdir / "summary.csv", results)
        write_csv(outdir / "target_fit_long.csv", target_rows)
        (outdir / "latest_completed_case.json").write_text(
            json.dumps(jsonable(result), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(
            f"completed {label}: S={result['population_scale']:.5f}, "
            f"p={result['entry_price']:.6f}, TFR={result['entry_tfr']:.5f}, "
            f"fiscal={result['scaled_fiscal_residual']:.2e}, "
            f"market={result['scaled_market_residual']:.2e}, "
            f"elapsed={result['elapsed_seconds']:.1f}s",
            flush=True,
        )

    baseline = results[0]
    for row in results:
        row["fixed_tfr_change_level"] = row["fixed_tfr"] - baseline["fixed_tfr"]
        row["fixed_tfr_change_percent"] = 100.0 * (
            row["fixed_tfr"] / baseline["fixed_tfr"] - 1.0
        )
        row["entry_tfr_change_level"] = row["entry_tfr"] - baseline["entry_tfr"]
        row["entry_tfr_change_percent"] = 100.0 * (
            row["entry_tfr"] / baseline["entry_tfr"] - 1.0
        )
        row["fixed_price_change_percent"] = 100.0 * (
            row["fixed_price"] / baseline["fixed_price"] - 1.0
        )
        row["entry_price_change_percent"] = 100.0 * (
            row["entry_price"] / baseline["entry_price"] - 1.0
        )
        row["total_births_change_percent"] = 100.0 * (
            row["scaled_total_births"] / baseline["scaled_total_births"] - 1.0
        )
    write_csv(outdir / "summary.csv", results)
    make_plot(outdir, results)

    lines = [
        "# Funded housing policy with entry adjustment",
        "",
        "Fixed current-M parameters. The government budget, house price, entry probability, "
        "and stationary population scale are solved consistently in every case.",
        "",
        "| Case | TFR | TFR change | Population scale | Total births | Price change |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in results:
        lines.append(
            f"| {row['case']} | {row['entry_tfr']:.4f} | "
            f"{row['entry_tfr_change_level']:+.4f} "
            f"({row['entry_tfr_change_percent']:+.2f}%) | "
            f"{row['population_scale']:.4f} | "
            f"{row['total_births_change_percent']:+.2f}% | "
            f"{row['entry_price_change_percent']:+.2f}% |"
        )
    lines.extend(
        [
            "",
            "The fixed-population columns in `summary.csv` are decomposition results only. "
            "The entry-adjusting columns are the headline general-equilibrium results.",
        ]
    )
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    metadata.update(
        status="complete",
        elapsed_seconds=time.perf_counter() - started,
        benchmark_outside_objects=outside,
        fixed_metadata_status=fixed_metadata.get("status"),
    )
    (outdir / "metadata.json").write_text(
        json.dumps(jsonable(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
