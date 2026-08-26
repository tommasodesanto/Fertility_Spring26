#!/usr/bin/env python3
"""Evaluate one fixed policy point under the coherent terminal person law.

The supplied price, transfer, tax rate, and fertility-preference intercept are
not solved by this driver. The driver solves only the inner stationary mapping:
household policy implies births per head; the annual person law pins the finite
population and head-age stocks; the household distribution is then iterated to
the economic-state composition consistent with those stocks.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
MODEL_ROOT = ROOT / "code/model"
TOOLS_ROOT = MODEL_ROOT / "tools"
for search_path in (MODEL_ROOT, TOOLS_ROOT):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

import run_dynamic_population_transition as calendar  # noqa: E402
import run_e5f_open_population_transition as transition  # noqa: E402
import run_e5f_perfect_foresight_person_demography as person_pf  # noqa: E402
import run_e5f_perfect_foresight_rebated_property_tax as funded  # noqa: E402
import run_e5f_perfect_foresight_transition as pf  # noqa: E402
import run_e5f_post2023_coven_property_tax_smoke as impact  # noqa: E402
import run_e5f_post2023_no_policy_continuations as continuation  # noqa: E402


DEFAULT_OUTPUT = ROOT / "output/model/e5f_person_demography_terminal_point_20260826a"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=sorted(impact.CASES), required=True)
    parser.add_argument("--asset-price", type=float, required=True)
    parser.add_argument("--equal-transfer", type=float, required=True)
    parser.add_argument("--psi-child", type=float, required=True)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument(
        "--source-dir", type=Path, default=person_pf.DEFAULT_SOURCE_DIR
    )
    parser.add_argument(
        "--headship-dir", type=Path, default=person_pf.DEFAULT_HEADSHIP_DIR
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--maximum-iterations", type=int, default=250)
    parser.add_argument("--damping", type=float, default=0.50)
    return parser.parse_args()


def write_terminal_figure(
    result: person_pf.TerminalHouseholdPersonFixedPoint,
    parameters: Any,
    output: Path,
) -> None:
    iterations = np.asarray([row["iteration"] for row in result.history])
    mapping = np.asarray(
        [row["distribution_mapping_relative_l1"] for row in result.history]
    )
    birth_gap = np.asarray(
        [row["annual_birth_rate_relative_gap"] for row in result.history]
    )
    ages = np.arange(result.persons.persons.shape[1])
    model_ages = float(parameters.age_start) + float(parameters.da) * np.arange(
        int(parameters.J)
    )
    model_heads = person_pf._age_masses(result.g_pre)

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), constrained_layout=True)
    axes[0, 0].semilogy(iterations, np.maximum(mapping, 1e-16), label="Distribution")
    axes[0, 0].semilogy(
        iterations,
        np.maximum(np.where(np.isfinite(birth_gap), birth_gap, np.nan), 1e-16),
        label="Birth rate",
    )
    axes[0, 0].set(title="Inner fixed-point residuals", xlabel="Iteration")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].plot(
        iterations,
        [row["annual_births_per_head"] for row in result.history],
        color="#B4433C",
    )
    axes[0, 1].set(title="Annual births per household head", xlabel="Iteration")
    axes[1, 0].plot(
        ages,
        np.sum(result.persons.persons, axis=0),
        label="Resident persons",
    )
    axes[1, 0].plot(
        ages,
        np.sum(result.persons.heads, axis=0),
        label="Household heads",
    )
    axes[1, 0].set(title="Terminal age profiles", xlabel="Single year of age")
    axes[1, 0].legend(frameon=False)
    axes[1, 1].plot(model_ages, model_heads, marker="o", label="Household model")
    person_heads = person_pf.aggregate_heads_to_model_age_cells(
        result.persons,
        age_start=int(round(float(parameters.age_start))),
        cell_width=int(round(float(parameters.da))),
        number_of_cells=int(parameters.J),
    )
    axes[1, 1].plot(
        model_ages,
        person_heads,
        marker="x",
        linestyle="none",
        label="Person law",
    )
    axes[1, 1].set(title="Four-year head-stock identity", xlabel="Age-cell lower bound")
    axes[1, 1].legend(frameon=False)
    for axis in axes.flat:
        axis.grid(alpha=0.18)
    figure.savefig(output, dpi=190)
    figure.savefig(output.with_suffix(".pdf"))
    plt.close(figure)


def main() -> None:
    args = parse_args()
    for label, value in (
        ("asset price", args.asset_price),
        ("equal transfer", args.equal_transfer),
        ("psi_child", args.psi_child),
    ):
        if not math.isfinite(float(value)):
            raise ValueError(f"{label} must be finite")
    if float(args.asset_price) <= 0.0 or float(args.equal_transfer) < 0.0:
        raise ValueError("Price must be positive and transfer nonnegative")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    chain, model = transition.configure_sequential_model()
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    contracts, _ = pf.load_diagnostic_contracts(
        args.selected_report.resolve(),
        args.selected_transition.resolve(),
        args.source.resolve(),
    )
    prepared = continuation.prepare_model(
        contracts, args.source.resolve(), chain, model
    )
    initial_state, _, history_gate, impact_price = pf.reconstruct_2023_state(
        prepared,
        contracts,
        args.selected_transition.resolve(),
        market_tolerance=2e-4,
    )
    primitives = person_pf.build_annual_demographic_primitives(
        initial_state.g_pre,
        prepared.parameters,
        source_dir=args.source_dir,
        headship_dir=args.headship_dir,
    )

    anchor_state = continuation.DynamicState(
        g_pre=initial_state.g_pre.copy(),
        scheduled_entries=list(initial_state.scheduled_entries),
        scheduled_raw_entries=list(initial_state.scheduled_raw_entries),
        price_guess=float(impact_price),
        initial_policy=None,
    )
    anchor_parameters = copy.deepcopy(prepared.parameters)
    anchor_parameters.psi_child = float(contracts["report_best"]["new_psi_child"])
    anchor = impact.solve_fixed_price_equal_rebate(
        state=anchor_state,
        P=anchor_parameters,
        b_grid=prepared.b_grid,
        price=float(impact_price),
        case=impact.CASES["rebated-tax1-baseline"],
        counter=calendar.SolveCounter(),
    )
    supply_rule, supply_contract = impact.reanchor_supply_rule(
        inherited_price=float(impact_price),
        baseline_demand=float(anchor.evaluation.demand_by_loc[0]),
        elasticity=1.75,
    )

    case = impact.CASES[args.case]
    stationary = funded._stationary_evaluation(
        prepared=prepared,
        case=case,
        asset_price=float(args.asset_price),
        transfer=float(args.equal_transfer),
        psi_child=float(args.psi_child),
    )
    result = person_pf.solve_terminal_household_person_fixed_point(
        policy=stationary.policy,
        parameters=stationary.parameters,
        b_grid=stationary.b_grid,
        initial_g_pre=stationary.unit_g_pre,
        demographic_primitives=primitives,
        supply_rule=supply_rule,
        maximum_iterations=int(args.maximum_iterations),
        damping=float(args.damping),
    )
    pf.write_csv(output_dir / "fixed_point_history.csv", result.history)
    write_terminal_figure(result, stationary.parameters, output_dir / "diagnostics.png")

    summary = {
        "status": (
            "complete_fixed_policy_terminal_mapping_not_market_fiscal_root"
            if result.converged
            else "failed_fixed_policy_terminal_mapping"
        ),
        "promotion_status": "not_promoted",
        "fixed_policy_inputs": {
            "case": case.name,
            "annual_property_tax_rate": float(case.annual_tax_rate),
            "asset_price": float(args.asset_price),
            "equal_transfer_period_units": float(args.equal_transfer),
            "psi_child": float(args.psi_child),
        },
        "interpretation": (
            "Only the inner demographic-household stationary mapping is solved. "
            "The reported housing and fiscal residuals diagnose the supplied point; "
            "price, transfer, and psi_child are not roots of the corrected model."
        ),
        "initial_history_reproduction_status": history_gate.get("status"),
        "supply_contract": supply_contract,
        "converged": result.converged,
        "iterations": result.iterations,
        "resident_persons_model_units": float(np.sum(result.persons.persons)),
        "resident_persons_actual_scale": float(np.sum(result.persons.persons))
        / primitives.scale_model_units_per_person,
        "household_heads_model_units": float(np.sum(result.g_pre)),
        "annual_births_per_head": result.annual_births_per_head,
        "model_period_births": result.model_period_births,
        "renewal_ratio": result.renewal_ratio,
        "housing_demand": result.housing_demand,
        "housing_supply": result.housing_supply,
        "housing_excess_demand_relative": result.housing_excess_demand_relative,
        "housing_market_relative_residual": result.housing_market_relative_residual,
        "government_budget_residual": result.government_budget_residual,
        "scaled_government_budget_residual": (
            result.scaled_government_budget_residual
        ),
        "implied_equal_transfer": result.implied_equal_transfer,
        "equal_transfer_gap": result.equal_transfer_gap,
        "distribution_mapping_relative_l1": result.distribution_mapping_relative_l1,
        "annual_birth_rate_relative_gap": result.annual_birth_rate_relative_gap,
        "person_one_step_relative_l1": result.person_one_step_relative_l1,
        "head_one_step_relative_l1": result.head_one_step_relative_l1,
        "age_head_one_step_max_abs": result.age_head_one_step_max_abs,
        "household_person_head_gap": result.household_person_head_gap,
        "headship_alignment_factors": (
            primitives.model_age_headship_alignment_factors.tolist()
        ),
        "source_hashes": {
            "driver": pf.file_sha256(Path(__file__).resolve()),
            "person_demography": pf.file_sha256(
                TOOLS_ROOT / "run_e5f_perfect_foresight_person_demography.py"
            ),
            "perfect_foresight": pf.file_sha256(
                TOOLS_ROOT / "run_e5f_perfect_foresight_transition.py"
            ),
            "funded_policy": pf.file_sha256(
                TOOLS_ROOT / "run_e5f_perfect_foresight_rebated_property_tax.py"
            ),
            "solver": pf.file_sha256(
                MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
            ),
            **{
                f"input_{name}": pf.file_sha256(path)
                for name, path in primitives.source_paths.items()
            },
        },
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not result.converged:
        raise RuntimeError("Corrected terminal household/person mapping did not converge")


if __name__ == "__main__":
    main()
