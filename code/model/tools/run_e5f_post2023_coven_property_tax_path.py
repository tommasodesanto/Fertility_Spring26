#!/usr/bin/env python3
"""Short-run rebated property-tax paths from the fitted 2023 state.

The baseline permanently rebates a 1 percent annual property tax; the reform
permanently rebates a 2 percent annual property tax. Both paths begin from the
same fitted 2023 pre-fertility distribution and birth-vintage queue. At every
four-year date the housing price and equal per-household transfer solve the
housing market and government budget jointly. Population is nationally closed
after 2023 (M=0, rho=1), and no purchase grant is available.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_coven_property_tax_smoke as impact
import run_e5f_post2023_no_policy_continuations as baseline
import run_e5f_post2023_policy_mechanisms as mechanism
import run_e5f_post2023_rebated_property_tax_smoke as tax


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = ROOT / "output/model/e5f_post2023_coven_property_tax_path"
DEFAULT_POST_2023_PERIODS = 3
DEFAULT_HOUSING_SUPPLY_ELASTICITY = 1.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, required=True)
    parser.add_argument("--selected-case-transition", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--expected-case-transition-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-scientific-contract-sha256", required=True)
    parser.add_argument("--expected-selection-sha256", required=True)
    parser.add_argument("--market-tol", type=float, default=tax.MARKET_TOLERANCE)
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument(
        "--post-2023-periods", type=int, default=DEFAULT_POST_2023_PERIODS
    )
    parser.add_argument(
        "--housing-supply-elasticity",
        type=float,
        default=DEFAULT_HOUSING_SUPPLY_ELASTICITY,
    )
    return parser.parse_args()


def copy_state(state: baseline.DynamicState) -> baseline.DynamicState:
    return baseline.DynamicState(
        g_pre=state.g_pre.copy(),
        scheduled_entries=list(state.scheduled_entries),
        scheduled_raw_entries=list(state.scheduled_raw_entries),
        price_guess=float(state.price_guess),
        initial_policy=None,
    )


def gate_period_row(
    row: dict[str, Any], next_state: baseline.DynamicState
) -> dict[str, Any]:
    gates = {
        "calendar_year": int(row["calendar_year"]),
        "market_residual": abs(float(row["relative_market_residual"])),
        "fiscal_residual": abs(float(row["government_budget_residual"])),
        "mass_residual": abs(float(row["mass_accounting_residual"])),
        "nonfinite_distribution_count": int(row["nonfinite_distribution_count"]),
        "minimum_distribution_mass": float(row["min_distribution_mass"]),
        "feasibility_projection_mass": float(
            row["feasibility_frontier_projection_mass"]
        ),
        "next_distribution_sha256": baseline.array_sha256(next_state.g_pre),
        "next_queue_sha256": baseline.canonical_json_sha256(
            next_state.scheduled_entries
        ),
    }
    if gates["market_residual"] > tax.MARKET_TOLERANCE:
        raise RuntimeError(f"Market gate failed: {gates}")
    if gates["fiscal_residual"] > tax.FISCAL_ABSOLUTE_TOLERANCE:
        raise RuntimeError(f"Fiscal gate failed: {gates}")
    if gates["mass_residual"] > tax.MASS_TOLERANCE:
        raise RuntimeError(f"Mass gate failed: {gates}")
    if gates["nonfinite_distribution_count"] != 0:
        raise RuntimeError(f"Nonfinite distribution: {gates}")
    if gates["minimum_distribution_mass"] < -1.0e-14:
        raise RuntimeError(f"Negative distribution mass: {gates}")
    if gates["feasibility_projection_mass"] > 1.0e-10:
        raise RuntimeError(f"Feasibility projection gate failed: {gates}")
    return gates


def run_case_path(
    *,
    case: tax.TaxCase,
    prepared: baseline.PreparedModel,
    state_2023: baseline.DynamicState,
    parameters_2023: SimpleNamespace,
    mass_2023: float,
    supply_rule: calendar.HousingSupplyRule,
    housing_supply_elasticity: float,
    initial_price: float,
    initial_transfer: float,
    post_2023_periods: int,
    progress_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    state = copy_state(state_2023)
    state.price_guess = float(initial_price)
    parameters = copy.deepcopy(parameters_2023)
    transfer_guess = float(initial_transfer)
    counter = calendar.SolveCounter()
    rows: list[dict[str, Any]] = []
    period_gates: list[dict[str, Any]] = []

    for post_index in range(post_2023_periods + 1):
        period = baseline.TRANSITION_PERIODS + post_index
        solved = tax.solve_joint_rebated_period(
            state=state,
            P=parameters,
            b_grid=prepared.b_grid,
            supply_rule=supply_rule,
            counter=counter,
            case=case,
            initial_price=float(state.price_guess),
            initial_transfer=transfer_guess,
        )
        row, next_state = baseline.advance_from_evaluation(
            label=case.name,
            period_from_2007=period,
            evaluation=solved.evaluation,
            state=state,
            P=parameters,
            b_grid=prepared.b_grid,
            shared=solved.shared,
            supply_rule=supply_rule,
            outside_flow=0.0,
            retention=1.0,
            initial_mass_2007=prepared.initial_mass_2007,
            mass_2023=mass_2023,
            next_bridge_year=None,
            grid_fallback=False,
        )
        row.update(solved.ledger)
        row.update(
            policy_case=case.name,
            policy_label=case.label,
            annual_property_tax_rate=case.annual_tax_rate,
            property_tax_rate_period_units=(
                case.annual_tax_rate * float(parameters.period_years)
            ),
            equal_rebate=True,
            housing_supply_elasticity=float(housing_supply_elasticity),
            housing_user_cost=float(parameters.user_cost_rate) * float(solved.price),
            user_cost_rate=float(parameters.user_cost_rate),
            joint_iterations=int(solved.joint_iterations),
            joint_model_evaluations=int(solved.joint_model_evaluations),
            cumulative_model_solve_count=int(counter.total),
        )
        gates = gate_period_row(row, next_state)
        rows.append(row)
        period_gates.append(gates)
        baseline.write_csv(progress_dir / "path_progress.csv", rows)
        baseline.write_json(
            progress_dir / "latest_completed_period.json",
            {
                "status": "running",
                "policy_case": case.name,
                "calendar_year": int(row["calendar_year"]),
                "completed_dates": len(rows),
                "requested_dates": post_2023_periods + 1,
                "latest": row,
                "latest_gates": gates,
            },
        )
        print(
            f"REBATED_TAX_PATH case={case.name} year={row['calendar_year']} "
            f"price={row['asset_price']:.8f} transfer={row['equal_transfer_period_units']:.8f} "
            f"births={row['topcode_adjusted_births_per_adult']:.8f}",
            flush=True,
        )
        state = next_state
        transfer_guess = float(solved.transfer)
    return rows, period_gates, int(counter.total)


def build_effect_rows(
    baseline_rows: list[dict[str, Any]], reform_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    if len(baseline_rows) != len(reform_rows):
        raise ValueError("Baseline and reform paths have different lengths")
    effects: list[dict[str, Any]] = []
    for base, reform in zip(baseline_rows, reform_rows):
        if int(base["calendar_year"]) != int(reform["calendar_year"]):
            raise ValueError("Baseline and reform calendars differ")
        row = {
            "calendar_year": int(base["calendar_year"]),
            "years_from_2023": int(base["years_from_2023"]),
            "asset_price_percent_change": 100.0
            * (float(reform["asset_price"]) / float(base["asset_price"]) - 1.0),
            "housing_user_cost_percent_change": 100.0
            * (
                float(reform["housing_user_cost"])
                / float(base["housing_user_cost"])
                - 1.0
            ),
            "equal_transfer_percent_change": 100.0
            * (
                float(reform["equal_transfer_period_units"])
                / float(base["equal_transfer_period_units"])
                - 1.0
            ),
            "owner_rate_pp_change": 100.0
            * (float(reform["owner_rate"]) - float(base["owner_rate"])),
            "dependent_child_owner_rate_pp_change": 100.0
            * (
                float(reform["dependent_child_owner_rate"])
                - float(base["dependent_child_owner_rate"])
            ),
            "births_per_adult_percent_change": 100.0
            * (
                float(reform["topcode_adjusted_births_per_adult"])
                / float(base["topcode_adjusted_births_per_adult"])
                - 1.0
            ),
            "adult_population_percent_change": 100.0
            * (
                float(reform["adult_population"])
                / float(base["adult_population"])
                - 1.0
            ),
        }
        effects.append(row)
    return effects


def make_figure(paths: dict[str, list[dict[str, Any]]], outdir: Path) -> None:
    fields = (
        ("asset_price", "House price"),
        ("owner_rate", "Ownership rate"),
        ("dependent_child_owner_rate", "Dependent-child ownership"),
        ("topcode_adjusted_births_per_adult", "Adjusted births per adult"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.3), constrained_layout=True)
    for axis, (field, title) in zip(axes.flat, fields):
        for case_name, rows in paths.items():
            label = "Rebated 1%" if case_name == "rebated-tax1-baseline" else "Rebated 2%"
            axis.plot(
                [int(row["calendar_year"]) for row in rows],
                [float(row[field]) for row in rows],
                marker="o",
                label=label,
            )
        axis.set_title(title)
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    axes.flat[0].legend(frameon=False)
    fig.savefig(outdir / "rebated_property_tax_path.png", dpi=220)
    fig.savefig(outdir / "rebated_property_tax_path.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    post_periods = int(args.post_2023_periods)
    if not 1 <= post_periods <= 10:
        raise ValueError("--post-2023-periods must lie in [1,10]")
    elasticity = float(args.housing_supply_elasticity)
    if not math.isfinite(elasticity) or elasticity <= 0.0:
        raise ValueError("--housing-supply-elasticity must be finite and positive")
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    chain, model = transition.configure_sequential_model()
    contracts = baseline.validate_input_contracts(
        report_path=args.selected_report,
        task_summary_path=None,
        case_dir=None,
        case_transition_path=args.selected_case_transition,
        source_path=args.source,
        expected_report_sha256=args.expected_report_sha256,
        expected_task_summary_sha256=None,
        expected_case_transition_sha256=args.expected_case_transition_sha256,
        expected_source_sha256=args.expected_source_sha256,
        expected_target_fingerprint=args.expected_target_fingerprint,
        expected_code_bundle_sha256=args.expected_code_bundle_sha256,
        expected_renewal_contract_sha256=args.expected_renewal_contract_sha256,
        expected_scientific_contract_sha256=args.expected_scientific_contract_sha256,
        expected_selection_sha256=args.expected_selection_sha256,
        chain=chain,
        model=model,
    )
    prepared = baseline.prepare_model(contracts, args.source.resolve(), chain, model)
    history, state_2023, parameters_2023, mass_2023, reconstruction = (
        mechanism.reconstruct_matched_2023(
            prepared,
            contracts,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
            progress_dir=outdir,
        )
    )
    inherited_price = impact.inherited_2023_price(
        args.selected_case_transition.resolve()
    )
    anchor_counter = calendar.SolveCounter()
    anchor = impact.solve_fixed_price_equal_rebate(
        state=state_2023,
        P=copy.deepcopy(parameters_2023),
        b_grid=prepared.b_grid,
        price=inherited_price,
        case=impact.CASES["rebated-tax1-baseline"],
        counter=anchor_counter,
    )
    baseline_demand = float(anchor.evaluation.demand_by_loc[0])
    supply_rule, reanchor = impact.reanchor_supply_rule(
        inherited_price=inherited_price,
        baseline_demand=baseline_demand,
        elasticity=elasticity,
    )
    reanchor.update(
        fixed_price_equal_transfer=float(anchor.transfer),
        fixed_price_fiscal_residual=float(anchor.ledger["government_budget_residual"]),
        fixed_price_model_evaluations=int(anchor.model_evaluations),
    )
    baseline.write_json(outdir / "supply_reanchor_gate.json", reanchor)

    common_distribution_sha = baseline.array_sha256(state_2023.g_pre)
    common_queue_sha = baseline.canonical_json_sha256(state_2023.scheduled_entries)
    paths: dict[str, list[dict[str, Any]]] = {}
    gates: dict[str, list[dict[str, Any]]] = {}
    solve_counts: dict[str, int] = {}
    for case in impact.CASES.values():
        case_dir = outdir / case.name
        case_dir.mkdir(parents=True, exist_ok=False)
        rows, case_gates, solve_count = run_case_path(
            case=case,
            prepared=prepared,
            state_2023=state_2023,
            parameters_2023=parameters_2023,
            mass_2023=mass_2023,
            supply_rule=supply_rule,
            housing_supply_elasticity=elasticity,
            initial_price=inherited_price,
            initial_transfer=(
                float(anchor.transfer) * (float(case.annual_tax_rate) / 0.01)
            ),
            post_2023_periods=post_periods,
            progress_dir=case_dir,
        )
        for row in rows:
            row["common_initial_distribution_sha256"] = common_distribution_sha
            row["common_initial_queue_sha256"] = common_queue_sha
        baseline.write_csv(case_dir / "path.csv", rows)
        baseline.write_json(
            case_dir / "summary.json",
            {
                "status": "complete",
                "case": case.name,
                "rows": len(rows),
                "solve_count": solve_count,
                "period_gates": case_gates,
            },
        )
        paths[case.name] = rows
        gates[case.name] = case_gates
        solve_counts[case.name] = solve_count

    effects = build_effect_rows(
        paths["rebated-tax1-baseline"], paths["rebated-tax2-reform"]
    )
    maximum_pre_2043_population_gap = max(
        abs(float(row["adult_population_percent_change"]))
        for row in effects
        if int(row["calendar_year"]) < 2043
    )
    if maximum_pre_2043_population_gap > 1.0e-10:
        raise RuntimeError(
            "Population differs before the twenty-year birth-to-entry lag: "
            f"{maximum_pre_2043_population_gap:.3e} percent"
        )

    all_rows = [row for case_rows in paths.values() for row in case_rows]
    reproduction_rows = reconstruction.pop("history_reproduction_rows")
    baseline.write_csv(outdir / "shared_2007_2019_history.csv", history)
    baseline.write_csv(outdir / "history_reproduction_audit.csv", reproduction_rows)
    baseline.write_csv(outdir / "policy_paths.csv", all_rows)
    baseline.write_csv(outdir / "effects_by_date.csv", effects)
    make_figure(paths, outdir)
    summary = {
        "status": "complete_exact_rebated_property_tax_path",
        "scope": (
            "permanent rebated 1% baseline versus permanent rebated 2% reform; "
            "fitted 2023 initial state; closed M=0,rho=1; no grant or welfare claim"
        ),
        "calendar_years": [int(row["calendar_year"]) for row in effects],
        "post_2023_periods": post_periods,
        "housing_supply_elasticity": elasticity,
        "housing_supply_elasticity_source": reanchor["supply_elasticity_source"],
        "common_initial_state": {
            "distribution_sha256": common_distribution_sha,
            "queue_sha256": common_queue_sha,
            "adult_population": mass_2023,
        },
        "supply_reanchor": reanchor,
        "population_lag": tax.renewal_lag_gate(),
        "maximum_pre_2043_population_percent_gap": maximum_pre_2043_population_gap,
        "case_solve_counts": solve_counts,
        "case_gates": gates,
        "effects_by_date": effects,
        "reconstruction": reconstruction,
        "input_hashes": contracts["hashes"],
        "elapsed_seconds": time.perf_counter() - started,
    }
    baseline.write_json(outdir / "summary.json", summary)
    print(
        f"REBATED_TAX_PATH_COMPLETE dates={len(effects)} "
        f"terminal_year={effects[-1]['calendar_year']} "
        f"seconds={summary['elapsed_seconds']:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
