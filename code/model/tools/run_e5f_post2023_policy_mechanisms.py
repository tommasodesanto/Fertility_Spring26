#!/usr/bin/env python3
"""Finite-horizon housing-policy mechanisms from the fitted 2023 state.

Each case reconstructs the certified 2007--2023 calibration, verifies the
reconstructed history, and then applies one transparent housing intervention
to the same 2023 pre-fertility distribution and birth-vintage queue.  The
post-2023 population law is the closed national benchmark, ``M=0, rho=1``.

This is a mechanism experiment, not a welfare calculation.  Property-tax
revenue is discarded rather than rebated, and no purchase grant is available.
The exercise makes no claim that the terminal date is a steady state.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as baseline


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = ROOT / "output/model/e5f_post2023_policy_mechanisms"
POLICY_START_YEAR = 2023
DEFAULT_POST_PERIODS = 10


@dataclass(frozen=True)
class PolicySpec:
    name: str
    label: str
    supply_multiplier: float
    dependent_child_ltv95: bool
    annual_property_tax_rate: float


POLICIES = {
    "baseline": PolicySpec(
        name="baseline",
        label="No policy",
        supply_multiplier=1.0,
        dependent_child_ltv95=False,
        annual_property_tax_rate=0.01,
    ),
    "supply-plus-20": PolicySpec(
        name="supply-plus-20",
        label="Housing supply +20%",
        supply_multiplier=1.2,
        dependent_child_ltv95=False,
        annual_property_tax_rate=0.01,
    ),
    "dependent-child-ltv95": PolicySpec(
        name="dependent-child-ltv95",
        label="5% down payment for dependent-child households",
        supply_multiplier=1.0,
        dependent_child_ltv95=True,
        annual_property_tax_rate=0.01,
    ),
    "combined": PolicySpec(
        name="combined",
        label="Supply +20% and family credit",
        supply_multiplier=1.2,
        dependent_child_ltv95=True,
        annual_property_tax_rate=0.01,
    ),
    "property-tax-2pct-no-rebate": PolicySpec(
        name="property-tax-2pct-no-rebate",
        label="Property tax: 1% to 2% (no rebate)",
        supply_multiplier=1.0,
        dependent_child_ltv95=False,
        annual_property_tax_rate=0.02,
    ),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, required=True)
    parser.add_argument("--selected-task-summary", type=Path)
    parser.add_argument("--selected-case-dir", type=Path)
    parser.add_argument("--selected-case-transition", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--policy-case", choices=tuple(POLICIES), required=True)
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--expected-task-summary-sha256")
    parser.add_argument("--expected-case-transition-sha256", required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-target-fingerprint", required=True)
    parser.add_argument("--expected-code-bundle-sha256", required=True)
    parser.add_argument("--expected-renewal-contract-sha256", required=True)
    parser.add_argument("--expected-scientific-contract-sha256", required=True)
    parser.add_argument("--expected-selection-sha256", required=True)
    parser.add_argument("--post-2023-periods", type=int, default=DEFAULT_POST_PERIODS)
    parser.add_argument("--market-tol", type=float, default=2e-4)
    parser.add_argument("--market-max-iter", type=int, default=30)
    return parser.parse_args()


def policy_supply_rule(
    rule: calendar.HousingSupplyRule, spec: PolicySpec
) -> calendar.HousingSupplyRule:
    if rule.mode != "static-elastic":
        raise RuntimeError(
            "The policy packet requires the calibrated static-elastic supply rule"
        )
    return calendar.HousingSupplyRule(
        mode=rule.mode,
        initial_price=float(rule.initial_price),
        initial_stock=float(rule.initial_stock) * float(spec.supply_multiplier),
        elasticity=float(rule.elasticity),
    )


def apply_policy(parameters: Any, spec: PolicySpec) -> None:
    parameters.tau_H = float(spec.annual_property_tax_rate) * float(
        parameters.period_years
    )
    parameters.user_cost_rate = (
        float(parameters.q) + float(parameters.delta) + float(parameters.tau_H)
    )
    parameters.property_tax_lump_sum_transfer = 0.0
    parameters.birth_entry_grant = False
    parameters.birth_entry_grant_amount = 0.0
    parameters.birth_entry_grant_locations = np.array([], dtype=int)
    parameters.birth_entry_grant_owner_rungs = np.array([], dtype=int)
    parameters.propagate_birth_entry_grant = True
    parameters.parent_dp_waiver = bool(spec.dependent_child_ltv95)
    parameters.parent_dp_waiver_phi = (
        0.95 if spec.dependent_child_ltv95 else 1.0
    )
    parameters.parent_dp_waiver_locations = np.array([], dtype=int)
    parameters.parent_dp_waiver_owner_rungs = np.array([], dtype=int)
    parameters.parent_dp_waiver_birth_state_only = False


def reconstruct_matched_2023(
    prepared: baseline.PreparedModel,
    contracts: dict[str, Any],
    *,
    market_tol: float,
    market_max_iter: int,
    progress_dir: Path,
) -> tuple[
    list[dict[str, Any]],
    baseline.DynamicState,
    Any,
    float,
    dict[str, Any],
]:
    """Reconstruct the certified history and return the unadvanced 2023 state."""
    best = contracts["report_best"]
    old_psi = float(best["old_psi_child"])
    new_psi = float(best["new_psi_child"])
    old = contracts["renewal_old_state"]
    outside_flow = float(old["outside_flow_M"])
    retention = float(old["retention_rho"])
    raw_B_old = (
        transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        * prepared.births_old_raw
    )
    state = baseline.DynamicState(
        g_pre=prepared.initial_g_pre_2007.copy(),
        scheduled_entries=[prepared.B_old] * baseline.QUEUE_WAITING_SLOTS,
        scheduled_raw_entries=[raw_B_old] * baseline.QUEUE_WAITING_SLOTS,
        price_guess=prepared.old_price,
        initial_policy=None,
    )
    parameters = copy.deepcopy(prepared.parameters)
    counter = calendar.SolveCounter()
    age_target_years = transition.census_household_age_target_years()
    history: list[dict[str, Any]] = []
    previous_psi: float | None = None
    for period in range(baseline.TRANSITION_PERIODS):
        dated_psi = transition.preference_shifter_at_date(
            period,
            old_psi,
            new_psi,
            baseline.TRANSITION_PERIODS,
        )
        parameters.psi_child = dated_psi
        apply_policy(parameters, POLICIES["baseline"])
        if previous_psi is not None and not math.isclose(
            dated_psi, previous_psi, rel_tol=0.0, abs_tol=1e-14
        ):
            state.initial_policy = None
        previous_psi = dated_psi
        evaluation, shared, fallback = baseline.evaluate_state(
            state,
            parameters,
            prepared.b_grid,
            prepared.supply_rule,
            counter,
            market_tol,
            market_max_iter,
        )
        row, state = baseline.advance_from_evaluation(
            label="shared_calibrated_history",
            period_from_2007=period,
            evaluation=evaluation,
            state=state,
            P=parameters,
            b_grid=prepared.b_grid,
            shared=shared,
            supply_rule=prepared.supply_rule,
            outside_flow=outside_flow,
            retention=retention,
            initial_mass_2007=prepared.initial_mass_2007,
            mass_2023=1.0,
            next_bridge_year=age_target_years[period + 1],
            grid_fallback=fallback,
        )
        history.append(row)
        baseline.write_csv(progress_dir / "history_progress.csv", history)

    state_2023 = baseline.DynamicState(
        g_pre=state.g_pre.copy(),
        scheduled_entries=list(state.scheduled_entries),
        scheduled_raw_entries=list(state.scheduled_raw_entries),
        price_guess=float(state.price_guess),
        initial_policy=None,
    )
    mass_2023 = float(np.sum(state_2023.g_pre))
    parameters.psi_child = new_psi
    apply_policy(parameters, POLICIES["baseline"])
    audit_eval, audit_shared, audit_fallback = baseline.evaluate_state(
        state_2023,
        parameters,
        prepared.b_grid,
        prepared.supply_rule,
        counter,
        market_tol,
        market_max_iter,
    )
    audit_row, _ = baseline.advance_from_evaluation(
        label="matched_2023_no_policy_audit",
        period_from_2007=baseline.TRANSITION_PERIODS,
        evaluation=audit_eval,
        state=state_2023,
        P=parameters,
        b_grid=prepared.b_grid,
        shared=audit_shared,
        supply_rule=prepared.supply_rule,
        outside_flow=outside_flow,
        retention=retention,
        initial_mass_2007=prepared.initial_mass_2007,
        mass_2023=mass_2023,
        next_bridge_year=None,
        grid_fallback=audit_fallback,
    )
    reproduction = baseline.validate_reconstructed_history(
        history,
        [audit_row],
        Path(contracts["paths"]["selected_case_transition"]),
    )
    return history, state_2023, parameters, mass_2023, {
        "history_reproduction": {
            key: value for key, value in reproduction.items() if key != "comparisons"
        },
        "history_reproduction_rows": reproduction["comparisons"],
        "matched_2023_distribution_sha256": baseline.array_sha256(state_2023.g_pre),
        "matched_2023_queue_sha256": baseline.canonical_json_sha256(
            state_2023.scheduled_entries
        ),
        "solve_count_through_2023_audit": counter.total,
    }


def run_policy_path(
    prepared: baseline.PreparedModel,
    state_2023: baseline.DynamicState,
    parameters_2023: Any,
    mass_2023: float,
    spec: PolicySpec,
    *,
    post_2023_periods: int,
    market_tol: float,
    market_max_iter: int,
    progress_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    state = baseline.DynamicState(
        g_pre=state_2023.g_pre.copy(),
        scheduled_entries=list(state_2023.scheduled_entries),
        scheduled_raw_entries=list(state_2023.scheduled_raw_entries),
        price_guess=float(state_2023.price_guess),
        initial_policy=None,
    )
    parameters = copy.deepcopy(parameters_2023)
    parameters.psi_child = float(parameters_2023.psi_child)
    calibrated_annual_tax = float(parameters.tau_H) / float(parameters.period_years)
    if not math.isclose(calibrated_annual_tax, 0.01, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(
            "Expected the fitted 2023 state to carry the calibrated 1 percent "
            f"property tax; found {calibrated_annual_tax}"
        )
    apply_policy(parameters, spec)
    supply_rule = policy_supply_rule(prepared.supply_rule, spec)
    counter = calendar.SolveCounter()
    rows: list[dict[str, Any]] = []
    for post_index in range(post_2023_periods + 1):
        period = baseline.TRANSITION_PERIODS + post_index
        evaluation, shared, fallback = baseline.evaluate_state(
            state,
            parameters,
            prepared.b_grid,
            supply_rule,
            counter,
            market_tol,
            market_max_iter,
        )
        property_tax_revenue = float(
            calendar.model.property_tax_revenue_from_distribution(
                evaluation.g_current,
                evaluation.policy.hR_pol,
                evaluation.policy.price,
                parameters,
            )
        )
        row, state = baseline.advance_from_evaluation(
            label=spec.name,
            period_from_2007=period,
            evaluation=evaluation,
            state=state,
            P=parameters,
            b_grid=prepared.b_grid,
            shared=shared,
            supply_rule=supply_rule,
            outside_flow=0.0,
            retention=1.0,
            initial_mass_2007=prepared.initial_mass_2007,
            mass_2023=mass_2023,
            next_bridge_year=None,
            grid_fallback=fallback,
        )
        row["policy_case"] = spec.name
        row["policy_label"] = spec.label
        row["supply_intercept_multiplier"] = spec.supply_multiplier
        row["dependent_child_ltv95"] = spec.dependent_child_ltv95
        row["annual_property_tax_rate"] = spec.annual_property_tax_rate
        row["housing_user_cost"] = (
            float(parameters.user_cost_rate) * float(row["asset_price"])
        )
        row["property_tax_revenue_discarded"] = property_tax_revenue
        row["property_tax_revenue_discarded_per_adult"] = (
            property_tax_revenue
            / max(float(np.sum(evaluation.g_current)), 1e-15)
        )
        row["property_tax_lump_sum_transfer"] = 0.0
        row["purchase_grant"] = False
        rows.append(row)
        baseline.write_csv(progress_dir / "policy_path_progress.csv", rows)
        baseline.write_json(
            progress_dir / "latest_completed_period.json",
            {
                "status": "running",
                "policy_case": spec.name,
                "calendar_year": int(row["calendar_year"]),
                "completed_dates": len(rows),
                "requested_dates": post_2023_periods + 1,
                "latest": row,
            },
        )
        print(
            f"POLICY_PATH case={spec.name} year={row['calendar_year']} "
            f"price={row['asset_price']:.6f} pop={row['population_index_2023']:.6f} "
            f"births_per_adult={row['topcode_adjusted_births_per_adult']:.6f}",
            flush=True,
        )
    gates = {
        "solve_count": counter.total,
        "maximum_market_residual": max(
            abs(float(row["relative_market_residual"])) for row in rows
        ),
        "maximum_mass_residual": max(
            abs(float(row["mass_accounting_residual"])) for row in rows
        ),
        "maximum_projection_mass": max(
            abs(float(row["feasibility_frontier_projection_mass"])) for row in rows
        ),
        "maximum_nonfinite_count": max(
            int(row["nonfinite_distribution_count"]) for row in rows
        ),
        "minimum_distribution_mass": min(
            float(row["min_distribution_mass"]) for row in rows
        ),
    }
    return rows, gates


def make_case_figure(rows: list[dict[str, Any]], outdir: Path) -> None:
    years = [int(row["calendar_year"]) for row in rows]
    metrics = (
        ("asset_price_index", "Housing cost (2007=1)"),
        ("owner_rate", "Ownership rate"),
        ("topcode_adjusted_births_per_adult", "Adjusted births per adult"),
        ("population_index_2023", "Adult households (2023=1)"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.5), constrained_layout=True)
    for axis, (field, label) in zip(axes.flat, metrics):
        axis.plot(years, [float(row[field]) for row in rows], marker="o")
        axis.set_title(label)
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    fig.suptitle(str(rows[0]["policy_label"]))
    fig.savefig(outdir / "policy_path.png", dpi=220)
    fig.savefig(outdir / "policy_path.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not 1 <= int(args.post_2023_periods) <= 40:
        raise ValueError("--post-2023-periods must lie in [1, 40]")
    if not 0.0 < float(args.market_tol) <= 2e-4:
        raise ValueError("--market-tol must lie in (0, 2e-4]")
    outdir = args.output_dir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    chain, model = transition.configure_sequential_model()
    contracts = baseline.validate_input_contracts(
        report_path=args.selected_report,
        task_summary_path=args.selected_task_summary,
        case_dir=args.selected_case_dir,
        case_transition_path=args.selected_case_transition,
        source_path=args.source,
        expected_report_sha256=args.expected_report_sha256,
        expected_task_summary_sha256=args.expected_task_summary_sha256,
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
        reconstruct_matched_2023(
            prepared,
            contracts,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
            progress_dir=outdir,
        )
    )
    spec = POLICIES[str(args.policy_case)]
    rows, gates = run_policy_path(
        prepared,
        state_2023,
        parameters_2023,
        mass_2023,
        spec,
        post_2023_periods=int(args.post_2023_periods),
        market_tol=float(args.market_tol),
        market_max_iter=int(args.market_max_iter),
        progress_dir=outdir,
    )
    baseline.write_csv(outdir / "shared_2007_2019_history.csv", history)
    baseline.write_csv(
        outdir / "history_reproduction_audit.csv",
        reconstruction.pop("history_reproduction_rows"),
    )
    baseline.write_csv(outdir / "policy_path.csv", rows)
    make_case_figure(rows, outdir)

    summary = {
        "status": "complete_post2023_policy_mechanism_case",
        "scope": (
            "finite-horizon closed-national mechanism experiment from the fitted "
            "2023 state; property-tax revenue is discarded; no welfare or "
            "steady-state claim"
        ),
        "policy_start_year": POLICY_START_YEAR,
        "policy": {
            "case": spec.name,
            "label": spec.label,
            "supply_intercept_multiplier": spec.supply_multiplier,
            "dependent_child_ltv95": spec.dependent_child_ltv95,
            "annual_property_tax_rate": spec.annual_property_tax_rate,
            "fiscal_treatment": (
                "property-tax revenue is discarded; no rebate and no purchase grant"
            ),
            "ltv_interpretation": (
                "When active, households with dependent children may finance 95% "
                "of an owned housing choice; this is not restricted to first-time buyers."
            ),
        },
        "population_closure": "closed national benchmark: M=0, rho=1",
        "initial_state": (
            "exact common fitted 2023 pre-fertility distribution and inherited "
            "four-slot birth-vintage queue"
        ),
        "renewal_lag": (
            "A policy-induced birth flow first changes adult-household entry twenty "
            "years later."
        ),
        "post_2023_periods": int(args.post_2023_periods),
        "terminal_year": int(rows[-1]["calendar_year"]),
        "input_hashes": contracts["hashes"],
        "driver_sha256": sha256_file(Path(__file__).resolve()),
        "reconstruction": reconstruction,
        "numerical_gates": gates,
        "impact_2023": rows[0],
        "year_2043": next(
            (row for row in rows if int(row["calendar_year"]) == 2043), None
        ),
        "terminal": rows[-1],
        "elapsed_seconds": time.perf_counter() - started,
    }
    baseline.write_json(outdir / "summary.json", summary)
    artifact_names = (
        "shared_2007_2019_history.csv",
        "history_reproduction_audit.csv",
        "policy_path.csv",
        "policy_path.png",
        "policy_path.pdf",
        "summary.json",
    )
    manifest = {
        "schema": "e5f_post2023_policy_mechanism_manifest_v1",
        "artifacts": {
            name: sha256_file(outdir / name) for name in artifact_names
        },
    }
    baseline.write_json(outdir / "manifest.json", manifest)
    baseline.write_json(
        outdir / "latest_completed_period.json",
        {
            "status": "complete",
            "policy_case": spec.name,
            "terminal_year": int(rows[-1]["calendar_year"]),
            "summary_sha256": sha256_file(outdir / "summary.json"),
            "manifest_sha256": sha256_file(outdir / "manifest.json"),
        },
    )
    print(
        f"POLICY_CASE_COMPLETE case={spec.name} dates={len(rows)} "
        f"elapsed={summary['elapsed_seconds']:.2f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
