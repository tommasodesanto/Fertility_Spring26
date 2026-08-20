#!/usr/bin/env python3
"""Exact one-date property-tax/rebate smoke from the fitted 2023 state.

The three cases share the certified fitted 2023 pre-fertility distribution and
four-slot birth-vintage queue.  The two rebated cases jointly solve the housing
price and equal per-household transfer.  Population is nationally closed after
2023 (``M=0, rho=1``).  Purchase grants are deliberately outside this driver.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as baseline
import run_e5f_post2023_policy_mechanisms as mechanism


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = ROOT / "output/model/e5f_post2023_rebated_property_tax_smoke"
MARKET_TOLERANCE = 2.0e-4
JOINT_ROOT_TOLERANCE = 1.0e-4
FISCAL_ABSOLUTE_TOLERANCE = 2.5e-5
MASS_TOLERANCE = 2.0e-8


@dataclass(frozen=True)
class TaxCase:
    name: str
    label: str
    annual_tax_rate: float
    rebate_revenue: bool


CASES = {
    "status-quo-tax1-unrebated": TaxCase(
        name="status-quo-tax1-unrebated",
        label="Status quo: 1% property tax, revenue not rebated",
        annual_tax_rate=0.01,
        rebate_revenue=False,
    ),
    "tax1-equal-rebate": TaxCase(
        name="tax1-equal-rebate",
        label="1% property tax, equal rebate",
        annual_tax_rate=0.01,
        rebate_revenue=True,
    ),
    "tax2-equal-rebate": TaxCase(
        name="tax2-equal-rebate",
        label="2% property tax, equal rebate",
        annual_tax_rate=0.02,
        rebate_revenue=True,
    ),
}


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
    parser.add_argument("--market-tol", type=float, default=MARKET_TOLERANCE)
    parser.add_argument("--market-max-iter", type=int, default=30)
    return parser.parse_args()


def set_tax_policy(P: SimpleNamespace, case: TaxCase, transfer: float) -> None:
    """Set only the property-tax and equal-transfer primitives."""
    P.tau_H = float(case.annual_tax_rate) * float(P.period_years)
    P.user_cost_rate = float(P.q) + float(P.delta) + float(P.tau_H)
    P.property_tax_lump_sum_transfer = float(transfer)
    P.birth_entry_grant = False
    P.birth_entry_grant_amount = 0.0
    P.birth_entry_grant_locations = np.array([], dtype=int)
    P.birth_entry_grant_owner_rungs = np.array([], dtype=int)
    P.propagate_birth_entry_grant = False


def fiscal_ledger(
    evaluation: calendar.PeriodEvaluation,
    P: SimpleNamespace,
) -> dict[str, float]:
    revenue = float(
        calendar.model.property_tax_revenue_from_distribution(
            evaluation.g_current,
            evaluation.policy.hR_pol,
            evaluation.policy.price,
            P,
        )
    )
    current_mass = float(np.sum(evaluation.g_current))
    transfer = float(P.property_tax_lump_sum_transfer)
    outlays = transfer * current_mass
    residual = revenue - outlays
    return {
        "property_tax_revenue": revenue,
        "equal_transfer_period_units": transfer,
        "equal_transfer_outlays": outlays,
        "government_budget_residual": residual,
        "scaled_government_budget_residual": residual
        / max(abs(revenue), abs(outlays), 1.0e-12),
        "implied_equal_transfer": revenue / max(current_mass, 1.0e-12),
    }


def solve_joint_rebated_period(
    *,
    state: baseline.DynamicState,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    supply_rule: calendar.HousingSupplyRule,
    counter: calendar.SolveCounter,
    case: TaxCase,
    tolerance: float = JOINT_ROOT_TOLERANCE,
    initial_price: float | None = None,
    initial_transfer: float | None = None,
) -> SimpleNamespace:
    """Jointly solve log price and the equal balanced-budget transfer."""
    if not case.rebate_revenue:
        raise ValueError("The joint fiscal root is only for rebated cases")
    cache: dict[tuple[float, float], SimpleNamespace] = {}

    def evaluate(point: np.ndarray) -> SimpleNamespace:
        log_price = float(point[0])
        transfer = float(point[1])
        price = float(math.exp(log_price))
        key = (round(log_price, 10), round(transfer, 10))
        if key not in cache:
            set_tax_policy(P, case, transfer)
            shared = calendar.model.precompute_shared(P, b_grid)
            period = calendar.evaluate_period(
                np.array([price]),
                state.g_pre,
                P,
                b_grid,
                shared,
                counter,
                supply_rule=supply_rule,
            )
            ledger = fiscal_ledger(period, P)
            supply = float(period.supply_by_loc[0])
            market_residual = float(period.demand_by_loc[0] - supply) / max(
                abs(supply), 1.0e-14
            )
            fiscal_scale = max(
                abs(float(ledger["property_tax_revenue"])),
                abs(float(ledger["equal_transfer_outlays"])),
                0.1,
            )
            cache[key] = SimpleNamespace(
                point=np.array([log_price, transfer]),
                price=price,
                transfer=transfer,
                evaluation=period,
                shared=shared,
                ledger=ledger,
                residual=np.array(
                    [
                        market_residual,
                        float(ledger["government_budget_residual"]) / fiscal_scale,
                    ]
                ),
            )
        return cache[key]

    price_guess = float(state.price_guess if initial_price is None else initial_price)
    approximate_transfer = (
        float(case.annual_tax_rate)
        * float(P.period_years)
        * price_guess
        * float(supply_rule.quantity(np.array([price_guess]))[0])
        / max(float(np.sum(state.g_pre)), 1.0e-12)
    )
    transfer_guess = (
        approximate_transfer if initial_transfer is None else float(initial_transfer)
    )
    point = np.array([math.log(price_guess), transfer_guess])
    result = evaluate(point)
    steps = np.array([0.01, 0.01])

    def finite_difference_jacobian(at: np.ndarray, base_result: SimpleNamespace) -> np.ndarray:
        return np.column_stack(
            [
                (evaluate(at + np.array([steps[0], 0.0])).residual - base_result.residual)
                / steps[0],
                (evaluate(at + np.array([0.0, steps[1]])).residual - base_result.residual)
                / steps[1],
            ]
        )

    jacobian = finite_difference_jacobian(point, result)
    iterations = 0
    for iteration in range(16):
        if not np.all(np.isfinite(result.residual)):
            raise RuntimeError(f"Nonfinite joint residual: {result.residual.tolist()}")
        normalized_residuals_pass = (
            float(np.max(np.abs(result.residual))) <= float(tolerance)
        )
        absolute_fiscal_residual_passes = (
            abs(float(result.ledger["government_budget_residual"]))
            <= FISCAL_ABSOLUTE_TOLERANCE
        )
        if normalized_residuals_pass and absolute_fiscal_residual_passes:
            break
        try:
            direction = np.linalg.solve(jacobian, -result.residual)
        except np.linalg.LinAlgError:
            direction = np.linalg.lstsq(jacobian, -result.residual, rcond=None)[0]
        direction = np.clip(direction, [-0.25, -0.15], [0.25, 0.15])
        current_norm = float(np.linalg.norm(result.residual))
        accepted: SimpleNamespace | None = None
        accepted_point: np.ndarray | None = None
        for damping in (1.0, 0.5, 0.25, 0.125):
            trial_point = point + damping * direction
            trial_point[0] = float(
                np.clip(
                    trial_point[0],
                    math.log(max(float(getattr(P, "p_min", 1.0e-4)), 1.0e-8)),
                    math.log(float(getattr(P, "p_max", 100.0))),
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
            jacobian = finite_difference_jacobian(point, result)
            continue
        displacement = accepted_point - point
        residual_change = accepted.residual - result.residual
        denominator = float(displacement @ displacement)
        if denominator > 1.0e-14:
            jacobian += np.outer(
                residual_change - jacobian @ displacement, displacement
            ) / denominator
        point = accepted_point
        result = accepted
        iterations = iteration + 1

    if float(np.max(np.abs(result.residual))) > float(tolerance):
        raise RuntimeError(
            "Housing and fiscal markets did not clear jointly: "
            f"residuals={result.residual.tolist()}, evaluations={len(cache)}"
        )
    if abs(float(result.ledger["government_budget_residual"])) > FISCAL_ABSOLUTE_TOLERANCE:
        raise RuntimeError(f"Absolute fiscal gate failed: {result.ledger}")
    result.joint_iterations = iterations
    result.joint_model_evaluations = len(cache)
    return result


def run_impact_case(
    *,
    case: TaxCase,
    prepared: baseline.PreparedModel,
    state_2023: baseline.DynamicState,
    parameters_2023: SimpleNamespace,
    mass_2023: float,
    market_tol: float,
    market_max_iter: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    state = baseline.DynamicState(
        g_pre=state_2023.g_pre.copy(),
        scheduled_entries=list(state_2023.scheduled_entries),
        scheduled_raw_entries=list(state_2023.scheduled_raw_entries),
        price_guess=float(state_2023.price_guess),
        initial_policy=None,
    )
    P = copy.deepcopy(parameters_2023)
    counter = calendar.SolveCounter()
    if case.rebate_revenue:
        solved = solve_joint_rebated_period(
            state=state,
            P=P,
            b_grid=prepared.b_grid,
            supply_rule=prepared.supply_rule,
            counter=counter,
            case=case,
        )
        evaluation = solved.evaluation
        shared = solved.shared
        ledger = solved.ledger
        fallback = False
        joint_iterations = int(solved.joint_iterations)
        joint_evaluations = int(solved.joint_model_evaluations)
    else:
        set_tax_policy(P, case, 0.0)
        evaluation, shared, fallback = baseline.evaluate_state(
            state,
            P,
            prepared.b_grid,
            prepared.supply_rule,
            counter,
            market_tol,
            market_max_iter,
        )
        ledger = fiscal_ledger(evaluation, P)
        joint_iterations = 0
        joint_evaluations = 0

    row, next_state = baseline.advance_from_evaluation(
        label=case.name,
        period_from_2007=baseline.TRANSITION_PERIODS,
        evaluation=evaluation,
        state=state,
        P=P,
        b_grid=prepared.b_grid,
        shared=shared,
        supply_rule=prepared.supply_rule,
        outside_flow=0.0,
        retention=1.0,
        initial_mass_2007=prepared.initial_mass_2007,
        mass_2023=mass_2023,
        next_bridge_year=None,
        grid_fallback=fallback,
    )
    row.update(ledger)
    row.update(
        policy_case=case.name,
        policy_label=case.label,
        policy_active=case.rebate_revenue,
        annual_property_tax_rate=case.annual_tax_rate,
        property_tax_rate_period_units=case.annual_tax_rate * float(P.period_years),
        revenue_rebated=case.rebate_revenue,
        revenue_disposition=("equal household rebate" if case.rebate_revenue else "not rebated"),
        fiscal_budget_gate_applicable=case.rebate_revenue,
        model_solve_count=int(counter.total),
        joint_iterations=joint_iterations,
        joint_model_evaluations=joint_evaluations,
    )
    gates = {
        "market_residual": abs(float(row["relative_market_residual"])),
        "mass_residual": abs(float(row["mass_accounting_residual"])),
        "fiscal_residual": (
            abs(float(row["government_budget_residual"]))
            if case.rebate_revenue
            else None
        ),
        "nonfinite_distribution_count": int(row["nonfinite_distribution_count"]),
        "minimum_distribution_mass": float(row["min_distribution_mass"]),
        "feasibility_projection_mass": float(row["feasibility_frontier_projection_mass"]),
        "next_queue_sha256": baseline.canonical_json_sha256(next_state.scheduled_entries),
    }
    if gates["market_residual"] > float(market_tol):
        raise RuntimeError(f"{case.name}: market gate failed: {gates}")
    if gates["mass_residual"] > MASS_TOLERANCE:
        raise RuntimeError(f"{case.name}: mass gate failed: {gates}")
    if case.rebate_revenue and float(gates["fiscal_residual"]) > FISCAL_ABSOLUTE_TOLERANCE:
        raise RuntimeError(f"{case.name}: fiscal gate failed: {gates}")
    if gates["nonfinite_distribution_count"] != 0:
        raise RuntimeError(f"{case.name}: nonfinite distribution")
    if gates["minimum_distribution_mass"] < -1.0e-14:
        raise RuntimeError(f"{case.name}: negative distribution mass")
    return row, gates


def renewal_lag_gate() -> dict[str, Any]:
    queue = [1.0, 2.0, 3.0, 4.0]
    scheduled, updated = transition.advance_birth_vintage_queue(queue, 9.0, 1.0 / 2.1)
    passed = (
        baseline.QUEUE_WAITING_SLOTS == 4
        and math.isclose(scheduled, 1.0)
        and np.allclose(updated, [2.0, 3.0, 4.0, 9.0 / 2.1])
    )
    if not passed:
        raise RuntimeError("The four-slot births/2.1 renewal lag contract failed")
    return {
        "passed": True,
        "waiting_slots": baseline.QUEUE_WAITING_SLOTS,
        "birth_to_population_effect_model_dates": baseline.QUEUE_WAITING_SLOTS + 1,
        "birth_to_population_effect_years": int(
            (baseline.QUEUE_WAITING_SLOTS + 1)
            * (baseline.TERMINAL_YEAR - baseline.START_YEAR)
            / baseline.TRANSITION_PERIODS
        ),
        "conversion": 1.0 / 2.1,
    }


def main() -> None:
    args = parse_args()
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
    initial_distribution_sha = baseline.array_sha256(state_2023.g_pre)
    initial_queue_sha = baseline.canonical_json_sha256(state_2023.scheduled_entries)

    rows: list[dict[str, Any]] = []
    case_gates: dict[str, Any] = {}
    for case in CASES.values():
        print(f"TAX_SMOKE_CASE_START case={case.name}", flush=True)
        row, gates = run_impact_case(
            case=case,
            prepared=prepared,
            state_2023=state_2023,
            parameters_2023=parameters_2023,
            mass_2023=mass_2023,
            market_tol=float(args.market_tol),
            market_max_iter=int(args.market_max_iter),
        )
        row["common_initial_distribution_sha256"] = initial_distribution_sha
        row["common_initial_queue_sha256"] = initial_queue_sha
        rows.append(row)
        case_gates[case.name] = gates
        baseline.write_csv(outdir / "impact_results_progress.csv", rows)
        print(
            f"TAX_SMOKE_CASE_COMPLETE case={case.name} price={row['asset_price']:.8f} "
            f"transfer={row['equal_transfer_period_units']:.8f} "
            f"fiscal={row['government_budget_residual']:.3e}",
            flush=True,
        )

    baseline_row = rows[0]
    for row in rows:
        row["asset_price_percent_change_vs_status_quo"] = 100.0 * (
            float(row["asset_price"]) / float(baseline_row["asset_price"]) - 1.0
        )
        row["births_per_adult_percent_change_vs_status_quo"] = 100.0 * (
            float(row["topcode_adjusted_births_per_adult"])
            / float(baseline_row["topcode_adjusted_births_per_adult"])
            - 1.0
        )
        row["owner_rate_pp_change_vs_status_quo"] = 100.0 * (
            float(row["owner_rate"]) - float(baseline_row["owner_rate"])
        )
        row["population_percent_change_vs_status_quo"] = 100.0 * (
            float(row["adult_population"]) / float(baseline_row["adult_population"]) - 1.0
        )

    impact_population_gap = max(
        abs(float(row["adult_population"]) - float(baseline_row["adult_population"]))
        for row in rows
    )
    if impact_population_gap > 1.0e-12:
        raise RuntimeError(
            "Policy changes contemporaneous adult population despite the renewal lag: "
            f"gap={impact_population_gap:.3e}"
        )
    lag = renewal_lag_gate()
    reconstruction_rows = reconstruction.pop("history_reproduction_rows")
    baseline.write_csv(outdir / "shared_2007_2019_history.csv", history)
    baseline.write_csv(outdir / "history_reproduction_audit.csv", reconstruction_rows)
    baseline.write_csv(outdir / "impact_results.csv", rows)
    summary = {
        "status": "complete_exact_rebated_property_tax_smoke",
        "scope": (
            "one-date 2023 impact comparison from the fitted state; closed national "
            "population M=0,rho=1; no purchase grant and no welfare claim"
        ),
        "cases": list(CASES),
        "current_status_quo": "1% annual property tax with revenue not rebated",
        "rebate_experiment": (
            "return contemporaneous property-tax revenue as an equal transfer to each "
            "current adult household"
        ),
        "common_initial_state": {
            "distribution_sha256": initial_distribution_sha,
            "queue_sha256": initial_queue_sha,
            "adult_population": mass_2023,
        },
        "population_closure": "closed national benchmark: M=0, rho=1",
        "renewal_contract": lag,
        "impact_population_identity_max_abs_gap": impact_population_gap,
        "reconstruction": reconstruction,
        "case_gates": case_gates,
        "impact_results": rows,
        "input_hashes": contracts["hashes"],
        "elapsed_seconds": time.perf_counter() - started,
    }
    baseline.write_json(outdir / "summary.json", summary)
    print(
        f"TAX_SMOKE_COMPLETE cases={len(rows)} seconds={summary['elapsed_seconds']:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
