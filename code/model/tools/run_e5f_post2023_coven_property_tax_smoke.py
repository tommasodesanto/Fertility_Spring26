#!/usr/bin/env python3
"""Exact 2023 Coven-style property-tax impact benchmark.

Both cases rebate all contemporaneous property-tax revenue equally across
current adult households.  They share the fitted 2023 pre-fertility state and
birth-vintage queue and a nationally closed population contract (M=0, rho=1).
The default housing-supply elasticity is the California estimate 0.232 from
Baum-Snow and Han (2024), as used by Coven et al. (2025).  The supply rule is
re-anchored so the rebated 1 percent baseline clears at the inherited fitted
2023 price.  Alternative elasticities are explicit diagnostic sensitivities.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as baseline
import run_e5f_post2023_policy_mechanisms as mechanism
import run_e5f_post2023_rebated_property_tax_smoke as tax


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = ROOT / "output/model/e5f_post2023_coven_property_tax_smoke"
COVEN_CALIFORNIA_SUPPLY_ELASTICITY = 0.232
CASES = {
    "rebated-tax1-baseline": tax.TaxCase(
        name="rebated-tax1-baseline",
        label="Rebated 1% property tax",
        annual_tax_rate=0.01,
        rebate_revenue=True,
    ),
    "rebated-tax2-reform": tax.TaxCase(
        name="rebated-tax2-reform",
        label="Rebated 2% property tax",
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
    parser.add_argument("--market-tol", type=float, default=tax.MARKET_TOLERANCE)
    parser.add_argument("--market-max-iter", type=int, default=30)
    parser.add_argument(
        "--housing-supply-elasticity",
        type=float,
        default=COVEN_CALIFORNIA_SUPPLY_ELASTICITY,
        help=(
            "Static housing-supply elasticity. The default 0.232 is the "
            "Baum-Snow--Han California estimate used by Coven et al.; any "
            "other value is recorded as a diagnostic sensitivity."
        ),
    )
    return parser.parse_args()


def inherited_2023_price(path: Path) -> float:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    matches = [
        row
        for row in rows
        if int(float(row["period"])) == baseline.TRANSITION_PERIODS
        and baseline.START_YEAR + int(float(row["years_from_start"]))
        == baseline.TERMINAL_YEAR
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one fitted 2023 row; found {len(matches)}")
    price = float(matches[0]["asset_price"])
    if not math.isfinite(price) or price <= 0.0:
        raise RuntimeError(f"Invalid inherited 2023 price: {price}")
    return price


def evaluate_fixed_price(
    *,
    state: baseline.DynamicState,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    price: float,
    transfer: float,
    case: tax.TaxCase,
    counter: calendar.SolveCounter,
) -> SimpleNamespace:
    tax.set_tax_policy(P, case, transfer)
    shared = calendar.model.precompute_shared(P, b_grid)
    evaluation = calendar.evaluate_period(
        np.array([price]),
        state.g_pre,
        P,
        b_grid,
        shared,
        counter,
        supply_rule=None,
    )
    return SimpleNamespace(
        transfer=float(transfer),
        evaluation=evaluation,
        shared=shared,
        ledger=tax.fiscal_ledger(evaluation, P),
    )


def solve_fixed_price_equal_rebate(
    *,
    state: baseline.DynamicState,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    price: float,
    case: tax.TaxCase,
    counter: calendar.SolveCounter,
) -> SimpleNamespace:
    """Solve only the baseline fiscal fixed point at the inherited price."""
    cache: dict[float, SimpleNamespace] = {}

    def evaluate(transfer: float) -> SimpleNamespace:
        key = round(float(transfer), 10)
        if key not in cache:
            cache[key] = evaluate_fixed_price(
                state=state,
                P=P,
                b_grid=b_grid,
                price=price,
                transfer=float(transfer),
                case=case,
                counter=counter,
            )
        return cache[key]

    left = evaluate(0.0)
    right_transfer = max(float(left.ledger["implied_equal_transfer"]) * 1.25, 0.05)
    right = evaluate(right_transfer)
    while float(right.ledger["government_budget_residual"]) > 0.0:
        right_transfer *= 2.0
        if right_transfer > 8.0:
            raise RuntimeError("Could not bracket the fixed-price fiscal root")
        right = evaluate(right_transfer)
    lo, hi = 0.0, right_transfer
    flo = float(left.ledger["government_budget_residual"])
    fhi = float(right.ledger["government_budget_residual"])
    result = right
    for iteration in range(24):
        if abs(float(result.ledger["government_budget_residual"])) <= tax.FISCAL_ABSOLUTE_TOLERANCE:
            result.iterations = iteration
            result.model_evaluations = len(cache)
            return result
        if abs(fhi - flo) > 1.0e-14:
            trial = hi - fhi * (hi - lo) / (fhi - flo)
        else:
            trial = 0.5 * (lo + hi)
        if not lo < trial < hi:
            trial = 0.5 * (lo + hi)
        result = evaluate(trial)
        residual = float(result.ledger["government_budget_residual"])
        if residual > 0.0:
            lo, flo = trial, residual
        else:
            hi, fhi = trial, residual
    raise RuntimeError(
        "Fixed-price fiscal root failed: "
        f"residual={result.ledger['government_budget_residual']:.3e}, "
        f"evaluations={len(cache)}"
    )


def reanchor_supply_rule(
    *, inherited_price: float, baseline_demand: float, elasticity: float
) -> tuple[calendar.HousingSupplyRule, dict[str, Any]]:
    if not math.isfinite(elasticity) or elasticity <= 0.0:
        raise ValueError(f"Housing-supply elasticity must be finite and positive: {elasticity}")
    rule = calendar.HousingSupplyRule(
        mode="static-elastic",
        initial_price=float(inherited_price),
        initial_stock=float(baseline_demand),
        elasticity=float(elasticity),
    )
    supply = float(rule.quantity(np.array([inherited_price]))[0])
    gap = abs(supply - baseline_demand) / max(abs(baseline_demand), 1.0e-15)
    if gap > 1.0e-12:
        raise RuntimeError(f"Coven baseline supply re-anchor failed: {gap:.3e}")
    return rule, {
        "passed": True,
        "inherited_2023_asset_price": float(inherited_price),
        "rebated_tax1_baseline_housing_demand": float(baseline_demand),
        "reanchored_housing_stock": float(supply),
        "relative_identity_gap": float(gap),
        "supply_elasticity": float(elasticity),
        "supply_elasticity_source": (
            "Baum-Snow and Han (2024), California, as used by Coven et al. (2025)"
            if abs(float(elasticity) - COVEN_CALIFORNIA_SUPPLY_ELASTICITY) <= 1.0e-15
            else "explicit diagnostic sensitivity"
        ),
    }


def run_case(
    *,
    case: tax.TaxCase,
    prepared: baseline.PreparedModel,
    state_2023: baseline.DynamicState,
    parameters_2023: SimpleNamespace,
    mass_2023: float,
    supply_rule: calendar.HousingSupplyRule,
    housing_supply_elasticity: float,
    inherited_price: float,
    initial_transfer: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    state = baseline.DynamicState(
        g_pre=state_2023.g_pre.copy(),
        scheduled_entries=list(state_2023.scheduled_entries),
        scheduled_raw_entries=list(state_2023.scheduled_raw_entries),
        price_guess=float(inherited_price),
        initial_policy=None,
    )
    P = copy.deepcopy(parameters_2023)
    counter = calendar.SolveCounter()
    solved = tax.solve_joint_rebated_period(
        state=state,
        P=P,
        b_grid=prepared.b_grid,
        supply_rule=supply_rule,
        counter=counter,
        case=case,
        initial_price=float(inherited_price),
        initial_transfer=float(initial_transfer),
    )
    row, next_state = baseline.advance_from_evaluation(
        label=case.name,
        period_from_2007=baseline.TRANSITION_PERIODS,
        evaluation=solved.evaluation,
        state=state,
        P=P,
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
        property_tax_rate_period_units=case.annual_tax_rate * float(P.period_years),
        equal_rebate=True,
        housing_supply_elasticity=float(housing_supply_elasticity),
        housing_user_cost=float(P.user_cost_rate) * float(solved.price),
        user_cost_rate=float(P.user_cost_rate),
        joint_iterations=int(solved.joint_iterations),
        joint_model_evaluations=int(solved.joint_model_evaluations),
        model_solve_count=int(counter.total),
    )
    gates = {
        "market_residual": abs(float(row["relative_market_residual"])),
        "fiscal_residual": abs(float(row["government_budget_residual"])),
        "mass_residual": abs(float(row["mass_accounting_residual"])),
        "nonfinite_distribution_count": int(row["nonfinite_distribution_count"]),
        "minimum_distribution_mass": float(row["min_distribution_mass"]),
        "feasibility_projection_mass": float(row["feasibility_frontier_projection_mass"]),
        "next_queue_sha256": baseline.canonical_json_sha256(next_state.scheduled_entries),
    }
    if gates["market_residual"] > tax.MARKET_TOLERANCE:
        raise RuntimeError(f"{case.name}: market gate failed: {gates}")
    if gates["fiscal_residual"] > tax.FISCAL_ABSOLUTE_TOLERANCE:
        raise RuntimeError(f"{case.name}: fiscal gate failed: {gates}")
    if gates["mass_residual"] > tax.MASS_TOLERANCE:
        raise RuntimeError(f"{case.name}: mass gate failed: {gates}")
    if gates["nonfinite_distribution_count"] != 0:
        raise RuntimeError(f"{case.name}: nonfinite distribution")
    return row, gates


def main() -> None:
    args = parse_args()
    if (
        not math.isfinite(float(args.housing_supply_elasticity))
        or float(args.housing_supply_elasticity) <= 0.0
    ):
        raise ValueError(
            "--housing-supply-elasticity must be finite and positive; got "
            f"{args.housing_supply_elasticity}"
        )
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
    inherited_price = inherited_2023_price(args.selected_case_transition.resolve())
    anchor_counter = calendar.SolveCounter()
    anchor_P = copy.deepcopy(parameters_2023)
    anchor = solve_fixed_price_equal_rebate(
        state=state_2023,
        P=anchor_P,
        b_grid=prepared.b_grid,
        price=inherited_price,
        case=CASES["rebated-tax1-baseline"],
        counter=anchor_counter,
    )
    baseline_demand = float(anchor.evaluation.demand_by_loc[0])
    supply_rule, reanchor = reanchor_supply_rule(
        inherited_price=inherited_price,
        baseline_demand=baseline_demand,
        elasticity=float(args.housing_supply_elasticity),
    )
    reanchor.update(
        fixed_price_equal_transfer=float(anchor.transfer),
        fixed_price_fiscal_residual=float(anchor.ledger["government_budget_residual"]),
        fixed_price_model_evaluations=int(anchor.model_evaluations),
    )
    baseline.write_json(outdir / "supply_reanchor_gate.json", reanchor)
    print(
        "COVEN_REANCHOR_PASS "
        f"price={inherited_price:.8f} stock={baseline_demand:.8f} "
        f"gap={reanchor['relative_identity_gap']:.3e} transfer={anchor.transfer:.8f}",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    gates: dict[str, Any] = {}
    transfer_guess = float(anchor.transfer)
    for case in CASES.values():
        print(f"COVEN_TAX_CASE_START case={case.name}", flush=True)
        row, case_gate = run_case(
            case=case,
            prepared=prepared,
            state_2023=state_2023,
            parameters_2023=parameters_2023,
            mass_2023=mass_2023,
            supply_rule=supply_rule,
            housing_supply_elasticity=float(args.housing_supply_elasticity),
            inherited_price=inherited_price,
            initial_transfer=transfer_guess * (case.annual_tax_rate / 0.01),
        )
        row["common_initial_distribution_sha256"] = baseline.array_sha256(state_2023.g_pre)
        row["common_initial_queue_sha256"] = baseline.canonical_json_sha256(
            state_2023.scheduled_entries
        )
        rows.append(row)
        gates[case.name] = case_gate
        baseline.write_csv(outdir / "impact_results_progress.csv", rows)
        print(
            f"COVEN_TAX_CASE_COMPLETE case={case.name} p={row['asset_price']:.8f} "
            f"uc={row['housing_user_cost']:.8f} T={row['equal_transfer_period_units']:.8f}",
            flush=True,
        )

    base, reform = rows
    for field in (
        "asset_price",
        "housing_user_cost",
        "equal_transfer_period_units",
        "topcode_adjusted_births_per_adult",
        "owner_rate",
        "dependent_child_owner_rate",
    ):
        reform[f"{field}_change_vs_baseline"] = float(reform[field]) - float(base[field])
        reform[f"{field}_percent_change_vs_baseline"] = 100.0 * (
            float(reform[field]) / max(abs(float(base[field])), 1.0e-15) - 1.0
        )
    population_gap = abs(float(reform["adult_population"]) - float(base["adult_population"]))
    if population_gap > 1.0e-12:
        raise RuntimeError(f"Contemporaneous population lag gate failed: {population_gap:.3e}")

    reproduction_rows = reconstruction.pop("history_reproduction_rows")
    baseline.write_csv(outdir / "shared_2007_2019_history.csv", history)
    baseline.write_csv(outdir / "history_reproduction_audit.csv", reproduction_rows)
    baseline.write_csv(outdir / "impact_results.csv", rows)
    summary = {
        "status": "complete_exact_coven_property_tax_impact_smoke",
        "scope": (
            "fitted 2023 pre-fertility state; rebated 1% baseline versus rebated 2% "
            "impact; explicit static housing-supply elasticity; closed M=0,rho=1; "
            "no grant"
        ),
        "housing_supply_elasticity": float(args.housing_supply_elasticity),
        "housing_supply_elasticity_source": reanchor["supply_elasticity_source"],
        "cases": list(CASES),
        "supply_reanchor": reanchor,
        "population_lag": tax.renewal_lag_gate(),
        "impact_population_identity_gap": population_gap,
        "common_initial_state": {
            "distribution_sha256": baseline.array_sha256(state_2023.g_pre),
            "queue_sha256": baseline.canonical_json_sha256(state_2023.scheduled_entries),
        },
        "anchor_model_solve_count": int(anchor_counter.total),
        "case_gates": gates,
        "impact_results": rows,
        "reconstruction": reconstruction,
        "input_hashes": contracts["hashes"],
        "elapsed_seconds": time.perf_counter() - started,
    }
    baseline.write_json(outdir / "summary.json", summary)
    print(
        f"COVEN_TAX_SMOKE_COMPLETE cases={len(rows)} seconds={summary['elapsed_seconds']:.2f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
