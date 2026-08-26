#!/usr/bin/env python3
"""Perfect-foresight rebated property-tax comparison for the active E5F model.

The experiment keeps the live policy contract: a permanent rebated 1 percent
annual property tax is the funded-policy reference and a permanent rebated
2 percent tax is the reform.  Both cases start from the same reconstructed 2023
state, use the same declared exogenous fertility-preference path, close national
population renewal with M=0 and rho=1, and use the same re-anchored housing-
supply schedule.  Households perfectly foresee prices, rents, transfers, and
the preference path.  Housing markets and the government budget clear at every
date.  Each fiscal regime has its own jointly solved terminal steady state.

This remains an unpromoted diagnostic.  It does not estimate the preference
path and it makes no welfare claim.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Sequence

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

import audit_closed_reproductive_closure as closure
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_perfect_foresight_transition as pf
import run_e5f_post2023_coven_property_tax_smoke as impact
import run_e5f_post2023_no_policy_continuations as continuation
import run_e5f_post2023_rebated_property_tax_smoke as tax

DEFAULT_OUTPUT = (
    ROOT
    / "output/model/e5f_perfect_foresight_rebated_property_tax_diagnostic_20260825"
)
DEFAULT_HOUSING_SUPPLY_ELASTICITY = 1.75
DEFAULT_FISCAL_TOLERANCE = tax.FISCAL_ABSOLUTE_TOLERANCE
DEFAULT_TRANSFER_TERMINAL_RELATIVE_TOLERANCE = 0.01
DEFAULT_POLICY_RENEWAL_TOLERANCE = 1.0e-5


@dataclass
class PolicyTerminalSteadyState:
    case: tax.TaxCase
    psi_child: float
    asset_price: float
    renter_price: float
    equal_transfer: float
    population_scale: float
    entry_flow: float
    raw_queue_flow: float
    policy: calendar.PolicyBundle
    parameters: SimpleNamespace
    state: pf.PFInitialState
    endpoint: dict[str, Any]
    reconstruction: dict[str, Any]
    fixed_point_gates: dict[str, Any]


@dataclass
class FundedPFSolution:
    case: str
    status: str
    converged: bool
    iterations: int
    prices: np.ndarray
    rents: np.ndarray
    transfers: np.ndarray
    psi_path: np.ndarray
    rows: list[dict[str, Any]]
    iteration_history: list[dict[str, Any]]
    maximum_market_residual: float
    maximum_fiscal_residual: float
    maximum_log_price_residual: float
    maximum_transfer_residual: float
    maximum_policy_reproduction_error: float
    maximum_mass_accounting_error: float
    maximum_feasibility_projection_mass: float
    terminal_convergence: dict[str, Any]
    bellman_solves: int
    elapsed_seconds: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, default=pf.DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=pf.DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=pf.DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--horizons", type=int, nargs="+", default=[32, 64, 160])
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=sorted(impact.CASES),
        default=sorted(impact.CASES),
    )
    parser.add_argument("--psi-persistence", type=float, default=pf.DEFAULT_PSI_PERSISTENCE)
    parser.add_argument("--terminal-psi-child", type=float)
    parser.add_argument(
        "--housing-supply-elasticity",
        type=float,
        default=DEFAULT_HOUSING_SUPPLY_ELASTICITY,
    )
    parser.add_argument("--endpoint-price-min-ratio", type=float, default=1.0e-4)
    parser.add_argument("--endpoint-price-max-ratio", type=float, default=3.0)
    parser.add_argument("--endpoint-price-grid-points", type=int, default=25)
    parser.add_argument(
        "--endpoint-root-tol",
        type=float,
        default=DEFAULT_POLICY_RENEWAL_TOLERANCE,
    )
    parser.add_argument(
        "--endpoint-transfer-tol",
        type=float,
        default=tax.FISCAL_ABSOLUTE_TOLERANCE,
    )
    parser.add_argument("--market-tol", type=float, default=2.0e-4)
    parser.add_argument(
        "--fiscal-tol", type=float, default=DEFAULT_FISCAL_TOLERANCE
    )
    parser.add_argument("--max-path-iterations", type=int, default=35)
    parser.add_argument("--price-damping", type=float, default=0.25)
    parser.add_argument("--transfer-damping", type=float, default=0.50)
    parser.add_argument("--maximum-log-price-step", type=float, default=0.12)
    parser.add_argument("--maximum-transfer-step", type=float, default=0.08)
    parser.add_argument(
        "--initial-path",
        action="append",
        default=[],
        metavar="CASE=CSV",
        help=(
            "Optional solved path used only as a numerical seed. Repeat once per "
            "case. The CSV needs period, asset_price, and "
            "equal_transfer_period_units columns."
        ),
    )
    return parser.parse_args()


def policy_overrides(
    prepared: continuation.PreparedModel,
    case: tax.TaxCase,
    transfer: float,
) -> dict[str, Any]:
    overrides = copy.deepcopy(prepared.base_overrides)
    period_tax = float(case.annual_tax_rate) * float(prepared.parameters.period_years)
    overrides.update(
        tau_H=period_tax,
        user_cost_rate=(
            float(prepared.parameters.q)
            + float(prepared.parameters.delta)
            + period_tax
        ),
        property_tax_lump_sum_transfer=float(transfer),
        birth_entry_grant=False,
        birth_entry_grant_amount=0.0,
        birth_entry_grant_locations=np.array([], dtype=int),
        birth_entry_grant_owner_rungs=np.array([], dtype=int),
        propagate_birth_entry_grant=False,
    )
    return overrides


def _stationary_evaluation(
    *,
    prepared: continuation.PreparedModel,
    case: tax.TaxCase,
    asset_price: float,
    transfer: float,
    psi_child: float,
) -> SimpleNamespace:
    overrides = policy_overrides(prepared, case, transfer)
    solution, parameters, price_array, elapsed = closure.solve_pe(
        prepared.chain,
        overrides,
        asset_price=float(asset_price),
        psi_child=float(psi_child),
    )
    solved_price = float(np.asarray(price_array, dtype=float).reshape(-1)[0])
    if not math.isclose(solved_price, float(asset_price), rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("Stationary policy solve did not reproduce the fixed price.")
    expected_tax = float(case.annual_tax_rate) * float(parameters.period_years)
    if not math.isclose(float(parameters.tau_H), expected_tax, rel_tol=0.0, abs_tol=1e-13):
        raise RuntimeError("Stationary policy solve changed the property-tax rate.")
    if not math.isclose(
        float(parameters.property_tax_lump_sum_transfer),
        float(transfer),
        rel_tol=0.0,
        abs_tol=1e-13,
    ):
        raise RuntimeError("Stationary policy solve changed the equal transfer.")
    b_grid = np.asarray(solution.b_grid, dtype=float)
    if b_grid.shape != prepared.b_grid.shape or not np.allclose(
        b_grid, prepared.b_grid, rtol=0.0, atol=0.0
    ):
        raise RuntimeError("Stationary policy solve changed the asset grid.")
    parameters._fert2_probs = np.asarray(solution.fert2_probs, dtype=float).copy()
    shared = calendar.model.precompute_shared(parameters, b_grid)
    policy = calendar.policy_from_solution(
        solution, price_array, parameters, b_grid, shared
    )
    unit_g_pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        solution, policy, parameters, b_grid, shared
    )
    unit_mass = float(np.sum(unit_g_pre))
    if not math.isclose(unit_mass, 1.0, rel_tol=0.0, abs_tol=2e-10):
        raise RuntimeError(f"Stationary distribution is not normalized: {unit_mass}")
    evaluation = calendar.evaluate_period(
        price_array,
        unit_g_pre,
        parameters,
        b_grid,
        shared,
        calendar.SolveCounter(),
        supply_rule=None,
        supplied_policy=policy,
    )
    ledger = tax.fiscal_ledger(evaluation, parameters)
    accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        parameters,
    )
    current_mass = float(np.sum(evaluation.g_current))
    entry = float(np.sum(unit_g_pre[:, :, :, 0, :, :, :])) / unit_mass
    adjusted_births = float(accounting["topcode_adjusted_birth_children"]) / current_mass
    mature = transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION * adjusted_births
    demand_per_adult = float(evaluation.demand_by_loc[0]) / current_mass
    return SimpleNamespace(
        solution=solution,
        parameters=parameters,
        price_array=price_array,
        b_grid=b_grid,
        shared=shared,
        policy=policy,
        unit_g_pre=unit_g_pre,
        reconstruction=reconstruction,
        evaluation=evaluation,
        ledger=ledger,
        renewal_ratio=mature / max(entry, 1e-15),
        entry_per_adult=entry,
        adjusted_births_per_adult=adjusted_births,
        mature_per_adult=mature,
        demand_per_adult=demand_per_adult,
        solve_seconds=float(elapsed),
    )


def solve_stationary_equal_transfer(
    *,
    prepared: continuation.PreparedModel,
    case: tax.TaxCase,
    asset_price: float,
    psi_child: float,
    tolerance: float,
) -> tuple[SimpleNamespace, int]:
    evaluations = 0

    def evaluate(transfer: float) -> SimpleNamespace:
        nonlocal evaluations
        evaluations += 1
        return _stationary_evaluation(
            prepared=prepared,
            case=case,
            asset_price=float(asset_price),
            transfer=float(transfer),
            psi_child=float(psi_child),
        )

    left_transfer = 0.0
    left = evaluate(left_transfer)
    left_gap = float(left.ledger["implied_equal_transfer"]) - left_transfer
    if left_gap < 0.0:
        raise RuntimeError("The equal-transfer fixed point has the wrong sign at zero.")
    right_transfer = max(1.25 * float(left.ledger["implied_equal_transfer"]), 0.05)
    del left
    right = evaluate(right_transfer)
    right_gap = float(right.ledger["implied_equal_transfer"]) - right_transfer
    while right_gap > 0.0:
        right_transfer *= 2.0
        if right_transfer > 8.0:
            raise RuntimeError("Could not bracket the stationary fiscal fixed point.")
        right = evaluate(right_transfer)
        right_gap = float(right.ledger["implied_equal_transfer"]) - right_transfer

    result = right
    del right
    for _ in range(32):
        result_gap = (
            float(result.ledger["implied_equal_transfer"])
            - float(result.parameters.property_tax_lump_sum_transfer)
        )
        if abs(result_gap) <= float(tolerance):
            return result, evaluations
        if abs(right_gap - left_gap) > 1e-14:
            trial = right_transfer - right_gap * (
                right_transfer - left_transfer
            ) / (right_gap - left_gap)
        else:
            trial = 0.5 * (left_transfer + right_transfer)
        if not left_transfer < trial < right_transfer:
            trial = 0.5 * (left_transfer + right_transfer)
        result = evaluate(trial)
        result_gap = float(result.ledger["implied_equal_transfer"]) - trial
        if result_gap > 0.0:
            left_transfer, left_gap = trial, result_gap
        else:
            right_transfer, right_gap = trial, result_gap
    raise RuntimeError(
        "Stationary fiscal fixed point failed: "
        f"case={case.name}, price={asset_price}, evaluations={evaluations}, "
        f"bracket=({left_transfer},{right_transfer}), "
        f"gaps=({left_gap},{right_gap}), last_gap={result_gap}"
    )


def _endpoint_row(
    result: SimpleNamespace,
    *,
    case: tax.TaxCase,
    supply_rule: calendar.HousingSupplyRule,
    model_evaluations: int,
) -> dict[str, Any]:
    price = float(result.price_array[0])
    supply = float(supply_rule.quantity(np.array([price]))[0])
    scale = supply / max(float(result.demand_per_adult), 1e-15)
    transfer = float(result.parameters.property_tax_lump_sum_transfer)
    return {
        "policy_case": case.name,
        "annual_property_tax_rate": float(case.annual_tax_rate),
        "asset_price": price,
        "renter_price": float(result.parameters.user_cost_rate) * price,
        "equal_transfer_period_units": transfer,
        "implied_equal_transfer": float(result.ledger["implied_equal_transfer"]),
        "government_budget_residual_per_unit_mass": float(
            result.ledger["government_budget_residual"]
        ),
        "scaled_government_budget_residual": float(
            result.ledger["scaled_government_budget_residual"]
        ),
        "stationary_entry_flow_E": float(result.entry_per_adult),
        "topcode_adjusted_birth_children_per_adult": float(
            result.adjusted_births_per_adult
        ),
        "owner_rate": pf._owner_rate(result.evaluation.g_current),
        "queue_mature_flow_B": float(result.mature_per_adult),
        "queue_B_over_E": float(result.renewal_ratio),
        "renewal_residual_ratio": float(result.renewal_ratio) - 1.0,
        "housing_demand_per_adult": float(result.demand_per_adult),
        "stationary_housing_supply": supply,
        "stationary_population_scale": scale,
        "scaled_housing_demand": scale * float(result.demand_per_adult),
        "housing_clearing_residual": scale * float(result.demand_per_adult) - supply,
        "stationary_fiscal_model_evaluations": int(model_evaluations),
        "stationary_policy_solve_seconds": float(result.solve_seconds),
    }


def solve_policy_terminal_steady_state(
    prepared: continuation.PreparedModel,
    *,
    case: tax.TaxCase,
    terminal_psi_child: float,
    supply_rule: calendar.HousingSupplyRule,
    price_min_ratio: float,
    price_max_ratio: float,
    price_grid_points: int,
    root_tolerance: float,
    transfer_tolerance: float,
    output_dir: Path,
) -> PolicyTerminalSteadyState:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows_by_price: dict[float, dict[str, Any]] = {}

    def evaluate(
        asset_price: float, *, enforce_scaled_fiscal_gate: bool = False
    ) -> tuple[dict[str, Any], SimpleNamespace]:
        key = round(float(asset_price), 14)
        result, count = solve_stationary_equal_transfer(
            prepared=prepared,
            case=case,
            asset_price=float(asset_price),
            psi_child=float(terminal_psi_child),
            tolerance=float(transfer_tolerance),
        )
        row = _endpoint_row(
            result,
            case=case,
            supply_rule=supply_rule,
            model_evaluations=count,
        )
        if enforce_scaled_fiscal_gate:
            # The equal transfer is per adult, while the live fiscal gate is an
            # absolute aggregate residual.  Tighten only candidate-root solves
            # enough to preserve that gate after the stationary distribution is
            # scaled to clear the housing market.  Broad schedule points retain
            # the declared diagnostic tolerance, which is robust to policy kinks.
            for _ in range(3):
                projected_absolute_residual = abs(
                    float(row["government_budget_residual_per_unit_mass"])
                ) * float(row["stationary_population_scale"])
                if projected_absolute_residual <= tax.FISCAL_ABSOLUTE_TOLERANCE:
                    break
                tighter_tolerance = min(
                    float(transfer_tolerance),
                    0.8
                    * tax.FISCAL_ABSOLUTE_TOLERANCE
                    / max(float(row["stationary_population_scale"]), 1.0),
                )
                result, added_count = solve_stationary_equal_transfer(
                    prepared=prepared,
                    case=case,
                    asset_price=float(asset_price),
                    psi_child=float(terminal_psi_child),
                    tolerance=tighter_tolerance,
                )
                count += added_count
                row = _endpoint_row(
                    result,
                    case=case,
                    supply_rule=supply_rule,
                    model_evaluations=count,
                )
            projected_absolute_residual = abs(
                float(row["government_budget_residual_per_unit_mass"])
            ) * float(row["stationary_population_scale"])
            if projected_absolute_residual > tax.FISCAL_ABSOLUTE_TOLERANCE:
                raise RuntimeError(
                    "Candidate terminal fiscal root cannot satisfy the scaled "
                    f"absolute budget gate: {row}"
                )
        row["projected_government_budget_absolute_residual"] = abs(
            float(row["government_budget_residual_per_unit_mass"])
        ) * float(row["stationary_population_scale"])
        row["scaled_fiscal_gate_enforced"] = bool(enforce_scaled_fiscal_gate)
        rows_by_price[key] = row
        rows = sorted(rows_by_price.values(), key=lambda x: x["asset_price"])
        pf.write_csv(output_dir / "closed_stationary_schedule_progress.csv", rows)
        pf.write_json(
            output_dir / "latest_endpoint_search.json",
            {"status": "running", "evaluated_prices": len(rows), "latest": row},
        )
        return row, result

    prices = float(prepared.old_price) * np.geomspace(
        float(price_min_ratio), float(price_max_ratio), int(price_grid_points)
    )
    schedule = []
    for price in prices:
        row, result = evaluate(float(price))
        schedule.append(row)
        del result
    brackets = continuation.closed_schedule_brackets(schedule)
    if len(brackets) != 1:
        pf.write_csv(output_dir / "closed_stationary_schedule.csv", schedule)
        raise RuntimeError(
            f"{case.name} needs one positive renewal root; found {len(brackets)}: {brackets}"
        )
    lower, upper = brackets[0]
    if lower == upper:
        root_row, root_result = evaluate(lower, enforce_scaled_fiscal_gate=True)
        final_bracket_rows = (root_row, root_row)
    else:
        low_row, low_result = evaluate(lower, enforce_scaled_fiscal_gate=True)
        del low_result
        high_row, high_result = evaluate(upper, enforce_scaled_fiscal_gate=True)
        del high_result
        low_gap = float(low_row["renewal_residual_ratio"])
        high_gap = float(high_row["renewal_residual_ratio"])
        if low_gap * high_gap > 0.0:
            raise RuntimeError(
                "The funded-policy renewal bracket did not survive the final "
                f"scaled fiscal gate: lower={low_row}, upper={high_row}"
            )
        root_row: dict[str, Any]
        root_result: SimpleNamespace
        for _ in range(50):
            midpoint = 0.5 * (lower + upper)
            root_row, root_result = evaluate(
                midpoint, enforce_scaled_fiscal_gate=True
            )
            gap = float(root_row["renewal_residual_ratio"])
            if abs(gap) <= float(root_tolerance):
                break
            if low_gap * gap <= 0.0:
                upper = midpoint
                high_row = root_row
            else:
                lower, low_gap = midpoint, gap
                low_row = root_row
        final_bracket_rows = (low_row, high_row)

    endpoint = dict(root_row)
    endpoint.update(
        status="complete_usable_positive_funded_policy_root",
        usable_closed_root=True,
        audited_root_bracket=brackets[0],
        audited_final_root_bracket=[
            float(final_bracket_rows[0]["asset_price"]),
            float(final_bracket_rows[1]["asset_price"]),
        ],
        audited_final_root_bracket_residuals=[
            float(final_bracket_rows[0]["renewal_residual_ratio"]),
            float(final_bracket_rows[1]["renewal_residual_ratio"]),
        ],
        renewal_root_absolute_residual=abs(float(root_row["renewal_residual_ratio"])),
        renewal_root_declared_tolerance=float(root_tolerance),
        renewal_root_interpretation=(
            "Closest funded-policy renewal point within the declared tolerance; "
            "the signed final bracket is retained to expose any discrete-policy kink."
        ),
        fiscal_root_absolute_transfer_gap=abs(
            float(root_row["equal_transfer_period_units"])
            - float(root_row["implied_equal_transfer"])
        ),
        fiscal_root_declared_tolerance=float(transfer_tolerance),
    )
    if endpoint["renewal_root_absolute_residual"] > float(root_tolerance):
        raise RuntimeError(f"Terminal renewal root failed: {endpoint}")
    if endpoint["fiscal_root_absolute_transfer_gap"] > float(transfer_tolerance):
        raise RuntimeError(f"Terminal fiscal root failed: {endpoint}")

    scale = float(endpoint["stationary_population_scale"])
    unit_g_pre = np.asarray(root_result.unit_g_pre, dtype=float)
    scaled_g_pre = scale * unit_g_pre
    parameters = root_result.parameters
    policy = root_result.policy
    entry_flow = float(np.sum(scaled_g_pre[:, :, :, 0, :, :, :]))
    evaluation = calendar.evaluate_period(
        root_result.price_array,
        scaled_g_pre,
        parameters,
        root_result.b_grid,
        root_result.shared,
        calendar.SolveCounter(),
        supply_rule=supply_rule,
        supplied_policy=policy,
    )
    accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        parameters,
    )
    conversion = transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
    implied_queue_flow = conversion * float(
        accounting["topcode_adjusted_birth_children"]
    )
    raw_queue_flow = conversion * float(evaluation.births)
    state = pf.PFInitialState(
        g_pre=scaled_g_pre,
        scheduled_entries=[entry_flow] * continuation.QUEUE_WAITING_SLOTS,
        scheduled_raw_entries=[raw_queue_flow] * continuation.QUEUE_WAITING_SLOTS,
    )
    transfer = float(parameters.property_tax_lump_sum_transfer)
    price = float(root_result.price_array[0])
    one_step = pf.evaluate_path_at_prices(
        prices=np.array([price]),
        psi_path=np.array([float(terminal_psi_child)]),
        terminal_price=price,
        terminal_V=policy.V,
        base_parameters=parameters,
        b_grid=root_result.b_grid,
        initial_state=state,
        supply_rule=supply_rule,
        birth_to_entry_conversion=conversion,
        transfer_path=np.array([transfer]),
    )
    next_mass = float(np.sum(one_step.terminal_state.g_pre))
    normalized_l1 = float(
        np.sum(
            np.abs(
                one_step.terminal_state.g_pre / max(next_mass, 1e-15)
                - scaled_g_pre / max(scale, 1e-15)
            )
        )
    )
    ledger = tax.fiscal_ledger(evaluation, parameters)
    gates = {
        "population_scale": scale,
        "stationary_entry_flow_from_distribution": entry_flow,
        "implied_birth_queue_flow": implied_queue_flow,
        "birth_queue_relative_gap": abs(implied_queue_flow - entry_flow)
        / max(entry_flow, 1e-15),
        "housing_market_relative_residual": float(evaluation.relative_market_residual),
        "government_budget_absolute_residual": abs(
            float(ledger["government_budget_residual"])
        ),
        "equal_transfer_absolute_gap": abs(
            transfer - float(ledger["implied_equal_transfer"])
        ),
        "one_step_population_absolute_gap": abs(next_mass - scale),
        "one_step_normalized_distribution_l1": normalized_l1,
        "one_step_policy_reproduction_max_abs": one_step.maximum_policy_reproduction_error,
        "one_step_mass_accounting_max_abs": one_step.maximum_mass_accounting_error,
        "renewal_root_absolute_residual": endpoint["renewal_root_absolute_residual"],
    }
    checks = {
        "birth_queue_relative_gap": gates["birth_queue_relative_gap"]
        <= max(5.0 * float(root_tolerance), 5e-8),
        "housing_market_relative_residual": abs(gates["housing_market_relative_residual"])
        <= 2.5e-5,
        "government_budget_absolute_residual": gates[
            "government_budget_absolute_residual"
        ]
        <= tax.FISCAL_ABSOLUTE_TOLERANCE,
        "equal_transfer_absolute_gap": gates["equal_transfer_absolute_gap"]
        <= tax.FISCAL_ABSOLUTE_TOLERANCE,
        "one_step_population_absolute_gap": gates["one_step_population_absolute_gap"]
        <= 2e-8,
        "one_step_normalized_distribution_l1": gates[
            "one_step_normalized_distribution_l1"
        ]
        <= 2e-8,
        "one_step_policy_reproduction_max_abs": gates[
            "one_step_policy_reproduction_max_abs"
        ]
        <= 2e-8,
        "one_step_mass_accounting_max_abs": gates[
            "one_step_mass_accounting_max_abs"
        ]
        <= 2e-8,
    }
    gates.update(status="passed" if all(checks.values()) else "failed", checks=checks)
    if gates["status"] != "passed":
        raise RuntimeError(f"Funded-policy terminal fixed-point gates failed: {gates}")

    schedule = sorted(rows_by_price.values(), key=lambda x: x["asset_price"])
    pf.write_csv(output_dir / "closed_stationary_schedule.csv", schedule)
    pf.write_json(output_dir / "closed_stationary_endpoint.json", endpoint)
    pf.write_json(
        output_dir / "latest_endpoint_search.json",
        {"status": "complete", "closed_stationary_endpoint": endpoint},
    )
    return PolicyTerminalSteadyState(
        case=case,
        psi_child=float(terminal_psi_child),
        asset_price=price,
        renter_price=float(parameters.user_cost_rate) * price,
        equal_transfer=transfer,
        population_scale=scale,
        entry_flow=entry_flow,
        raw_queue_flow=raw_queue_flow,
        policy=policy,
        parameters=parameters,
        state=state,
        endpoint=endpoint,
        reconstruction=root_result.reconstruction,
        fixed_point_gates=gates,
    )


def solve_funded_perfect_foresight(
    *,
    case: tax.TaxCase,
    initial_prices: Sequence[float],
    initial_transfers: Sequence[float],
    psi_path: Sequence[float],
    terminal: PolicyTerminalSteadyState,
    b_grid: np.ndarray,
    initial_state: pf.PFInitialState,
    supply_rule: calendar.HousingSupplyRule,
    market_tolerance: float,
    fiscal_tolerance: float,
    max_iterations: int,
    price_damping: float,
    transfer_damping: float,
    maximum_log_price_step: float,
    maximum_transfer_step: float,
    progress_callback: Callable[[dict[str, Any], pf.PathEvaluation], None] | None = None,
) -> FundedPFSolution:
    started = time.perf_counter()
    prices = np.asarray(initial_prices, dtype=float).reshape(-1).copy()
    transfers = np.asarray(initial_transfers, dtype=float).reshape(-1).copy()
    psi_values = np.asarray(psi_path, dtype=float).reshape(-1)
    if not (prices.shape == transfers.shape == psi_values.shape):
        raise ValueError("Price, transfer, and preference paths must have equal length.")
    if np.any(prices <= 0.0) or np.any(transfers < 0.0):
        raise ValueError("Initial prices must be positive and transfers nonnegative.")
    current_price_damping = float(price_damping)
    current_transfer_damping = float(transfer_damping)
    history: list[dict[str, Any]] = []
    best: tuple[
        float,
        float,
        float,
        np.ndarray,
        np.ndarray,
        int,
    ] | None = None
    total_bellman = 0
    previous_score = math.inf
    converged = False
    for iteration in range(1, int(max_iterations) + 1):
        evaluation = pf.evaluate_path_at_prices(
            prices=prices,
            psi_path=psi_values,
            terminal_price=terminal.asset_price,
            terminal_V=terminal.policy.V,
            base_parameters=terminal.parameters,
            b_grid=b_grid,
            initial_state=initial_state,
            supply_rule=supply_rule,
            birth_to_entry_conversion=transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION,
            transfer_path=transfers,
        )
        total_bellman += evaluation.bellman_solves
        demand = np.array([float(row["housing_demand"]) for row in evaluation.rows])
        target_prices = pf.inverse_supply_price(demand, supply_rule)
        target_transfers = np.array(
            [float(row["implied_equal_transfer"]) for row in evaluation.rows]
        )
        log_residual = np.log(target_prices) - np.log(prices)
        transfer_residual = target_transfers - transfers
        market_error = float(evaluation.maximum_market_residual)
        fiscal_error = max(
            abs(float(row["government_budget_residual"]))
            for row in evaluation.rows
        )
        max_log = float(np.max(np.abs(log_residual)))
        max_transfer = float(np.max(np.abs(transfer_residual)))
        score = max(
            market_error / float(market_tolerance),
            fiscal_error / float(fiscal_tolerance),
        )
        if best is None or (score, market_error, fiscal_error) < best[:3]:
            best = (
                score,
                market_error,
                fiscal_error,
                evaluation.prices.copy(),
                transfers.copy(),
                iteration,
            )
        record = {
            "iteration": iteration,
            "maximum_market_residual": market_error,
            "maximum_fiscal_residual": fiscal_error,
            "maximum_log_price_residual": max_log,
            "maximum_transfer_residual": max_transfer,
            "price_damping": current_price_damping,
            "transfer_damping": current_transfer_damping,
            "minimum_asset_price": float(np.min(prices)),
            "maximum_asset_price": float(np.max(prices)),
            "minimum_renter_price": float(np.min(evaluation.rents)),
            "maximum_renter_price": float(np.max(evaluation.rents)),
            "minimum_equal_transfer": float(np.min(transfers)),
            "maximum_equal_transfer": float(np.max(transfers)),
            "bellman_solves_this_iteration": evaluation.bellman_solves,
            "bellman_solves_cumulative": total_bellman,
            "elapsed_seconds": time.perf_counter() - started,
        }
        history.append(record)
        if market_error <= float(market_tolerance) and fiscal_error <= float(
            fiscal_tolerance
        ):
            converged = True
            best = (
                score,
                market_error,
                fiscal_error,
                evaluation.prices.copy(),
                transfers.copy(),
                iteration,
            )
            if progress_callback is not None:
                progress_callback(record, evaluation)
            break
        if score > previous_score * 1.02:
            current_price_damping = max(0.04, 0.5 * current_price_damping)
            current_transfer_damping = max(0.08, 0.5 * current_transfer_damping)
        elif score < previous_score * 0.80:
            current_price_damping = min(float(price_damping), 1.08 * current_price_damping)
            current_transfer_damping = min(
                float(transfer_damping), 1.08 * current_transfer_damping
            )
        price_step = np.clip(
            log_residual,
            -float(maximum_log_price_step),
            float(maximum_log_price_step),
        )
        transfer_step = np.clip(
            transfer_residual,
            -float(maximum_transfer_step),
            float(maximum_transfer_step),
        )
        candidate_prices = np.exp(
            np.log(prices) + current_price_damping * price_step
        )
        candidate_transfers = np.maximum(
            transfers + current_transfer_damping * transfer_step, 0.0
        )
        prices, feasibility = project_price_path_to_positive_rents(
            candidate_prices,
            terminal=terminal,
            minimum_rent_share=1.0e-6,
        )
        transfers = candidate_transfers
        record.update(
            next_price_feasibility_adjusted_period_count=feasibility[
                "adjusted_period_count"
            ],
            next_price_feasibility_maximum_relative_change=feasibility[
                "maximum_price_relative_change"
            ],
            next_price_minimum_implied_renter_price=feasibility[
                "minimum_implied_renter_price"
            ],
        )
        if progress_callback is not None:
            progress_callback(record, evaluation)
        previous_score = score
        if iteration < int(max_iterations):
            # Release the full value-function stack before constructing the
            # next path evaluation.  Assignment evaluates its right-hand side
            # before replacing the old object, which otherwise doubles peak
            # memory at long horizons.
            del evaluation

    if best is None:
        raise RuntimeError("Joint policy-path iteration produced no evaluation.")
    (
        _,
        best_market,
        best_fiscal,
        best_prices,
        best_transfers,
        best_iteration,
    ) = best
    if best_iteration == len(history):
        best_evaluation = evaluation
    else:
        best_evaluation = pf.evaluate_path_at_prices(
            prices=best_prices,
            psi_path=psi_values,
            terminal_price=terminal.asset_price,
            terminal_V=terminal.policy.V,
            base_parameters=terminal.parameters,
            b_grid=b_grid,
            initial_state=initial_state,
            supply_rule=supply_rule,
            birth_to_entry_conversion=transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION,
            transfer_path=best_transfers,
        )
        total_bellman += best_evaluation.bellman_solves
    diagnostics = pf.terminal_convergence_diagnostics(
        evaluation=best_evaluation,
        psi_path=psi_values,
        reference_state=terminal.state,
        reference_entry_flow=terminal.entry_flow,
        reference_price=terminal.asset_price,
        reference_psi=terminal.psi_child,
        base_parameters=terminal.parameters,
    )
    last_transfer = float(best_transfers[-1])
    transfer_relative_gap = abs(last_transfer - terminal.equal_transfer) / max(
        abs(terminal.equal_transfer), 1e-15
    )
    diagnostics["terminal_equal_transfer_last_endogenous_date"] = last_transfer
    diagnostics["stationary_equal_transfer"] = terminal.equal_transfer
    diagnostics["metrics"]["equal_transfer_relative_gap"] = transfer_relative_gap
    diagnostics["tolerances"]["equal_transfer_relative_gap"] = (
        DEFAULT_TRANSFER_TERMINAL_RELATIVE_TOLERANCE
    )
    diagnostics["checks"]["equal_transfer_relative_gap"] = (
        transfer_relative_gap <= DEFAULT_TRANSFER_TERMINAL_RELATIVE_TOLERANCE
    )
    diagnostics["all_checks_pass"] = all(diagnostics["checks"].values())
    diagnostics["status"] = (
        "passed" if diagnostics["all_checks_pass"] else "not_converged"
    )
    for row in best_evaluation.rows:
        row.update(
            policy_case=case.name,
            policy_label=case.label,
            annual_property_tax_rate=float(case.annual_tax_rate),
            property_tax_rate_period_units=float(terminal.parameters.tau_H),
            equal_rebate=True,
            housing_supply_elasticity=float(supply_rule.elasticity),
        )
    return FundedPFSolution(
        case=case.name,
        status="converged" if converged else "maximum_iterations_reached",
        converged=converged,
        iterations=len(history),
        prices=best_evaluation.prices,
        rents=best_evaluation.rents,
        transfers=best_transfers,
        psi_path=psi_values,
        rows=best_evaluation.rows,
        iteration_history=history,
        maximum_market_residual=best_market,
        maximum_fiscal_residual=best_fiscal,
        maximum_log_price_residual=float(
            max(abs(math.log(float(row["housing_demand"]) / float(row["housing_supply"]))) / float(supply_rule.elasticity) for row in best_evaluation.rows)
        ),
        maximum_transfer_residual=max(
            abs(float(row["implied_equal_transfer"]) - float(row["equal_transfer_period_units"]))
            for row in best_evaluation.rows
        ),
        maximum_policy_reproduction_error=best_evaluation.maximum_policy_reproduction_error,
        maximum_mass_accounting_error=best_evaluation.maximum_mass_accounting_error,
        maximum_feasibility_projection_mass=best_evaluation.maximum_feasibility_projection_mass,
        terminal_convergence=diagnostics,
        bellman_solves=total_bellman,
        elapsed_seconds=time.perf_counter() - started,
    )


def parse_initial_paths(specifications: Sequence[str]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for specification in specifications:
        if "=" not in specification:
            raise ValueError(f"Initial path must have CASE=CSV form: {specification}")
        case, raw_path = specification.split("=", 1)
        if case not in impact.CASES:
            raise ValueError(f"Unknown initial-path case: {case}")
        if case in paths:
            raise ValueError(f"Duplicate initial path for {case}")
        paths[case] = Path(raw_path).resolve()
    return paths


def project_price_path_to_positive_rents(
    prices: Sequence[float],
    *,
    terminal: PolicyTerminalSteadyState,
    minimum_rent_share: float = 0.01,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Project a numerical price path onto the positive-rent domain."""
    projected = np.asarray(prices, dtype=float).reshape(-1).copy()
    if projected.size < 1 or np.any(~np.isfinite(projected)) or np.any(projected <= 0.0):
        raise ValueError("The initial asset-price seed must be finite and positive.")
    margin = float(minimum_rent_share)
    if not 0.0 < margin < 1.0:
        raise ValueError("The seed minimum-rent share must lie in (0,1).")
    carrying_factor = (
        float(terminal.parameters.R_gross)
        + float(terminal.parameters.delta)
        + float(terminal.parameters.tau_H)
    )
    if not math.isfinite(carrying_factor) or carrying_factor <= 0.0:
        raise ValueError("The owner carrying factor must be finite and positive.")

    original = projected.copy()
    next_price = float(terminal.asset_price)
    adjusted_periods: list[int] = []
    for period in range(projected.size - 1, -1, -1):
        minimum_price = (1.0 + margin) * next_price / carrying_factor
        if projected[period] < minimum_price:
            projected[period] = minimum_price
            adjusted_periods.append(period)
        next_price = float(projected[period])
    rents = pf.rents_from_asset_prices(
        projected, terminal.asset_price, terminal.parameters
    )
    relative_changes = np.abs(projected / original - 1.0)
    return projected, {
        "minimum_rent_share_of_next_asset_price": margin,
        "adjusted_period_count": len(adjusted_periods),
        "adjusted_periods": sorted(adjusted_periods),
        "maximum_price_relative_change": float(np.max(relative_changes)),
        "minimum_implied_renter_price": float(np.min(rents)),
        "equilibrium_contract_changed": False,
    }


def project_price_seed_to_positive_rents(
    prices: Sequence[float],
    *,
    terminal: PolicyTerminalSteadyState,
    minimum_rent_share: float = 0.01,
) -> tuple[np.ndarray, dict[str, Any]]:
    projected, contract = project_price_path_to_positive_rents(
        prices,
        terminal=terminal,
        minimum_rent_share=minimum_rent_share,
    )
    contract.update(
        purpose="numerical_seed_feasibility_only",
        maximum_seed_price_relative_change=contract[
            "maximum_price_relative_change"
        ],
    )
    return projected, contract


def load_initial_path(
    path: Path,
    *,
    horizon: int,
    terminal: PolicyTerminalSteadyState,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"period", "asset_price", "equal_transfer_period_units"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Initial policy path must contain {sorted(required)}: {path}")
    if len(rows) > int(horizon):
        raise ValueError("Initial policy path is longer than the requested horizon.")
    periods = [int(row["period"]) for row in rows]
    if periods != list(range(len(rows))):
        raise ValueError("Initial policy-path periods must be contiguous from zero.")
    prices = np.full(int(horizon), terminal.asset_price, dtype=float)
    transfers = np.full(int(horizon), terminal.equal_transfer, dtype=float)
    prices[: len(rows)] = [float(row["asset_price"]) for row in rows]
    transfers[: len(rows)] = [
        float(row["equal_transfer_period_units"]) for row in rows
    ]
    prices, feasibility = project_price_seed_to_positive_rents(
        prices, terminal=terminal
    )
    return prices, transfers, {
        "mode": "recorded_policy_path_then_terminal_extension",
        "source": str(path),
        "source_sha256": pf.file_sha256(path),
        "source_rows": len(rows),
        "requested_horizon": int(horizon),
        "terminal_extension_periods": int(horizon) - len(rows),
        "price_seed_feasibility_projection": feasibility,
    }


def analytic_initial_paths(
    *,
    impact_price: float,
    impact_transfer: float,
    terminal: PolicyTerminalSteadyState,
    persistence: float,
    horizon: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    prices = pf.initial_price_path(
        impact_price, terminal.asset_price, persistence, horizon
    )
    weights = float(persistence) ** np.arange(int(horizon), dtype=float)
    transfers = terminal.equal_transfer + weights * (
        float(impact_transfer) - terminal.equal_transfer
    )
    prices, feasibility = project_price_seed_to_positive_rents(
        prices, terminal=terminal
    )
    return prices, transfers, {
        "mode": "analytic_impact_to_terminal_decay",
        "impact_asset_price": float(impact_price),
        "impact_equal_transfer": float(impact_transfer),
        "terminal_asset_price": terminal.asset_price,
        "terminal_equal_transfer": terminal.equal_transfer,
        "persistence": float(persistence),
        "requested_horizon": int(horizon),
        "price_seed_feasibility_projection": feasibility,
    }


def effect_rows(
    baseline_rows: Sequence[dict[str, Any]], reform_rows: Sequence[dict[str, Any]]
) -> list[dict[str, Any]]:
    if len(baseline_rows) != len(reform_rows):
        raise ValueError("Funded-policy paths have different lengths.")
    rows: list[dict[str, Any]] = []
    for baseline, reform in zip(baseline_rows, reform_rows, strict=True):
        base_births = float(baseline["birth_children_topcode_adjusted"]) / float(
            baseline["adult_population"]
        )
        reform_births = float(reform["birth_children_topcode_adjusted"]) / float(
            reform["adult_population"]
        )
        rows.append(
            {
                "period": int(baseline["period"]),
                "calendar_year": int(baseline["calendar_year"]),
                "asset_price_percent_change": 100.0
                * (float(reform["asset_price"]) / float(baseline["asset_price"]) - 1.0),
                "renter_price_percent_change": 100.0
                * (float(reform["renter_price"]) / float(baseline["renter_price"]) - 1.0),
                "equal_transfer_percent_change": 100.0
                * (
                    float(reform["equal_transfer_period_units"])
                    / float(baseline["equal_transfer_period_units"])
                    - 1.0
                ),
                "owner_rate_pp_change": 100.0
                * (float(reform["owner_rate"]) - float(baseline["owner_rate"])),
                "births_per_adult_percent_change": 100.0
                * (reform_births / base_births - 1.0),
                "adult_population_percent_change": 100.0
                * (
                    float(reform["adult_population"])
                    / float(baseline["adult_population"])
                    - 1.0
                ),
            }
        )
    return rows


def make_comparison_figure(
    solutions: dict[str, FundedPFSolution], output_dir: Path
) -> None:
    fields = (
        ("asset_price", "House price"),
        ("renter_price", "Rent"),
        ("equal_transfer_period_units", "Equal transfer"),
        ("owner_rate", "Ownership rate"),
        ("birth_children_topcode_adjusted", "Adjusted births"),
        ("adult_population", "Adult population"),
    )
    figure, axes = plt.subplots(3, 2, figsize=(10.0, 9.0), constrained_layout=True)
    for axis, (field, title) in zip(axes.flat, fields, strict=True):
        for case_name, solution in solutions.items():
            axis.plot(
                [int(row["calendar_year"]) for row in solution.rows],
                [float(row[field]) for row in solution.rows],
                label=impact.CASES[case_name].label,
            )
        axis.set_title(title)
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    axes.flat[0].legend(frameon=False)
    figure.savefig(output_dir / "perfect_foresight_rebated_property_tax.png", dpi=200)
    figure.savefig(output_dir / "perfect_foresight_rebated_property_tax.pdf")
    plt.close(figure)


def make_case_figure(
    solution: FundedPFSolution,
    terminal: PolicyTerminalSteadyState,
    output_dir: Path,
) -> None:
    years = np.array([int(row["calendar_year"]) for row in solution.rows])
    figure, axes = plt.subplots(3, 2, figsize=(10.2, 9.0), constrained_layout=True)
    axes[0, 0].plot(years, solution.prices, label="House price")
    axes[0, 0].plot(years, solution.rents, label="Rent")
    axes[0, 0].axhline(terminal.asset_price, color="0.4", linestyle="--", linewidth=1)
    axes[0, 0].set_title("Housing prices")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].plot(years, solution.transfers, color="#7A4E2D")
    axes[0, 1].axhline(
        terminal.equal_transfer, color="0.4", linestyle="--", linewidth=1
    )
    axes[0, 1].set_title("Equal property-tax rebate")
    births = np.array(
        [float(row["birth_children_topcode_adjusted"]) for row in solution.rows]
    )
    entries = np.array([float(row["entry_flow_E"]) for row in solution.rows])
    axes[1, 0].plot(years, births, label="Adjusted births")
    axes[1, 0].plot(years, entries, label="Adult entrants")
    axes[1, 0].set_title("Demographic flows")
    axes[1, 0].legend(frameon=False)
    population = np.array([float(row["adult_population"]) for row in solution.rows])
    axes[1, 1].plot(years, population, color="#21476B")
    axes[1, 1].axhline(
        terminal.population_scale, color="0.4", linestyle="--", linewidth=1
    )
    axes[1, 1].set_title("Adult population")
    renewal = np.array([float(row["queue_B_over_current_E"]) for row in solution.rows])
    axes[2, 0].plot(years, renewal, color="#4D7C57")
    axes[2, 0].axhline(1.0, color="0.4", linestyle="--", linewidth=1)
    axes[2, 0].set_title(r"Renewal ratio $B/E$")
    market = np.array(
        [abs(float(row["relative_market_residual"])) for row in solution.rows]
    )
    fiscal = np.array(
        [abs(float(row["government_budget_residual"])) for row in solution.rows]
    )
    axes[2, 1].semilogy(years, np.maximum(market, 1e-12), label="Housing")
    axes[2, 1].semilogy(years, np.maximum(fiscal, 1e-12), label="Fiscal")
    axes[2, 1].axhline(2e-4, color="0.55", linestyle="--", linewidth=1)
    axes[2, 1].axhline(2.5e-5, color="0.70", linestyle=":", linewidth=1)
    axes[2, 1].set_title("Equilibrium residuals")
    axes[2, 1].legend(frameon=False)
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    figure.suptitle(terminal.case.label)
    figure.savefig(output_dir / "diagnostics.png", dpi=200)
    figure.savefig(output_dir / "diagnostics.pdf")
    plt.close(figure)


def solution_stability(
    solutions: dict[int, FundedPFSolution], comparison_periods: int = 5
) -> dict[str, Any]:
    compatible = {horizon: solution for horizon, solution in solutions.items()}
    return pf.horizon_stability(compatible, comparison_periods=comparison_periods)


def main() -> None:
    args = parse_args()
    horizons = sorted(set(int(value) for value in args.horizons))
    if any(horizon < 2 for horizon in horizons):
        raise ValueError("Every horizon must contain at least two dates.")
    if not 0.0 < float(args.market_tol) <= 2e-4:
        raise ValueError("--market-tol must lie in (0,2e-4].")
    if not 0.0 < float(args.fiscal_tol) <= DEFAULT_FISCAL_TOLERANCE:
        raise ValueError(
            f"--fiscal-tol must lie in (0,{DEFAULT_FISCAL_TOLERANCE}]."
        )
    initial_path_sources = parse_initial_paths(args.initial_path)
    if initial_path_sources and len(horizons) != 1:
        raise ValueError("Recorded initial paths require exactly one requested horizon.")
    output_dir = args.output_dir.resolve()
    launcher_files = {"launch_contract.json", "driver.log", "heartbeat.json"}
    unexpected = (
        sorted(path.name for path in output_dir.iterdir() if path.name not in launcher_files)
        if output_dir.exists()
        else []
    )
    if unexpected:
        raise RuntimeError(f"Refusing to overwrite non-launcher output: {unexpected}")
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    chain, model = transition.configure_sequential_model()
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    contracts, provenance = pf.load_diagnostic_contracts(
        args.selected_report, args.selected_transition, args.source
    )
    prepared = continuation.prepare_model(
        contracts, args.source.resolve(), chain, model
    )
    best = contracts["report_best"]
    old_psi = float(best["old_psi_child"])
    new_psi = float(best["new_psi_child"])
    terminal_psi = old_psi if args.terminal_psi_child is None else float(args.terminal_psi_child)

    initial_state, history, history_gate, impact_price = pf.reconstruct_2023_state(
        prepared,
        contracts,
        args.selected_transition.resolve(),
        market_tolerance=float(args.market_tol),
    )
    pf.write_csv(output_dir / "reconstructed_2007_2023_history.csv", history)
    pf.write_json(
        output_dir / "history_reproduction_gate.json",
        {key: value for key, value in history_gate.items() if key != "comparisons"},
    )
    inherited_price = impact.inherited_2023_price(args.selected_transition.resolve())
    if not math.isclose(impact_price, inherited_price, rel_tol=0.0, abs_tol=2e-4):
        raise RuntimeError(
            f"The reconstructed and inherited 2023 prices differ: {impact_price}, {inherited_price}"
        )
    anchor_state = continuation.DynamicState(
        g_pre=initial_state.g_pre.copy(),
        scheduled_entries=list(initial_state.scheduled_entries),
        scheduled_raw_entries=list(initial_state.scheduled_raw_entries),
        price_guess=inherited_price,
        initial_policy=None,
    )
    anchor_parameters = copy.deepcopy(prepared.parameters)
    anchor_parameters.psi_child = new_psi
    anchor = impact.solve_fixed_price_equal_rebate(
        state=anchor_state,
        P=anchor_parameters,
        b_grid=prepared.b_grid,
        price=inherited_price,
        case=impact.CASES["rebated-tax1-baseline"],
        counter=calendar.SolveCounter(),
    )
    supply_rule, reanchor = impact.reanchor_supply_rule(
        inherited_price=inherited_price,
        baseline_demand=float(anchor.evaluation.demand_by_loc[0]),
        elasticity=float(args.housing_supply_elasticity),
    )
    reanchor.update(
        fixed_price_equal_transfer=float(anchor.transfer),
        fixed_price_fiscal_residual=float(anchor.ledger["government_budget_residual"]),
        fixed_price_model_evaluations=int(anchor.model_evaluations),
        normalization_interpretation=(
            "The inherited static 2023 funded-policy normalization is held fixed "
            "when expectations are changed to perfect foresight."
        ),
    )
    pf.write_json(output_dir / "supply_reanchor_gate.json", reanchor)

    terminals: dict[str, PolicyTerminalSteadyState] = {}
    all_solutions: dict[str, dict[int, FundedPFSolution]] = {}
    seed_contracts: dict[str, dict[str, Any]] = {}
    for case_name in args.cases:
        case = impact.CASES[case_name]
        case_dir = output_dir / case_name
        case_dir.mkdir(parents=True, exist_ok=False)
        terminal = solve_policy_terminal_steady_state(
            prepared,
            case=case,
            terminal_psi_child=terminal_psi,
            supply_rule=supply_rule,
            price_min_ratio=float(args.endpoint_price_min_ratio),
            price_max_ratio=float(args.endpoint_price_max_ratio),
            price_grid_points=int(args.endpoint_price_grid_points),
            root_tolerance=float(args.endpoint_root_tol),
            transfer_tolerance=float(args.endpoint_transfer_tol),
            output_dir=case_dir / "terminal_stationary_endpoint",
        )
        terminals[case_name] = terminal
        pf.write_json(
            case_dir / "terminal_fixed_point_gate.json",
            {
                "endpoint": terminal.endpoint,
                "stationary_reconstruction": terminal.reconstruction,
                "fixed_point_gates": terminal.fixed_point_gates,
            },
        )
        case_solutions: dict[int, FundedPFSolution] = {}
        seed_contracts[case_name] = {}
        for horizon in horizons:
            horizon_dir = case_dir / f"horizon_{horizon:03d}"
            horizon_dir.mkdir(parents=True, exist_ok=False)
            psi_values = pf.preference_path(
                new_psi, terminal.psi_child, float(args.psi_persistence), horizon
            )
            if case_name in initial_path_sources:
                initial_prices, initial_transfers, seed = load_initial_path(
                    initial_path_sources[case_name],
                    horizon=horizon,
                    terminal=terminal,
                )
            else:
                initial_prices, initial_transfers, seed = analytic_initial_paths(
                    impact_price=inherited_price,
                    impact_transfer=float(anchor.transfer)
                    * (float(case.annual_tax_rate) / 0.01),
                    terminal=terminal,
                    persistence=float(args.psi_persistence),
                    horizon=horizon,
                )
            seed_contracts[case_name][str(horizon)] = seed
            pf.write_json(output_dir / "initial_path_seed_contract.json", seed_contracts)

            def progress(record: dict[str, Any], evaluation: pf.PathEvaluation) -> None:
                pf.write_json(
                    horizon_dir / "latest_iteration.json",
                    {
                        "status": "running",
                        **record,
                        "price_path": evaluation.prices,
                        "rent_path": evaluation.rents,
                        "transfer_path": [
                            float(row["equal_transfer_period_units"])
                            for row in evaluation.rows
                        ],
                    },
                )
                pf.write_csv(horizon_dir / "latest_path.csv", evaluation.rows)

            solution = solve_funded_perfect_foresight(
                case=case,
                initial_prices=initial_prices,
                initial_transfers=initial_transfers,
                psi_path=psi_values,
                terminal=terminal,
                b_grid=prepared.b_grid,
                initial_state=initial_state,
                supply_rule=supply_rule,
                market_tolerance=float(args.market_tol),
                fiscal_tolerance=float(args.fiscal_tol),
                max_iterations=int(args.max_path_iterations),
                price_damping=float(args.price_damping),
                transfer_damping=float(args.transfer_damping),
                maximum_log_price_step=float(args.maximum_log_price_step),
                maximum_transfer_step=float(args.maximum_transfer_step),
                progress_callback=progress,
            )
            case_solutions[horizon] = solution
            pf.write_csv(horizon_dir / "transition_path.csv", solution.rows)
            pf.write_csv(horizon_dir / "joint_iteration_history.csv", solution.iteration_history)
            pf.write_json(
                horizon_dir / "summary.json",
                {
                    key: value
                    for key, value in vars(solution).items()
                    if key not in {"rows", "iteration_history"}
                },
            )
            make_case_figure(solution, terminal, horizon_dir)
        all_solutions[case_name] = case_solutions
        stability = solution_stability(case_solutions)
        if stability.get("rows"):
            pf.write_csv(case_dir / "horizon_stability.csv", stability["rows"])
        pf.write_json(
            case_dir / "horizon_stability.json",
            {key: value for key, value in stability.items() if key != "rows"},
        )

    comparison: dict[str, Any] = {"status": "not_available_single_case_run"}
    if set(impact.CASES).issubset(all_solutions):
        common_horizons = sorted(
            set(all_solutions["rebated-tax1-baseline"])
            & set(all_solutions["rebated-tax2-reform"])
        )
        maximum_horizon = common_horizons[-1]
        chosen = {
            case_name: all_solutions[case_name][maximum_horizon]
            for case_name in impact.CASES
        }
        effects = effect_rows(
            chosen["rebated-tax1-baseline"].rows,
            chosen["rebated-tax2-reform"].rows,
        )
        pf.write_csv(output_dir / "policy_effects_longest_horizon.csv", effects)
        make_comparison_figure(chosen, output_dir)
        comparison = {
            "status": "computed",
            "horizon": maximum_horizon,
            "impact_effects_2023": effects[0],
            "terminal_steady_state_effects": {
                "asset_price_percent_change": 100.0
                * (
                    terminals["rebated-tax2-reform"].asset_price
                    / terminals["rebated-tax1-baseline"].asset_price
                    - 1.0
                ),
                "renter_price_percent_change": 100.0
                * (
                    terminals["rebated-tax2-reform"].renter_price
                    / terminals["rebated-tax1-baseline"].renter_price
                    - 1.0
                ),
                "equal_transfer_percent_change": 100.0
                * (
                    terminals["rebated-tax2-reform"].equal_transfer
                    / terminals["rebated-tax1-baseline"].equal_transfer
                    - 1.0
                ),
                "adult_population_percent_change": 100.0
                * (
                    terminals["rebated-tax2-reform"].population_scale
                    / terminals["rebated-tax1-baseline"].population_scale
                    - 1.0
                ),
                "births_per_adult_percent_change": 100.0
                * (
                    float(
                        terminals["rebated-tax2-reform"].endpoint[
                            "topcode_adjusted_birth_children_per_adult"
                        ]
                    )
                    / float(
                        terminals["rebated-tax1-baseline"].endpoint[
                            "topcode_adjusted_birth_children_per_adult"
                        ]
                    )
                    - 1.0
                ),
                "owner_rate_pp_change": 100.0
                * (
                    float(terminals["rebated-tax2-reform"].endpoint["owner_rate"])
                    - float(
                        terminals["rebated-tax1-baseline"].endpoint["owner_rate"]
                    )
                ),
            },
        }

    summary = {
        "status": "complete_unpromoted_perfect_foresight_funded_policy_diagnostic",
        "promotion_status": "not_promoted",
        "economic_scope": (
            "Permanent rebated 1% funded-policy reference versus permanent "
            "rebated 2% reform; common fitted 2023 state; common exogenous psi path; "
            "perfect foresight over prices, rents, transfers, and psi; date-by-date "
            "housing and fiscal clearing; closed national renewal M=0,rho=1; no grant."
        ),
        "welfare_claim": "none",
        "shock_interpretation": "declared exogenous path, not estimated or filtered",
        "cases": list(args.cases),
        "horizons": horizons,
        "psi_persistence_per_four_year_period": float(args.psi_persistence),
        "terminal_psi_child": terminal_psi,
        "housing_supply": reanchor,
        "common_initial_state": {
            "distribution_sha256": continuation.array_sha256(initial_state.g_pre),
            "queue_sha256": continuation.canonical_json_sha256(
                initial_state.scheduled_entries
            ),
            "adult_population": float(np.sum(initial_state.g_pre)),
        },
        "terminal_steady_states": {
            case_name: {
                "asset_price": terminal.asset_price,
                "renter_price": terminal.renter_price,
                "equal_transfer": terminal.equal_transfer,
                "adult_population": terminal.population_scale,
                "entry_flow": terminal.entry_flow,
                "births_per_adult": float(
                    terminal.endpoint["topcode_adjusted_birth_children_per_adult"]
                ),
                "owner_rate": float(terminal.endpoint["owner_rate"]),
                "fixed_point_gates": terminal.fixed_point_gates,
            }
            for case_name, terminal in terminals.items()
        },
        "solutions": {
            case_name: {
                str(horizon): {
                    "status": solution.status,
                    "converged": solution.converged,
                    "market_residual": solution.maximum_market_residual,
                    "fiscal_residual": solution.maximum_fiscal_residual,
                    "terminal_convergence": solution.terminal_convergence,
                }
                for horizon, solution in case_solutions.items()
            }
            for case_name, case_solutions in all_solutions.items()
        },
        "comparison": comparison,
        "provenance": {
            **provenance,
            "driver": str(Path(__file__).resolve()),
            "driver_sha256": pf.file_sha256(Path(__file__).resolve()),
            "perfect_foresight_driver": str(Path(pf.__file__).resolve()),
            "perfect_foresight_driver_sha256": pf.file_sha256(Path(pf.__file__).resolve()),
        },
        "elapsed_seconds": time.perf_counter() - started,
    }
    pf.write_json(output_dir / "summary.json", summary)
    print(json.dumps(pf.jsonable(summary), indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
