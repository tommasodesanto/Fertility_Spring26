#!/usr/bin/env python3
"""Isolated perfect-foresight transition diagnostic for the active E5F model.

The experiment starts from the reconstructed 2023 pre-fertility distribution
and closes national population renewal after 2023.  For a declared permanent
fertility-preference regime, it first solves the terminal stationary household
problem jointly with closed renewal, housing clearing, and the population
scale.  Given an asset-price path to that endpoint, the household problem is
solved backward in calendar time.  The household distribution and birth-vintage
queue are then advanced forward.  The complete price path is updated until
housing demand clears the date-0-normalized static supply schedule at every
date.

This is a diagnostic solution-method exercise.  The preference path is not an
estimate, and no result produced here is promoted into the live calibration.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
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
import run_e5f_post2023_no_policy_continuations as continuation


DEFAULT_REPORT = (
    ROOT
    / "output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/summary.json"
)
DEFAULT_SELECTED_TRANSITION = DEFAULT_REPORT.parent / "selected_transition_path.csv"
DEFAULT_SOURCE = closure.E5F_FLOOR_SOURCE
DEFAULT_OUTPUT = (
    ROOT / "output/model/e5f_perfect_foresight_transition_diagnostic_20260823"
)
CALENDAR_START_YEAR = 2023
DEFAULT_PSI_PERSISTENCE = 0.60
DEFAULT_TERMINAL_TOLERANCES = {
    "population_relative_gap": 0.01,
    "normalized_distribution_l1": 0.02,
    "birth_queue_maximum_relative_gap": 0.01,
    "asset_price_relative_gap": 0.01,
    "renter_price_relative_gap": 0.01,
    "psi_absolute_gap": 1e-3,
}


@dataclass
class PFInitialState:
    g_pre: np.ndarray
    scheduled_entries: list[float]
    scheduled_raw_entries: list[float]


@dataclass
class TerminalSteadyState:
    psi_child: float
    asset_price: float
    renter_price: float
    population_scale: float
    entry_flow: float
    raw_queue_flow: float
    policy: calendar.PolicyBundle
    parameters: SimpleNamespace
    state: PFInitialState
    endpoint: dict[str, Any]
    reconstruction: dict[str, Any]
    fixed_point_gates: dict[str, Any]


@dataclass
class PathEvaluation:
    prices: np.ndarray
    rents: np.ndarray
    values: list[np.ndarray]
    rows: list[dict[str, Any]]
    terminal_state: PFInitialState
    maximum_market_residual: float
    maximum_policy_reproduction_error: float
    maximum_mass_accounting_error: float
    maximum_feasibility_projection_mass: float
    bellman_solves: int
    elapsed_seconds: float


@dataclass
class PFSolution:
    status: str
    converged: bool
    iterations: int
    prices: np.ndarray
    rents: np.ndarray
    psi_path: np.ndarray
    rows: list[dict[str, Any]]
    iteration_history: list[dict[str, Any]]
    maximum_market_residual: float
    maximum_log_price_residual: float
    maximum_policy_reproduction_error: float
    maximum_mass_accounting_error: float
    maximum_feasibility_projection_mass: float
    terminal_convergence: dict[str, Any]
    bellman_solves: int
    elapsed_seconds: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--selected-transition", type=Path, default=DEFAULT_SELECTED_TRANSITION
    )
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--horizons", type=int, nargs="+", default=[8, 12])
    parser.add_argument(
        "--psi-persistence", type=float, default=DEFAULT_PSI_PERSISTENCE
    )
    parser.add_argument(
        "--terminal-psi-child",
        type=float,
        help=(
            "Permanent terminal fertility-preference intercept. The selected old "
            "steady-state value is used when omitted."
        ),
    )
    parser.add_argument("--endpoint-price-min-ratio", type=float, default=1.0e-4)
    parser.add_argument("--endpoint-price-max-ratio", type=float, default=3.0)
    parser.add_argument("--endpoint-price-grid-points", type=int, default=25)
    parser.add_argument("--endpoint-root-tol", type=float, default=1.0e-8)
    parser.add_argument("--market-tol", type=float, default=2e-4)
    parser.add_argument("--max-price-iterations", type=int, default=30)
    parser.add_argument("--price-damping", type=float, default=0.25)
    parser.add_argument("--maximum-log-price-step", type=float, default=0.12)
    parser.add_argument(
        "--initial-price-path",
        type=Path,
        help=(
            "Optional transition-path CSV used only as the numerical initial "
            "guess for a single requested horizon. It must contain contiguous "
            "period and asset_price columns starting at period zero. If it is "
            "shorter than the requested horizon, the remaining dates are seeded "
            "at the solved terminal asset price."
        ),
    )
    return parser.parse_args()


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(jsonable(payload), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(jsonable(list(rows)))
    temporary.replace(path)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preference_path(
    initial_psi: float,
    terminal_psi: float,
    persistence: float,
    periods: int,
) -> np.ndarray:
    if periods < 1:
        raise ValueError("The perfect-foresight horizon must be positive.")
    if not 0.0 <= float(persistence) < 1.0:
        raise ValueError("Preference persistence must lie in [0,1).")
    dates = np.arange(int(periods), dtype=float)
    return float(terminal_psi) + float(persistence) ** dates * (
        float(initial_psi) - float(terminal_psi)
    )


def rents_from_asset_prices(
    prices: Sequence[float], terminal_price: float, P: SimpleNamespace
) -> np.ndarray:
    """Apply the one-period owner/renter no-arbitrage identity.

    Owners earn next period's asset price and pay depreciation and property
    tax.  Hence r_t + p_{t+1} = (R + delta + tau_H) p_t.  At a constant
    price this reduces exactly to the stationary user-cost identity used by
    the existing model.
    """
    current = np.asarray(prices, dtype=float).reshape(-1)
    if current.size < 1 or np.any(~np.isfinite(current)) or np.any(current <= 0.0):
        raise ValueError("Asset prices must be finite and strictly positive.")
    terminal = float(terminal_price)
    if not math.isfinite(terminal) or terminal <= 0.0:
        raise ValueError("The terminal asset price must be finite and positive.")
    next_prices = np.r_[current[1:], terminal]
    carrying_factor = float(P.R_gross) + float(P.delta) + float(P.tau_H)
    rents = carrying_factor * current - next_prices
    if np.any(~np.isfinite(rents)) or np.any(rents <= 0.0):
        raise ValueError(
            "The candidate asset-price path implies a nonpositive renter price."
        )
    return rents


def inverse_supply_price(
    demand: Sequence[float], supply_rule: calendar.HousingSupplyRule
) -> np.ndarray:
    quantities = np.asarray(demand, dtype=float).reshape(-1)
    if np.any(~np.isfinite(quantities)) or np.any(quantities <= 0.0):
        raise ValueError("Housing demand must be finite and strictly positive.")
    if supply_rule.mode != "static-elastic":
        raise ValueError(
            "Perfect-foresight price iteration requires the static-elastic supply rule."
        )
    elasticity = float(supply_rule.elasticity)
    if elasticity <= 0.0:
        raise ValueError("Housing-supply elasticity must be positive.")
    return float(supply_rule.initial_price) * (
        quantities / float(supply_rule.initial_stock)
    ) ** (1.0 / elasticity)


def terminal_convergence_diagnostics(
    *,
    evaluation: PathEvaluation,
    psi_path: Sequence[float],
    reference_state: PFInitialState,
    reference_entry_flow: float,
    reference_price: float,
    reference_psi: float,
    base_parameters: SimpleNamespace,
    tolerances: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Check that the endogenous endpoint is genuinely near the terminal SS.

    The terminal value function and next-date asset price are boundary
    conditions, so they do not by themselves establish convergence.  This gate
    instead compares the forward-simulated distribution and birth-entry queue
    at T+1, plus the last endogenous price, rent, and preference at T, with the
    stationary objects.
    """
    active_tolerances = dict(DEFAULT_TERMINAL_TOLERANCES)
    if tolerances is not None:
        unknown = sorted(set(tolerances) - set(active_tolerances))
        if unknown:
            raise ValueError(f"Unknown terminal-convergence tolerances: {unknown}")
        active_tolerances.update(
            {key: float(value) for key, value in tolerances.items()}
        )
    if any(value <= 0.0 for value in active_tolerances.values()):
        raise ValueError("Every terminal-convergence tolerance must be positive.")

    candidate_g = np.asarray(evaluation.terminal_state.g_pre, dtype=float)
    reference_g = np.asarray(reference_state.g_pre, dtype=float)
    if candidate_g.shape != reference_g.shape:
        raise ValueError(
            "Terminal and stationary distributions have different shapes: "
            f"{candidate_g.shape} versus {reference_g.shape}."
        )
    candidate_mass = float(np.sum(candidate_g))
    reference_mass = float(np.sum(reference_g))
    if candidate_mass <= 0.0 or reference_mass <= 0.0:
        raise ValueError("Terminal and stationary distribution mass must be positive.")
    distribution_l1 = float(np.sum(np.abs(candidate_g - reference_g)))
    normalized_distribution_l1 = float(
        np.sum(
            np.abs(candidate_g / candidate_mass - reference_g / reference_mass)
        )
    )
    population_relative_gap = abs(candidate_mass - reference_mass) / reference_mass

    candidate_queue = np.asarray(
        evaluation.terminal_state.scheduled_entries, dtype=float
    ).reshape(-1)
    reference_queue = np.asarray(
        reference_state.scheduled_entries, dtype=float
    ).reshape(-1)
    if candidate_queue.shape != reference_queue.shape:
        raise ValueError(
            "Terminal and stationary birth-entry queues have different shapes."
        )
    entry_scale = max(abs(float(reference_entry_flow)), 1e-15)
    queue_gap = float(np.max(np.abs(candidate_queue - reference_queue)))
    queue_relative_gap = queue_gap / entry_scale

    last_price = float(evaluation.prices[-1])
    last_rent = float(evaluation.rents[-1])
    stationary_price = float(reference_price)
    stationary_rent = float(base_parameters.user_cost_rate) * stationary_price
    last_psi = float(np.asarray(psi_path, dtype=float).reshape(-1)[-1])
    price_relative_gap = abs(last_price - stationary_price) / stationary_price
    rent_relative_gap = abs(last_rent - stationary_rent) / stationary_rent
    psi_absolute_gap = abs(last_psi - float(reference_psi))

    metrics = {
        "population_relative_gap": population_relative_gap,
        "normalized_distribution_l1": normalized_distribution_l1,
        "birth_queue_maximum_relative_gap": queue_relative_gap,
        "asset_price_relative_gap": price_relative_gap,
        "renter_price_relative_gap": rent_relative_gap,
        "psi_absolute_gap": psi_absolute_gap,
    }
    checks = {
        key: metrics[key] <= active_tolerances[key] for key in active_tolerances
    }
    return {
        "status": "passed" if all(checks.values()) else "not_converged",
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "metrics": metrics,
        "tolerances": active_tolerances,
        "terminal_distribution_l1_from_stationary": distribution_l1,
        "terminal_population": candidate_mass,
        "stationary_population": reference_mass,
        "terminal_birth_queue": candidate_queue,
        "stationary_birth_queue": reference_queue,
        "terminal_asset_price_last_endogenous_date": last_price,
        "stationary_asset_price": stationary_price,
        "terminal_renter_price_last_endogenous_date": last_rent,
        "stationary_renter_price": stationary_rent,
        "terminal_psi_last_endogenous_date": last_psi,
        "stationary_psi": float(reference_psi),
    }


def policy_from_objects(
    objects: tuple[Any, ...],
    price: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
) -> calendar.PolicyBundle:
    (
        V,
        c_pol,
        hR_pol,
        bp_pol,
        tenure_choice,
        tenure_probs,
        loc_probs,
        fert_probs,
        fert_value,
        _,
    ) = objects
    price_array = np.array([float(price)])
    return calendar.PolicyBundle(
        V=V,
        c_pol=c_pol,
        hR_pol=hR_pol,
        bp_pol=bp_pol,
        tenure_choice=tenure_choice,
        tenure_probs=tenure_probs,
        loc_probs=loc_probs,
        fert_probs=fert_probs,
        fert_value=fert_value,
        price=price_array,
        maps=calendar.build_transition_maps(price_array, P, b_grid, shared),
    )


def solve_date_policy(
    *,
    price: float,
    rent: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    shared: SimpleNamespace,
    continuation_V: np.ndarray,
) -> calendar.PolicyBundle:
    objects = calendar.model.solve_bellman_full_markov_income(
        np.array([float(rent)]),
        np.array([float(price)]),
        P,
        b_grid,
        shared,
        continuation_V=continuation_V,
    )
    return policy_from_objects(objects, price, P, b_grid, shared)


def backward_value_path(
    *,
    prices: np.ndarray,
    rents: np.ndarray,
    psi_path: np.ndarray,
    terminal_V: np.ndarray,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    transfer_path: Sequence[float] | None = None,
) -> tuple[list[np.ndarray], int]:
    periods = int(len(prices))
    if rents.shape != prices.shape or psi_path.shape != prices.shape:
        raise ValueError("Price, rent, and preference paths must have the same length.")
    transfers = (
        np.zeros(periods, dtype=float)
        if transfer_path is None
        else np.asarray(transfer_path, dtype=float).reshape(-1)
    )
    if transfers.shape != prices.shape:
        raise ValueError("The equal-transfer path must match the price path.")
    if np.any(~np.isfinite(transfers)) or np.any(transfers < 0.0):
        raise ValueError("Equal transfers must be finite and nonnegative.")
    values: list[np.ndarray | None] = [None] * (periods + 1)
    values[-1] = np.asarray(terminal_V, dtype=float)
    solves = 0
    for period in range(periods - 1, -1, -1):
        parameters = copy.deepcopy(base_parameters)
        parameters.psi_child = float(psi_path[period])
        parameters.property_tax_lump_sum_transfer = float(transfers[period])
        shared = calendar.model.precompute_shared(parameters, b_grid)
        policy = solve_date_policy(
            price=float(prices[period]),
            rent=float(rents[period]),
            P=parameters,
            b_grid=b_grid,
            shared=shared,
            continuation_V=np.asarray(values[period + 1], dtype=float),
        )
        values[period] = policy.V
        solves += 1
    return [np.asarray(value, dtype=float) for value in values], solves


def _owner_rate(distribution: np.ndarray) -> float:
    total = float(np.sum(distribution))
    return float(np.sum(distribution[:, 1:, ...])) / max(total, 1e-15)


def evaluate_path_at_prices(
    *,
    prices: Sequence[float],
    psi_path: Sequence[float],
    terminal_price: float,
    terminal_V: np.ndarray,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    initial_state: PFInitialState,
    supply_rule: calendar.HousingSupplyRule,
    birth_to_entry_conversion: float,
    transfer_path: Sequence[float] | None = None,
) -> PathEvaluation:
    started = time.perf_counter()
    price_path = np.asarray(prices, dtype=float).reshape(-1)
    psi_values = np.asarray(psi_path, dtype=float).reshape(-1)
    transfer_values = (
        np.zeros_like(price_path)
        if transfer_path is None
        else np.asarray(transfer_path, dtype=float).reshape(-1)
    )
    if transfer_values.shape != price_path.shape:
        raise ValueError("The equal-transfer path must match the price path.")
    if np.any(~np.isfinite(transfer_values)) or np.any(transfer_values < 0.0):
        raise ValueError("Equal transfers must be finite and nonnegative.")
    rents = rents_from_asset_prices(price_path, terminal_price, base_parameters)
    values, backward_solves = backward_value_path(
        prices=price_path,
        rents=rents,
        psi_path=psi_values,
        terminal_V=terminal_V,
        base_parameters=base_parameters,
        b_grid=b_grid,
        transfer_path=transfer_values,
    )

    state = PFInitialState(
        g_pre=np.asarray(initial_state.g_pre, dtype=float).copy(),
        scheduled_entries=list(initial_state.scheduled_entries),
        scheduled_raw_entries=list(initial_state.scheduled_raw_entries),
    )
    rows: list[dict[str, Any]] = []
    reproduction_error = 0.0
    maximum_mass_error = 0.0
    maximum_projection = 0.0
    forward_solves = 0
    for period, (price, rent, psi, transfer_value) in enumerate(
        zip(price_path, rents, psi_values, transfer_values, strict=True)
    ):
        parameters = copy.deepcopy(base_parameters)
        parameters.psi_child = float(psi)
        parameters.property_tax_lump_sum_transfer = float(transfer_value)
        shared = calendar.model.precompute_shared(parameters, b_grid)
        policy = solve_date_policy(
            price=float(price),
            rent=float(rent),
            P=parameters,
            b_grid=b_grid,
            shared=shared,
            continuation_V=values[period + 1],
        )
        forward_solves += 1
        current_reproduction = float(np.max(np.abs(policy.V - values[period])))
        reproduction_error = max(reproduction_error, current_reproduction)
        evaluation = calendar.evaluate_period(
            np.array([float(price)]),
            state.g_pre,
            parameters,
            b_grid,
            shared,
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
        adjusted_births = float(accounting["topcode_adjusted_birth_children"])
        due_entry, next_queue = transition.advance_birth_vintage_queue(
            state.scheduled_entries,
            adjusted_births,
            birth_to_entry_conversion,
        )
        due_raw, next_raw_queue = transition.advance_birth_vintage_queue(
            state.scheduled_raw_entries,
            float(evaluation.births),
            birth_to_entry_conversion,
        )
        empty_next, model_mature_by_loc, deaths, _ = (
            transition.advance_sequential_calendar_distribution(
                evaluation,
                np.zeros(int(parameters.I)),
                parameters,
                b_grid,
                shared,
            )
        )
        entry_shares = np.asarray(parameters.entry_shares, dtype=float).reshape(-1)
        entry_shares /= float(np.sum(entry_shares))
        entrants_next = float(due_entry) * entry_shares
        empty_next[:, :, :, 0, :, :, :] = calendar.entrant_cohort(
            entrants_next, parameters, b_grid
        )
        expected_mass = (
            float(np.sum(evaluation.g_post_fertility))
            - float(deaths)
            + float(np.sum(entrants_next))
        )
        mass_error = float(np.sum(empty_next)) - expected_mass
        maximum_mass_error = max(maximum_mass_error, abs(mass_error))
        maximum_projection = max(
            maximum_projection, float(evaluation.feasibility_projection_mass)
        )
        health = calendar.distribution_health(
            {
                "pre": evaluation.g_pre,
                "post_fertility": evaluation.g_post_fertility,
                "current": evaluation.g_current,
                "next_pre": empty_next,
            }
        )
        if (
            int(health["nonfinite_distribution_count"]) != 0
            or health["min_distribution_mass"] is None
            or float(health["min_distribution_mass"]) < -1e-13
        ):
            raise RuntimeError(f"Distribution-health gate failed: {health}")
        current_mass = float(np.sum(evaluation.g_current))
        property_tax_revenue = float(
            calendar.model.property_tax_revenue_from_distribution(
                evaluation.g_current,
                evaluation.policy.hR_pol,
                evaluation.policy.price,
                parameters,
            )
        )
        equal_transfer_outlays = float(transfer_value) * current_mass
        government_budget_residual = (
            property_tax_revenue - equal_transfer_outlays
        )
        entry_flow = float(
            np.sum(evaluation.g_pre[:, :, :, 0, :, :, :])
        )
        rows.append(
            {
                "period": period,
                "calendar_year": CALENDAR_START_YEAR
                + period * int(parameters.period_years),
                "psi_child": float(psi),
                "asset_price": float(price),
                "renter_price": float(rent),
                "housing_demand": float(evaluation.demand_by_loc[0]),
                "housing_supply": float(evaluation.supply_by_loc[0]),
                "relative_market_residual": float(
                    evaluation.relative_market_residual
                ),
                "adult_population": current_mass,
                "property_tax_revenue": property_tax_revenue,
                "equal_transfer_period_units": float(transfer_value),
                "equal_transfer_outlays": equal_transfer_outlays,
                "government_budget_residual": government_budget_residual,
                "scaled_government_budget_residual": (
                    government_budget_residual
                    / max(
                        abs(property_tax_revenue),
                        abs(equal_transfer_outlays),
                        1e-12,
                    )
                ),
                "implied_equal_transfer": property_tax_revenue
                / max(current_mass, 1e-12),
                "entry_flow_E": entry_flow,
                "birth_children": float(evaluation.births),
                "birth_children_topcode_adjusted": adjusted_births,
                "effective_mature_entrant_flow_B": float(due_entry),
                "raw_state_scheduled_mature_entrant_flow_B": float(due_raw),
                "entrant_flow_next": float(np.sum(entrants_next)),
                "queue_B_over_current_E": float(due_entry)
                / max(entry_flow, 1e-15),
                "owner_rate": _owner_rate(evaluation.g_current),
                "adult_deaths": float(deaths),
                "model_state_same_period_mature_flow_B": float(
                    parameters.entrant_conversion_factor
                )
                * float(np.sum(model_mature_by_loc)),
                "mass_accounting_residual": mass_error,
                "policy_reproduction_max_abs": current_reproduction,
                "feasibility_frontier_projection_mass": float(
                    evaluation.feasibility_projection_mass
                ),
                "minimum_distribution_mass": float(
                    health["min_distribution_mass"]
                ),
                "nonfinite_distribution_count": int(
                    health["nonfinite_distribution_count"]
                ),
            }
        )
        state = PFInitialState(
            g_pre=empty_next,
            scheduled_entries=list(next_queue),
            scheduled_raw_entries=list(next_raw_queue),
        )

    return PathEvaluation(
        prices=price_path,
        rents=rents,
        values=values,
        rows=rows,
        terminal_state=state,
        maximum_market_residual=max(
            abs(float(row["relative_market_residual"])) for row in rows
        ),
        maximum_policy_reproduction_error=reproduction_error,
        maximum_mass_accounting_error=maximum_mass_error,
        maximum_feasibility_projection_mass=maximum_projection,
        bellman_solves=backward_solves + forward_solves,
        elapsed_seconds=time.perf_counter() - started,
    )


def solve_perfect_foresight(
    *,
    initial_prices: Sequence[float],
    psi_path: Sequence[float],
    terminal_price: float,
    terminal_V: np.ndarray,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    initial_state: PFInitialState,
    supply_rule: calendar.HousingSupplyRule,
    birth_to_entry_conversion: float,
    market_tolerance: float,
    max_iterations: int,
    damping: float,
    maximum_log_price_step: float,
    terminal_reference_state: PFInitialState | None = None,
    terminal_reference_entry_flow: float | None = None,
    terminal_reference_psi: float | None = None,
    terminal_tolerances: dict[str, float] | None = None,
    progress_callback: Callable[[dict[str, Any], PathEvaluation], None] | None = None,
) -> PFSolution:
    started = time.perf_counter()
    prices = np.asarray(initial_prices, dtype=float).reshape(-1).copy()
    psi_values = np.asarray(psi_path, dtype=float).reshape(-1)
    if prices.shape != psi_values.shape:
        raise ValueError("Initial prices and preferences must have the same length.")
    if max_iterations < 1 or not 0.0 < damping <= 1.0:
        raise ValueError("Invalid price-iteration controls.")
    iteration_history: list[dict[str, Any]] = []
    best: tuple[float, float, PathEvaluation] | None = None
    total_bellman = 0
    current_damping = float(damping)
    previous_market_error = math.inf
    converged = False
    last_log_residual = math.inf
    for iteration in range(1, int(max_iterations) + 1):
        evaluation = evaluate_path_at_prices(
            prices=prices,
            psi_path=psi_values,
            terminal_price=terminal_price,
            terminal_V=terminal_V,
            base_parameters=base_parameters,
            b_grid=b_grid,
            initial_state=initial_state,
            supply_rule=supply_rule,
            birth_to_entry_conversion=birth_to_entry_conversion,
        )
        total_bellman += evaluation.bellman_solves
        demand = np.array(
            [float(row["housing_demand"]) for row in evaluation.rows]
        )
        target_prices = inverse_supply_price(demand, supply_rule)
        log_residual = np.log(target_prices) - np.log(prices)
        maximum_log_residual = float(np.max(np.abs(log_residual)))
        last_log_residual = maximum_log_residual
        market_error = float(evaluation.maximum_market_residual)
        if best is None or (market_error, maximum_log_residual) < (best[0], best[1]):
            best = (market_error, maximum_log_residual, evaluation)
        record = {
            "iteration": iteration,
            "maximum_market_residual": market_error,
            "maximum_log_price_residual": maximum_log_residual,
            "price_damping": current_damping,
            "minimum_asset_price": float(np.min(prices)),
            "maximum_asset_price": float(np.max(prices)),
            "minimum_renter_price": float(np.min(evaluation.rents)),
            "maximum_renter_price": float(np.max(evaluation.rents)),
            "bellman_solves_this_iteration": evaluation.bellman_solves,
            "bellman_solves_cumulative": total_bellman,
            "elapsed_seconds": time.perf_counter() - started,
        }
        iteration_history.append(record)
        if progress_callback is not None:
            progress_callback(record, evaluation)
        if market_error <= float(market_tolerance):
            converged = True
            best = (market_error, maximum_log_residual, evaluation)
            break
        if market_error > previous_market_error * 1.02:
            current_damping = max(0.05, 0.5 * current_damping)
        elif market_error < previous_market_error * 0.80:
            current_damping = min(float(damping), 1.08 * current_damping)
        step = np.clip(
            log_residual,
            -float(maximum_log_price_step),
            float(maximum_log_price_step),
        )
        candidate = np.exp(np.log(prices) + current_damping * step)
        # The next-date asset payoff cannot exceed the gross carrying value or
        # the implied rent would be nonpositive.  Shrink the update until the
        # complete candidate path satisfies the household budget identity.
        shrink = 1.0
        while True:
            proposed = np.exp(
                np.log(prices)
                + shrink * (np.log(candidate) - np.log(prices))
            )
            try:
                rents_from_asset_prices(proposed, terminal_price, base_parameters)
                prices = proposed
                break
            except ValueError:
                shrink *= 0.5
                if shrink < 1e-4:
                    raise RuntimeError(
                        "Price iteration could not maintain a positive renter-price path."
                    )
        previous_market_error = market_error

    if best is None:
        raise RuntimeError("Perfect-foresight price iteration produced no evaluation.")
    best_market_error, best_log_residual, best_evaluation = best
    if terminal_reference_state is None:
        terminal_convergence = {
            "status": "not_evaluated",
            "all_checks_pass": False,
        }
    else:
        if terminal_reference_entry_flow is None or terminal_reference_psi is None:
            raise ValueError(
                "Terminal convergence requires the stationary entry flow and psi."
            )
        terminal_convergence = terminal_convergence_diagnostics(
            evaluation=best_evaluation,
            psi_path=psi_values,
            reference_state=terminal_reference_state,
            reference_entry_flow=float(terminal_reference_entry_flow),
            reference_price=float(terminal_price),
            reference_psi=float(terminal_reference_psi),
            base_parameters=base_parameters,
            tolerances=terminal_tolerances,
        )
    return PFSolution(
        status=("converged" if converged else "maximum_iterations_reached"),
        converged=converged,
        iterations=len(iteration_history),
        prices=best_evaluation.prices,
        rents=best_evaluation.rents,
        psi_path=psi_values,
        rows=best_evaluation.rows,
        iteration_history=iteration_history,
        maximum_market_residual=best_market_error,
        maximum_log_price_residual=(
            best_log_residual if best is not None else last_log_residual
        ),
        maximum_policy_reproduction_error=(
            best_evaluation.maximum_policy_reproduction_error
        ),
        maximum_mass_accounting_error=(
            best_evaluation.maximum_mass_accounting_error
        ),
        maximum_feasibility_projection_mass=(
            best_evaluation.maximum_feasibility_projection_mass
        ),
        terminal_convergence=terminal_convergence,
        bellman_solves=total_bellman,
        elapsed_seconds=time.perf_counter() - started,
    )


def load_diagnostic_contracts(
    report_path: Path, selected_transition_path: Path, source_path: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Load the selected calibration while recording the intentional code fork."""
    report_path = report_path.resolve()
    selected_transition_path = selected_transition_path.resolve()
    source_path = source_path.resolve()
    for path in (report_path, selected_transition_path, source_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    raw = json.loads(report_path.read_text(encoding="utf-8"))
    if raw.get("schema") != "e5f_transition_ridge_refinement_report_v1":
        raise RuntimeError("The selected report has the wrong schema.")
    if raw.get("status") != "complete_refinement_with_two_independent_identity_repeats":
        raise RuntimeError("The selected report is not a completed refinement.")
    if not bool(raw.get("promotion_eligible", False)):
        raise RuntimeError("The selected report is not promotion eligible.")
    scientific = dict(raw.get("scientific_contract") or {})
    best = dict(raw.get("best_candidate") or {})
    if not scientific or not best:
        raise RuntimeError("The selected report lacks its scientific contract or winner.")
    source_hash = file_sha256(source_path)
    if source_hash != str(scientific.get("source_sha256")):
        raise RuntimeError("The stationary source does not match the selected report.")
    report = dict(scientific)
    report.update(
        best_candidate=best,
        old_psi_child=best.get("old_psi_child"),
        renewal_accounting_old_state=best.get("renewal_accounting_old_state"),
    )
    contracts = {
        "report": report,
        "raw_report": raw,
        "task": copy.deepcopy(report),
        "report_best": best,
        "task_best": copy.deepcopy(best),
        "renewal_contract": scientific.get("renewal_accounting_contract"),
        "renewal_old_state": best.get("renewal_accounting_old_state"),
        "outside_origin_entry_share": scientific.get(
            "outside_origin_entry_share"
        ),
    }
    provenance = {
        "selected_report": str(report_path),
        "selected_report_sha256": file_sha256(report_path),
        "selected_transition": str(selected_transition_path),
        "selected_transition_sha256": file_sha256(selected_transition_path),
        "stationary_source": str(source_path),
        "stationary_source_sha256": source_hash,
        "selected_code_bundle_sha256": (
            scientific.get("code_fingerprints") or {}
        ).get("bundle_sha256"),
        "perfect_foresight_solver_file": str(Path(__file__).resolve()),
        "perfect_foresight_solver_sha256": file_sha256(Path(__file__).resolve()),
        "calendar_bellman_file": str(
            MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
        ),
        "calendar_bellman_sha256": file_sha256(
            MODEL_ROOT / "intergen_eqscale_seq_optimized/solver.py"
        ),
        "closed_endpoint_solver_file": str(
            TOOLS_ROOT / "run_e5f_post2023_no_policy_continuations.py"
        ),
        "closed_endpoint_solver_sha256": file_sha256(
            TOOLS_ROOT / "run_e5f_post2023_no_policy_continuations.py"
        ),
        "code_contract_interpretation": (
            "The selected economic parameters and stationary source are frozen. "
            "The code-bundle hash is intentionally different because this isolated "
            "driver adds the calendar-time continuation argument; the default-path "
            "bitwise nesting test is reported separately."
        ),
    }
    return contracts, provenance


def recover_stationary_pre_fertility(
    prepared: continuation.PreparedModel,
) -> tuple[np.ndarray, dict[str, Any]]:
    parameters = prepared.parameters
    policy = prepared.old_policy
    b_grid = prepared.b_grid
    shared = calendar.model.precompute_shared(parameters, b_grid)
    _, stats = calendar.model.forward_distribution_markov_income(
        policy.bp_pol,
        policy.hR_pol,
        policy.tenure_choice,
        policy.loc_probs,
        policy.fert_probs,
        policy.V,
        parameters.user_cost_rate * policy.price,
        policy.price,
        parameters,
        b_grid,
        shared,
        fast_stats=False,
        tenure_probs=policy.tenure_probs,
    )
    stationary_solution = SimpleNamespace(
        g_beginning_distribution=np.asarray(
            stats.g_beginning_distribution, dtype=float
        ),
        entry_by_loc=np.asarray(stats.entry_by_loc, dtype=float),
    )
    stationary_g_pre, diagnostics = calendar.reconstruct_stationary_pre_fertility(
        stationary_solution,
        policy,
        parameters,
        b_grid,
        shared,
    )
    return stationary_g_pre, diagnostics


def solve_terminal_steady_state(
    prepared: continuation.PreparedModel,
    *,
    terminal_psi_child: float,
    price_min_ratio: float,
    price_max_ratio: float,
    price_grid_points: int,
    root_tolerance: float,
    output_dir: Path,
) -> TerminalSteadyState:
    """Solve and verify the terminal distribution, queue, price, and scale."""
    endpoint, _ = continuation.solve_closed_stationary_endpoint(
        chain=prepared.chain,
        base_overrides=prepared.base_overrides,
        old_price=prepared.old_price,
        new_psi_child=float(terminal_psi_child),
        supply_rule=prepared.supply_rule,
        price_min_ratio=float(price_min_ratio),
        price_max_ratio=float(price_max_ratio),
        grid_points=int(price_grid_points),
        outdir=output_dir,
        root_residual_tolerance=float(root_tolerance),
    )
    if not bool(endpoint.get("usable_closed_root", False)):
        raise RuntimeError(
            "The declared terminal preference regime has no unique usable closed "
            f"stationary endpoint: {endpoint.get('status')}"
        )

    endpoint_price = float(endpoint["asset_price"])
    solution, parameters, price_array, _ = closure.solve_pe(
        prepared.chain,
        prepared.base_overrides,
        asset_price=endpoint_price,
        psi_child=float(terminal_psi_child),
    )
    solved_price = float(np.asarray(price_array, dtype=float).reshape(-1)[0])
    if not math.isclose(solved_price, endpoint_price, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(
            "The terminal policy solve did not reproduce the endpoint price: "
            f"{solved_price} versus {endpoint_price}."
        )
    b_grid = np.asarray(solution.b_grid, dtype=float)
    if b_grid.shape != prepared.b_grid.shape or not np.allclose(
        b_grid, prepared.b_grid, rtol=0.0, atol=0.0
    ):
        raise RuntimeError("The terminal steady state changed the asset grid.")
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
        raise RuntimeError(
            "The terminal stationary distribution is not unit normalized: "
            f"mass={unit_mass:.12g}."
        )
    scale = float(endpoint["stationary_population_scale"])
    scaled_g_pre = scale * unit_g_pre
    entry_flow = float(np.sum(scaled_g_pre[:, :, :, 0, :, :, :]))
    endpoint_entry_flow = scale * float(endpoint["stationary_entry_flow_E"])
    entry_flow_gap = abs(entry_flow - endpoint_entry_flow)

    evaluation = calendar.evaluate_period(
        price_array,
        scaled_g_pre,
        parameters,
        b_grid,
        shared,
        calendar.SolveCounter(),
        supply_rule=prepared.supply_rule,
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
    queue_relative_gap = abs(implied_queue_flow - entry_flow) / max(
        entry_flow, 1e-15
    )
    state = PFInitialState(
        g_pre=scaled_g_pre,
        scheduled_entries=[entry_flow] * continuation.QUEUE_WAITING_SLOTS,
        scheduled_raw_entries=[raw_queue_flow] * continuation.QUEUE_WAITING_SLOTS,
    )
    one_step = evaluate_path_at_prices(
        prices=np.array([endpoint_price]),
        psi_path=np.array([float(terminal_psi_child)]),
        terminal_price=endpoint_price,
        terminal_V=policy.V,
        base_parameters=parameters,
        b_grid=b_grid,
        initial_state=state,
        supply_rule=prepared.supply_rule,
        birth_to_entry_conversion=conversion,
    )
    next_mass = float(np.sum(one_step.terminal_state.g_pre))
    distribution_l1 = float(
        np.sum(np.abs(one_step.terminal_state.g_pre - scaled_g_pre))
    )
    normalized_distribution_l1 = float(
        np.sum(
            np.abs(
                one_step.terminal_state.g_pre / max(next_mass, 1e-15)
                - scaled_g_pre / scale
            )
        )
    )
    gates = {
        "status": "passed",
        "population_scale": scale,
        "unit_stationary_distribution_mass": unit_mass,
        "scaled_stationary_distribution_mass": float(np.sum(scaled_g_pre)),
        "stationary_entry_flow_from_distribution": entry_flow,
        "stationary_entry_flow_from_endpoint": endpoint_entry_flow,
        "entry_flow_absolute_gap": entry_flow_gap,
        "implied_birth_queue_flow": implied_queue_flow,
        "birth_queue_relative_gap": queue_relative_gap,
        "housing_market_relative_residual": float(
            evaluation.relative_market_residual
        ),
        "one_step_population_absolute_gap": abs(next_mass - scale),
        "one_step_distribution_l1": distribution_l1,
        "one_step_normalized_distribution_l1": normalized_distribution_l1,
        "one_step_policy_reproduction_max_abs": (
            one_step.maximum_policy_reproduction_error
        ),
        "one_step_mass_accounting_max_abs": (
            one_step.maximum_mass_accounting_error
        ),
        "renewal_root_absolute_residual": float(
            endpoint["renewal_root_absolute_residual"]
        ),
        "renewal_root_declared_tolerance": float(root_tolerance),
    }
    acceptance = {
        "entry_flow_absolute_gap": entry_flow_gap <= 2e-10,
        "birth_queue_relative_gap": queue_relative_gap
        <= max(5.0 * float(root_tolerance), 5e-8),
        "housing_market_relative_residual": float(
            evaluation.relative_market_residual
        )
        <= 2.5e-5,
        "one_step_population_absolute_gap": abs(next_mass - scale) <= 2e-8,
        "one_step_normalized_distribution_l1": normalized_distribution_l1
        <= 2e-8,
        "one_step_policy_reproduction_max_abs": (
            one_step.maximum_policy_reproduction_error <= 2e-8
        ),
        "one_step_mass_accounting_max_abs": (
            one_step.maximum_mass_accounting_error <= 2e-8
        ),
    }
    gates["checks"] = acceptance
    gates["status"] = "passed" if all(acceptance.values()) else "failed"
    if gates["status"] != "passed":
        raise RuntimeError(f"Terminal fixed-point gates failed: {gates}")
    return TerminalSteadyState(
        psi_child=float(terminal_psi_child),
        asset_price=endpoint_price,
        renter_price=float(parameters.user_cost_rate) * endpoint_price,
        population_scale=scale,
        entry_flow=entry_flow,
        raw_queue_flow=raw_queue_flow,
        policy=policy,
        parameters=parameters,
        state=state,
        endpoint=endpoint,
        reconstruction=reconstruction,
        fixed_point_gates=gates,
    )


def reconstruct_2023_state(
    prepared: continuation.PreparedModel,
    contracts: dict[str, Any],
    selected_transition_path: Path,
    *,
    market_tolerance: float,
) -> tuple[PFInitialState, list[dict[str, Any]], dict[str, Any], float]:
    best = contracts["report_best"]
    old_psi = float(best["old_psi_child"])
    new_psi = float(best["new_psi_child"])
    report_old = contracts["renewal_old_state"]
    outside_flow = float(report_old["outside_flow_M"])
    retention = float(report_old["retention_rho"])
    raw_B = (
        transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        * prepared.births_old_raw
    )
    state = continuation.DynamicState(
        g_pre=prepared.initial_g_pre_2007.copy(),
        scheduled_entries=[prepared.B_old] * continuation.QUEUE_WAITING_SLOTS,
        scheduled_raw_entries=[raw_B] * continuation.QUEUE_WAITING_SLOTS,
        price_guess=prepared.old_price,
        initial_policy=None,
    )
    parameters = copy.deepcopy(prepared.parameters)
    counter = calendar.SolveCounter()
    target_years = transition.census_household_age_target_years()
    common_history: list[dict[str, Any]] = []
    previous_psi: float | None = None
    for period in range(continuation.TRANSITION_PERIODS):
        dated_psi = transition.preference_shifter_at_date(
            period, old_psi, new_psi, continuation.TRANSITION_PERIODS
        )
        parameters.psi_child = dated_psi
        if previous_psi is not None and not math.isclose(
            dated_psi, previous_psi, rel_tol=0.0, abs_tol=1e-14
        ):
            state.initial_policy = None
        previous_psi = dated_psi
        evaluation, shared, fallback = continuation.evaluate_state(
            state,
            parameters,
            prepared.b_grid,
            prepared.supply_rule,
            counter,
            market_tolerance,
            30,
        )
        row, state = continuation.advance_from_evaluation(
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
            next_bridge_year=target_years[period + 1],
            grid_fallback=fallback,
        )
        common_history.append(row)

    state_2023 = continuation.DynamicState(
        g_pre=state.g_pre.copy(),
        scheduled_entries=list(state.scheduled_entries),
        scheduled_raw_entries=list(state.scheduled_raw_entries),
        price_guess=state.price_guess,
        initial_policy=None,
    )
    parameters.psi_child = new_psi
    evaluation_2023, shared_2023, fallback_2023 = continuation.evaluate_state(
        state_2023,
        parameters,
        prepared.b_grid,
        prepared.supply_rule,
        counter,
        market_tolerance,
        30,
    )
    row_2023, _ = continuation.advance_from_evaluation(
        label="open_cbsa_sensitivity",
        period_from_2007=continuation.TRANSITION_PERIODS,
        evaluation=evaluation_2023,
        state=copy.deepcopy(state_2023),
        P=parameters,
        b_grid=prepared.b_grid,
        shared=shared_2023,
        supply_rule=prepared.supply_rule,
        outside_flow=outside_flow,
        retention=retention,
        initial_mass_2007=prepared.initial_mass_2007,
        mass_2023=float(np.sum(state_2023.g_pre)),
        next_bridge_year=None,
        grid_fallback=fallback_2023,
    )
    reproduction = continuation.validate_reconstructed_history(
        common_history,
        [row_2023],
        selected_transition_path,
    )
    return (
        PFInitialState(
            g_pre=state_2023.g_pre.copy(),
            scheduled_entries=list(state_2023.scheduled_entries),
            scheduled_raw_entries=list(state_2023.scheduled_raw_entries),
        ),
        common_history + [row_2023],
        reproduction,
        float(evaluation_2023.policy.price[0]),
    )


def zero_shock_test(
    prepared: continuation.PreparedModel,
    stationary_g_pre: np.ndarray,
    *,
    periods: int = 3,
) -> dict[str, Any]:
    parameters = copy.deepcopy(prepared.parameters)
    old_psi = float(parameters.psi_child)
    shared = calendar.model.precompute_shared(parameters, prepared.b_grid)
    stationary_rule, normalization = calendar.normalize_date0_housing_supply(
        stationary_g_pre,
        prepared.old_policy,
        parameters,
        prepared.b_grid,
        shared,
        "static-elastic",
    )
    evaluation = calendar.evaluate_period(
        prepared.old_policy.price,
        stationary_g_pre,
        parameters,
        prepared.b_grid,
        shared,
        calendar.SolveCounter(),
        supply_rule=stationary_rule,
        supplied_policy=prepared.old_policy,
    )
    accounting = transition.calendar_topcode_birth_accounting(
        evaluation.g_pre,
        evaluation.g_post_fertility,
        float(evaluation.births),
        parameters,
    )
    adjusted_births = float(accounting["topcode_adjusted_birth_children"])
    exact_conversion = prepared.E_old / adjusted_births
    raw_conversion = prepared.E_old / max(float(evaluation.births), 1e-15)
    initial = PFInitialState(
        g_pre=stationary_g_pre.copy(),
        scheduled_entries=[prepared.E_old] * continuation.QUEUE_WAITING_SLOTS,
        scheduled_raw_entries=[prepared.E_old] * continuation.QUEUE_WAITING_SLOTS,
    )
    path = evaluate_path_at_prices(
        prices=np.full(periods, prepared.old_price),
        psi_path=np.full(periods, old_psi),
        terminal_price=prepared.old_price,
        terminal_V=prepared.old_policy.V,
        base_parameters=parameters,
        b_grid=prepared.b_grid,
        initial_state=initial,
        supply_rule=stationary_rule,
        birth_to_entry_conversion=exact_conversion,
    )
    distribution_l1 = float(
        np.sum(np.abs(path.terminal_state.g_pre - stationary_g_pre))
    )
    stationary_rent = float(parameters.user_cost_rate * prepared.old_price)
    rent_identity_error = float(np.max(np.abs(path.rents - stationary_rent)))
    return {
        "status": (
            "passed"
            if max(
                path.maximum_market_residual,
                path.maximum_policy_reproduction_error,
                path.maximum_mass_accounting_error,
                distribution_l1,
                rent_identity_error,
            )
            <= 5e-8
            else "failed"
        ),
        "periods": periods,
        "maximum_market_residual": path.maximum_market_residual,
        "maximum_policy_reproduction_error": (
            path.maximum_policy_reproduction_error
        ),
        "maximum_mass_accounting_error": path.maximum_mass_accounting_error,
        "terminal_distribution_l1_from_stationary": distribution_l1,
        "stationary_rent_identity_error": rent_identity_error,
        "operator_exact_birth_to_entry_conversion": exact_conversion,
        "raw_flow_exact_birth_to_entry_conversion": raw_conversion,
        "maintained_birth_to_entry_conversion": (
            transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
        ),
        "maintained_B_over_E": prepared.B_old / prepared.E_old,
        "supply_normalization": normalization,
    }


def initial_price_path(
    impact_price: float,
    terminal_price: float,
    persistence: float,
    periods: int,
) -> np.ndarray:
    weights = float(persistence) ** np.arange(periods, dtype=float)
    return np.exp(
        np.log(float(terminal_price))
        + weights * (np.log(float(impact_price)) - np.log(float(terminal_price)))
    )


def load_initial_price_path(
    path: Path,
    *,
    horizon: int,
    terminal_price: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load a recorded price-path seed and extend it at the terminal price.

    The seed affects only the starting point of time-path iteration. Market
    clearing and every terminal-convergence gate are evaluated on the solved
    path, so accepting a seed cannot relax the equilibrium contract.
    """
    source = path.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    requested_horizon = int(horizon)
    if requested_horizon < 1:
        raise ValueError("The seeded perfect-foresight horizon must be positive.")
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Initial price seed is empty: {source}")
    missing = {"period", "asset_price"} - set(rows[0])
    if missing:
        raise ValueError(
            f"Initial price seed is missing required columns {sorted(missing)}: "
            f"{source}"
        )
    if len(rows) > requested_horizon:
        raise ValueError(
            f"Initial price seed has {len(rows)} rows but the requested horizon "
            f"has only {requested_horizon}: {source}"
        )
    periods: list[int] = []
    prices: list[float] = []
    for row_number, row in enumerate(rows, start=2):
        try:
            period = int(row["period"])
            price = float(row["asset_price"])
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Invalid period or asset_price in {source} at CSV row "
                f"{row_number}."
            ) from error
        periods.append(period)
        prices.append(price)
    expected_periods = list(range(len(rows)))
    if periods != expected_periods:
        raise ValueError(
            "Initial price seed periods must be contiguous and start at zero: "
            f"observed={periods}; expected={expected_periods}."
        )
    price_array = np.asarray(prices, dtype=float)
    if np.any(~np.isfinite(price_array)) or np.any(price_array <= 0.0):
        raise ValueError("Every seeded asset price must be finite and positive.")
    seeded = np.full(requested_horizon, float(terminal_price), dtype=float)
    seeded[: len(price_array)] = price_array
    return seeded, {
        "mode": "recorded_transition_path_then_terminal_extension",
        "source": source,
        "source_sha256": file_sha256(source),
        "source_rows": len(rows),
        "requested_horizon": requested_horizon,
        "terminal_extension_periods": requested_horizon - len(rows),
        "extension_asset_price": float(terminal_price),
    }


def make_figures(solution: PFSolution, output_dir: Path, label: str) -> None:
    years = np.array([int(row["calendar_year"]) for row in solution.rows])
    figure, axes = plt.subplots(3, 2, figsize=(10.5, 9.2), constrained_layout=True)
    axes[0, 0].plot(years, solution.prices, marker="o", label="Asset price")
    axes[0, 0].plot(years, solution.rents, marker="s", label="Renter price")
    axes[0, 0].set_title("Housing prices")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].plot(years, solution.psi_path, marker="o", color="#B4433C")
    axes[0, 1].set_title("Diagnostic fertility-preference path")
    births = np.array(
        [float(row["birth_children_topcode_adjusted"]) for row in solution.rows]
    )
    entries = np.array([float(row["entry_flow_E"]) for row in solution.rows])
    axes[1, 0].plot(years, births, marker="o", label="Adjusted births")
    axes[1, 0].plot(years, entries, marker="s", label="Entrants")
    axes[1, 0].set_title("Demographic flows")
    axes[1, 0].legend(frameon=False)
    population = np.array(
        [float(row["adult_population"]) for row in solution.rows]
    )
    terminal_population = float(
        solution.terminal_convergence["stationary_population"]
    )
    axes[1, 1].plot(years, population, marker="o", color="#21476B")
    axes[1, 1].axhline(
        terminal_population,
        color="0.35",
        linestyle="--",
        linewidth=1,
        label="Terminal steady state",
    )
    axes[1, 1].set_title("Adult population")
    axes[1, 1].legend(frameon=False)
    renewal = np.array(
        [float(row["queue_B_over_current_E"]) for row in solution.rows]
    )
    axes[2, 0].plot(years, renewal, marker="o", color="#4D7C57")
    axes[2, 0].axhline(1.0, color="0.35", linestyle="--", linewidth=1)
    axes[2, 0].set_title(r"Renewal ratio $B/E$")
    residual = np.array(
        [abs(float(row["relative_market_residual"])) for row in solution.rows]
    )
    axes[2, 1].semilogy(years, np.maximum(residual, 1e-12), marker="o")
    axes[2, 1].axhline(2e-4, color="0.5", linestyle="--", linewidth=1)
    axes[2, 1].set_title("Housing-market residual")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.set_xlabel("Year")
    for suffix in ("png", "pdf"):
        figure.savefig(output_dir / f"{label}_diagnostics.{suffix}", dpi=180)
    plt.close(figure)


def horizon_stability(
    solutions: dict[int, PFSolution], comparison_periods: int = 5
) -> dict[str, Any]:
    horizons = sorted(solutions)
    if len(horizons) < 2:
        return {"status": "not_applicable", "horizons": horizons}
    reference_horizon = horizons[-1]
    reference = solutions[reference_horizon]
    rows: list[dict[str, Any]] = []
    for horizon in horizons[:-1]:
        candidate = solutions[horizon]
        common = min(comparison_periods, horizon, reference_horizon)
        for name, left, right in (
            ("asset_price", candidate.prices, reference.prices),
            ("renter_price", candidate.rents, reference.rents),
            (
                "adjusted_births",
                np.array(
                    [
                        float(row["birth_children_topcode_adjusted"])
                        for row in candidate.rows
                    ]
                ),
                np.array(
                    [
                        float(row["birth_children_topcode_adjusted"])
                        for row in reference.rows
                    ]
                ),
            ),
            (
                "adult_population",
                np.array(
                    [float(row["adult_population"]) for row in candidate.rows]
                ),
                np.array(
                    [float(row["adult_population"]) for row in reference.rows]
                ),
            ),
        ):
            scale = np.maximum(np.abs(right[:common]), 1e-12)
            rows.append(
                {
                    "short_horizon": horizon,
                    "reference_horizon": reference_horizon,
                    "comparison_periods": common,
                    "variable": name,
                    "maximum_absolute_gap": float(
                        np.max(np.abs(left[:common] - right[:common]))
                    ),
                    "maximum_relative_gap": float(
                        np.max(np.abs(left[:common] - right[:common]) / scale)
                    ),
                }
            )
    return {
        "status": "computed",
        "horizons": horizons,
        "reference_horizon": reference_horizon,
        "comparison_periods": comparison_periods,
        "maximum_relative_gap": max(
            float(row["maximum_relative_gap"]) for row in rows
        ),
        "rows": rows,
    }


def main() -> None:
    args = parse_args()
    requested_horizons = sorted(set(int(value) for value in args.horizons))
    if any(horizon < 2 for horizon in requested_horizons):
        raise ValueError("Every requested horizon must contain at least two dates.")
    if args.initial_price_path is not None and len(requested_horizons) != 1:
        raise ValueError(
            "--initial-price-path requires exactly one requested horizon."
        )
    if not 0.0 < float(args.market_tol) <= 2e-4:
        raise ValueError("--market-tol must lie in (0,2e-4].")
    if args.terminal_psi_child is not None and not math.isfinite(
        float(args.terminal_psi_child)
    ):
        raise ValueError("--terminal-psi-child must be finite.")
    started = time.perf_counter()
    output_dir = args.output_dir.resolve()
    launcher_files = {"launch_contract.json", "driver.log", "heartbeat.json"}
    unexpected_existing = (
        sorted(path.name for path in output_dir.iterdir() if path.name not in launcher_files)
        if output_dir.exists()
        else []
    )
    if unexpected_existing:
        raise RuntimeError(
            "Refusing to overwrite output directory containing non-launcher files: "
            f"{output_dir}; found={unexpected_existing}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    chain, model = transition.configure_sequential_model()
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = (
        transition.advance_sequential_calendar_distribution
    )
    calendar.distribution_rows = transition.independent_child_distribution_rows
    contracts, provenance = load_diagnostic_contracts(
        args.selected_report, args.selected_transition, args.source
    )
    prepared = continuation.prepare_model(
        contracts, args.source.resolve(), chain, model
    )
    stationary_g_pre, stationary_reconstruction = recover_stationary_pre_fertility(
        prepared
    )
    zero_shock = zero_shock_test(prepared, stationary_g_pre)
    write_json(output_dir / "zero_shock_test.json", zero_shock)
    if zero_shock["status"] != "passed":
        raise RuntimeError(f"Zero-shock perfect-foresight nesting failed: {zero_shock}")

    best = contracts["report_best"]
    old_psi = float(best["old_psi_child"])
    new_psi = float(best["new_psi_child"])
    terminal_psi = (
        old_psi
        if args.terminal_psi_child is None
        else float(args.terminal_psi_child)
    )
    terminal = solve_terminal_steady_state(
        prepared,
        terminal_psi_child=terminal_psi,
        price_min_ratio=float(args.endpoint_price_min_ratio),
        price_max_ratio=float(args.endpoint_price_max_ratio),
        price_grid_points=int(args.endpoint_price_grid_points),
        root_tolerance=float(args.endpoint_root_tol),
        output_dir=output_dir / "terminal_stationary_endpoint",
    )
    write_json(
        output_dir / "terminal_fixed_point_gate.json",
        {
            "endpoint": terminal.endpoint,
            "stationary_reconstruction": terminal.reconstruction,
            "fixed_point_gates": terminal.fixed_point_gates,
        },
    )

    initial_state, history, history_gate, impact_price = reconstruct_2023_state(
        prepared,
        contracts,
        args.selected_transition.resolve(),
        market_tolerance=float(args.market_tol),
    )
    write_csv(output_dir / "reconstructed_2007_2023_history.csv", history)
    write_json(
        output_dir / "history_reproduction_gate.json",
        {key: value for key, value in history_gate.items() if key != "comparisons"},
    )
    solutions: dict[int, PFSolution] = {}
    initial_price_contracts: dict[str, dict[str, Any]] = {}
    for horizon in requested_horizons:
        label = f"horizon_{horizon:02d}"
        horizon_dir = output_dir / label
        horizon_dir.mkdir(parents=True, exist_ok=False)
        psi_values = preference_path(
            new_psi, terminal.psi_child, float(args.psi_persistence), horizon
        )
        if args.initial_price_path is None:
            initial_prices = initial_price_path(
                impact_price,
                terminal.asset_price,
                float(args.psi_persistence),
                horizon,
            )
            price_contract = {
                "mode": "analytic_impact_to_terminal_decay",
                "requested_horizon": horizon,
                "impact_asset_price": impact_price,
                "terminal_asset_price": terminal.asset_price,
                "persistence": float(args.psi_persistence),
            }
        else:
            initial_prices, price_contract = load_initial_price_path(
                args.initial_price_path,
                horizon=horizon,
                terminal_price=terminal.asset_price,
            )
        initial_price_contracts[str(horizon)] = price_contract
        write_json(
            output_dir / "initial_price_seed_contract.json",
            initial_price_contracts,
        )

        def progress(record: dict[str, Any], evaluation: PathEvaluation) -> None:
            write_json(
                horizon_dir / "latest_iteration.json",
                {
                    "status": "running",
                    **record,
                    "price_path": evaluation.prices,
                    "rent_path": evaluation.rents,
                },
            )
            write_csv(horizon_dir / "latest_path.csv", evaluation.rows)

        solution = solve_perfect_foresight(
            initial_prices=initial_prices,
            psi_path=psi_values,
            terminal_price=terminal.asset_price,
            terminal_V=terminal.policy.V,
            base_parameters=terminal.parameters,
            b_grid=prepared.b_grid,
            initial_state=initial_state,
            supply_rule=prepared.supply_rule,
            birth_to_entry_conversion=(
                transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION
            ),
            market_tolerance=float(args.market_tol),
            max_iterations=int(args.max_price_iterations),
            damping=float(args.price_damping),
            maximum_log_price_step=float(args.maximum_log_price_step),
            terminal_reference_state=terminal.state,
            terminal_reference_entry_flow=terminal.entry_flow,
            terminal_reference_psi=terminal.psi_child,
            progress_callback=progress,
        )
        solutions[horizon] = solution
        write_csv(horizon_dir / "transition_path.csv", solution.rows)
        write_csv(horizon_dir / "price_iteration_history.csv", solution.iteration_history)
        write_json(
            horizon_dir / "summary.json",
            {
                key: value
                for key, value in vars(solution).items()
                if key not in {"rows", "iteration_history"}
            },
        )
        make_figures(solution, horizon_dir, label)

    stability = horizon_stability(solutions)
    terminally_converged_horizons = [
        horizon
        for horizon, solution in sorted(solutions.items())
        if bool(solution.terminal_convergence.get("all_checks_pass", False))
    ]
    if stability.get("rows"):
        write_csv(output_dir / "horizon_stability.csv", stability["rows"])
    write_json(
        output_dir / "horizon_stability.json",
        {key: value for key, value in stability.items() if key != "rows"},
    )
    summary = {
        "status": "complete_isolated_perfect_foresight_diagnostic",
        "promotion_status": "not_promoted",
        "economic_scope": (
            "Closed national renewal after 2023; observed 2023 household state; "
            "exogenous diagnostic convergence of psi_child to a declared permanent "
            "regime whose distribution, queue, price, and population scale are "
            "jointly solved; households have perfect foresight over prices and psi."
        ),
        "shock_interpretation": "diagnostic path, not estimated or filtered",
        "psi_persistence_per_four_year_period": float(args.psi_persistence),
        "terminal_steady_state": {
            "psi_child": terminal.psi_child,
            "asset_price": terminal.asset_price,
            "renter_price": terminal.renter_price,
            "adult_population": terminal.population_scale,
            "entry_flow_E": terminal.entry_flow,
            "completed_fertility": float(terminal.endpoint["tfr_topcoded"]),
            "B_over_E": float(terminal.endpoint["queue_B_over_E"]),
            "fixed_point_gates": terminal.fixed_point_gates,
        },
        "initial_2023": {
            "psi_child": new_psi,
            "temporary_equilibrium_price_used_only_as_initial_guess": impact_price,
            "adult_population": float(np.sum(initial_state.g_pre)),
            "terminal_population_ratio": terminal.population_scale
            / float(np.sum(initial_state.g_pre)),
        },
        "no_arbitrage_identity": (
            "r_t + p_{t+1} = (R_gross + delta + tau_H) p_t"
        ),
        "population_law": "M=0, rho=1, four-slot birth-vintage queue",
        "terminal_convergence_status": (
            "passed"
            if terminally_converged_horizons
            else "no_requested_horizon_passed"
        ),
        "terminally_converged_horizons": terminally_converged_horizons,
        "terminal_convergence_rule": (
            "A transition is terminally converged only if the forward-simulated "
            "population level, normalized household distribution, inherited "
            "birth-entry queue, last endogenous asset price and renter price, "
            "and psi are all within their recorded tolerances of the terminal "
            "steady state."
        ),
        "zero_shock_test": zero_shock,
        "history_reproduction_gate": {
            key: value for key, value in history_gate.items() if key != "comparisons"
        },
        "stationary_reconstruction": stationary_reconstruction,
        "initial_price_seeds": initial_price_contracts,
        "horizons": {
            str(horizon): {
                "status": solution.status,
                "converged": solution.converged,
                "iterations": solution.iterations,
                "maximum_market_residual": solution.maximum_market_residual,
                "maximum_log_price_residual": solution.maximum_log_price_residual,
                "maximum_policy_reproduction_error": (
                    solution.maximum_policy_reproduction_error
                ),
                "maximum_mass_accounting_error": (
                    solution.maximum_mass_accounting_error
                ),
                "maximum_feasibility_projection_mass": (
                    solution.maximum_feasibility_projection_mass
                ),
                "terminal_convergence": solution.terminal_convergence,
                "bellman_solves": solution.bellman_solves,
                "elapsed_seconds": solution.elapsed_seconds,
            }
            for horizon, solution in solutions.items()
        },
        "horizon_stability": {
            key: value for key, value in stability.items() if key != "rows"
        },
        "provenance": provenance,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json(output_dir / "summary.json", summary)
    write_json(
        output_dir / "latest_completed_case.json",
        {
            "status": "complete",
            "output_dir": output_dir,
            "horizons": sorted(solutions),
            "elapsed_seconds": summary["elapsed_seconds"],
        },
    )
    readme = f"""# Isolated perfect-foresight transition diagnostic

This packet is not part of the live calibration and has not been promoted.
It starts from the reconstructed 2023 pre-fertility household distribution,
sets post-2023 outside entry to zero, retains every inherited birth-vintage
queue entry, and solves household choices backward against the complete price
and preference paths.  Distributions and the birth queue move forward.

The fertility-preference path is diagnostic: its four-year persistence is
`{float(args.psi_persistence):.6g}` and its terminal value is the old closed
steady-state value.  It is not estimated or filtered from data.

The owner/renter identity is
`r_t + p_(t+1) = (R_gross + delta + tau_H) p_t`.  Housing supply remains the
existing date-0-normalized static-elastic schedule in the asset price.

Inspect `summary.json`, `zero_shock_test.json`, `horizon_stability.json`, and
each `horizon_XX/transition_path.csv` plus the paired diagnostic figure.

`initial_price_seed_contract.json` records the numerical starting path and its
hash when a previous converged transition is used to accelerate time-path
iteration.  A seed never changes the equilibrium equations or acceptance gates.

The packet now treats the terminal value function and next-date price only as
boundary conditions.  A horizon is labeled terminally converged only when its
forward distribution, population level, birth-entry queue, last endogenous
prices, and preference shifter all pass the explicit steady-state proximity
checks recorded in that horizon's `summary.json`.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")


if __name__ == "__main__":
    main()
