"""Dated financing-share path, copied from the pinned perfect-foresight functions.

Only the added phi_path input and per-date P.phi assignments differ.  The
existing household, population, fiscal and diagnostic gates remain in force.
The review fixture compares the non-policy AST with the upstream source.
"""
from __future__ import annotations
import copy
import time
from types import SimpleNamespace
from typing import Any, Callable, Sequence
import numpy as np
import run_e5f_perfect_foresight_transition as base
from run_e5f_perfect_foresight_transition import (
    PFInitialState, PathEvaluation, HistoricalConditioning, calendar,
    social_security, transition, copy_birth_queue, entry_clock_timing,
    rents_from_asset_prices, _owner_rate, CALENDAR_START_YEAR,
    validate_entry_queues,
)

def backward_value_path(
    *,
    prices: np.ndarray,
    rents: np.ndarray,
    psi_path: np.ndarray,
    phi_path: np.ndarray,
    terminal_V: np.ndarray,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    transfer_path: Sequence[float] | None = None,
    pension_path: Sequence[float] | None = None,
    payroll_tax_path: Sequence[float] | None = None,
) -> tuple[list[np.ndarray], int]:
    periods = int(len(prices))
    pensions, payroll_taxes = social_security.validated_fiscal_paths(
        periods, pension_path, payroll_tax_path)
    phi_values = np.asarray(phi_path, dtype=float).reshape(-1)
    if (rents.shape != prices.shape or psi_path.shape != prices.shape
            or phi_values.shape != prices.shape or not np.isfinite(phi_values).all()
            or np.any((phi_values < 0.0) | (phi_values > 1.0))):
        raise ValueError("Price, rent, preference, and financing paths must align; phi must lie in [0,1].")
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
        parameters.phi = np.full_like(np.asarray(parameters.phi, dtype=float), float(phi_values[period]))
        parameters.property_tax_lump_sum_transfer = float(transfers[period])
        social_security.apply_fiscal_date(parameters, period, pensions, payroll_taxes)
        shared = calendar.model.precompute_shared(parameters, b_grid)
        policy = base.solve_date_policy(
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

def evaluate_path_at_prices(
    *,
    prices: Sequence[float],
    psi_path: Sequence[float],
    phi_path: Sequence[float],
    terminal_price: float,
    terminal_V: np.ndarray,
    base_parameters: SimpleNamespace,
    b_grid: np.ndarray,
    initial_state: PFInitialState,
    supply_rule: calendar.HousingSupplyRule,
    birth_to_entry_conversion: float,
    transfer_path: Sequence[float] | None = None,
    historical_conditioning: HistoricalConditioning | None = None,
    pension_path: Sequence[float] | None = None,
    payroll_tax_path: Sequence[float] | None = None,
    dated_observer: Callable[[int, calendar.PeriodEvaluation, SimpleNamespace,
                              np.ndarray, SimpleNamespace, np.ndarray], None] | None = None,
) -> PathEvaluation:
    # Observer receives date-t evaluation and the actual date-(t+1) entrant cohort.
    # It may audit or raise; it must not mutate the model objects supplied to it.
    validate_entry_queues(initial_state, base_parameters)
    if dated_observer is not None and not callable(dated_observer):
        raise TypeError("Dated observer must be callable")
    started = time.perf_counter()
    price_path = np.asarray(prices, dtype=float).reshape(-1)
    psi_values = np.asarray(psi_path, dtype=float).reshape(-1)
    phi_values = np.asarray(phi_path, dtype=float).reshape(-1)
    if (phi_values.shape != price_path.shape or not np.isfinite(phi_values).all()
            or np.any((phi_values < 0.0) | (phi_values > 1.0))):
        raise ValueError("Financing path must match prices and stay in [0,1].")
    pensions, payroll_taxes = social_security.validated_fiscal_paths(
        len(price_path), pension_path, payroll_tax_path)
    if historical_conditioning is not None:
        historical_conditioning.validate(base_parameters, len(price_path), initial_state,
                                         birth_to_entry_conversion)
        if psi_values.shape != price_path.shape or not np.isfinite(psi_values).all():
            raise ValueError("Historical preference path must match the finite price path")
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
        phi_path=phi_values,
        terminal_V=terminal_V,
        base_parameters=base_parameters,
        b_grid=b_grid,
        transfer_path=transfer_values,
        pension_path=pensions,
        payroll_tax_path=payroll_taxes,
    )

    state = PFInitialState(
        g_pre=np.asarray(initial_state.g_pre, dtype=float).copy(),
        scheduled_entries=copy_birth_queue(initial_state.scheduled_entries),
        scheduled_raw_entries=copy_birth_queue(initial_state.scheduled_raw_entries),
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
        parameters.phi = np.full_like(np.asarray(parameters.phi, dtype=float), float(phi_values[period]))
        parameters.property_tax_lump_sum_transfer = float(transfer_value)
        social_security.apply_fiscal_date(parameters, period, pensions, payroll_taxes)
        shared = calendar.model.precompute_shared(parameters, b_grid)
        policy = base.solve_date_policy(
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
        if historical_conditioning is not None and historical_conditioning.observer is not None:
            historical_conditioning.observer(period, evaluation, parameters, b_grid, shared)
        accounting = transition.calendar_topcode_birth_accounting(
            evaluation.g_pre,
            evaluation.g_post_fertility,
            float(evaluation.births),
            parameters,
        )
        adjusted_births = float(accounting["topcode_adjusted_birth_children"])
        due_entry, next_queue = transition.advance_adult_entry_clock(
            state.scheduled_entries,
            adjusted_births,
            birth_to_entry_conversion,
            entry_clock_timing(parameters),
        )
        due_raw, next_raw_queue = transition.advance_adult_entry_clock(
            state.scheduled_raw_entries,
            float(evaluation.births),
            birth_to_entry_conversion,
            entry_clock_timing(parameters),
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
        bridge_audit = None
        if historical_conditioning is not None:
            entrants_next = historical_conditioning.outside_flow * entry_shares
            entrants_next[0] += historical_conditioning.retention * float(due_entry)
        next_entrant_cohort = calendar.entrant_cohort(entrants_next, parameters, b_grid)
        empty_next[:, :, :, 0, :, :, :] = next_entrant_cohort
        if dated_observer is not None:
            dated_observer(period, evaluation, parameters, b_grid, shared, next_entrant_cohort)
        if historical_conditioning is not None:
            next_year = historical_conditioning.next_age_targets.get(period + 1)
            if next_year is not None:
                ages = float(parameters.age_start) + np.arange(int(parameters.J)) * float(parameters.da)
                empty_next, bridge_audit = transition.reweight_distribution_to_observed_age_path(
                    empty_next, ages, year=next_year, initial_mass=historical_conditioning.initial_mass,
                )
        expected_mass = (
            float(np.sum(evaluation.g_post_fertility))
            - float(deaths)
            + float(np.sum(entrants_next))
        )
        if bridge_audit is not None:
            expected_mass += float(bridge_audit["net_residual"])
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
                "calendar_year": (CALENDAR_START_YEAR if historical_conditioning is None
                                  else historical_conditioning.start_year)
                + period * int(parameters.period_years),
                "psi_child": float(psi),
                "phi": float(phi_values[period]),
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
        if pensions is not None or payroll_taxes is not None:
            rows[-1].update(social_security.fiscal_accounts(evaluation.g_current, parameters))
        if historical_conditioning is not None:
            rows[-1].update(
                historical_conditioning_scope="Conditional historical PF evaluation; supplied terminal boundary",
                historical_outside_flow=float(historical_conditioning.outside_flow),
                historical_retained_due_entry=historical_conditioning.retention * float(due_entry),
                historical_next_bridge_year=historical_conditioning.next_age_targets.get(period + 1),
                historical_bridge_audit=bridge_audit,
                historical_bridge_net_residual=(float(bridge_audit["net_residual"])
                                               if bridge_audit is not None else 0.0),
            )
        state = PFInitialState(
            g_pre=empty_next,
            scheduled_entries=copy_birth_queue(next_queue),
            scheduled_raw_entries=copy_birth_queue(next_raw_queue),
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
