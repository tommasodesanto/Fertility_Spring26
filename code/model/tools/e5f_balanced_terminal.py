"""Actual-population terminal housing/PAYGO adapter for the unrebated baseline.

The caller supplies the selected parameters, original supply law, empirical
demographic primitives, numerical budgets and process watchdog. No inputs are
loaded or calibrated here, and no files or model jobs are created. A successful
receipt certifies this stationary endpoint only, not the historical path,
horizon, empirical target contract or production promotion.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
import math
import time
from typing import Any

import numpy as np

from e5f_matched_pf_endpoint import EndpointControls, evaluate_endpoint
from e5f_social_security import bind_social_security_income, fiscal_accounts
from e5f_social_security_root import solve_social_security_path


@dataclass(frozen=True)
class TerminalAuditControls:
    """Explicit tolerances for the retained household and reconstruction checks."""

    reconstruction_tolerance: float
    feasibility_projection_tolerance: float
    probability_tolerance: float
    occupied_mass_tolerance: float
    value_drop_tolerance: float


@dataclass
class BalancedTerminalEndpoint:
    parameters: Any
    b_grid: np.ndarray
    policy: Any
    endpoint: Any
    social_security: dict[str, Any]
    diagnostics: dict[str, Any]
    gates: dict[str, bool]

    @property
    def fixed_point(self):
        return self.endpoint.fixed_point

    @property
    def mapping_valid(self):
        return bool(self.endpoint.mapping_valid and all(self.gates.values()))


@dataclass
class BalancedTerminalResult:
    """A large endpoint object plus a separate, tensor-free numerical receipt."""

    endpoint: BalancedTerminalEndpoint | None
    root_receipt: dict[str, Any]

    @property
    def production_eligible(self):
        # External source, target, historical and horizon certificates remain
        # the caller's responsibility; this is endpoint numerical eligibility.
        return bool(self.root_receipt['endpoint_production_eligible'])


def _runtime():
    """Lazy model imports: importing/testing the adapter never starts Numba."""
    import run_e5f_matched_pf_smoke as primitive
    transition, calendar = primitive.transition, primitive.calendar
    transition.configure_sequential_model()
    calendar.model = primitive.model
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows
    return primitive.model, calendar, primitive


def _validate_inputs(parameters, b_grid, demographic_primitives, supply_rule,
                     controls, audit_controls, fiscal_tolerance):
    P = parameters
    if (int(P.I) != 1 or float(P.period_years) != 4.
            or not bool(P.scale_flows_to_period)
            or not math.isclose(float(P.tau_H), .04, rel_tol=0., abs_tol=1e-15)
            or float(P.property_tax_lump_sum_transfer) != 0.
            or not math.isclose(float(P.tau_pay), .179, rel_tol=0., abs_tol=1e-15)):
        raise ValueError('Require one market, four-year flows, 1% annual property tax, zero rebate and payroll tax .179')
    if any(bool(getattr(P, name, False)) for name in
           ('joint_nested_choice', 'fertility_nest_choice', 'two_shock_choice')):
        raise ValueError('Balanced baseline adapter requires the approved sequential architecture')
    if not bool(getattr(P, 'exhaustive_saving_control', False)):
        raise ValueError('Approved terminal requires exhaustive_saving_control=True')
    # Defer this import so simply importing the root adapter remains lightweight.
    # The common validator checks the complete approved utility/lifecycle contract.
    from e5f_parenthood_utility import validate_parenthood_utility
    validate_parenthood_utility(P)
    if not math.isfinite(float(P.psi_child)):
        raise ValueError('An explicit finite terminal preference is required')
    grid = np.asarray(b_grid, dtype=float)
    if (grid.ndim != 1 or len(grid) < 2 or int(P.Nb) != len(grid)
            or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0)):
        raise ValueError('Explicit grid must match parameters and be finite and increasing')
    if (supply_rule.mode != 'static-elastic'
            or any(not math.isfinite(float(getattr(supply_rule, name)))
                   or float(getattr(supply_rule, name)) <= 0.
                   for name in ('initial_price', 'initial_stock', 'elasticity'))):
        raise ValueError('The inherited static-elastic supply rule is required')
    if not math.isclose(float(supply_rule.elasticity), .63, rel_tol=0., abs_tol=1e-15):
        raise ValueError('Approved terminal supply elasticity must remain .63')
    parameter_supply = [np.asarray(getattr(P, name), dtype=float)
                        for name in ('H0', 'r_bar', 'xi_supply')]
    if (any(a.shape != (1,) or not np.isfinite(a).all() or np.any(a <= 0.)
            for a in parameter_supply)
            or not np.isfinite(P.user_cost_rate) or P.user_cost_rate <= 0.
            or not math.isclose(float(parameter_supply[2][0]), .63, rel_tol=0., abs_tol=1e-15)
            or not callable(getattr(supply_rule, 'quantity', None))):
        raise ValueError('Approved one-market parameter supply must be finite, positive and have elasticity .63')
    H0, r_bar, elasticity = (float(a[0]) for a in parameter_supply)
    # evaluate_endpoint deliberately passes supply_rule to every calendar
    # evaluation, overriding the P.H0 fallback. Verify that this explicit law is
    # the same calibrated initial curve, not a silently reanchored terminal law.
    for price in (float(supply_rule.initial_price), float(supply_rule.initial_price) * 1.1):
        original = H0 * (float(P.user_cost_rate) * price / r_bar) ** elasticity
        explicit = np.asarray(supply_rule.quantity(np.array([price])), dtype=float)
        if (explicit.shape != (1,) or not np.isfinite(explicit).all()
                or not math.isfinite(original)
                or not math.isclose(float(explicit[0]), original, rel_tol=1e-12, abs_tol=0.)):
            raise ValueError('Explicit terminal supply law differs from the calibrated initial parameter curve')
    if (int(demographic_primitives.start_year) != 2023
            or int(demographic_primitives.last_empirical_year) < 2023):
        raise ValueError('Explicit frozen terminal demographics must retain the 2023 anchor')
    if not isinstance(controls, EndpointControls) or not isinstance(audit_controls, TerminalAuditControls):
        raise ValueError('Explicit endpoint and household audit controls are required')
    if (isinstance(controls.maximum_inner_iterations, bool)
            or int(controls.maximum_inner_iterations) != controls.maximum_inner_iterations
            or controls.maximum_inner_iterations < 1
            or not math.isfinite(float(controls.inner_damping))
            or not 0. < controls.inner_damping <= 1.):
        raise ValueError('Endpoint iteration budget and damping must be admissible')
    for name in ('distribution_tolerance', 'birth_rate_tolerance',
                 'one_step_tolerance', 'fiscal_absolute_tolerance',
                 'accounting_absolute_tolerance'):
        value = float(getattr(controls, name))
        if not math.isfinite(value) or value <= 0.:
            raise ValueError(f'Explicit positive endpoint {name} required')
    if (not math.isfinite(float(controls.market_tolerance))
            or not 0 < controls.market_tolerance <= 2e-4
            or not math.isfinite(float(fiscal_tolerance))
            or not 0 < fiscal_tolerance <= 1e-6):
        raise ValueError('Housing/fiscal tolerances cannot exceed 2e-4/1e-6')
    ceilings = dict(reconstruction_tolerance=5e-9,
        feasibility_projection_tolerance=1e-6, probability_tolerance=1e-12,
        occupied_mass_tolerance=1e-12, value_drop_tolerance=1e-7)
    for name, ceiling in ceilings.items():
        value = float(getattr(audit_controls, name))
        if not math.isfinite(value) or not 0 <= value <= ceiling:
            raise ValueError(f'{name} must be explicit, nonnegative and no looser than {ceiling}')
    return grid


def _household_checks(evaluation, parameters, shared, grid, rent, primitive,
                      audit_controls):
    """Retained dated budget and numerical-audit conventions, without file I/O.

    The wealth screen matches policy_array_audit: an occupied lower node and a
    decrease in pre-choice V between adjacent wealth nodes. Probability bounds
    match the matched-PF primitive check. Calendar evaluation already executes
    the core feasibility/dead-mass gates; its projected mass is checked below.
    """
    budget = primitive.dated_budget(evaluation, parameters, shared, grid, rent)
    arrays = primitive.policy_arrays(evaluation.policy)
    finite = all(np.isfinite(a).all() for a in arrays.values())
    probabilities = {}
    for name, values in arrays.items():
        if 'prob' in name:
            a = np.asarray(values)
            probabilities[name] = dict(minimum=float(np.min(a)),
                maximum=float(np.max(a)), nonfinite=int(np.count_nonzero(~np.isfinite(a))))
    tol = audit_controls.probability_tolerance
    probability_gate = all(row['nonfinite'] == 0 and row['minimum'] >= -tol
        and row['maximum'] <= 1. + tol for row in probabilities.values())
    V = np.asarray(evaluation.policy.V)
    mass = np.asarray(evaluation.g_pre)
    if V.shape != mass.shape:
        raise ValueError('Terminal value and household distribution shapes differ')
    occupied = mass[:-1] > audit_controls.occupied_mass_tolerance
    steps = np.diff(V, axis=0)
    drops = occupied & (steps < -audit_controls.value_drop_tolerance)
    projection = float(evaluation.feasibility_projection_mass)
    distribution_gate = all(np.isfinite(g).all() and np.all(g >= 0.) for g in
        (evaluation.g_pre, evaluation.g_post_fertility, evaluation.g_current))
    diagnostics = dict(budget=budget, probabilities=probabilities,
        occupied_negative_steps=int(drops.sum()),
        maximum_occupied_value_drop=float(max(0., -np.min(steps[occupied]))) if occupied.any() else 0.,
        feasibility_projection_mass=projection)
    gates = dict(household_budget=math.isfinite(float(budget['budget_excess_mass']))
            and budget['budget_excess_mass'] <= 2e-10,
        finite_policy_arrays=bool(finite), probability_bounds=bool(probability_gate),
        occupied_value_monotonicity=not bool(drops.any()),
        finite_nonnegative_household_mass=bool(distribution_gate),
        feasibility_projection=math.isfinite(projection)
            and 0. <= projection <= audit_controls.feasibility_projection_tolerance)
    return diagnostics, gates


def _evaluate_terminal_trial(*, parameters, b_grid, demographic_primitives,
                             supply_rule, controls, audit_controls,
                             asset_price, pension_period):
    """Re-solve households and the actual terminal population for one trial."""
    model, calendar, primitive = _runtime()
    P = copy.deepcopy(parameters)
    # This must precede shared household inputs AND the Bellman solve. Pensions
    # are already four-year units; never call the annual income resolver here.
    bind_social_security_income(P, pension_period=pension_period, payroll_tax=.179)
    shared = model.precompute_shared(P, b_grid)
    price = np.array([asset_price], dtype=float)
    solution = model.solve_markov_income_at_prices(price, P, b_grid, SD=shared)
    policy = calendar.policy_from_solution(solution, price, P, b_grid, shared)
    seed, reconstruction = calendar.reconstruct_stationary_pre_fertility(
        solution, policy, P, b_grid, shared)
    del solution
    endpoint = evaluate_endpoint(parameters=P, b_grid=b_grid, policy=policy,
        asset_price=asset_price, transfer=0., psi_child=float(P.psi_child),
        demographic_primitives=demographic_primitives, initial_g_pre=seed,
        supply_rule=supply_rule, fiscal_regime='fixed_transfer', controls=controls)
    del seed
    current = calendar.evaluate_period(price, endpoint.fixed_point.g_pre, P,
        b_grid, shared, calendar.SolveCounter(), supply_rule=supply_rule,
        supplied_policy=policy)
    diagnostics, gates = _household_checks(current, P, shared, b_grid,
        float(P.user_cost_rate) * asset_price, primitive, audit_controls)
    # Fiscal exposure belongs to actual terminal households, never the analytic
    # seed age distribution and never total resident persons (which include nonheads).
    accounts = fiscal_accounts(current.g_current, P)
    original = np.asarray(endpoint.fixed_point.g_pre)
    before = original.sum(axis=(0, 1, 5, 6))
    after = np.asarray(current.g_current).sum(axis=(0, 1, 5, 6))
    exposure_mass_gap = float(np.max(np.abs(before - after)))
    housing_gap = max(abs(float(current.demand_by_loc[0]) - endpoint.residuals['housing_demand']),
                      abs(float(current.supply_by_loc[0]) - endpoint.residuals['housing_supply']))
    diagnostics.update(reconstruction=dict(reconstruction),
        fiscal_age_income_mass_gap=exposure_mass_gap,
        terminal_housing_recheck_gap=housing_gap,
        property_tax_revenue=endpoint.residuals['tax_revenue'],
        property_tax_transfer_outlays=0.,
        property_tax_disposition='retained unrebated baseline revenue; outside PAYGO budget')
    accounting_tol = controls.accounting_absolute_tolerance
    gates.update(
        seed_reconstruction=all(math.isfinite(float(reconstruction[name]))
            and abs(float(reconstruction[name])) <= audit_controls.reconstruction_tolerance
            for name in ('stationary_post_fertility_nesting_l1', 'stationary_post_fertility_nesting_max_abs')),
        seed_feasibility=math.isfinite(float(reconstruction['stationary_feasibility_projection_mass']))
            and 0. <= float(reconstruction['stationary_feasibility_projection_mass'])
            <= audit_controls.feasibility_projection_tolerance,
        fiscal_population_identity=math.isfinite(exposure_mass_gap)
            and exposure_mass_gap <= accounting_tol * max(1., float(original.sum())),
        terminal_housing_recheck=math.isfinite(housing_gap)
            and housing_gap <= accounting_tol * max(1., abs(endpoint.residuals['housing_supply'])),
        finite_renewal=math.isfinite(endpoint.residuals['renewal_ratio'])
            and 0. <= endpoint.residuals['renewal_ratio'] < 1.)
    return BalancedTerminalEndpoint(P, b_grid, policy, endpoint, accounts, diagnostics, gates)


def solve_balanced_terminal(*, parameters, b_grid, demographic_primitives,
                           supply_rule, controls, audit_controls,
                           start_price, start_pension_period, price_bounds,
                           pension_bounds, fiscal_tolerance, market_slope,
                           fiscal_slope, max_log_step, damping, max_evaluations,
                           deadline_monotonic, max_condition_number,
                           worsening_factor, final_reproduction_tolerance,
                           callback, initial_jacobian=None):
    """Two unknowns: asset price and period pension; fixed payroll tax .179.

    All controls and bounds are required; callback may explicitly be None.
    The generic root reserves a fresh final household/population evaluation
    inside max_evaluations. Its deadline checks occur between evaluations; the
    caller must bound a running evaluation with its own process watchdog.
    A failed/unfinished root may return the best endpoint as a diagnostic,
    always with endpoint_production_eligible=False. No model tensor is placed
    in root payloads, histories or callbacks, nor cached across all trials.
    """
    grid = _validate_inputs(parameters, b_grid, demographic_primitives,
        supply_rule, controls, audit_controls, fiscal_tolerance)
    lo, hi = map(float, price_bounds)
    if (not np.isfinite([lo, hi, start_price]).all()
            or not 0 < lo < hi or not lo <= start_price <= hi):
        raise ValueError('Explicit finite positive price bounds must contain the start')
    if callback is not None and not callable(callback):
        raise ValueError('callback must be callable or explicitly None')
    last, best = None, None
    trial_count = 0

    def evaluate(prices, pensions):
        nonlocal last, trial_count
        if time.monotonic() >= deadline_monotonic:
            raise TimeoutError('Terminal evaluation deadline reached')
        trial_count += 1
        last = None  # Release the rejected last tensor before allocating another.
        last = _evaluate_terminal_trial(parameters=parameters, b_grid=grid,
            demographic_primitives=demographic_primitives, supply_rule=supply_rule,
            controls=controls, audit_controls=audit_controls,
            asset_price=float(prices[0]), pension_period=float(pensions[0]))
        accounts = last.social_security
        # Zero-budget rule: both zero implies residual zero. A positive outlay
        # with zero revenue has residual -1 and cannot pass the fiscal gate.
        revenue, outlays = accounts['payroll_tax_revenue'], accounts['pension_outlays']
        scale = max(abs(revenue), abs(outlays), 1e-12)
        residual = (revenue - outlays) / scale
        if revenue == 0. and outlays != 0.:
            residual = -1.
        return dict(market_residual=[last.endpoint.residuals['housing_relative']],
            fiscal_residual=[residual], mapping_valid=last.mapping_valid,
            payload=dict(trial=trial_count, asset_price=float(prices[0]),
                pension_period=float(pensions[0]), social_security=dict(accounts),
                fiscal_residual_scale=scale, scaled_pension_budget_residual=residual,
                endpoint_residuals=dict(last.endpoint.residuals),
                endpoint_gates=dict(last.endpoint.gates), household_gates=dict(last.gates),
                diagnostics=last.diagnostics))

    def progress(record):
        nonlocal best
        if record.get('new_best'):
            best = last
        if callback is not None:
            callback(record)

    receipt = solve_social_security_path(closure='fixed_tax',
        initial_prices=[start_price], initial_fiscal_values=[start_pension_period],
        evaluate=evaluate, project_prices=lambda prices: np.clip(prices, lo, hi),
        fiscal_bounds=pension_bounds, market_tolerance=controls.market_tolerance,
        fiscal_tolerance=fiscal_tolerance, market_slope=market_slope,
        fiscal_slope=fiscal_slope, max_log_step=max_log_step, damping=damping,
        max_evaluations=max_evaluations, deadline_monotonic=deadline_monotonic,
        max_condition_number=max_condition_number, worsening_factor=worsening_factor,
        final_reproduction_tolerance=final_reproduction_tolerance,
        callback=progress, initial_jacobian=initial_jacobian)
    final = receipt['final']
    matched_final = bool(final is not None and last is not None
        and final['payload']['trial'] == trial_count
        and float(final['prices'][0]) == float(last.policy.price[0])
        and float(final['fiscal_values'][0]) == float(last.parameters.pension))
    eligible = bool(receipt['converged'] and matched_final and last.mapping_valid
                    and all(receipt['gates'].values()))
    receipt.update(schema='e5f_balanced_terminal_v1',
        endpoint_production_eligible=eligible, fresh_endpoint_matches_final=matched_final,
        returned_endpoint='fresh_final' if matched_final else 'best_diagnostic',
        endpoint_scope='stationary terminal only; external contracts and path/horizon certification required',
        fiscal_accounting_population='actual terminal post-choice household heads',
        fiscal_residual_definition='(payroll revenue - pension outlays)/max(revenue,outlays,1e-12)',
        fiscal_zero_budget_rule='both zero: zero residual; zero revenue with positive outlays: residual -1',
        payroll_tax=.179, annual_property_tax=.01, period_property_tax=.04,
        property_tax_rebate=0., price_bounds=(lo, hi),
        endpoint_controls=vars(controls).copy(), audit_controls=vars(audit_controls).copy())
    return BalancedTerminalResult(last if matched_final else best, receipt)
