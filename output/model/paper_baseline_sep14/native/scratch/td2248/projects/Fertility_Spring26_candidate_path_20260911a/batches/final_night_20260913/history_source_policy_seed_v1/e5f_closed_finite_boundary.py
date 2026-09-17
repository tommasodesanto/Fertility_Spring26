"""Provisional finite-horizon household boundary at actual carried population.

This is a dated snapshot, NOT a stationary population equilibrium. At the
boundary, households value their finite remaining lifetimes at constant prices,
pensions and equal rebates. The caller jointly solves the boundary coordinates
and the interior path so the boundary's actual distribution clears housing and
both fiscal budgets. Constant conditions after this snapshot are a numerical
truncation: their future fiscal feasibility is NOT established by this adapter.

No population transition, entry, migration, stationary KFE, population scaling,
or demographic primitive is supplied here. The caller must carry its zero-
migration state to this boundary and certify that demographic path separately.
The native finite-life Bellman recursion retains its original death/bequest
condition. Never pass a fabricated stationary endpoint to the old validators.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
import time
from typing import Any

import numpy as np


@dataclass
class FiniteBoundaryEvaluation:
    parameters: Any
    b_grid: np.ndarray
    g_pre: np.ndarray
    policy: Any
    residuals: dict
    actual_accounts: dict
    diagnostics: dict
    gates: dict
    boundary_kind: str = 'closed_finite_constant_conditions_snapshot_v1'
    horizon_verified: bool = False
    production_eligible: bool = False

    @property
    def V(self):
        return self.policy.V

    @property
    def mapping_valid(self):
        # Valid evaluation does not mean its three equilibrium residuals clear.
        return bool(all(self.gates.values()))


@dataclass
class BoundaryPolicy:
    """Lifetime value at boundary conditions; no population or equilibrium claim."""
    parameters: Any
    b_grid: np.ndarray
    policy: Any


def _runtime():
    import e5f_balanced_terminal as balanced
    import e5f_social_security as social
    import run_e5f_perfect_foresight_transition as pf
    model, calendar, primitive = balanced._runtime()
    return model, calendar, primitive, balanced, social, pf


def _scaled_difference(revenue, outlays):
    magnitude = max(abs(revenue), abs(outlays))
    return (revenue - outlays) / magnitude if magnitude else 0.


def boundary_policy(*, parameters, grid, price, pension, transfer,
                    deadline_monotonic=None, callback=None):
    """Construct a lifetime policy without evaluating an arbitrary population.

    Every actual population must subsequently pass its dated household and
    fiscal/market checks. This object supplies only the backward value boundary.
    """
    def progress(phase):
        if deadline_monotonic is not None and (
                not np.isfinite(deadline_monotonic) or time.monotonic() >= deadline_monotonic):
            raise TimeoutError('Finite-boundary policy deadline exhausted')
        if callback is not None:callback(dict(phase=phase, boundary_kind='lifetime_policy_only'))
    progress('boundary_policy_validate')
    P=parameters;b_grid=np.asarray(grid,dtype=float)
    if (b_grid.ndim!=1 or len(b_grid)<2 or not np.isfinite(b_grid).all()
            or np.any(np.diff(b_grid)<=0) or int(P.Nb)!=len(b_grid) or int(P.I)!=1):
        raise ValueError('Boundary policy requires the retained one-market wealth grid')
    if not np.isfinite([price,pension,transfer]).all() or price<=0 or pension<0 or transfer<0:
        raise ValueError('Positive finite price and nonnegative fiscal values required')
    if (not bool(P.exhaustive_saving_control) or any(bool(getattr(P,n,False))
            for n in ('joint_nested_choice','fertility_nest_choice','two_shock_choice'))):
        raise ValueError('Retained exhaustive sequential household specification required')
    model,_,_,_,social,pf=_runtime()
    P=copy.deepcopy(P);P.property_tax_lump_sum_transfer=float(transfer)
    social.bind_social_security_income(P,pension_period=float(pension),payroll_tax=float(P.tau_pay))
    prices=np.array([float(price)])
    rent=float(pf.rents_from_asset_prices(prices,float(price),P)[0])
    shared=model.precompute_shared(P,b_grid)
    progress('boundary_bellman')
    objects=model.solve_bellman_full_markov_income(np.array([rent]),prices,P,b_grid,shared,continuation_V=None)
    policy=pf.policy_from_objects(objects,float(price),P,b_grid,shared)
    progress('boundary_policy_ready')
    return BoundaryPolicy(P,b_grid.copy(),policy)


def boundary_evaluation(*, parameters, g_pre, grid, supply_rule, price,
                        pension, transfer, audit_controls=None,
                        deadline_monotonic=None, callback=None):
    """Evaluate one boundary price/pension/rebate trial on the actual state.

    All fiscal values use model-period units. The caller owns the joint root,
    exact final replay, checkpoints and a process watchdog (a deadline cannot
    interrupt a native Bellman call). ``callback`` receives small phase records.
    This function never clears markets, iterates population or writes files.
    """
    def progress(phase):
        if deadline_monotonic is not None:
            if not np.isfinite(deadline_monotonic) or time.monotonic() >= deadline_monotonic:
                raise TimeoutError('Finite-boundary evaluation deadline exhausted')
        if callback is not None:
            callback(dict(phase=phase, boundary_kind='closed_finite_constant_conditions_snapshot_v1'))

    progress('boundary_validate')
    g = np.asarray(g_pre, dtype=float)
    b_grid = np.asarray(grid, dtype=float)
    P = parameters
    if (g.ndim != 7 or not np.isfinite(g).all() or np.any(g < 0.)
            or not np.isfinite(g.sum()) or g.sum() <= 0.):
        raise ValueError('Actual positive finite nonnegative carried household state required')
    if (b_grid.ndim != 1 or len(b_grid) < 2 or not np.isfinite(b_grid).all()
            or np.any(np.diff(b_grid) <= 0.) or g.shape[0] != len(b_grid)
            or int(P.Nb) != len(b_grid) or g.shape[2] != int(P.I)
            or int(P.I) != 1 or g.shape[3] != int(P.J)
            or g.shape[4] != len(P.z_grid)):
        raise ValueError('Boundary grid and age/income/state dimensions must match')
    if not np.isfinite([price, pension, transfer]).all() or price <= 0. or pension < 0. or transfer < 0.:
        raise ValueError('Positive finite price and nonnegative period fiscal values required')
    if supply_rule is None:
        raise ValueError('Explicit inherited supply rule required; no boundary reanchoring')
    if (not bool(P.exhaustive_saving_control) or any(bool(getattr(P, n, False))
            for n in ('joint_nested_choice', 'fertility_nest_choice', 'two_shock_choice'))):
        raise ValueError('Retained exhaustive sequential household specification required')
    model, calendar, primitive, balanced, social, pf = _runtime()
    if audit_controls is None:
        audit_controls = balanced.TerminalAuditControls(5e-9, 1e-6, 1e-12, 1e-12, 1e-7)
    for name, cap in dict(reconstruction_tolerance=5e-9,
            feasibility_projection_tolerance=1e-6, probability_tolerance=1e-12,
            occupied_mass_tolerance=1e-12, value_drop_tolerance=1e-7).items():
        value = float(getattr(audit_controls, name))
        if not np.isfinite(value) or not 0. <= value <= cap:
            raise ValueError('Boundary audit tolerance may not be relaxed: ' + name)
    def policy_progress(record):
        # Keep the population-evaluation callback contract unchanged.
        if record['phase'] == 'boundary_bellman':
            progress('boundary_bellman')
    template=boundary_policy(parameters=P,grid=b_grid,price=price,pension=pension,
        transfer=transfer,deadline_monotonic=deadline_monotonic,callback=policy_progress)
    P,policy=template.parameters,template.policy
    prices = np.array([float(price)])
    # This also checks tax/user-cost consistency and positive rents.
    rent = float(pf.rents_from_asset_prices(prices, float(price), P)[0])
    shared = model.precompute_shared(P, b_grid)
    progress('boundary_actual_population')
    actual = calendar.evaluate_period(prices, g.copy(), P, b_grid, shared,
        calendar.SolveCounter(), supply_rule=supply_rule, supplied_policy=policy)
    diagnostics, gates = balanced._household_checks(actual, P, shared, b_grid,
        rent, primitive, audit_controls)
    accounts = social.fiscal_accounts(actual.g_current, P)
    heads = float(np.sum(actual.g_current))
    demand = float(actual.demand_by_loc[0])
    supply = float(actual.supply_by_loc[0])
    revenue = float(model.property_tax_revenue_from_distribution(
        actual.g_current, actual.policy.hR_pol, actual.policy.price, P))
    outlays = float(transfer) * heads
    if not np.isfinite([heads, demand, supply, revenue, outlays]).all() or heads <= 0. or supply <= 0.:
        raise ValueError('Boundary actual accounting requires finite positive heads and supply')
    mass_gap = abs(heads - float(g.sum()))
    gates = dict(gates, actual_head_mass_preserved=mass_gap <= 2e-9 * max(1., float(g.sum())))
    accounts = dict(accounts, household_heads=heads, property_tax_revenue=revenue,
        equal_transfer_period_units=float(transfer), equal_transfer_outlays=outlays,
        property_tax_budget_residual=revenue-outlays,
        implied_equal_transfer_period=revenue / heads)
    residuals = dict(housing_relative=(demand-supply)/supply,
        pension_relative=float(accounts['scaled_pension_budget_residual']),
        rebate_relative=_scaled_difference(revenue, outlays),
        housing_demand=demand, housing_supply=supply,
        housing_absolute=demand-supply,
        pension_absolute=float(accounts['pension_budget_residual']),
        rebate_absolute=revenue-outlays)
    diagnostics = dict(diagnostics, actual_head_mass_gap=mass_gap,
        household_lifetime_age_cells=int(P.J), constant_conditions_rent=rent,
        stationary_population_computed=False, demographic_transition_computed=False,
        future_fiscal_consistency_verified=False,
        boundary_interpretation='Constant-conditions finite-lifetime household value; actual dated boundary budgets only',
        horizon_status='unverified_finite_truncation')
    progress('boundary_complete')
    return FiniteBoundaryEvaluation(P, b_grid.copy(), g.copy(), actual.policy,
        residuals, accounts, diagnostics, gates)
