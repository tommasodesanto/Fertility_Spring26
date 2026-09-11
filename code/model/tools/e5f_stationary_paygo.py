"""Balanced initial pensions for the one-market exogenous earnings process.

This shortcut applies only to the pre-announcement stationary economy. Dated
person/household demographic equilibria must use their actual distributions.
Every solved stationary result is checked against the proposed age/earnings
marginal, so a change in the forward operator cannot silently invalidate it.
"""
from __future__ import annotations

import copy
import numpy as np

from e5f_social_security import bind_social_security_income, fiscal_accounts


def stationary_age_income_mass(P):
    """Return normalized m[location, age, earnings] from entrant transitions.

The source forward operator enters households with z_weights and propagates
the row-stochastic Pi_z, multiplied by age survival. Fertility, tenure, wealth
and child aging redistribute mass within these cells. One location rules out
endogenous spatial sorting of the earnings tax base.
"""
    if int(P.I) != 1 or int(P.J) < 2:
        raise ValueError('Analytic initial pension requires one market and at least two ages')
    if str(getattr(P, 'income_type_transition', '')).lower() not in {'markov', 'stochastic', 'persistent'}:
        raise ValueError('Explicit exogenous Markov income process required')
    z = np.asarray(P.z_grid, dtype=float)
    weights = np.asarray(P.z_weights, dtype=float)
    transition = np.asarray(P.Pi_z, dtype=float)
    if (z.ndim != 1 or not len(z) or weights.shape != z.shape
            or transition.shape != (len(z), len(z))
            or not all(np.isfinite(a).all() for a in (z, weights, transition))
            or np.any(z <= 0) or np.any(weights < 0) or weights.sum() <= 0
            or np.any(transition < 0)
            or not np.allclose(transition.sum(axis=1), 1., rtol=0, atol=1e-12)):
        raise ValueError('Invalid exogenous entrant/earnings transition inputs')
    survival = (np.asarray(P.survival_probs, dtype=float)[:int(P.J)-1]
                if bool(getattr(P, 'use_age_survival', False)) else np.ones(int(P.J)-1))
    if (survival.shape != (int(P.J)-1,) or not np.isfinite(survival).all()
            or np.any(survival < 0) or np.any(survival > 1)):
        raise ValueError('Invalid age-only survival schedule')
    mass = np.zeros((1, int(P.J), len(z)))
    mass[0, 0] = weights / weights.sum()
    transition = transition / transition.sum(axis=1)[:, None]
    for age in range(1, int(P.J)):
        mass[0, age] = survival[age-1] * (mass[0, age-1] @ transition)
    return mass / mass.sum()


def bind_initial_balanced_pension(P, *, payroll_tax):
    """Return a copied parameter object with the stationary balanced benefit."""
    out = copy.deepcopy(P)
    mass = stationary_age_income_mass(out)
    # singleton wealth, tenure, parity and dependent-child axes preserve the
    # actual-ledger API without inventing reference income-state weights.
    g = mass[None, None, :, :, :, None, None]
    bind_social_security_income(out, payroll_tax=payroll_tax)
    accounts = fiscal_accounts(g, out)
    pension = accounts['implied_balanced_pension_period']
    if pension is None or not np.isfinite(pension) or pension <= 0:
        raise ValueError('Initial stationary budget requires positive retiree exposure and revenue')
    bind_social_security_income(out, pension_period=pension, payroll_tax=payroll_tax)
    return out, dict(method='exogenous_one_market_age_income_recursion',
                     predicted_accounts=fiscal_accounts(g, out))


def certify_initial_pension(g, P, *, marginal_tolerance, fiscal_tolerance):
    """Verify the shortcut on the actual solved household distribution."""
    if (not np.isfinite([marginal_tolerance, fiscal_tolerance]).all()
            or not 0 < marginal_tolerance <= 1e-8 or not 0 < fiscal_tolerance <= 1e-6):
        raise ValueError('Explicit tight marginal and fiscal gates required')
    accounts = fiscal_accounts(g, P)
    actual = np.asarray(g, dtype=float).sum(axis=(0, 1, 5, 6))
    if actual.sum() <= 0:
        raise ValueError('Actual stationary household mass must be positive')
    gap = float(np.max(np.abs(actual / actual.sum() - stationary_age_income_mass(P))))
    fiscal_gap = abs(accounts['scaled_pension_budget_residual'])
    if gap > marginal_tolerance or fiscal_gap > fiscal_tolerance:
        raise RuntimeError(f'Initial pension certification failed: marginal={gap}, fiscal={fiscal_gap}')
    return dict(actual_accounts=accounts, normalized_age_income_max_gap=gap,
                marginal_tolerance=marginal_tolerance, fiscal_tolerance=fiscal_tolerance,
                marginal_gate=True, fiscal_gate=True)


def solve_balanced_initial_equilibrium(*, model, parameters, b_grid, initial_prices,
                                       payroll_tax, marginal_tolerance, fiscal_tolerance):
    """Solve households/prices with the balancing pension already anticipated."""
    P, receipt = bind_initial_balanced_pension(parameters, payroll_tax=payroll_tax)
    sol, P, prices = model.solve_markov_income_equilibrium(
        np.asarray(initial_prices, dtype=float), P, b_grid, verbose=False)
    if not bool(sol.converged) or not bool(sol.timings.get('strict_converged', False)):
        raise RuntimeError('Initial housing equilibrium failed its unchanged strict gate')
    receipt.update(certify_initial_pension(sol.g, P,
        marginal_tolerance=marginal_tolerance, fiscal_tolerance=fiscal_tolerance))
    return sol, P, prices, receipt


def rebase_initial_supply(P, *, asset_prices, elasticity):
    """Change elasticity while preserving the supplied old price/stock point.

The model uses H=H0*(user_cost_rate*p/r_bar)**xi_supply. This is a
one-time seed mapping; candidate H0 subsequently varies in calibration.
"""
    out = copy.deepcopy(P)
    price = np.asarray(asset_prices, dtype=float)
    old_level = np.asarray(P.H0, dtype=float)
    old_elasticity = np.asarray(P.xi_supply, dtype=float)
    reference = np.asarray(P.r_bar, dtype=float)
    if (any(a.shape != (int(P.I),) for a in (price, old_level, old_elasticity, reference))
            or not all(np.isfinite(a).all() for a in (price, old_level, old_elasticity, reference))
            or any(np.any(a <= 0) for a in (price, old_level, old_elasticity, reference))
            or not np.isfinite(elasticity) or elasticity <= 0
            or not np.isfinite(P.user_cost_rate) or P.user_cost_rate <= 0):
        raise ValueError('Invalid supply-rebase inputs')
    ratio = float(P.user_cost_rate) * price / reference
    old_stock = old_level * ratio ** old_elasticity
    out.xi_supply = np.full(int(P.I), float(elasticity))
    out.H0 = old_stock / ratio ** out.xi_supply
    if not np.isfinite(out.H0).all() or np.any(out.H0 <= 0):
        raise ValueError('Supply rebase produced an invalid level')
    return out, dict(old_elasticity=old_elasticity.tolist(), new_elasticity=out.xi_supply.tolist(),
                     old_H0=old_level.tolist(), rebased_H0=out.H0.tolist(),
                     anchor_asset_prices=price.tolist(), anchor_housing_stock=old_stock.tolist())
