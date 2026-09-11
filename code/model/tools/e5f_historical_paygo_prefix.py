"""Diagnostic PAYGO predictions for the five observed-age stocks, 2007--2023.

No empirical files, household solves or root changes. Economically, within-age
earnings composition follows the exogenous Markov matrix; observed age totals
replace the cohort totals afterward. Stored probability rounding, pruning and
cohort mass corrections can prevent bitwise equality with the actual forward
operator. Actual dated ledgers and marginal checks remain necessary.
"""
from __future__ import annotations

import copy
import numpy as np

from e5f_social_security import bind_social_security_income, fiscal_accounts
from e5f_stationary_paygo import stationary_age_income_mass


YEARS = (2007, 2011, 2015, 2019, 2023)
SCHEMA = 'e5f_historical_paygo_prefix_v1'


def _years(years):
    if not np.array_equal(np.asarray(years), np.asarray(YEARS)):
        raise ValueError('Supply exactly the five observed dates 2007 through 2023')


def _ledger(marginal, parameters, pension):
    # Binding replaces attributes, not shared primitive arrays. Avoid copying
    # large cached policy arrays that a solved parameter object may carry.
    dated = copy.copy(parameters)
    bind_social_security_income(dated, pension_period=pension, payroll_tax=.179)
    # wealth, tenure, location, age, income, parity, dependent-child axes
    mass = np.asarray(marginal)[None, None, None, :, :, None, None]
    return fiscal_accounts(mass, dated)


def predict_historical_paygo_prefix(*, initial_age_income_mass, parameters,
                                   observed_age_masses, years):
    """Return five diagnostic marginals and period pensions.

    `initial_age_income_mass` has shape (J, Nz), after ACS2007 reweighting.
    `observed_age_masses` has shape (5, J) and includes its 2007 age totals.
    Synthetic J/Nz are allowed; production callers own full-grid/source pins.
    All nonnegative mass rows must be positive. No post-2023 extrapolation.
    """
    _years(years)
    P = parameters
    if (int(P.I) != 1 or float(P.period_years) != 4.
            or not bool(P.scale_flows_to_period) or float(P.tau_pay) != .179
            or not bool(getattr(P, 'sequential_births', False))
            or any(bool(getattr(P, name, False)) for name in
                   ('joint_nested_choice', 'fertility_nest_choice', 'two_shock_choice'))):
        raise ValueError('Retain one-market sequential four-year income and fixed payroll tax .179')
    # Reuse the stationary helper's source-consistent Markov/entry/survival
    # validation; its stationary distribution is NOT used as a dated marginal.
    stationary_age_income_mass(P)
    if bool(getattr(P, 'use_age_survival', False)) and np.any(
            np.asarray(P.survival_probs, dtype=float)[:int(P.J)-1] <= 0.):
        raise ValueError('Positive observed age stocks require positive source survival')
    z = np.asarray(P.z_grid, dtype=float)
    transition = np.asarray(P.Pi_z, dtype=float)
    transition = transition / transition.sum(axis=1)[:, None]
    entry = np.asarray(P.z_weights, dtype=float)
    entry = entry / entry.sum()
    initial = np.asarray(initial_age_income_mass, dtype=float)
    age_totals = np.asarray(observed_age_masses, dtype=float)
    if (initial.shape != (int(P.J), len(z))
            or age_totals.shape != (5, int(P.J))
            or not np.isfinite(initial).all() or np.any(initial < 0.)
            or not np.isfinite(age_totals).all() or np.any(age_totals <= 0.)):
        raise ValueError('Finite nonnegative (J,Nz) marginal and positive (5,J) age totals required')
    initial_totals = initial.sum(axis=1)
    if (not np.isfinite(initial_totals).all() or np.any(initial_totals <= 0.)
            or not np.isfinite(age_totals.sum(axis=1)).all()
            or not np.allclose(initial_totals, age_totals[0], rtol=1e-10, atol=0.)):
        raise ValueError('Initial age/income marginal must match the supplied 2007 age totals')
    marginals = np.empty((5, int(P.J), len(z)))
    marginals[0] = initial.copy()
    conditional = initial / initial_totals[:, None]
    for date in range(1, 5):
        following = np.empty_like(conditional)
        following[0] = entry
        following[1:] = conditional[:-1] @ transition
        marginals[date] = age_totals[date, :, None] * following
        conditional = following
    pensions, ledgers = [], []
    for marginal in marginals:
        trial = _ledger(marginal, P, 0.)
        benefit = trial['implied_balanced_pension_period']
        if benefit is None or not np.isfinite(benefit) or benefit <= 0.:
            raise ValueError('Every date requires positive payroll revenue and retiree exposure')
        pensions.append(float(benefit))
        ledgers.append(_ledger(marginal, P, benefit))
    return dict(schema=SCHEMA, years=list(YEARS),
        age_income_marginals=marginals, pensions_period=np.asarray(pensions),
        predicted_ledgers=ledgers,
        ledger_factor_receipt=dict(payroll_tax=.179, period_years=4.,
            first_retired_age_index=int(P.J_R), wage=np.asarray(P.w_hat).tolist(),
            income_age_profile=np.asarray(P.income_age_profile).tolist(),
            income_states=z.tolist(), normalized_entry_weights=entry.tolist(),
            row_stochastic_income_transition=transition.tolist(),
            retirement_income_z_scale=float(getattr(P, 'retirement_income_z_scale', 0.)),
            observed_age_masses=age_totals.tolist(),
            fiscal_definition='Existing fiscal_accounts: gross period payroll and retirement benefit exposure'),
        status='diagnostic_prediction_not_compared_with_actual_path',
        actual_marginals_verified=False, historical_equilibrium_certified=False,
        root_coordinates_reduced=False,
        numerical_caveat='Economic recursion; float32 tenure probabilities, pruning and cohort corrections require actual dated marginal and fiscal verification')


def compare_historical_paygo_prefix(*, prediction, actual_ledgers, actual_years,
                                   fiscal_tolerance=1e-6):
    """Compare with the first five actual path ledgers, without promoting a path.

    Compare implied balanced benefits and both ledger factors, even when the
    actual path used different trial pensions. Ledger-only agreement does not
    verify the full age/income marginal or establish fiscal path equilibrium.
    """
    _years(actual_years)
    if (prediction.get('schema') != SCHEMA
            or not np.array_equal(prediction.get('years'), np.asarray(YEARS))
            or not np.isfinite(fiscal_tolerance) or not 0. < fiscal_tolerance <= 1e-6):
        raise ValueError('Recognized five-date prediction and fiscal tolerance no looser than 1e-6 required')
    benefits = np.asarray(prediction['pensions_period'], dtype=float)
    ledgers = list(actual_ledgers)
    predicted = prediction['predicted_ledgers']
    if (benefits.shape != (5,) or not np.isfinite(benefits).all()
            or np.any(benefits <= 0.) or len(ledgers) != 5 or len(predicted) != 5):
        raise ValueError('Exactly five positive predictions and actual dated ledgers required')
    rows = []
    for year, benefit, actual, expected in zip(YEARS, benefits, ledgers, predicted):
        names = ('payroll_tax_base_period', 'retiree_benefit_exposure',
                 'payroll_tax_revenue', 'pension_outlays', 'pension_period_units', 'payroll_tax_rate')
        values = [float(actual[name]) for name in names]
        base, exposure, revenue, outlays, actual_benefit, tax = values
        if (not np.isfinite(values).all() or min(base, exposure, revenue) <= 0.
                or min(outlays, actual_benefit) < 0. or tax != .179
                or not np.isclose(revenue, tax*base, rtol=1e-12, atol=0.)
                or not np.isclose(outlays, actual_benefit*exposure, rtol=1e-12, atol=0.)):
            raise ValueError('Actual ledger has invalid fixed-tax factors or internally inconsistent accounts')
        predicted_outlays = float(benefit) * exposure
        residual = (revenue-predicted_outlays)/max(revenue, predicted_outlays)
        factor_gaps = {}
        for name, value in (('payroll_tax_base_period', base), ('retiree_benefit_exposure', exposure)):
            target = float(expected[name])
            if not np.isfinite(target) or target <= 0.:
                raise ValueError('Prediction has invalid positive ledger factors')
            factor_gaps[name] = abs(value-target)/max(value,target)
        actual_residual = (revenue-outlays)/max(revenue,outlays)
        rows.append(dict(year=year, predicted_pension_period=float(benefit),
            actual_implied_balanced_pension_period=revenue/exposure,
            actual_trial_pension_period=actual_benefit,
            predicted_pension_scaled_budget_residual_on_actual=residual,
            actual_path_scaled_budget_residual=actual_residual,
            relative_ledger_factor_gaps=factor_gaps,
            pension_prediction_pass=abs(residual) <= fiscal_tolerance,
            ledger_factors_pass=max(factor_gaps.values()) <= fiscal_tolerance,
            actual_trial_budget_pass=abs(actual_residual) <= fiscal_tolerance))
    passed = all(row['pension_prediction_pass'] and row['ledger_factors_pass'] for row in rows)
    return dict(schema=SCHEMA+'_comparison', years=list(YEARS), rows=rows,
        fiscal_tolerance=float(fiscal_tolerance), all_comparisons_pass=passed,
        actual_trial_budgets_pass=all(row['actual_trial_budget_pass'] for row in rows),
        status='passed_diagnostic_ledger_comparison' if passed else 'failed_diagnostic_ledger_comparison',
        actual_marginals_verified=False, historical_equilibrium_certified=False,
        root_coordinates_reduced=False)
