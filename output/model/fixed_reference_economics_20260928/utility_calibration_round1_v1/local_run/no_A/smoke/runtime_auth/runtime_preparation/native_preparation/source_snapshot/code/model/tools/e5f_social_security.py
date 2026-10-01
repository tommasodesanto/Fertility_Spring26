"""Explicit period-unit Social Security income and household-budget accounting.

The fiscal rule is chosen by the caller. These functions neither infer an
equilibrium nor use the reference stationary age distribution. Property-tax
rebates and resident non-head persons are outside this budget.
"""
from __future__ import annotations

import numpy as np


def _income_primitives(P):
    I, J, JR = int(P.I), int(P.J), int(P.J_R)
    if I < 1 or J < 1 or not 0 <= JR <= J:
        raise ValueError('Invalid location, age or retirement dimensions')
    wage = np.asarray(P.w_hat, dtype=float)
    profile = np.asarray(P.income_age_profile, dtype=float)
    if (wage.shape != (I,) or profile.shape != (J,)
            or not np.isfinite(wage).all() or np.any(wage <= 0)
            or not np.isfinite(profile).all() or np.any(profile < 0)):
        raise ValueError('Invalid gross wage or income age profile')
    scale = (float(getattr(P, 'period_years', getattr(P, 'da', 1.)))
             if bool(getattr(P, 'scale_flows_to_period', False)) else 1.)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Income period scale must be finite and positive')
    return I, J, JR, wage, profile, scale


def bind_social_security_income(P, *, pension_period=None, payroll_tax=None):
    """Bind dated income in place, once in period units; return the same object.

    Call after parameter overrides and before precomputing household inputs.
    No annual-income resolver is called, so repeated binding is idempotent.
    """
    I, J, JR, wage, profile, scale = _income_primitives(P)
    pension = float(P.pension if pension_period is None else pension_period)
    tax = float(P.tau_pay if payroll_tax is None else payroll_tax)
    if not np.isfinite(pension) or pension < 0:
        raise ValueError('Period pension must be finite and nonnegative')
    if not np.isfinite(tax) or not 0 <= tax < 1:
        raise ValueError('Payroll tax must be finite and in [0, 1)')
    income = np.full((I, J), pension, dtype=float)
    income[:, :JR] = scale * (1. - tax) * wage[:, None] * profile[None, :JR]
    P.tau_pay = tax
    P.pension = pension
    P.pension_by_loc = np.full(I, pension)
    P.income = income
    P.social_security_income_units = 'period'
    return P


def validated_fiscal_paths(periods, pension_path=None, payroll_tax_path=None):
    """Validate explicit paths; None preserves the supplied parameter income."""
    if type(periods) is not int or periods < 1:
        raise ValueError('Fiscal paths require a positive integer date count')
    paths = []
    for name, raw in (('pension', pension_path), ('payroll tax', payroll_tax_path)):
        if raw is None:
            paths.append(None)
            continue
        path = np.asarray(raw, dtype=float)
        if (path.shape != (periods,) or not np.isfinite(path).all()
                or np.any(path < 0) or (name == 'payroll tax' and np.any(path >= 1))):
            raise ValueError(f'Invalid {name} path: require one admissible value per date')
        paths.append(path.copy())
    return tuple(paths)


def apply_fiscal_date(P, period, pensions, taxes):
    """Use the same binding in backward values and forward household choices."""
    if pensions is not None or taxes is not None:
        bind_social_security_income(P,
            pension_period=None if pensions is None else pensions[period],
            payroll_tax=None if taxes is None else taxes[period])
    return P


def fiscal_accounts(g, P):
    """Integrate gross payroll and pension exposure over actual current heads.

    Axes are wealth, tenure, location, age, income, parity and dependent children.
    Undefined implied benefit/tax ratios are None, never epsilon-denominated.
    The residual may be nonzero in a trial; acceptance is the caller's job.
    """
    I, J, JR, wage, profile, scale = _income_primitives(P)
    mass = np.asarray(g, dtype=float)
    z = np.asarray(P.z_grid, dtype=float)
    if z.ndim != 1 or not len(z) or not np.isfinite(z).all() or np.any(z <= 0):
        raise ValueError('Income states must be finite and positive')
    if (mass.ndim != 7 or mass.shape[2:5] != (I, J, len(z))
            or any(n == 0 for n in mass.shape)
            or not np.isfinite(mass).all() or np.any(mass < 0)):
        raise ValueError('Expected finite nonnegative seven-dimensional household mass')
    tax, pension = float(P.tau_pay), float(P.pension)
    if not np.isfinite(tax) or not 0 <= tax < 1:
        raise ValueError('Payroll tax must be finite and in [0, 1)')
    if not np.isfinite(pension) or pension < 0:
        raise ValueError('Period pension must be finite and nonnegative')
    retirement_scale = float(getattr(P, 'retirement_income_z_scale', 0.))
    if not np.isfinite(retirement_scale) or retirement_scale < 0:
        raise ValueError('Retirement income scale must be finite and nonnegative')
    multiplier = 1. + retirement_scale * (z - 1.)
    if not np.isfinite(multiplier).all() or np.any(multiplier <= 0):
        raise ValueError('Retirement income multipliers must be finite and positive')
    age_income_mass = mass.sum(axis=(0, 1, 5, 6))
    payroll_base = float(scale * np.sum(age_income_mass[:, :JR, :]
        * wage[:, None, None] * profile[None, :JR, None] * z[None, None, :]))
    exposure = float(np.sum(age_income_mass[:, JR:, :] * multiplier[None, None, :]))
    revenue, outlays = tax * payroll_base, pension * exposure
    if not np.isfinite([payroll_base, exposure, revenue, outlays]).all():
        raise ValueError('Nonfinite aggregate payroll or pension accounting')
    residual = revenue - outlays
    magnitude = max(abs(revenue), abs(outlays))
    return dict(
        payroll_tax_base_period=payroll_base, retiree_benefit_exposure=exposure,
        worker_household_mass=float(age_income_mass[:, :JR, :].sum()),
        retiree_household_mass=float(age_income_mass[:, JR:, :].sum()),
        payroll_tax_rate=tax, pension_period_units=pension,
        payroll_tax_revenue=revenue, pension_outlays=outlays,
        pension_budget_residual=residual,
        scaled_pension_budget_residual=residual / magnitude if magnitude > 0 else 0.,
        implied_balanced_pension_period=revenue / exposure if exposure > 0 else None,
        implied_balanced_payroll_tax=outlays / payroll_base if payroll_base > 0 else None,
    )
