"""Two-block acceleration helpers for the four-successive-surprise transition.

The transition's root coordinates are log house prices and log period pensions.
Its residuals remain physical housing imbalance and unscaled PAYGO imbalance;
the Social Security root is solely responsible for its one-time internal fiscal
row scaling and for certifying both physical gates on its fresh final replay.
"""
from __future__ import annotations

from contextlib import contextmanager

import numpy as np

import e5f_social_security_root as social_security_root
from e5f_ssj_scaled_step_root import solve_price_path_scaled
from e5f_ssj_toeplitz_jacobian import toeplitz_block


UNKNOWN_BLOCKS = ("log_house_price", "log_period_pension")
RESIDUAL_BLOCKS = ("housing_imbalance", "pension_imbalance")


def _column(column, horizon):
    values = np.asarray(column, dtype=float)
    if values.shape != (2 * horizon,) or not np.isfinite(values).all():
        raise ValueError("each physical derivative column must be finite with shape (2T,)")
    return values


def assemble_jacobian(columns, horizon, perturbed_date):
    """Build a physical ``(2T, 2T)`` block-Toeplitz Jacobian.

    ``columns`` contains the two central-log-difference derivative columns in
    ``UNKNOWN_BLOCKS`` order.  Their rows are ordered by ``RESIDUAL_BLOCKS``.
    The supplied columns are already derivatives, not residual evaluations;
    this helper deliberately constructs no fake-news derivatives.
    """
    if type(horizon) is not int or horizon <= 0:
        raise ValueError("horizon must be a positive integer")
    if type(perturbed_date) is not int or not 0 <= perturbed_date < horizon:
        raise ValueError("perturbed date must lie inside the horizon")
    if len(columns) != 2:
        raise ValueError("exactly two derivative columns are required for the two-block root")

    lags = np.arange(-perturbed_date, horizon - perturbed_date)
    jacobian = np.zeros((2 * horizon, 2 * horizon))
    profiles = {}
    for unknown_index, supplied in enumerate(columns):
        blocks = _column(supplied, horizon).reshape(2, horizon)
        for residual_index, profile in enumerate(blocks):
            jacobian[residual_index * horizon:(residual_index + 1) * horizon,
                     unknown_index * horizon:(unknown_index + 1) * horizon] = toeplitz_block(
                         lags, profile, horizon)
            profiles[f"{RESIDUAL_BLOCKS[residual_index]}<-{UNKNOWN_BLOCKS[unknown_index]}"] = profile.tolist()

    measured = set(lags.tolist())
    zero_filled = sum(
        1 for date in range(horizon) for shock_date in range(horizon)
        if date - shock_date not in measured
    )
    receipt = dict(
        horizon=horizon,
        perturbed_date=perturbed_date,
        measured_lags=lags.tolist(),
        lag_profiles=profiles,
        zero_filled_entries_per_block=zero_filled,
        unknown_blocks=list(UNKNOWN_BLOCKS),
        residual_blocks=list(RESIDUAL_BLOCKS),
        coordinate_convention=("d physical residual / d log(coordinate); rows are "
                               "[housing imbalance, unscaled pension imbalance]"),
        fake_news_derivatives_constructed=False,
    )
    return jacobian, receipt


def _validate_initial_jacobian(initial_jacobian, horizon):
    if initial_jacobian is None:
        return
    matrix = np.asarray(initial_jacobian, dtype=float)
    if matrix.shape == (3 * horizon, 3 * horizon):
        raise ValueError("old three-block scaled-200 Jacobian is incompatible with the two-block physical root")
    if matrix.shape != (2 * horizon, 2 * horizon) or not np.isfinite(matrix).all():
        raise ValueError("initial_jacobian must be a finite 2T by 2T physical-log-coordinate matrix")


def extend_measured_jacobian(receipt, horizon):
    """Reuse measured lag responses as an approximate seed at another horizon.

    Unmeasured lags are zero, never extrapolated economic derivatives. The
    nonlinear root must still pass its original physical and replay gates.
    """
    if type(horizon) is not int or horizon<=0:
        raise ValueError('Positive integer forecast horizon required')
    if (receipt.get('unknown_blocks')!=list(UNKNOWN_BLOCKS) or
            receipt.get('residual_blocks')!=list(RESIDUAL_BLOCKS) or
            receipt.get('residual_units')!='physical_unscaled' or
            receipt.get('coordinate_order')!=list(UNKNOWN_BLOCKS)):
        raise ValueError('Measured seed must use the two-block physical/log convention')
    measured=receipt['horizon'];date=receipt['perturbed_date']
    if type(measured) is not int or measured<=0 or type(date) is not int or not 0<=date<measured:
        raise ValueError('Invalid measured horizon/date')
    lags=list(range(-date,measured-date))
    if receipt['measured_lags']!=lags:raise ValueError('Measured lag support changed')
    names={r+'<-'+u for r in RESIDUAL_BLOCKS for u in UNKNOWN_BLOCKS}
    if set(receipt['lag_profiles'])!=names:raise ValueError('Complete two-block lag profiles required')
    matrix=np.zeros((2*horizon,2*horizon))
    for i,residual in enumerate(RESIDUAL_BLOCKS):
        for j,unknown in enumerate(UNKNOWN_BLOCKS):
            profile=np.asarray(receipt['lag_profiles'][residual+'<-'+unknown],float)
            if profile.shape!=(measured,) or not np.isfinite(profile).all():
                raise ValueError('Finite measured lag profile required')
            matrix[i*horizon:(i+1)*horizon,j*horizon:(j+1)*horizon]=toeplitz_block(lags,profile,horizon)
    return matrix


@contextmanager
def _scaled_generic_root():
    """Install the scaled-step generic root for one serial root invocation."""
    original = social_security_root.solve_price_path
    social_security_root.solve_price_path = solve_price_path_scaled
    try:
        yield
    finally:
        social_security_root.solve_price_path = original


def solve_joint_with_acceleration(**kwargs):
    """Run the existing two-block joint root with a uniform log-step direction.

    The wrapper is intentionally serial: it temporarily replaces the imported
    generic root binding in :mod:`e5f_social_security_root`, restores it even
    after an exception, and otherwise delegates all economics and certification
    to ``solve_social_security_path`` unchanged.
    """
    if "initial_prices" not in kwargs:
        raise ValueError("initial_prices is required to validate the two-block Jacobian")
    prices = np.asarray(kwargs["initial_prices"], dtype=float)
    if prices.ndim != 1 or not prices.size:
        raise ValueError("initial_prices must be a nonempty one-dimensional path")
    _validate_initial_jacobian(kwargs.get("initial_jacobian"), prices.size)
    with _scaled_generic_root():
        return social_security_root.solve_social_security_path(**kwargs)
