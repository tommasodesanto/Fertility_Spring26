"""Stationary block-Toeplitz initial Jacobian for the dated native root.

Pure numpy; no model imports.  This is the cheapest object in the sequence-
space family: it does *not* build fake-news derivatives of household policies,
distribution operators or queue laws.  It measures, by central finite
differences of the exact native mapping, the response of every dated residual
to one perturbation date ``s`` in each unknown block, and then assumes the
response depends only on the lag ``t - s`` (stationary Toeplitz structure).
Lags outside the measured window are set to zero and reported.

Coordinate convention (matches ``e5f_matched_pf_path_root.solve_price_path``):
every root coordinate is differentiated in *logs*, including pensions and
rebates, because the retained Broyden solver works on ``x = log(prices)`` for
the whole stacked vector.  Residual ordering follows
``e5f_rebated_surprises.stack_dated_residuals``: ``[housing_0..T-1,
paygo_0..T-1, rebate_0..T-1]`` with both fiscal rows already scaled by 200.
"""
from __future__ import annotations

import numpy as np

UNKNOWN_BLOCKS = ("log_house_price", "log_pension", "log_rebate")
RESIDUAL_BLOCKS = ("housing_relative_imbalance",
                   "paygo_relative_imbalance_scaled_200",
                   "rebate_relative_imbalance_scaled_200")


def _check_residual(residual, horizon):
    residual = np.asarray(residual, dtype=float)
    if residual.shape != (3 * horizon,) or not np.isfinite(residual).all():
        raise ValueError("stacked residual must be finite with shape (3T,)")
    return residual


def central_column(residual_plus, residual_minus, step, horizon):
    """d residual / d log u_{j,s} as a stacked ``(3T,)`` vector."""
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("finite-difference step must be positive")
    plus = _check_residual(residual_plus, horizon)
    minus = _check_residual(residual_minus, horizon)
    return (plus - minus) / (2.0 * step)


def lag_profiles(column, horizon, perturbed_date):
    """Split a stacked column into three lag profiles keyed by ``t - s``.

    Returns ``(lags, profiles)`` with ``lags = arange(-s, T - s)`` and
    ``profiles`` of shape ``(3, T)`` where row ``i`` is the response of
    residual block ``i`` at lag ``lags[k]``.
    """
    column = _check_residual(column, horizon)
    if not 0 <= perturbed_date < horizon:
        raise ValueError("perturbed date must lie inside the horizon")
    lags = np.arange(-perturbed_date, horizon - perturbed_date)
    profiles = column.reshape(3, horizon)
    return lags, profiles.copy()


def toeplitz_block(lags, profile, horizon):
    """Dense ``(T, T)`` block with ``B[t, s] = profile[t - s]`` (zero outside)."""
    lags = np.asarray(lags, dtype=int)
    profile = np.asarray(profile, dtype=float)
    if lags.shape != profile.shape or lags.ndim != 1 or not np.isfinite(profile).all():
        raise ValueError("lags and profile must be equally sized finite vectors")
    table = dict(zip(lags.tolist(), profile.tolist()))
    block = np.zeros((horizon, horizon))
    for t in range(horizon):
        for s in range(horizon):
            block[t, s] = table.get(t - s, 0.0)
    return block


def assemble_jacobian(columns, horizon, perturbed_date):
    """Build the ``(3T, 3T)`` block-Toeplitz Jacobian.

    ``columns`` is a length-3 sequence of stacked derivative columns, one per
    unknown block in ``UNKNOWN_BLOCKS`` order, all measured at the same date.
    Returns ``(jacobian, receipt)``; the receipt lists the measured lag window,
    every lag profile, and the count of zero-filled unmeasured entries.
    """
    if len(columns) != 3:
        raise ValueError("one measured column per unknown block is required")
    jacobian = np.zeros((3 * horizon, 3 * horizon))
    profiles = {}
    lags = None
    for j, column in enumerate(columns):
        lags, block_profiles = lag_profiles(column, horizon, perturbed_date)
        for i in range(3):
            block = toeplitz_block(lags, block_profiles[i], horizon)
            jacobian[i * horizon:(i + 1) * horizon, j * horizon:(j + 1) * horizon] = block
            profiles[f"{RESIDUAL_BLOCKS[i]}<-{UNKNOWN_BLOCKS[j]}"] = block_profiles[i].tolist()
    measured = set(lags.tolist())
    unmeasured = sum(1 for t in range(horizon) for s in range(horizon) if (t - s) not in measured)
    receipt = dict(horizon=horizon, perturbed_date=perturbed_date,
                   measured_lags=lags.tolist(), lag_profiles=profiles,
                   zero_filled_entries_per_block=unmeasured,
                   unknown_blocks=list(UNKNOWN_BLOCKS), residual_blocks=list(RESIDUAL_BLOCKS),
                   coordinate_convention="d residual / d log(coordinate) for all three blocks",
                   fake_news_derivatives_constructed=False)
    if not np.isfinite(jacobian).all():
        raise ValueError("assembled Jacobian is not finite")
    return jacobian, receipt


def diagonal_default(horizon, slope, fiscal_scale=200.0):
    """The retained solver's default reset Jacobian for comparison."""
    return np.diag(np.concatenate([np.full(horizon, -float(slope)),
                                   np.full(2 * horizon, -float(fiscal_scale))]))


def condition_report(jacobian, default):
    """Scalar diagnostics comparing the measured Jacobian with the diagonal default."""
    jacobian = np.asarray(jacobian, dtype=float)
    default = np.asarray(default, dtype=float)
    return dict(condition_number=float(np.linalg.cond(jacobian)),
                default_condition_number=float(np.linalg.cond(default)),
                relative_frobenius_gap_to_default=float(np.linalg.norm(jacobian - default)
                                                        / np.linalg.norm(default)),
                min_abs_diagonal=float(np.min(np.abs(np.diag(jacobian)))))
