"""Construct a generic four-year persistent-plus-transitory earnings process.

The constructor works at the model period frequency.  It deliberately does not
map annual estimates into four-year parameters: that source-specific decision
belongs outside this reusable numerical adapter.  The persistent component is
a Rouwenhorst approximation to a stationary log AR(1), and the transitory
component is an independent mean-one Gauss--Hermite lognormal rule.

Only override dictionaries and auditable moment metadata are returned.  No
household model is imported or solved.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

import build_literature_period_income as _literature


_PERIOD_YEARS = 4
_DIAGNOSTIC_PERSISTENT_COUNTS = (5, 9, 15)
_DIAGNOSTIC_IID_COUNTS = (3, 5)


def _positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return value


def _integer_count(name: str, value: int) -> int:
    numeric = float(value)
    if not np.isfinite(numeric) or not numeric.is_integer() or numeric < 2:
        raise ValueError(f"{name} must be an integer at least two")
    return int(numeric)


def _validate_inputs(
    rho_period: float,
    persistent_innovation_sd_period: float,
    transitory_sd_period: float,
    n_persistent: int,
    n_iid: int,
) -> tuple[float, float, float, int, int]:
    rho_period = float(rho_period)
    if not np.isfinite(rho_period) or not 0.0 < rho_period < 1.0:
        raise ValueError("rho_period must be finite and strictly between zero and one")
    persistent_sd = _positive_finite(
        "persistent_innovation_sd_period", persistent_innovation_sd_period
    )
    transitory_sd = float(transitory_sd_period)
    if not np.isfinite(transitory_sd) or transitory_sd < 0.0:
        raise ValueError("transitory_sd_period must be finite and nonnegative")
    n_persistent = _integer_count("n_persistent", n_persistent)
    n_iid_numeric = float(n_iid)
    if not np.isfinite(n_iid_numeric) or not n_iid_numeric.is_integer():
        raise ValueError("n_iid must be an integer")
    n_iid = int(n_iid_numeric)
    if transitory_sd == 0.0:
        if n_iid != 1:
            raise ValueError("zero transitory SD requires n_iid=1")
    elif n_iid < 2:
        raise ValueError("positive transitory SD requires n_iid at least 2")
    else:
        n_iid = _integer_count("n_iid", n_iid)
    return rho_period, persistent_sd, transitory_sd, n_persistent, n_iid


def _components(
    rho_period: float,
    persistent_sd: float,
    transitory_sd: float,
    n_persistent: int,
    n_iid: int,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    persistent_variance = persistent_sd**2 / (1.0 - rho_period**2)
    persistent, persistent_weights, persistent_transition = _literature._rouwenhorst_period(
        n_persistent, rho_period, persistent_variance
    )
    if transitory_sd == 0.0:
        iid = np.ones(1, dtype=float)
        iid_weights = np.ones(1, dtype=float)
    else:
        iid, iid_weights = _literature._iid_lognormal_rule(transitory_sd, n_iid)
    return (
        {
            "persistent_levels": persistent,
            "persistent_weights": persistent_weights,
            "persistent_transition": persistent_transition,
            "iid_levels": iid,
            "iid_weights": iid_weights,
        },
        {
            "persistent_log_variance": np.asarray(persistent_variance),
            "transitory_log_variance": np.asarray(transitory_sd**2),
        },
    )


def _flatten_components(
    components: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    persistent = components["persistent_levels"]
    persistent_weights = components["persistent_weights"]
    persistent_transition = components["persistent_transition"]
    iid = components["iid_levels"]
    iid_weights = components["iid_weights"]

    z_grid = np.multiply.outer(persistent, iid).reshape(-1)
    z_weights = np.multiply.outer(persistent_weights, iid_weights).reshape(-1)
    iid_transition = np.broadcast_to(iid_weights, (iid_weights.size, iid_weights.size)).copy()
    pi_z = np.kron(persistent_transition, iid_transition)
    z_grid /= float(z_weights @ z_grid)
    z_weights /= float(z_weights.sum())
    pi_z /= pi_z.sum(axis=1, keepdims=True)
    return z_grid, z_weights, pi_z


def _continuous_covariances(
    rho_period: float, persistent_sd: float, transitory_sd: float, max_lag: int = 4
) -> tuple[np.ndarray, np.ndarray]:
    vp = persistent_sd**2 / (1.0 - rho_period**2)
    ve = transitory_sd**2
    log_cov = np.array(
        [vp + ve] + [vp * rho_period**lag for lag in range(1, max_lag + 1)],
        dtype=float,
    )
    return log_cov, np.expm1(log_cov)


def _json_components(components: dict[str, np.ndarray]) -> dict[str, list[Any]]:
    return {key: np.asarray(value).tolist() for key, value in components.items()}


def _resolution_table(
    rho_period: float,
    persistent_sd: float,
    transitory_sd: float,
    continuous_level_covariances: np.ndarray,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    iid_counts = _DIAGNOSTIC_IID_COUNTS if transitory_sd > 0.0 else (1,)
    for n_persistent in _DIAGNOSTIC_PERSISTENT_COUNTS:
        for n_iid in iid_counts:
            raw_components, _ = _components(
                rho_period,
                persistent_sd,
                transitory_sd,
                n_persistent,
                n_iid,
            )
            z_grid, z_weights, pi_z = _flatten_components(raw_components)
            discrete = _literature._chain_level_covariances(z_grid, z_weights, pi_z, 4)
            rows.append(
                {
                    "n_persistent": n_persistent,
                    "n_iid": n_iid,
                    "joint_states": int(z_grid.size),
                    "level_covariances": discrete.tolist(),
                    "errors_vs_continuous": (discrete - continuous_level_covariances).tolist(),
                    "maximum_absolute_error": float(np.max(np.abs(discrete - continuous_level_covariances))),
                    "mean": float(z_weights @ z_grid),
                }
            )
    return rows


def build_period_earnings_process(
    *,
    rho_period: float,
    persistent_innovation_sd_period: float,
    transitory_sd_period: float,
    n_persistent: int = 5,
    n_iid: int = 3,
) -> tuple[dict[str, Any], dict[str, Any]]:
    r"""Return native income overrides and a four-year-frequency receipt.

    ``rho_period`` is persistence per four-year model period.  Both standard
    deviations are standard deviations of log earnings shocks at that same
    frequency.  The stationary persistent log variance is
    \(\sigma_p^2/(1-\rho^2)\), while the independent transitory log variance is
    \(\sigma_e^2\).  A zero transitory standard deviation requires the explicit
    single iid node ``n_iid=1``, giving a pure persistent process.
    """
    rho_period, persistent_sd, transitory_sd, n_persistent, n_iid = _validate_inputs(
        rho_period,
        persistent_innovation_sd_period,
        transitory_sd_period,
        n_persistent,
        n_iid,
    )
    components, component_variances = _components(
        rho_period, persistent_sd, transitory_sd, n_persistent, n_iid
    )
    z_grid, z_weights, pi_z = _flatten_components(components)
    continuous_log, continuous_level = _continuous_covariances(
        rho_period, persistent_sd, transitory_sd, 4
    )
    discrete_log = _literature._chain_log_covariances(z_grid, z_weights, pi_z, 4)
    discrete_level = _literature._chain_level_covariances(z_grid, z_weights, pi_z, 4)
    drift_row = z_weights @ pi_z
    stationarity_error = drift_row - z_weights

    overrides: dict[str, Any] = {
        "use_income_types": True,
        "income_type_transition": "markov",
        "income_shock_persistence": rho_period,
        "z_grid": z_grid,
        "z_weights": z_weights,
        "Pi_z": pi_z,
        "permanent_income_levels_enabled": False,
        "permanent_income_log_variance": 0.0,
    }
    metadata: dict[str, Any] = {
        "frequency": "four_year",
        "period_years": _PERIOD_YEARS,
        "n_persistent": n_persistent,
        "n_iid": n_iid,
        "joint_states": int(z_grid.size),
        "rho_period": rho_period,
        "persistent_innovation_sd_period": persistent_sd,
        "persistent_innovation_variance_period": persistent_sd**2,
        "transitory_sd_period": transitory_sd,
        "transitory_log_variance_period": transitory_sd**2,
        "stationary_persistent_log_variance": float(component_variances["persistent_log_variance"]),
        "stationary_log_variance": float(continuous_log[0]),
        "continuous_period_log_covariances": continuous_log.tolist(),
        "continuous_period_level_covariances": continuous_level.tolist(),
        "discrete_period_log_covariances": discrete_log.tolist(),
        "discrete_period_level_covariances": discrete_level.tolist(),
        # Short aliases match the existing endpoint adapter's consumer-facing
        # receipt while the period-qualified names above remain unambiguous.
        "continuous_log_covariances": continuous_log.tolist(),
        "continuous_level_covariances": continuous_level.tolist(),
        "discrete_log_covariances": discrete_log.tolist(),
        "discrete_level_covariances": discrete_level.tolist(),
        "discrete_level_errors_vs_continuous": (discrete_level - continuous_level).tolist(),
        "discrete_log_errors_vs_continuous": (discrete_log - continuous_log).tolist(),
        "drift_row": drift_row.tolist(),
        "driftrow": drift_row.tolist(),
        "stationary_weights": z_weights.tolist(),
        "stationarity_error": stationarity_error.tolist(),
        "stationarity_error_max_abs": float(np.max(np.abs(stationarity_error))),
        "stationarity": bool(np.allclose(drift_row, z_weights, atol=1e-12, rtol=0.0)),
        "stationary": bool(np.allclose(drift_row, z_weights, atol=1e-12, rtol=0.0)),
        "mean": float(z_weights @ z_grid),
        "row_sum_error_max_abs": float(np.max(np.abs(pi_z.sum(axis=1) - 1.0))),
        "iid_transition_independent": bool(
            np.allclose(
                pi_z.reshape(n_persistent, n_iid, n_persistent, n_iid)[:, 0, :, :],
                pi_z.reshape(n_persistent, n_iid, n_persistent, n_iid)[:, -1, :, :],
                atol=1e-14,
                rtol=0.0,
            )
        ),
        "sd_variance_units": {
            "rho_period": "persistence per four-year period",
            "persistent_innovation_sd_period": "standard deviation of period log innovation",
            "persistent_innovation_variance_period": "variance of period log innovation",
            "transitory_sd_period": "standard deviation of period log iid shock",
            "transitory_log_variance_period": "variance of period log iid shock",
        },
        "components": _json_components(components),
        "persistent_component": _json_components(components),
        "resolution_table": _resolution_table(
            rho_period, persistent_sd, transitory_sd, continuous_level
        ),
        "permanent_income_levels_enabled": False,
    }
    return overrides, metadata


__all__ = ["build_period_earnings_process"]
