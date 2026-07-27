"""Externally measured permanent income levels for the E6b variant."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from .externals import flhsv_income_overrides

E6B_PROFILE_NAME = "eqscale_seq_e6b_permanent_income_levels_20260727"
E6B_FIXED_LOG_VARIANCE = 0.3930530
E6B_FIXED_LOG_VARIANCE_SE = 0.0332553
E6B_FIXED_LOG_VARIANCE_P025 = 0.3299568
E6B_FIXED_LOG_VARIANCE_P975 = 0.4528349
E6B_LEVEL_COUNT = 3


def permanent_income_levels(
    log_variance: float = E6B_FIXED_LOG_VARIANCE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Three-node Gauss-Hermite rule for a mean-one log income multiplier.

    The unshifted log nodes ``{-sqrt(3)sigma, 0, sqrt(3)sigma}`` with weights
    ``{1/6, 2/3, 1/6}`` match a centered normal's first five moments. A common
    log shift normalizes the arithmetic mean multiplier to one without
    changing the log variance.
    """

    variance = float(log_variance)
    if not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("permanent income log variance must be positive")
    sigma = math.sqrt(variance)
    weights = np.array([1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0])
    raw_log_nodes = sigma * np.array([-math.sqrt(3.0), 0.0, math.sqrt(3.0)])
    log_nodes = raw_log_nodes - math.log(float(weights @ np.exp(raw_log_nodes)))
    levels = np.exp(log_nodes)
    return levels, weights, log_nodes


def e6b_overrides() -> dict[str, Any]:
    """Multiply the five-state FL-HSV chain by permanent entry levels."""

    base = flhsv_income_overrides()
    base_grid = np.asarray(base["z_grid"], dtype=float)
    base_weights = np.asarray(base["z_weights"], dtype=float)
    base_transition = np.asarray(base["Pi_z"], dtype=float)
    levels, level_weights, _ = permanent_income_levels()
    combined_grid = np.kron(levels, base_grid)
    combined_weights = np.kron(level_weights, base_weights)
    combined_transition = np.kron(np.eye(E6B_LEVEL_COUNT), base_transition)
    level_index = np.repeat(np.arange(E6B_LEVEL_COUNT), base_grid.size)
    base_state_index = np.tile(np.arange(base_grid.size), E6B_LEVEL_COUNT)
    return {
        **base,
        "z_grid": combined_grid,
        "z_weights": combined_weights,
        "Pi_z": combined_transition,
        "permanent_income_levels_enabled": True,
        "permanent_income_log_variance": E6B_FIXED_LOG_VARIANCE,
        "permanent_income_level_values": levels,
        "permanent_income_level_weights": level_weights,
        "permanent_income_group_index": level_index,
        "permanent_income_base_state_index": base_state_index,
    }


def e6b_metadata() -> dict[str, Any]:
    """Serializable empirical restriction and discretization record."""

    levels, weights, log_nodes = permanent_income_levels()
    return {
        "profile": E6B_PROFILE_NAME,
        "source": (
            "PSID EARNINDRRC, reference persons ages 25-60, 1984-2019; "
            "fixed effect + AR(1) + transitory long-lag minimum distance"
        ),
        "fixed_log_variance": E6B_FIXED_LOG_VARIANCE,
        "bootstrap_se": E6B_FIXED_LOG_VARIANCE_SE,
        "bootstrap_p025": E6B_FIXED_LOG_VARIANCE_P025,
        "bootstrap_p975": E6B_FIXED_LOG_VARIANCE_P975,
        "level_values": levels,
        "level_weights": weights,
        "log_level_nodes": log_nodes,
        "level_count": E6B_LEVEL_COUNT,
        "base_income_state_count": 5,
        "combined_income_state_count": 5 * E6B_LEVEL_COUNT,
        "identification": (
            "Externally imposed permanent log-income variance; no added free "
            "calibration parameter."
        ),
    }
