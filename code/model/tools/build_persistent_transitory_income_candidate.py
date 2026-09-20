"""Build a persistent-plus-iid-transitory income-process candidate.

This module only constructs override dictionaries.  It does not mutate model
parameters or run a solve.  The state is flattened so the existing optimized
solver can consume it through its single ``z`` Markov-state axis.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from intergen_eqscale_seq_optimized.local_panel import income_process_overrides


def _iid_lognormal_rule(sd_period: float) -> tuple[np.ndarray, np.ndarray]:
    """Return a three-node mean-one Gauss-Hermite lognormal rule."""
    if not np.isfinite(sd_period) or sd_period < 0.0:
        raise ValueError("transitory_log_sd_period must be finite and nonnegative")
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    weights = weights / math.sqrt(math.pi)
    log_var = float(sd_period) ** 2
    raw = np.exp(-0.5 * log_var + math.sqrt(2.0 * log_var) * nodes)
    raw /= float(weights @ raw)
    return raw, weights


def build_persistent_transitory_income_candidate(
    *,
    rho_annual: float,
    persistent_innovation_sd_annual: float,
    transitory_log_sd_period: float,
    period_years: float,
    persistent_states: int = 5,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return ``(overrides, metadata)`` for a no-permanent-type candidate.

    ``transitory_log_sd_period`` is the period-level standard deviation of
    the log iid shock. Annual-to-period aggregation must be done by the caller.
    """
    from intergen_eqscale_seq_optimized.local_panel import PERIOD_YEARS
    if not np.isclose(period_years, PERIOD_YEARS):
        raise ValueError(f"period_years must equal frozen PERIOD_YEARS={PERIOD_YEARS}")
    if not np.isfinite(rho_annual) or not 0.0 < rho_annual < 1.0:
        raise ValueError("rho_annual must lie strictly between zero and one")
    if not np.isfinite(persistent_innovation_sd_annual) or persistent_innovation_sd_annual <= 0.0:
        raise ValueError("persistent_innovation_sd_annual must be positive")
    persistent = income_process_overrides(
        int(persistent_states),
        "rouwenhorst",
        float(persistent_innovation_sd_annual),
        float(rho_annual),
    )
    eps_grid, eps_weights = _iid_lognormal_rule(float(transitory_log_sd_period))
    z_p = np.asarray(persistent["z_grid"], dtype=float)
    w_p = np.asarray(persistent["z_weights"], dtype=float)
    pi_p = np.asarray(persistent["Pi_z"], dtype=float)
    z_grid = np.multiply.outer(z_p, eps_grid).reshape(-1)
    z_weights = np.multiply.outer(w_p, eps_weights).reshape(-1)
    pi_eps = np.broadcast_to(eps_weights, (eps_weights.size, eps_weights.size)).copy()
    pi_z = np.kron(pi_p, pi_eps)
    overrides = {
        "use_income_types": True,
        "income_type_transition": "markov",
        "income_shock_persistence": float(persistent["income_shock_persistence"]),
        "z_grid": z_grid,
        "z_weights": z_weights,
        "Pi_z": pi_z,
        "permanent_income_levels_enabled": False,
        "permanent_income_log_variance": 0.0,
    }
    metadata = {
        "persistent_states": int(persistent_states),
        "transitory_states": 3,
        "joint_states": int(z_grid.size),
        "rho_annual": float(rho_annual),
        "rho_period": float(rho_annual) ** float(period_years),
        "persistent_innovation_sd_annual": float(persistent_innovation_sd_annual),
        "transitory_log_sd_period": float(transitory_log_sd_period),
        "transitory_log_variance_period": float(transitory_log_sd_period) ** 2,
        "period_years": float(period_years),
        "transitory_rule": "3-point Gauss-Hermite lognormal, mean-one",
        "permanent_types": False,
    }
    return overrides, metadata


__all__ = ["build_persistent_transitory_income_candidate"]
