"""Canonical parameter-vector handling for the Python DT model."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Sequence

import numpy as np

LEGACY_PARAM_NAMES = [
    "beta",
    "b_entry_fixed",
    "psi_child",
    "h_bar_jump",
    "h_bar_n",
    "c_bar_n",
    "kappa_fert",
    "chi",
    "kappa_loc",
    "mu_move",
    "theta0",
    "theta_n",
]

PARAM_NAMES = LEGACY_PARAM_NAMES + ["h_bar_0"]

LB_FULL = np.array([0.90, 0.00, 0.02, 0.20, 0.10, 0.05, 1.25, 1.01, 0.50, 0.00, 0.10, 0.00, 1.50])
UB_FULL = np.array([0.99, 1.00, 0.12, 2.50, 1.00, 0.25, 8.00, 1.10, 4.00, 10.00, 1.50, 0.75, 4.50])


def theta_names(theta: Sequence[float], names: Sequence[str] | None = None) -> list[str]:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    if names is not None:
        if len(theta) != len(names):
            raise ValueError(f"theta length {len(theta)} does not match names length {len(names)}")
        return list(names)
    if len(theta) == len(PARAM_NAMES):
        return PARAM_NAMES
    if len(theta) == len(LEGACY_PARAM_NAMES):
        return LEGACY_PARAM_NAMES
    raise ValueError(f"theta has {len(theta)} entries; expected {len(LEGACY_PARAM_NAMES)} or {len(PARAM_NAMES)}")


def theta_bounds(theta: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    if len(theta) == len(PARAM_NAMES):
        return LB_FULL, UB_FULL
    if len(theta) == len(LEGACY_PARAM_NAMES):
        return LB_FULL[: len(LEGACY_PARAM_NAMES)], UB_FULL[: len(LEGACY_PARAM_NAMES)]
    raise ValueError(f"theta has {len(theta)} entries; expected {len(LEGACY_PARAM_NAMES)} or {len(PARAM_NAMES)}")


def validate_theta(theta: Sequence[float]) -> bool:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    lb, ub = theta_bounds(theta)
    return bool(np.all(theta >= lb) and np.all(theta <= ub) and np.all(np.isfinite(theta)))


def apply_theta(P: SimpleNamespace, theta: Sequence[float], names: Sequence[str] | None = None) -> SimpleNamespace:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    active_names = theta_names(theta, names)
    for name, value in zip(active_names, theta):
        setattr(P, name, float(value))
    P.eps_fert = P.kappa_fert
    P.eps_loc = P.kappa_loc
    return P

