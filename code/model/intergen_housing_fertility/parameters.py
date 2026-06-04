"""Parameter setup for the no-geography intergenerational housing model."""

from __future__ import annotations

import copy
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np


def clone_ns(obj: SimpleNamespace) -> SimpleNamespace:
    return copy.deepcopy(obj)


def apply_overrides(P: SimpleNamespace, overrides: Mapping[str, Any] | SimpleNamespace | None = None) -> SimpleNamespace:
    if overrides is None:
        return P
    items = vars(overrides).items() if isinstance(overrides, SimpleNamespace) else dict(overrides).items()
    for key, value in items:
        setattr(P, key, value)
    return finalize_parameters(P)


def setup_parameters(mode: str = "benchmark") -> SimpleNamespace:
    """Return baseline parameters.

    The model is intentionally not calibrated. `mode="smoke"` makes the grids
    smaller for fast verification.
    """

    mode = (mode or "benchmark").lower()
    P = SimpleNamespace()

    P.mode = mode
    P.age_start = 18
    P.J = 42
    P.retire_age = 66
    P.fertility_choice_age = 28
    P.beta = 0.96
    P.r = 0.03
    P.R = 1.0 + P.r
    P.sigma = 1.0
    P.alpha_h = 0.36
    P.beta_n = 0.515
    P.c_min = 1e-5
    P.h_need_base = 1.05
    P.h_need_child = 0.70
    P.child_cost_0 = 0.025
    P.child_cost_income = 0.010

    P.Nb = 55
    P.b_min = -2.0
    P.b_mid = 8.0
    P.b_max = 45.0
    P.b_grid_power = 1.6
    P.b_entry = 0.0

    P.z_grid = np.array([0.70, 1.00, 1.35])
    P.z_dist = np.array([0.25, 0.50, 0.25])
    P.Pi_z = np.array(
        [
            [0.82, 0.16, 0.02],
            [0.10, 0.80, 0.10],
            [0.02, 0.16, 0.82],
        ],
        dtype=float,
    )
    P.wage = 1.0
    P.pension_replacement = 0.55
    P.income_age_breaks = np.array([18.0, 25.0, 35.0, 45.0, 55.0, 66.0])
    P.income_age_values = np.array([0.58, 0.84, 1.00, 1.04, 0.98, 0.80])

    P.n_child_options = np.array([0, 1, 2])

    # Tenure index 0 is renter. Owner indices 1..K correspond to owner sizes.
    P.owner_h = np.array([2.2, 3.8, 5.8])
    # Service units per normalized adult in the lifecycle cross-section.
    P.owner_supply = np.array([0.45, 1.15, 0.55])
    P.renter_h = 2.1
    P.rent_user_cost = 0.078

    P.owner_user_cost = np.array([0.055, 0.062, 0.070])
    P.delta = 0.02
    P.tau_property = 0.012
    P.phi_ltv = 0.80
    P.psi_pti = 0.28
    P.mortgage_rate = 0.045
    P.mortgage_maturity = 30

    P.buyer_transaction_cost = 0.00
    P.owner_move_cost = 0.035
    P.old_retention_age = 62
    P.old_retention_wedge = 0.12

    P.max_iter_eq = 60
    P.tol_eq = 5e-4
    P.price_damping = 0.15
    P.price_min = 0.025
    P.price_max = 0.25

    if mode == "smoke":
        P.J = 18
        P.retire_age = 55
        P.fertility_choice_age = 24
        P.Nb = 24
        P.b_min = -1.0
        P.b_mid = 5.0
        P.b_max = 20.0
        P.z_grid = np.array([0.80, 1.20])
        P.z_dist = np.array([0.55, 0.45])
        P.Pi_z = np.array([[0.86, 0.14], [0.14, 0.86]], dtype=float)
        P.max_iter_eq = 8
        P.tol_eq = 2e-3
        P.price_damping = 0.10

    return finalize_parameters(P)


def finalize_parameters(P: SimpleNamespace) -> SimpleNamespace:
    P.z_grid = np.asarray(P.z_grid, dtype=float)
    P.z_dist = np.asarray(P.z_dist, dtype=float)
    P.z_dist = P.z_dist / P.z_dist.sum()
    P.Pi_z = np.asarray(P.Pi_z, dtype=float)
    P.Pi_z = P.Pi_z / P.Pi_z.sum(axis=1, keepdims=True)
    P.owner_h = np.asarray(P.owner_h, dtype=float)
    P.owner_supply = np.asarray(P.owner_supply, dtype=float)
    P.owner_user_cost = np.asarray(P.owner_user_cost, dtype=float)
    P.n_child_options = np.asarray(P.n_child_options, dtype=int)
    P.K = len(P.owner_h)
    P.Nt = P.K + 1
    P.Nn = len(P.n_child_options)
    P.rho_property = P.r + P.delta + P.tau_property
    P.fertility_choice_index = int(np.clip(P.fertility_choice_age - P.age_start, 0, P.J - 1))
    P.old_retention_index = int(np.clip(P.old_retention_age - P.age_start, 0, P.J - 1))
    return P
