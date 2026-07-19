"""Immutable case and tolerance contract for the optimized-package promotion battery."""

from __future__ import annotations

from typing import Any

import numpy as np

from .m5_profile import M5_PRICE, M5_THETA


# The exact live M5 search domain. beta_annual is transformed to the model's
# four-year beta before being passed to the solver.
FREE_PARAMETER_BOUNDS: tuple[tuple[str, float, float], ...] = (
    ("beta_annual", 0.80, 0.9995),
    ("alpha_cons", 0.02, 0.98),
    ("c_bar_0", 0.0, 2.0),
    ("c_bar_n", 0.0, 3.0),
    ("h_bar_0", 0.05, 5.80),
    ("h_bar_jump", 0.0, 8.0),
    ("h_bar_n", 0.0, 5.0),
    ("psi_child", -3.0, 3.0),
    ("kappa_fert", 0.02, 50.0),
    ("chi", 0.10, 5.0),
    ("H0", 0.20, 80.0),
    ("theta0", 0.0, 8.0),
    ("theta1", 0.02, 16.0),
    ("tenure_choice_kappa", 0.0, 0.12),
)
ARRAY_NAMES: tuple[str, ...] = (
    "V",
    "c_pol",
    "hR_pol",
    "bp_pol",
    "tenure_choice",
    "tenure_probs",
    "loc_probs",
    "fert_probs",
    "fert_value",
    "g",
)

PARITY_TOLERANCES = {
    "array_max_abs": 5.0e-12,
    "target_moment_max_abs": 5.0e-10,
    "mass_abs": 5.0e-12,
    "exception_dead_mass_abs": 5.0e-12,
    "exception_dead_mass_rel": 5.0e-8,
}

ROOT_TOLERANCES = {
    "strict_residual": 2.5e-5,
    "repeat_price_abs": 0.0,
    "repeat_moment_abs": 0.0,
    "fixed_map_excess_abs": 5.0e-10,
}


def solver_override(name: str, value: float) -> dict[str, Any]:
    """Translate one search coordinate into the runtime parameterization."""

    if name == "beta_annual":
        return {"beta": float(value) ** 4.0}
    if name == "H0":
        return {"H0": np.array([float(value)])}
    return {name: float(value)}


def bound_cases() -> tuple[dict[str, Any], ...]:
    cases: list[dict[str, Any]] = [{"name": "baseline", "changes": {}}]
    for name, lower, upper in FREE_PARAMETER_BOUNDS:
        cases.append(
            {
                "name": f"{name}__lower",
                "changes": solver_override(name, lower),
                "search_coordinate": {name: lower},
            }
        )
        cases.append(
            {
                "name": f"{name}__upper",
                "changes": solver_override(name, upper),
                "search_coordinate": {name: upper},
            }
        )
    return tuple(cases)


DEMAND_PRICES: tuple[float, ...] = (
    0.01,
    0.02,
    0.05,
    0.10,
    0.20,
    0.35,
    0.50,
    0.65,
    0.72,
    M5_PRICE,
    0.80,
    0.90,
    1.20,
    2.00,
    5.00,
    10.00,
    30.00,
)


# These cases exercise smooth, deterministic-choice, high-logit-temperature,
# supply, and housing-demand perturbations without redefining the M5 target
# system. Exact-bound feasibility is covered separately by bound_cases().
ROOT_CASES: tuple[dict[str, Any], ...] = (
    {"name": "baseline", "changes": {}},
    {"name": "deterministic_tenure", "changes": {"tenure_choice_kappa": 0.0}},
    {"name": "high_tenure_smoothing", "changes": {"tenure_choice_kappa": 0.12}},
    {"name": "supply_minus20", "changes": {"H0": np.array([0.8 * M5_THETA["H0"]])}},
    {"name": "supply_plus20", "changes": {"H0": np.array([1.2 * M5_THETA["H0"]])}},
    {"name": "larger_child_space", "changes": {"h_bar_n": 0.35}},
)


# A deterministic, deliberately small calibration-throughput surface. It is a
# numerical smoke, not a paper calibration or a new identifying specification.
CALIBRATION_CASES: tuple[dict[str, Any], ...] = (
    {"name": "baseline", "changes": {}},
    {"name": "beta_down", "changes": {"beta": 0.955}},
    {"name": "beta_up", "changes": {"beta": 0.975}},
    {"name": "cbar0_down", "changes": {"c_bar_0": 1.15}},
    {"name": "psi_child_down", "changes": {"psi_child": 0.10}},
    {"name": "hbar_n_up", "changes": {"h_bar_n": 0.35}},
    {"name": "tenure_deterministic", "changes": {"tenure_choice_kappa": 0.0}},
    {"name": "theta0_down", "changes": {"theta0": 0.28}},
    {"name": "theta1_up", "changes": {"theta1": 0.44}},
    {"name": "supply_plus5", "changes": {"H0": np.array([1.05 * M5_THETA["H0"]])}},
)
