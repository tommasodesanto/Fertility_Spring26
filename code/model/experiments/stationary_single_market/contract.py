"""Supported structural choices; numerical parameters and grids stay configurable."""
from __future__ import annotations
import numpy as np


def validate_contract(P):
    required = {
        "I": 1,
        "child_state_mode": "independent_count",
        "sequential_births": True,
        "joint_nested_choice": False,
        "two_shock_choice": False,
        "fertility_nest_choice": False,
        "permanent_income_levels_enabled": False,
        "income_type_transition": "markov",
        "native_exact_allocation_output": True,
        "n_parity": 4,
        "n_child_states": 4,
    }
    for name, expected in required.items():
        if not hasattr(P, name) or getattr(P, name) != expected:
            raise ValueError(f"Unsupported single-market contract: {name} must equal {expected!r}")
    if not bool(getattr(P, "use_loc_kernel", True)):
        raise ValueError("Single-market engine requires the compiled location stage")
    z = np.asarray(P.z_grid)
    pi = np.asarray(P.Pi_z)
    weights = np.asarray(P.z_weights)
    if int(P.Nz) != z.size or weights.shape != z.shape:
        raise ValueError("Nz and z_weights must match the Markov income grid")
    if (not np.all(np.isfinite(weights)) or np.any(weights < 0)
            or not np.isclose(weights.sum(), 1.0, rtol=0.0, atol=1e-12)):
        raise ValueError("Income-state weights must be finite nonnegative probabilities")
    if z.ndim != 1 or z.size < 1 or pi.shape != (z.size, z.size):
        raise ValueError("Single Markov income process requires z_grid and square Pi_z")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(pi)) or np.any(pi < 0):
        raise ValueError("Markov income inputs must be finite with nonnegative transition entries")
    if not np.allclose(pi.sum(axis=1), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("Markov transition rows must sum to one")
    return P
