"""Outside-option population closure prototype.

This file is an experiment copy. It is not imported by the validated Python
port. The purpose is to isolate the entry-mass calculation needed to replace
the current fixed-mass normalization once we are ready to alter the solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np


@dataclass(frozen=True)
class EntryClosureResult:
    """Entry masses implied by the outside-option city-size closure."""

    entry_by_loc: np.ndarray
    entry_shares: np.ndarray
    city_entry_total: float
    outside_entry: float
    potential_entrant_mass: float
    city_entry_prob: np.ndarray
    outside_prob: float
    entry_values: np.ndarray


def _logit_with_outside(entry_values: np.ndarray, outside_value: float, kappa_entry: float) -> tuple[np.ndarray, float]:
    """Return city and outside probabilities for a logit entry problem."""

    values = np.concatenate([np.asarray(entry_values, dtype=float).reshape(-1), [float(outside_value)]])
    kappa = max(float(kappa_entry), 1e-10)
    z = values / kappa
    z = z - np.max(z)
    probs = np.exp(z)
    probs = probs / np.sum(probs)
    return probs[:-1], float(probs[-1])


def interpolate_entry_values(
    value_function: np.ndarray,
    b_grid: np.ndarray,
    b_entry: float,
    *,
    renter_tenure: int = 0,
    first_age: int = 0,
    childless_parity: int = 0,
    no_child_state: int = 0,
) -> np.ndarray:
    """Compute entry values V_1(b_entry, R, i, 0, 0) by location.

    The Python port uses value-function arrays ordered as
    `(wealth, tenure, location, age, parity, child_state)`.
    """

    V = np.asarray(value_function, dtype=float)
    bg = np.asarray(b_grid, dtype=float).reshape(-1)
    n_locations = V.shape[2]
    out = np.empty(n_locations)
    b0 = float(np.clip(b_entry, bg[0], bg[-1]))
    for loc in range(n_locations):
        out[loc] = np.interp(
            b0,
            bg,
            V[:, renter_tenure, loc, first_age, childless_parity, no_child_state],
        )
    return out


def compute_outside_option_entry(
    entry_values: np.ndarray,
    *,
    mature_cityborn_flow: float,
    outside_entry_flow: float,
    outside_value: float,
    kappa_entry: float,
    local_birth_entry_weight: float = 1.0,
) -> EntryClosureResult:
    """Map entry values and mature local births into city entry masses.

    The potential entrant mass is

        M = outside_entry_flow + local_birth_entry_weight * mature_cityborn_flow.

    Potential entrants choose between the modeled city locations and the outside
    option. The outside option absorbs the people who do not enter the city, so
    fertility affects the stationary city size without mechanically creating
    growth inside the city.
    """

    mature_flow = max(float(mature_cityborn_flow), 0.0)
    outside_flow = max(float(outside_entry_flow), 0.0)
    local_weight = max(float(local_birth_entry_weight), 0.0)
    potential_mass = outside_flow + local_weight * mature_flow

    city_prob, outside_prob = _logit_with_outside(entry_values, outside_value, kappa_entry)
    entry_by_loc = potential_mass * city_prob
    city_entry_total = float(np.sum(entry_by_loc))
    if city_entry_total > 1e-14:
        entry_shares = entry_by_loc / city_entry_total
    else:
        entry_shares = np.ones_like(entry_by_loc) / len(entry_by_loc)

    return EntryClosureResult(
        entry_by_loc=entry_by_loc,
        entry_shares=entry_shares,
        city_entry_total=city_entry_total,
        outside_entry=float(potential_mass * outside_prob),
        potential_entrant_mass=float(potential_mass),
        city_entry_prob=city_prob,
        outside_prob=float(outside_prob),
        entry_values=np.asarray(entry_values, dtype=float).reshape(-1),
    )


def compute_entry_closure_from_value_function(
    value_function: np.ndarray,
    b_grid: np.ndarray,
    P: SimpleNamespace,
    *,
    mature_cityborn_flow: float,
) -> EntryClosureResult:
    """Convenience wrapper using the parameter names proposed for the solver."""

    entry_values = interpolate_entry_values(
        value_function,
        b_grid,
        float(getattr(P, "b_entry_fixed", 0.0)),
    )
    return compute_outside_option_entry(
        entry_values,
        mature_cityborn_flow=mature_cityborn_flow,
        outside_entry_flow=float(getattr(P, "outside_entry_flow", getattr(P, "E_total", 0.0))),
        outside_value=float(getattr(P, "outside_value", 0.0)),
        kappa_entry=float(getattr(P, "kappa_entry", getattr(P, "kappa_loc", 1.0))),
        local_birth_entry_weight=float(getattr(P, "local_birth_entry_weight", 1.0)),
    )

