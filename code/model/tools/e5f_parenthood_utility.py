"""Strict parameter adapter for the approved parenthood-only housing floor.

Use initialize_parenthood_utility once on the inherited jump-plus-slope P.
Use validate_parenthood_utility on reload and bind_parenthood_utility for
candidates. Binding avoids apply_overrides, which also rebuilds pensions and
income outside this adapter's scope. No existing solver or transform is changed.
"""
from __future__ import annotations

import copy
import math
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np
from intergen_eqscale_seq_optimized.e5f_income_entry_profile import E5F_INCOME_ENTRY_DOMAIN

UTILITY_CONTRACT = "e5f_parenthood_only_power_scale_20260911"
FIXED_ALPHA = 0.733
_FIXED_UTILITY = {
    "preference_spec": "eqscale", "eqscale_form": "power",
    "child_room_floor": True, "sigma": 2.0,
    "delta_alpha": 0.0, "delta_alpha_jump": 0.0,
    "c_bar_0": 0.0, "c_bar_n": 0.0,
}
# Retain all unchanged intervals and transforms. h_P uses the positive floor's
# log transform and the approved sum-of-inherited-intervals numerical bounds.
PARENTHOOD_SEARCH_DOMAIN = tuple(
    row for row in E5F_INCOME_ENTRY_DOMAIN
    if row[0] not in {"psi_child", "hbar_child_rooms"}
) + (("h_P", 0.10, 2.30, "log"),)
PARENTHOOD_SEARCH_NAMES = tuple(row[0] for row in PARENTHOOD_SEARCH_DOMAIN)


def _finite_scalar(value: Any, name: str) -> float:
    raw = np.asarray(value)
    if raw.size != 1 or raw.dtype.kind == "b":
        raise ValueError(f"{name} must be one finite numeric value")
    try:
        out = float(raw.reshape(-1)[0])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be one finite numeric value") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _validate_lifecycle(parameters: SimpleNamespace) -> None:
    if (float(parameters.period_years) != 4.0
            or str(parameters.child_state_mode) != "independent_count"
            or int(parameters.n_parity) != 4 or int(parameters.n_child_states) != 4
            or not bool(parameters.use_stochastic_aging)):
        raise ValueError("Parenthood utility requires the retained four-year independent child-count lifecycle")
    if _finite_scalar(parameters.alpha_cons, "alpha_cons") != FIXED_ALPHA:
        raise ValueError(f"Retain the externally fixed constant alpha_cons={FIXED_ALPHA}")


def initial_parenthood_requirement(parameters: SimpleNamespace) -> float:
    """Old first-child floor; this is a migration rule, not a reload rule."""
    jump = _finite_scalar(parameters.hbar_first_child_jump, "hbar_first_child_jump")
    slope = _finite_scalar(parameters.hbar_child_rooms, "hbar_child_rooms")
    if min(jump, slope) < 0.0:
        raise ValueError("Inherited child-room jump and slope must be nonnegative")
    return jump + slope


def validate_parenthood_utility(parameters: SimpleNamespace) -> None:
    """Fail on a stale/changed loaded P; never repair a reintroduced slope."""
    _validate_lifecycle(parameters)
    if _finite_scalar(parameters.hbar_child_rooms, "hbar_child_rooms") != 0.0:
        raise ValueError("Parenthood-only utility requires hbar_child_rooms == 0 exactly")
    for name, expected in _FIXED_UTILITY.items():
        value = getattr(parameters, name, None)
        if isinstance(expected, str):
            valid = value == expected
        elif isinstance(expected, bool):
            valid = value is expected
        else:
            valid = _finite_scalar(value, name) == expected
        if not valid:
            raise ValueError(f"Parenthood-only fixed utility requires {name}={expected!r}")
    floor = _finite_scalar(parameters.hbar_first_child_jump, "h_P")
    cap = _finite_scalar(parameters.hR_max, "hR_max")
    if not 0.0 <= floor < cap:
        raise ValueError("Parenthood floor must be nonnegative and strictly below hR_max")
    if not 0.0 < _finite_scalar(parameters.beta, "beta") < 1.0:
        raise ValueError("Period beta must lie strictly between zero and one")
    _finite_scalar(parameters.psi_child, "psi_child")


def initialize_parenthood_utility(parameters: SimpleNamespace) -> SimpleNamespace:
    """Copy P and map h_P=old jump+slope, preserving the first-child floor.

    This explicit one-time migration binds the author's fixed utility fields.
    Every other primitive, including fiscal income arrays, discounting, child
    transitions and the child preference intercept, remains unchanged.
    """
    _validate_lifecycle(parameters)
    floor = initial_parenthood_requirement(parameters)
    target = copy.deepcopy(parameters)
    for name, value in _FIXED_UTILITY.items():
        setattr(target, name, value)
    target.hbar_first_child_jump = floor
    target.hbar_child_rooms = 0.0
    validate_parenthood_utility(target)
    return target


def validate_parenthood_candidate(
    candidate: Mapping[str, Any] | None, *, require_complete: bool = False,
) -> dict[str, float]:
    """Validate physical coordinates; old jump/slope keys are always invalid.

    Partial candidates support bounded coordinate probes. A production launcher
    can require all nine coordinates explicitly. psi_child is normalized outside
    this structural vector.
    """
    values = {} if candidate is None else dict(candidate)
    unknown = set(values) - set(PARENTHOOD_SEARCH_NAMES)
    if unknown:
        raise ValueError(f"Unknown or fixed parenthood-only coordinate(s): {sorted(unknown)}")
    if require_complete and set(values) != set(PARENTHOOD_SEARCH_NAMES):
        raise ValueError("A complete parenthood candidate requires exactly nine coordinates")
    out = {}
    for name, lower, upper, _ in PARENTHOOD_SEARCH_DOMAIN:
        if name in values:
            value = _finite_scalar(values[name], name)
            if not lower <= value <= upper:
                raise ValueError(f"{name} must lie in [{lower}, {upper}]")
            out[name] = value
    return out


def parenthood_utility_overrides(
    parameters: SimpleNamespace, candidate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Translate validated coordinates into existing solver fields only."""
    validate_parenthood_utility(parameters)
    checked = validate_parenthood_candidate(candidate)
    overrides: dict[str, Any] = {}
    for name, value in checked.items():
        if name == "h_P":
            overrides["hbar_first_child_jump"] = value
        elif name == "beta_annual":
            overrides["beta"] = value ** 4
        elif name == "H0":
            shape = np.asarray(parameters.H0).shape
            if np.asarray(parameters.H0).size != 1:
                raise ValueError("A scalar H0 coordinate requires the retained one-market supply")
            overrides[name] = np.full(shape, value) if shape else value
        else:
            overrides[name] = value
    # Preserve the aliases maintained by the existing parameter construction.
    if "beta" in overrides:
        overrides["rho"] = 1.0 / overrides["beta"] - 1.0
        overrides["rho_hat"] = overrides["rho"]
    if "kappa_fert" in overrides:
        overrides["eps_fert"] = overrides["kappa_fert"]
    return overrides


def bind_parenthood_utility(
    parameters: SimpleNamespace, candidate: Mapping[str, Any] | None = None,
    *, copy_parameters: bool = True,
) -> SimpleNamespace:
    """Apply a candidate to validated new-utility P without fiscal rebuilding.

    A changed loaded slope fails before copying, even if the candidate omits it.
    Default deep copying leaves the caller's arrays and primitives independent.
    """
    overrides = parenthood_utility_overrides(parameters, candidate)
    # Validate an independent object even for explicit in-place binding, so an
    # infeasible candidate cannot leave its caller partially mutated.
    target = copy.deepcopy(parameters)
    for name, value in overrides.items():
        setattr(target, name, value)
    validate_parenthood_utility(target)
    if copy_parameters:
        return target
    for name in overrides:
        setattr(parameters, name, getattr(target, name))
    return parameters


def parenthood_utility_metadata() -> dict[str, Any]:
    """Serializable utility/search restrictions for a pinned launch receipt."""
    return {
        "contract": UTILITY_CONTRACT,
        "floor_formula": "hbar(m) = h_P * 1{m > 0}; m is current dependent children",
        "flow_utility": "-((2+0.7*m)/2)**0.7 / (c**alpha * (s-hbar(m))**(1-alpha)) + psi_child*m",
        "fixed_utility": {**_FIXED_UTILITY, "alpha_cons": FIXED_ALPHA},
        "fixed_slope": 0.0,
        "search_domain": PARENTHOOD_SEARCH_DOMAIN,
        "free_parameter_count": len(PARENTHOOD_SEARCH_DOMAIN),
        "h_P_bounds_status": "approved numerical sum-of-inherited-intervals bounds",
        "discount_conversion": "beta = beta_annual**4",
        "psi_child_status": "separately normalized, outside the nine structural coordinates",
        "scale_coefficients": {"adults": 2.0, "child_weight": 0.7, "power": 0.7},
    }
