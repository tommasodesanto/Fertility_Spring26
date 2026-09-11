"""Opt-in approximate initial CPS/NCHS fertility observation; no model solve.

This module is not connected to any production observer or target system.
Call ``observe_initial_fertility(evaluation, P, age_projection=...)`` explicitly
on a stationary sequential evaluation. Its metadata never certifies stationarity,
female exposure, empirical weights, or the separate 2.1 normalization.

The CPS window [40,45) overlaps the [38,42) and [42,46) model cells for two
and three years. Uniform birth-time interpolation averages the pre/post parity
stock using post shares .75 and .375, respectively. Constant-post-cell stock
is an explicitly different diagnostic approximation. Both use model age mass.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np


AGE_PROJECTIONS = ("uniform_birth_time", "constant_post_cell")
MASS_ATOL = 2.0e-10
FLOW_ATOL = 2.0e-10
MOMENT_ATOL = 2.0e-12


def _finite_vector(value: Any, name: str, size: int) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (size,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite length-{size} vector")
    return array


def _maximum_error(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.max(np.abs(left - right)))


def observe_initial_fertility(
    evaluation: Any,
    parameters: Any,
    *,
    age_projection: str,
) -> dict[str, Any]:
    """Return four diagnostic moments plus parity shares and accounting evidence.

    ``evaluation.g_pre`` and ``evaluation.g_post_fertility`` must be the same
    population immediately before/after the sequential fertility step, with
    axes wealth/tenure/location/age/income/parity/child-state. No distribution,
    parameter, normalization, target, or weight is changed. Birth-flow helpers
    are imported lazily from the existing transition measurement module.

    The caller must have configured that module's calendar model consistently
    with the supplied policy, as for its other sequential measurement calls.
    An accepted packet is only an internally consistent approximate observer;
    it is not an equilibrium or empirical-contract certification.
    """
    if age_projection not in AGE_PROJECTIONS:
        raise ValueError(f"Explicit age_projection must be one of {AGE_PROJECTIONS}")
    P = parameters
    geometry = np.asarray([P.J, P.age_start, P.da, P.period_years, P.n_parity], dtype=float)
    if not np.isfinite(geometry).all():
        raise ValueError("Age geometry must be finite")
    J = int(P.J)
    if (float(J) != float(P.J) or J < 7 or float(P.age_start) != 18.0
            or float(P.da) != 4.0 or float(P.period_years) != 4.0
            or float(P.n_parity) != 4.0):
        raise ValueError("Observer requires four-year cells from age 18 and literal parity 0/1/2/3+")
    if (not bool(getattr(P, "sequential_births", False))
            or str(getattr(P, "fertility_units", "")) != "literal_topcode"
            or str(getattr(P, "child_state_mode", "")) != "independent_count"
            or any(bool(getattr(P, name, False)) for name in
                   ("joint_nested_choice", "fertility_nest_choice", "two_shock_choice"))):
        raise ValueError("Observer requires the maintained sequential independent-child-count model")
    if float(P.A_f_start) != 1.0 or float(P.A_f_end) != 7.0:
        raise ValueError("First-birth timing requires fertility cells starting at ages 18 through 42")
    top_weight = float(P.tfr_top_bin_weight)
    if not math.isfinite(top_weight) or top_weight < 3.0:
        raise ValueError("Literal top-bin weight must be finite and at least three")

    pre = np.asarray(evaluation.g_pre, dtype=float)
    post = np.asarray(evaluation.g_post_fertility, dtype=float)
    if (pre.ndim != 7 or post.shape != pre.shape or pre.shape[3] != J
            or pre.shape[5] != 4 or pre.shape[6] < 4
            or any(size == 0 for size in pre.shape)):
        raise ValueError("Pre/post fertility distributions must share the seven-axis state shape")
    if (not np.isfinite(pre).all() or not np.isfinite(post).all()
            or np.any(pre < 0.0) or np.any(post < 0.0)):
        raise ValueError("Pre/post fertility mass must be finite and nonnegative")
    # Keep actual age masses, rather than equalizing ages or loading CPS weights.
    pre_parity = np.sum(pre, axis=(0, 1, 2, 4, 6))
    post_parity = np.sum(post, axis=(0, 1, 2, 4, 6))
    pre_age = pre_parity.sum(axis=1)
    post_age = post_parity.sum(axis=1)
    if not np.isfinite(pre_age).all() or not np.isfinite(post_age).all():
        raise ValueError("Aggregated fertility mass is nonfinite")
    mass_error = _maximum_error(pre_age, post_age)
    if mass_error > MASS_ATOL:
        raise RuntimeError(f"Fertility step changes age mass: {mass_error:.3e}")
    if float(pre_age.sum()) <= 0.0:
        raise ValueError("Initial population mass must be positive")

    from run_e5f_transition_calibration import (
        first_birth_accounting_by_age,
        period_fertility_diagnostics,
    )

    accounting = first_birth_accounting_by_age(evaluation, P)
    period = period_fertility_diagnostics(evaluation, P)
    first_flow = _finite_vector(accounting["flow"], "first-birth flow", J)
    at_risk = _finite_vector(accounting["at_risk"], "first-birth risk set", J)
    hazard = _finite_vector(accounting["hazard"], "first-birth hazard", J)
    flow_fields = ("birth_flow_first", "birth_flow_second", "birth_flow_third_bin_entry")
    flows = np.column_stack([_finite_vector(period[key], key, J) for key in flow_fields])
    if (np.any(first_flow < 0.0) or np.any(at_risk < 0.0) or np.any(flows < 0.0)
            or np.any(hazard < 0.0) or np.any(hazard > 1.0)):
        raise RuntimeError("Birth flows, risk sets, and hazards must be nonnegative probabilities/masses")
    if (np.any(at_risk > pre_parity[:, 0] + FLOW_ATOL)
            or np.any(first_flow > at_risk + FLOW_ATOL)
            or np.any(flows > pre_parity[:, :3] + FLOW_ATOL)):
        raise RuntimeError("A parity transition exceeds its pre-fertility risk-set mass")
    if _maximum_error(first_flow, at_risk * hazard) > FLOW_ATOL:
        raise RuntimeError("First-birth flow differs from risk set times hazard")
    if _maximum_error(first_flow, flows[:, 0]) > FLOW_ATOL:
        raise RuntimeError("Existing first-birth flow primitives disagree")
    inferred_flows = np.cumsum(pre_parity - post_parity, axis=1)[:, :3]
    parity_error = _maximum_error(inferred_flows, flows)
    if parity_error > FLOW_ATOL:
        raise RuntimeError(f"Pre/post parity stocks disagree with dated birth flows: {parity_error:.3e}")
    ages = float(P.age_start) + np.arange(J, dtype=float) * float(P.da)
    midpoints = ages + 0.5 * float(P.da)
    if (_maximum_error(_finite_vector(period["age_cell_start"], "period age starts", J), ages) > MOMENT_ATOL
            or _maximum_error(_finite_vector(period["age_cell_midpoint"], "period age midpoints", J), midpoints) > MOMENT_ATOL
            or _maximum_error(_finite_vector(period["age_mass"], "period age mass", J), pre_age) > MASS_ATOL):
        raise RuntimeError("Period primitive uses different age geometry or population mass")
    if np.any(flows[7:, :] > FLOW_ATOL):
        raise RuntimeError("Birth flow occurs outside the declared fertile age cells")
    total_births = float(flows.sum())
    recorded_births = float(evaluation.births)
    if (not math.isfinite(recorded_births) or recorded_births < 0.0
            or abs(total_births - recorded_births) > FLOW_ATOL):
        raise RuntimeError("Explicit parity flows do not reproduce evaluation births")
    first_total = float(first_flow.sum())
    if first_total <= 0.0:
        raise ValueError("First-birth timing denominator must be positive")
    mean_age = float(np.dot(midpoints, first_flow) / first_total)
    share30 = float(first_flow[ages >= 30.0].sum() / first_total)
    helper_mean = float(period["period_first_birth_mean_age"])
    helper_share = float(period["period_first_birth_share_age30plus"])
    if (not math.isfinite(helper_mean) or not math.isfinite(helper_share)
            or abs(helper_mean - mean_age) > MOMENT_ATOL
            or abs(helper_share - share30) > MOMENT_ATOL):
        raise RuntimeError("Period timing scalars disagree with midpoint-weighted first-birth flows")

    overlap_left = np.maximum(ages, 40.0)
    overlap_right = np.minimum(ages + float(P.da), 45.0)
    lengths = np.maximum(overlap_right - overlap_left, 0.0)
    overlap_weights = lengths / float(P.da)
    selected = lengths > 0.0
    if not math.isclose(float(lengths.sum()), 5.0, rel_tol=0.0, abs_tol=MOMENT_ATOL):
        raise ValueError("Age cells do not cover the entire CPS [40,45) interval")
    post_shares = np.ones(J, dtype=float)
    if age_projection == "uniform_birth_time":
        post_shares[selected] = (
            (overlap_left[selected] + overlap_right[selected]) / 2.0 - ages[selected]
        ) / float(P.da)
    projected = (1.0 - post_shares[:, None]) * pre_parity + post_shares[:, None] * post_parity
    window_parity_mass = (overlap_weights[:, None] * projected).sum(axis=0)
    window_mass = float(window_parity_mass.sum())
    mother_mass = float(window_parity_mass[1:].sum())
    if window_mass <= 0.0 or mother_mass <= 0.0:
        raise ValueError("CPS projected total and mother denominators must be positive")
    parity = window_parity_mass / window_mass
    if not np.isfinite(parity).all() or np.any(parity < 0.0) or np.any(parity > 1.0):
        raise RuntimeError("Projected parity shares are invalid")
    return {
        "moments": {
            "childless_rate_40_44": float(parity[0]),
            "exactly_one_among_mothers_40_44": float(window_parity_mass[1] / mother_mass),
            "period_mean_age_first_birth": helper_mean,
            "period_share_first_births_age30plus": helper_share,
        },
        "parity_shares_40_44": {label: float(value) for label, value in
                                zip(("0", "1", "2", "3plus"), parity, strict=True)},
        "accounting": {
            "window_parity_mass": window_parity_mass.tolist(),
            "window_population_mass": window_mass,
            "window_mother_mass": mother_mass,
            "age_cell_start": ages.tolist(),
            "age_cell_midpoint": midpoints.tolist(),
            "pre_age_mass": pre_age.tolist(),
            "pre_parity_mass_by_age": pre_parity.tolist(),
            "post_parity_mass_by_age": post_parity.tolist(),
            "overlap_weights": overlap_weights.tolist(),
            "post_parity_interpolation_share": [float(post_shares[j]) if selected[j] else None for j in range(J)],
            "parity_birth_flows_by_age": flows.tolist(),
            "first_birth_flow": first_total,
            "explicit_birth_flow": total_births,
            "maximum_age_mass_error": mass_error,
            "maximum_parity_flow_error": parity_error,
        },
        "metadata": {
            "observer_id": "initial_cps_nchs_fertility_diagnostic_v1",
            "diagnostic_only": True,
            "production_smm_eligible": False,
            "age_projection": age_projection,
            "cps_age_window": "integer ages 40-44, interpreted as [40,45)",
            "age_weights": "model age masses times interval overlap; no empirical age reweighting",
            "stock_phase": "pre/post fertility of the same evaluation; no aging or cohort splicing",
            "birth_time_assumption": ("one parity transition per cell, uniformly distributed within the four-year age interval"
                                      if age_projection == "uniform_birth_time" else
                                      "post-fertility parity stock held constant throughout each age cell"),
            "timing_operator": "first-birth flow weighted midpoint age; cell start >=30 for the 30+ share",
            "population_approximation": "model household reproductive member proxies for the maternal population; not certified female exposure",
            "empirical_boundary_rule": "NCHS ages below18 and above45 are collapsed into the first/last fertile cell; no new maternal-age states",
            "stationarity_certified_by_observer": False,
            "initial_normalization": "2.1 remains separate; no fertility normalization is performed or replaced",
            "weights_or_standard_errors_adopted": False,
            "mass_atol": MASS_ATOL,
            "flow_atol": FLOW_ATOL,
            "moment_atol": MOMENT_ATOL,
        },
    }
