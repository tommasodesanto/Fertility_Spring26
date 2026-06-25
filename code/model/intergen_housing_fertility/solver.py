"""Solver for the intergenerational housing/fertility lifecycle model."""

from __future__ import annotations

import copy
import math
import time
from types import SimpleNamespace
from typing import Any

import numpy as np

from .parameters import apply_overrides, finalize_location_choice_spec, setup_parameters
from .kernels import (
    NUMBA_AVAILABLE,
    eval_owner_block_kernel,
    eval_renter_block_kernel,
    full_owner_block_kernel,
    full_renter_block_kernel,
    golden_owner_kernel,
    golden_renter_kernel,
    forward_distribution_fast_kernel,
    interp_cols_kernel,
    interp_cols_preidx_kernel,
    location_logit_kernel,
    scatter_cols_kernel,
    scatter_cols_sameidx_kernel,
    scatter_vec_kernel,
    tenure_choice_kernel,
    tenure_logit_kernel,
)
from .utils import (
    flat_nc,
    interp_cols,
    interp_indices,
    interp_on_grid,
    interp_vector,
    logsumexp,
    make_grid,
    scatter_redistribute,
    scatter_redistribute_cols,
    scatter_redistribute_cols_sameidx,
    unflat_nc,
    weighted_median_from_cells,
)


OUTSIDE_OPTION_CLOSURES = {"outside_option", "outside_option_local_births", "local_births_outside", "open_city"}
ACCOUNTING_SCALE_PRICE_CLOSURES = {"accounting_scale_prices", "scaled_housing", "scaled_housing_accounting"}
BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES = {
    "outside_option_benchmark_normalized",
    "benchmark_outside_option_normalized",
}
RENEWAL_VALVE_FIXED_CLOSURES = {"renewal_valve", "renewal_scale", "demographic_valve"}
RENEWAL_VALVE_CALIBRATED_CLOSURES = {
    "renewal_valve_calibrated",
    "renewal_calibrated",
    "demographic_valve_calibrated",
    "benchmark_renewal_valve",
}
RENEWAL_VALVE_CLOSURES = RENEWAL_VALVE_FIXED_CLOSURES | RENEWAL_VALVE_CALIBRATED_CLOSURES


def uses_outside_option_closure(P: SimpleNamespace) -> bool:
    return str(getattr(P, "population_closure", "normalized")).lower() in OUTSIDE_OPTION_CLOSURES


def uses_accounting_scale_price_closure(P: SimpleNamespace) -> bool:
    mode = str(getattr(P, "population_closure", "normalized")).lower()
    return mode in ACCOUNTING_SCALE_PRICE_CLOSURES or mode in BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES


def uses_benchmark_normalized_outside_closure(P: SimpleNamespace) -> bool:
    return str(getattr(P, "population_closure", "normalized")).lower() in BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES


def uses_renewal_valve_closure(P: SimpleNamespace) -> bool:
    return str(getattr(P, "population_closure", "normalized")).lower() in RENEWAL_VALVE_CLOSURES


def uses_calibrated_renewal_valve_closure(P: SimpleNamespace) -> bool:
    mode = str(getattr(P, "population_closure", "normalized")).lower()
    return mode in RENEWAL_VALVE_CALIBRATED_CLOSURES or bool(getattr(P, "renewal_calibrate_outside_flow", False))


def normalize_population_mass(P: SimpleNamespace) -> bool:
    if uses_outside_option_closure(P):
        return False
    return bool(getattr(P, "normalize_population_mass", True))


def uses_income_types(P: SimpleNamespace) -> bool:
    return bool(getattr(P, "use_income_types", False)) and len(np.atleast_1d(getattr(P, "z_grid", [1.0]))) > 1


def uses_markov_income(P: SimpleNamespace) -> bool:
    transition = str(getattr(P, "income_type_transition", "permanent")).lower()
    return uses_income_types(P) and transition in {"markov", "stochastic", "persistent"}


def income_type_values(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(getattr(P, "z_grid", [1.0]), dtype=float).reshape(-1)
    weights = np.asarray(getattr(P, "z_weights", np.ones(len(z))), dtype=float).reshape(-1)
    if weights.size != z.size:
        weights = np.ones(z.size)
    weights = np.maximum(weights, 0.0)
    weights = weights / weights.sum() if weights.sum() > 0 else np.ones(z.size) / max(z.size, 1)
    return z, weights


def income_transition_values(P: SimpleNamespace) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z, weights = income_type_values(P)
    Pi = np.asarray(getattr(P, "Pi_z", np.eye(len(z))), dtype=float)
    if Pi.shape != (len(z), len(z)):
        Pi = np.eye(len(z))
    Pi = np.maximum(Pi, 0.0)
    row_sum = Pi.sum(axis=1)
    for row in range(Pi.shape[0]):
        if row_sum[row] > 0:
            Pi[row, :] /= row_sum[row]
        else:
            Pi[row, :] = weights
    return z, weights, Pi


def income_at_state(P: SimpleNamespace, i: int, j: int, z_value: float) -> float:
    y = float(P.income[i, j])
    if j < int(getattr(P, "J_R", P.J)):
        return y * float(z_value)
    return y


def housing_demand_normalizer(P: SimpleNamespace) -> float:
    if normalize_population_mass(P):
        return max(float(getattr(P, "N_target", 1.0)), 1e-12)
    return 1.0


def entry_values_by_location(V: np.ndarray, b_grid: np.ndarray, P: SimpleNamespace) -> np.ndarray:
    bg = np.asarray(b_grid, dtype=float).reshape(-1)
    b_entry = float(np.clip(getattr(P, "b_entry_fixed", 0.0), bg[0], bg[-1]))
    entry_values = np.empty(P.I)
    for i in range(P.I):
        entry_values[i] = np.interp(b_entry, bg, V[:, 0, i, 0, 0, 0])
    return entry_values


def outside_option_entry_target(
    V: np.ndarray, b_grid: np.ndarray, stats: SimpleNamespace, P: SimpleNamespace
) -> tuple[np.ndarray, dict[str, Any]]:
    entry_values = entry_values_by_location(V, b_grid, P)
    outside_value = float(getattr(P, "outside_value", 0.0))
    kappa = max(float(getattr(P, "kappa_entry", getattr(P, "kappa_loc", 1.0))), 1e-10)
    values = np.concatenate([entry_values, [outside_value]])
    z = values / kappa
    z = z - np.max(z)
    probs = np.exp(z)
    probs = probs / np.sum(probs)
    city_probs = probs[:-1]
    outside_prob = float(probs[-1])

    mature_flow = max(float(getattr(stats, "entrants_mature_total", 0.0)), 0.0)
    outside_flow = max(float(getattr(P, "outside_entry_flow", getattr(P, "E_total", 0.0))), 0.0)
    local_weight = max(float(getattr(P, "local_birth_entry_weight", 1.0)), 0.0)
    potential_mass = outside_flow + local_weight * mature_flow
    entry_by_loc = potential_mass * city_probs
    info = {
        "entry_values": entry_values.copy(),
        "city_entry_prob": city_probs.copy(),
        "outside_entry_prob": outside_prob,
        "potential_entrant_mass": float(potential_mass),
        "outside_entry_mass": float(potential_mass * outside_prob),
        "mature_cityborn_flow": mature_flow,
    }
    return entry_by_loc, info


def city_entry_probabilities(entry_values: np.ndarray, outside_value: float, kappa_entry: float) -> tuple[np.ndarray, float]:
    kappa = max(float(kappa_entry), 1e-10)
    values = np.concatenate([np.asarray(entry_values, dtype=float).reshape(-1), [float(outside_value)]])
    z = values / kappa
    z = z - np.max(z)
    probs = np.exp(z)
    probs = probs / np.sum(probs)
    return probs[:-1], float(probs[-1])


def calibrate_outside_value_for_entry_total(
    entry_values: np.ndarray,
    potential_entrant_mass: float,
    target_city_entry_total: float,
    kappa_entry: float,
) -> float:
    """Outside value that matches a target city entry mass at fixed entry values."""

    potential_mass = max(float(potential_entrant_mass), 1e-14)
    target_total = float(np.clip(target_city_entry_total, 1e-14, potential_mass - 1e-14))
    kappa = max(float(kappa_entry), 1e-10)
    z = np.asarray(entry_values, dtype=float).reshape(-1) / kappa
    zmax = float(np.max(z))
    log_city_sum = zmax + math.log(float(np.sum(np.exp(z - zmax))))
    target_city_prob = target_total / potential_mass
    log_outside = log_city_sum + math.log1p(-target_city_prob) - math.log(target_city_prob)
    return float(kappa * log_outside)


def accounting_population_scale(
    sol: SimpleNamespace,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    *,
    outside_entry_flow: float | None = None,
    local_birth_entry_weight: float | None = None,
    outside_value: float | None = None,
    kappa_entry: float | None = None,
    target_city_entry_total: float | None = None,
    calibrate_outside_value: bool = False,
) -> SimpleNamespace:
    """Fast stationary population accounting at fixed prices and policies.

    With fixed policies, mature city-born children are linear in total city
    entry. If the outside option keeps city-entry probabilities fixed, the
    stationary entry level solves a scalar fixed point:

        E = q [M^O + lambda_B b E],

    where q is the probability a potential entrant chooses one of the modeled
    city locations and b is mature city-born children per unit of city entry.
    """

    entry_values = entry_values_by_location(sol.V, b_grid, P)
    ref_entry_total = float(getattr(sol, "entry_rate", getattr(P, "E_total", 0.0)))
    ref_entry_total = max(ref_entry_total, 1e-14)
    ref_mature_flow = max(float(getattr(sol, "entrants_mature_total", 0.0)), 0.0)
    mature_per_entry = ref_mature_flow / ref_entry_total

    outside_flow = (
        float(outside_entry_flow)
        if outside_entry_flow is not None
        else float(getattr(P, "outside_entry_flow", getattr(P, "E_total", ref_entry_total)))
    )
    local_weight = (
        float(local_birth_entry_weight)
        if local_birth_entry_weight is not None
        else float(getattr(P, "local_birth_entry_weight", 1.0))
    )
    kappa = (
        float(kappa_entry)
        if kappa_entry is not None
        else float(getattr(P, "kappa_entry", getattr(P, "kappa_loc", 1.0)))
    )

    outside_value_supplied = outside_value is not None
    if outside_value is None and not calibrate_outside_value and hasattr(P, "outside_value"):
        outside_value = float(P.outside_value)

    calibrated_outside_value = bool(
        calibrate_outside_value or outside_value_supplied or getattr(P, "outside_value_is_calibrated", False)
    )
    if outside_value is None:
        target_entry = ref_entry_total if target_city_entry_total is None else float(target_city_entry_total)
        potential_at_ref = outside_flow + local_weight * ref_mature_flow
        outside_value = calibrate_outside_value_for_entry_total(entry_values, potential_at_ref, target_entry, kappa)
        calibrated_outside_value = True

    city_probs, outside_prob = city_entry_probabilities(entry_values, float(outside_value), kappa)
    city_prob_total = float(np.sum(city_probs))
    denominator = 1.0 - city_prob_total * local_weight * mature_per_entry
    finite_stationary_scale = denominator > 1e-10
    if finite_stationary_scale:
        entry_total = city_prob_total * outside_flow / denominator
        potential_mass = outside_flow + local_weight * mature_per_entry * entry_total
        entry_by_loc = potential_mass * city_probs
        outside_entry_mass = potential_mass * outside_prob
        scale_factor = entry_total / ref_entry_total
    else:
        entry_total = np.inf
        potential_mass = np.inf
        entry_by_loc = np.full_like(city_probs, np.inf)
        outside_entry_mass = np.inf
        scale_factor = np.inf

    if finite_stationary_scale:
        entry_residual = entry_total - city_prob_total * potential_mass
        scale_residual = scale_factor * ref_entry_total - city_prob_total * (
            outside_flow + local_weight * scale_factor * ref_mature_flow
        )
        relative_residual = abs(entry_residual) / max(abs(entry_total), abs(city_prob_total * potential_mass), 1e-14)
        conditional_entry_shares = city_probs / max(city_prob_total, 1e-14)
    else:
        entry_residual = np.nan
        scale_residual = np.nan
        relative_residual = np.nan
        conditional_entry_shares = np.full_like(city_probs, np.nan)

    base_hd = np.asarray(getattr(sol, "housing_demand", np.zeros(P.I)), dtype=float) * housing_demand_normalizer(P)
    return SimpleNamespace(
        finite_stationary_scale=bool(finite_stationary_scale),
        outside_value=float(outside_value),
        calibrated_outside_value=bool(calibrated_outside_value),
        target_city_entry_total=(
            float(target_city_entry_total) if target_city_entry_total is not None else float(ref_entry_total)
        ),
        outside_entry_flow=float(outside_flow),
        local_birth_entry_weight=float(local_weight),
        kappa_entry=float(kappa),
        entry_values=entry_values,
        city_entry_prob=city_probs,
        conditional_entry_shares=conditional_entry_shares,
        outside_entry_prob=float(outside_prob),
        city_entry_prob_total=city_prob_total,
        reference_entry_total=ref_entry_total,
        reference_mature_cityborn_flow=ref_mature_flow,
        mature_cityborn_per_entry=float(mature_per_entry),
        denominator=float(denominator),
        stationary_entry_residual=float(entry_residual),
        stationary_scale_residual=float(scale_residual),
        stationary_entry_relative_residual=float(relative_residual),
        implied_entry_total=float(entry_total),
        implied_entry_by_loc=entry_by_loc,
        implied_potential_entrant_mass=float(potential_mass),
        implied_outside_entry_mass=float(outside_entry_mass),
        scale_factor=float(scale_factor),
        reference_total_population=float(getattr(sol, "total_mass", np.nan)),
        implied_total_population=float(getattr(sol, "total_mass", np.nan) * scale_factor),
        reference_housing_demand=base_hd,
        implied_housing_demand=base_hd * scale_factor,
    )


def benchmark_normalized_outside_population_scale(
    sol: SimpleNamespace,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    *,
    target_total_population: float | None = None,
    target_city_entry_prob: float | None = None,
    calibrate_outside_value_to_q: bool | None = None,
) -> SimpleNamespace:
    """Benchmark outside-option closure with the baseline scale normalized.

    The benchmark solves a normalized stationary composition and treats city
    scale as the accounting object. Given

        S E_0(p) = q^E(p) [M + S B_0(p)],

    the benchmark normalization sets ``S`` and computes the residual
    outside-origin potential entrant mass

        M = S * (E_0(p) / q^E(p) - B_0(p)).

    Counterfactuals should then hold ``M`` and the outside value fixed.
    """

    ref_pop = max(float(getattr(sol, "total_mass", getattr(P, "N_target", 1.0))), 1e-14)
    entry_flow = max(float(getattr(sol, "entry_rate", getattr(P, "E_total", 0.0))), 0.0)
    mature_flow = max(float(getattr(sol, "entrants_mature_total", 0.0)), 0.0)
    entry_per_scale = entry_flow / ref_pop
    mature_per_scale = mature_flow / ref_pop
    local_weight = max(float(getattr(P, "local_birth_entry_weight", 1.0)), 0.0)
    target_pop = (
        float(target_total_population)
        if target_total_population is not None
        else float(getattr(P, "outside_benchmark_target_total_population", getattr(P, "N_target", ref_pop)))
    )
    target_pop = max(target_pop, 1e-14)

    entry_values = entry_values_by_location(sol.V, b_grid, P)
    kappa = max(float(getattr(P, "kappa_entry", getattr(P, "kappa_loc", 1.0))), 1e-10)
    do_q_calibration = (
        bool(calibrate_outside_value_to_q)
        if calibrate_outside_value_to_q is not None
        else bool(getattr(P, "calibrate_outside_value_to_entry_prob", False))
    )
    target_q = (
        float(target_city_entry_prob)
        if target_city_entry_prob is not None
        else float(getattr(P, "target_city_entry_prob", getattr(P, "outside_target_city_entry_prob", np.nan)))
    )
    if do_q_calibration:
        if not np.isfinite(target_q):
            raise ValueError("calibrate_outside_value_to_entry_prob=True requires target_city_entry_prob")
        target_q = float(np.clip(target_q, 1e-10, 1.0 - 1e-10))
        outside_value = calibrate_outside_value_for_entry_total(entry_values, 1.0, target_q, kappa)
        calibrated_outside_value = True
    else:
        outside_value = float(getattr(P, "outside_value", 0.0))
        calibrated_outside_value = bool(getattr(P, "outside_value_is_calibrated", False))

    city_probs, outside_prob = city_entry_probabilities(entry_values, outside_value, kappa)
    city_prob_total = float(np.sum(city_probs))
    denominator = entry_per_scale - city_prob_total * local_weight * mature_per_scale
    outside_flow = target_pop * (entry_per_scale / max(city_prob_total, 1e-14) - local_weight * mature_per_scale)
    finite_stationary_scale = city_prob_total > 1e-12 and denominator > 1e-12 and outside_flow >= 0.0

    ref_entry_by_loc = np.asarray(
        getattr(sol, "entry_by_loc", getattr(P, "entry_by_loc", np.zeros(P.I))),
        dtype=float,
    ).reshape(-1)
    if ref_entry_by_loc.size != P.I or float(np.sum(ref_entry_by_loc)) <= 0:
        ref_entry_by_loc = entry_flow * np.ones(P.I) / P.I
    entry_by_loc_per_scale = ref_entry_by_loc / ref_pop

    if finite_stationary_scale:
        implied_total_population = target_pop
        scale_factor = implied_total_population / ref_pop
        implied_entry_total = implied_total_population * entry_per_scale
        implied_mature_cityborn_flow = implied_total_population * mature_per_scale
        potential_mass = outside_flow + local_weight * implied_mature_cityborn_flow
        implied_entry_by_loc = potential_mass * city_probs
        implied_outside_entry_mass = potential_mass * outside_prob
        stationary_residual = implied_entry_total - city_prob_total * potential_mass
        location_required_entry = implied_total_population * entry_by_loc_per_scale
        location_residual = location_required_entry - implied_entry_by_loc
        location_relative_residual = np.abs(location_residual) / np.maximum(
            np.maximum(np.abs(location_required_entry), np.abs(implied_entry_by_loc)),
            1e-14,
        )
        relative_residual = abs(stationary_residual) / max(
            abs(implied_entry_total), abs(city_prob_total * potential_mass), 1e-14
        )
        conditional_entry_shares = city_probs / max(city_prob_total, 1e-14)
    else:
        implied_total_population = np.inf
        scale_factor = np.inf
        implied_entry_total = np.inf
        implied_mature_cityborn_flow = np.inf
        potential_mass = np.inf
        implied_entry_by_loc = np.full_like(city_probs, np.inf)
        implied_outside_entry_mass = np.inf
        stationary_residual = np.nan
        location_required_entry = np.full_like(city_probs, np.nan)
        location_residual = np.full_like(city_probs, np.nan)
        location_relative_residual = np.full_like(city_probs, np.nan)
        relative_residual = np.nan
        conditional_entry_shares = np.full_like(city_probs, np.nan)

    base_hd = np.asarray(getattr(sol, "housing_demand", np.zeros(P.I)), dtype=float) * housing_demand_normalizer(P)
    return SimpleNamespace(
        finite_stationary_scale=bool(finite_stationary_scale),
        benchmark_normalized=True,
        outside_value=float(outside_value),
        calibrated_outside_value=bool(calibrated_outside_value),
        target_city_entry_prob=(float(target_q) if np.isfinite(target_q) else np.nan),
        outside_entry_flow=float(outside_flow),
        local_birth_entry_weight=float(local_weight),
        kappa_entry=float(kappa),
        entry_values=entry_values,
        city_entry_prob=city_probs,
        conditional_entry_shares=conditional_entry_shares,
        outside_entry_prob=float(outside_prob),
        city_entry_prob_total=city_prob_total,
        reference_total_population=float(ref_pop),
        target_total_population=float(target_pop),
        reference_entry_total=float(entry_flow),
        reference_entry_by_loc=ref_entry_by_loc,
        reference_mature_cityborn_flow=float(mature_flow),
        entry_per_unit_scale=float(entry_per_scale),
        entry_by_loc_per_unit_scale=entry_by_loc_per_scale,
        mature_cityborn_per_unit_scale=float(mature_per_scale),
        mature_cityborn_per_entry=float(mature_flow / max(entry_flow, 1e-14)),
        denominator=float(denominator),
        stationary_entry_residual=float(stationary_residual),
        stationary_scale_residual=float(stationary_residual),
        stationary_entry_relative_residual=float(relative_residual),
        stationary_entry_by_loc_residual=location_residual,
        stationary_entry_by_loc_relative_residual=location_relative_residual,
        implied_entry_total=float(implied_entry_total),
        implied_entry_by_loc=implied_entry_by_loc,
        implied_potential_entrant_mass=float(potential_mass),
        implied_outside_entry_mass=float(implied_outside_entry_mass),
        implied_mature_cityborn_flow=float(implied_mature_cityborn_flow),
        scale_factor=float(scale_factor),
        implied_total_population=float(implied_total_population),
        reference_housing_demand=base_hd,
        implied_housing_demand=base_hd * scale_factor,
    )


def renewal_population_scale(
    sol: SimpleNamespace,
    P: SimpleNamespace,
    b_grid: np.ndarray | None = None,
    *,
    outside_entry_flow: float | None = None,
    renewal_retention: float | None = None,
    calibrate_outside_flow: bool | None = None,
    target_total_population: float | None = None,
) -> SimpleNamespace:
    """Stationary renewal-valve city scale at fixed policies.

    The normalized distribution supplies two per-unit-scale flows:

        E0: young-adult entry flow required to sustain one unit of city scale.
        B0: mature city-born children per unit of city scale.

    A reduced-form valve supplies an exogenous outside-born flow M and retains a
    fraction rho of mature city-born children. Stationarity requires

        S * E0 = M + rho * S * B0,

    so the finite stationary scale is S = M / (E0 - rho * B0).
    """

    ref_pop = max(float(getattr(sol, "total_mass", getattr(P, "N_target", 1.0))), 1e-14)
    entry_flow = max(float(getattr(sol, "entry_rate", getattr(P, "E_total", 0.0))), 0.0)
    mature_flow = max(float(getattr(sol, "entrants_mature_total", 0.0)), 0.0)
    entry_per_scale = entry_flow / ref_pop
    mature_per_scale = mature_flow / ref_pop

    rho = (
        float(renewal_retention)
        if renewal_retention is not None
        else float(getattr(P, "renewal_retention", 1.0))
    )
    rho = max(rho, 0.0)
    denominator = entry_per_scale - rho * mature_per_scale
    target_pop = (
        float(target_total_population)
        if target_total_population is not None
        else float(getattr(P, "renewal_target_total_population", getattr(P, "N_target", ref_pop)))
    )
    do_calibrate = (
        bool(calibrate_outside_flow)
        if calibrate_outside_flow is not None
        else uses_calibrated_renewal_valve_closure(P)
    )
    if do_calibrate:
        outside_flow = target_pop * denominator
    else:
        outside_flow = (
            float(outside_entry_flow)
            if outside_entry_flow is not None
            else float(getattr(P, "outside_entry_flow", getattr(P, "E_total", entry_flow)))
        )
    finite_stationary_scale = denominator > 1e-12 and outside_flow >= 0.0

    if finite_stationary_scale:
        implied_total_population = target_pop if do_calibrate else outside_flow / denominator
        scale_factor = implied_total_population / ref_pop
        implied_entry_total = implied_total_population * entry_per_scale
        implied_mature_cityborn_flow = implied_total_population * mature_per_scale
    else:
        implied_total_population = np.inf
        scale_factor = np.inf
        implied_entry_total = np.inf
        implied_mature_cityborn_flow = np.inf

    mature_by_loc = np.asarray(getattr(sol, "entrants_mature_by_loc", np.zeros(P.I)), dtype=float).reshape(-1)
    if mature_by_loc.size != P.I:
        mature_by_loc = np.zeros(P.I)
    mature_by_loc_per_scale = mature_by_loc / ref_pop

    outside_shares = np.asarray(
        getattr(P, "outside_entry_shares", getattr(P, "entry_shares", np.ones(P.I) / P.I)),
        dtype=float,
    ).reshape(-1)
    if outside_shares.size != P.I or float(np.sum(np.maximum(outside_shares, 0.0))) <= 0:
        outside_shares = np.ones(P.I) / P.I
    outside_shares = np.maximum(outside_shares, 0.0)
    outside_shares = outside_shares / np.sum(outside_shares)

    if finite_stationary_scale:
        outside_entry_by_loc = outside_flow * outside_shares
        retained_cityborn_by_loc = rho * implied_total_population * mature_by_loc_per_scale
        implied_entry_by_loc = outside_entry_by_loc + retained_cityborn_by_loc
        entry_residual = implied_entry_total - float(np.sum(implied_entry_by_loc))
        scale_residual = (
            implied_total_population * entry_per_scale
            - (outside_flow + rho * implied_total_population * mature_per_scale)
        )
        relative_residual = abs(scale_residual) / max(
            abs(implied_total_population * entry_per_scale),
            abs(outside_flow + rho * implied_total_population * mature_per_scale),
            1e-14,
        )
        conditional_entry_shares = implied_entry_by_loc / max(float(np.sum(implied_entry_by_loc)), 1e-14)
    else:
        outside_entry_by_loc = np.full(P.I, np.nan)
        retained_cityborn_by_loc = np.full(P.I, np.nan)
        implied_entry_by_loc = np.full(P.I, np.inf)
        entry_residual = np.nan
        scale_residual = np.nan
        relative_residual = np.nan
        conditional_entry_shares = np.full(P.I, np.nan)

    base_hd = np.asarray(getattr(sol, "housing_demand", np.zeros(P.I)), dtype=float) * housing_demand_normalizer(P)
    return SimpleNamespace(
        finite_stationary_scale=bool(finite_stationary_scale),
        outside_entry_flow=float(outside_flow),
        calibrated_outside_entry_flow=bool(do_calibrate),
        target_total_population=float(target_pop),
        renewal_retention=float(rho),
        reference_total_population=float(ref_pop),
        reference_entry_total=float(entry_flow),
        reference_mature_cityborn_flow=float(mature_flow),
        entry_per_unit_scale=float(entry_per_scale),
        mature_cityborn_per_unit_scale=float(mature_per_scale),
        mature_cityborn_per_entry=float(mature_flow / max(entry_flow, 1e-14)),
        denominator=float(denominator),
        implied_total_population=float(implied_total_population),
        scale_factor=float(scale_factor),
        implied_entry_total=float(implied_entry_total),
        implied_mature_cityborn_flow=float(implied_mature_cityborn_flow),
        outside_entry_shares=outside_shares,
        outside_entry_by_loc=outside_entry_by_loc,
        retained_cityborn_by_loc=retained_cityborn_by_loc,
        implied_entry_by_loc=implied_entry_by_loc,
        conditional_entry_shares=conditional_entry_shares,
        stationary_entry_residual=float(entry_residual),
        stationary_scale_residual=float(scale_residual),
        stationary_entry_relative_residual=float(relative_residual),
        reference_housing_demand=base_hd,
        implied_housing_demand=base_hd * scale_factor,
    )


def set_entry_masses(P: SimpleNamespace, entry_by_loc: np.ndarray) -> tuple[np.ndarray, float]:
    entry = np.maximum(np.asarray(entry_by_loc, dtype=float).reshape(-1), 0.0)
    total = float(np.sum(entry))
    if total > 1e-14:
        shares = entry / total
    else:
        shares = np.ones(P.I) / P.I
    P.entry_by_loc = entry
    P.E_total = total
    P.entry_shares = shares
    return shares, total


def run_model_cp_dt(P_override: Any | None = None, verbose: bool = True) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    t_start = time.perf_counter()
    P = setup_parameters()
    if P_override is not None:
        P = apply_overrides(P, P_override)
    if not hasattr(P, "beta") or P.beta is None:
        P.beta = 1 / (1 + P.rho) if hasattr(P, "rho") else 0.96
    P.rho = 1 / P.beta - 1
    P.rho_hat = P.rho
    P.user_cost_rate = P.q + P.delta + P.tau_H
    P.R_gross = 1 + P.q
    P.phi = np.asarray(P.phi, dtype=float).reshape(-1)
    if P.phi.size == 1:
        P.phi = P.phi.item() * np.ones(P.n_parity)
    if hasattr(P, "entry_init_override") and P.entry_init_override is not None:
        e0 = np.maximum(np.asarray(P.entry_init_override, dtype=float).reshape(-1), 0)
        P.entry_shares = e0 / e0.sum()
        P.entry_by_loc = (1 / P.J) * P.entry_shares
    if not hasattr(P, "use_stochastic_aging"):
        P.use_stochastic_aging = False
    P = finalize_location_choice_spec(P)
    if uses_outside_option_closure(P):
        P.normalize_population_mass = False

    if verbose:
        print("=" * 60)
        print("  INTERGENERATIONAL HOUSING/FERTILITY MODEL")
        print(f"  J={P.J}, beta={P.beta:.4f}, R={P.R_gross:.4f}, Nb={P.Nb}")
        print(f"  alpha={P.alpha_cons:.2f}, sigma={P.sigma:.2f}, chi={P.chi:.2f}, hR_max={P.hR_max:.1f}")
        print(f"  markets={P.I}, kappa_fert={P.kappa_fert:.2f}, PTI={bool(getattr(P, 'use_pti_constraint', False))}")
        print("=" * 60)

    b_grid = make_grid(P)
    p_init = P.r_bar / P.user_cost_rate
    if hasattr(P, "p_init_override") and P.p_init_override is not None:
        p_init = np.asarray(P.p_init_override, dtype=float).reshape(-1)

    solve_mode = str(getattr(P, "solve_mode", "ge")).lower()
    do_pe = solve_mode in ("pe", "partial", "partial_equilibrium", "fixed")
    if uses_markov_income(P):
        if do_pe:
            sol, P, p_eq = solve_markov_income_partial_equilibrium(
                np.asarray(P.p_fixed, dtype=float).reshape(-1),
                P,
                b_grid,
                verbose=verbose,
            )
        else:
            sol, P, p_eq = solve_markov_income_equilibrium(p_init, P, b_grid, verbose=verbose)
    elif uses_income_types(P):
        if do_pe:
            sol, P, p_eq = solve_income_type_partial_equilibrium(
                np.asarray(P.p_fixed, dtype=float).reshape(-1),
                P,
                b_grid,
                verbose=verbose,
            )
        else:
            sol, P, p_eq = solve_income_type_equilibrium(p_init, P, b_grid, verbose=verbose)
    elif do_pe:
        sol, P, p_eq = solve_partial_equilibrium_dt(
            np.asarray(P.p_fixed, dtype=float).reshape(-1),
            np.asarray(P.w_fixed, dtype=float).reshape(-1),
            np.asarray(P.entry_shares_fixed, dtype=float).reshape(-1),
            P,
            b_grid,
            verbose=verbose,
        )
    else:
        sol, P, p_eq = solve_equilibrium(p_init, P, b_grid, verbose=verbose)

    if verbose:
        elapsed = time.perf_counter() - t_start
        price_text = ",".join(f"{x:.3f}" for x in np.asarray(p_eq).reshape(-1))
        print("\n" + "=" * 20 + " RESULTS " + "=" * 20)
        print(f"TFR={2 * sol.mean_parity:.2f}, Own={100 * sol.own_rate:.1f}%, Prices=[{price_text}]")
        print(f"Pop mass={getattr(sol, 'total_mass', np.nan):.3f}, Housing demand={np.sum(sol.housing_demand):.3f}")
        print(f"Total: {elapsed:.1f} sec")
    return sol, P, p_eq


def solve_markov_income_partial_equilibrium(
    p_fixed: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    p_eq = np.asarray(p_fixed, dtype=float).reshape(-1).copy()
    P.eq_iter = 1
    sol = solve_markov_income_at_prices(p_eq, P, b_grid, verbose=False, fast_stats=False)
    if verbose:
        print(
            f"  Markov income PE: own={100 * sol.own_rate:.1f}% "
            f"TFR={2 * sol.mean_parity:.2f} p={','.join(f'{x:.3f}' for x in p_eq)}"
        )
    return sol, P, p_eq


def solve_markov_income_equilibrium(
    p_init: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    p = np.asarray(p_init, dtype=float).reshape(-1)
    z_grid, z_weights, Pi_z = income_transition_values(P)
    lam = float(getattr(P, "lambda_eq", 0.30))
    best_err = np.inf
    best_p = p.copy()
    best_sol: SimpleNamespace | None = None
    t_solve = 0.0
    final_err = np.nan
    iterations_completed = 0

    if verbose:
        z_text = ", ".join(f"{z:.2f}:{w:.2f}" for z, w in zip(z_grid, z_weights))
        print(f"  Markov income states: {z_text}")
        print(f"  Markov transition rows: {np.array2string(Pi_z, precision=3)}")

    SD_shared = precompute_shared(P, b_grid)
    for it in range(1, int(P.max_iter_eq) + 1):
        P.eq_iter = it
        t0 = time.perf_counter()
        sol_it = solve_markov_income_at_prices(p, P, b_grid, verbose=False, fast_stats=True, SD=SD_shared)
        t_solve += time.perf_counter() - t0
        Hd = np.asarray(sol_it.housing_demand, dtype=float).reshape(-1)
        p_target = np.zeros(P.I)
        for i in range(P.I):
            hd = max(float(Hd[i]), float(P.housing_demand_floor_for_supply))
            p_target[i] = P.r_bar[i] * (hd / P.H0[i]) ** (1.0 / P.xi_supply[i]) / P.user_cost_rate
        err = float(np.max(np.abs(p_target - p) / np.maximum(np.abs(p), 1e-6)))
        final_err = err
        iterations_completed = it
        if err < best_err:
            best_err = err
            best_p = p.copy()
            best_sol = sol_it
        if verbose:
            print(
                f"  ZM{it:3d}: ep={err:.4f} own={100 * sol_it.own_rate:.1f}% "
                f"TFR={2 * sol_it.mean_parity:.2f} p={','.join(f'{x:.3f}' for x in p)}"
            )
        if err < P.tol_eq:
            break
        p = p + lam * (p_target - p)
        if bool(getattr(P, "enforce_price_bounds", True)):
            p = np.clip(p, P.p_min, P.p_max)

    if best_sol is None:
        best_sol = solve_markov_income_at_prices(best_p, P, b_grid, verbose=False, fast_stats=True, SD=SD_shared)
    scalar_refine_info: dict[str, Any] = {"used": False}
    if P.I == 1 and bool(getattr(P, "scalar_market_refine", True)):
        best_sol, best_p, best_err, scalar_refine_info = refine_one_market_markov_income(
            best_p,
            best_sol,
            best_err,
            P,
            b_grid,
            verbose=verbose,
            SD=SD_shared,
        )
    best_sol = solve_markov_income_at_prices(best_p, P, b_grid, verbose=False, fast_stats=False, SD=SD_shared)
    best_sol.timings = {
        **getattr(best_sol, "timings", {}),
        "income_process": "markov",
        "iterations_completed": int(iterations_completed),
        "best_eq_error": float(best_err),
        "damped_final_eq_error": float(final_err),
        "final_eq_error": float(best_err),
        "accepted": bool(best_err < P.tol_eq),
        "strict_converged": bool(best_err < P.tol_eq),
        "convergence_reason": "strict_tol" if best_err < P.tol_eq else "max_iter",
        "markov_income_solve_time": float(t_solve),
        "scalar_market_refine": scalar_refine_info,
    }
    best_sol.converged = bool(best_err < P.tol_eq)
    return best_sol, P, best_p


def refine_one_market_markov_income(
    best_p: np.ndarray,
    best_sol: SimpleNamespace,
    best_err: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
    SD: SimpleNamespace | None = None,
) -> tuple[SimpleNamespace, np.ndarray, float, dict[str, Any]]:
    p0 = float(np.asarray(best_p, dtype=float).reshape(-1)[0])
    p_min = float(getattr(P, "p_min", 1e-4))
    p_max = float(getattr(P, "p_max", 100.0))
    expand = max(float(getattr(P, "scalar_market_refine_expand", 1.5)), 1.05)
    max_expand = max(0, int(getattr(P, "scalar_market_refine_max_expand", 8)))
    max_iter = max(1, int(getattr(P, "scalar_market_refine_iter", 24)))
    tol = float(getattr(P, "tol_eq", 1e-4))
    eval_count = 0
    eval_time = 0.0

    def eval_price(price: float) -> tuple[float, float, SimpleNamespace]:
        nonlocal eval_count, eval_time
        t_eval = time.perf_counter()
        sol = solve_markov_income_at_prices(np.array([price]), P, b_grid, verbose=False, fast_stats=True, SD=SD)
        eval_time += time.perf_counter() - t_eval
        eval_count += 1
        demand = float(np.asarray(sol.housing_demand).reshape(-1)[0])
        supply = float(np.asarray(sol.housing_supply).reshape(-1)[0])
        excess = demand - supply
        metric = abs(excess) / max(abs(supply), 1e-12)
        return excess, metric, sol

    best_excess = float(np.asarray(best_sol.housing_demand).reshape(-1)[0] - np.asarray(best_sol.housing_supply).reshape(-1)[0])
    best_metric = float(best_err)
    best_price = p0
    sol_best = best_sol

    lo = max(p_min, p0 / expand)
    hi = min(p_max, p0 * expand)
    ex_lo, metric_lo, sol_lo = eval_price(lo)
    ex_hi, metric_hi, sol_hi = eval_price(hi)
    for price, excess, metric, sol in ((lo, ex_lo, metric_lo, sol_lo), (hi, ex_hi, metric_hi, sol_hi)):
        if metric < best_metric:
            best_price, best_excess, best_metric, sol_best = price, excess, metric, sol

    expansions = 0
    while ex_lo * ex_hi > 0 and expansions < max_expand:
        expansions += 1
        lo = max(p_min, lo / expand)
        hi = min(p_max, hi * expand)
        ex_lo, metric_lo, sol_lo = eval_price(lo)
        ex_hi, metric_hi, sol_hi = eval_price(hi)
        for price, excess, metric, sol in ((lo, ex_lo, metric_lo, sol_lo), (hi, ex_hi, metric_hi, sol_hi)):
            if metric < best_metric:
                best_price, best_excess, best_metric, sol_best = price, excess, metric, sol
        if lo <= p_min and hi >= p_max:
            break

    info: dict[str, Any] = {
        "used": True,
        "method": str(getattr(P, "scalar_market_refine_method", "brent")).lower(),
        "bracket_found": bool(ex_lo * ex_hi <= 0),
        "initial_price": p0,
        "best_price": float(best_price),
        "best_metric": float(best_metric),
        "best_excess": float(best_excess),
        "expansions": int(expansions),
        "iterations": 0,
        "price_evaluations": int(eval_count),
        "price_evaluation_time_sec": float(eval_time),
    }

    if ex_lo * ex_hi > 0:
        if verbose:
            print(f"  Scalar refine: no bracket; best residual={best_metric:.3e}")
        return sol_best, np.array([best_price]), best_metric, info

    method = info["method"]
    if method == "brent":
        a, b_pt, fa, fb = lo, hi, ex_lo, ex_hi
        if abs(fa) < abs(fb):
            a, b_pt, fa, fb = b_pt, a, fb, fa
        c_pt, fc = a, fa
        mflag = True
        d_old = 0.0
        for k in range(1, max_iter + 1):
            if fa != fc and fb != fc:
                mid = (
                    a * fb * fc / ((fa - fb) * (fa - fc))
                    + b_pt * fa * fc / ((fb - fa) * (fb - fc))
                    + c_pt * fa * fb / ((fc - fa) * (fc - fb))
                )
            else:
                mid = b_pt - fb * (b_pt - a) / (fb - fa)
            lo_gate = min((3.0 * a + b_pt) / 4.0, b_pt)
            hi_gate = max((3.0 * a + b_pt) / 4.0, b_pt)
            use_bisect = (
                not (lo_gate < mid < hi_gate)
                or (mflag and abs(mid - b_pt) >= abs(b_pt - c_pt) / 2.0)
                or (not mflag and abs(mid - b_pt) >= abs(c_pt - d_old) / 2.0)
                or (mflag and abs(b_pt - c_pt) < 1e-12)
                or (not mflag and abs(c_pt - d_old) < 1e-12)
            )
            if use_bisect:
                mid = 0.5 * (a + b_pt)
                mflag = True
            else:
                mflag = False
            ex_mid, metric_mid, sol_mid = eval_price(mid)
            if metric_mid < best_metric:
                best_price, best_excess, best_metric, sol_best = mid, ex_mid, metric_mid, sol_mid
            info["iterations"] = int(k)
            info["price_evaluations"] = int(eval_count)
            info["price_evaluation_time_sec"] = float(eval_time)
            if metric_mid < tol:
                break
            d_old = c_pt
            c_pt, fc = b_pt, fb
            if fa * ex_mid < 0:
                b_pt, fb = mid, ex_mid
            else:
                a, fa = mid, ex_mid
            if abs(fa) < abs(fb):
                a, b_pt, fa, fb = b_pt, a, fb, fa
    else:
        left, right = lo, hi
        f_left, f_right = ex_lo, ex_hi
        side = 0
        for k in range(1, max_iter + 1):
            if method == "illinois":
                denom = f_right - f_left
                if abs(denom) > 1e-300:
                    mid = right - f_right * (right - left) / denom
                else:
                    mid = 0.5 * (left + right)
                if not (left < mid < right):
                    mid = 0.5 * (left + right)
            else:
                mid = 0.5 * (left + right)
            ex_mid, metric_mid, sol_mid = eval_price(mid)
            if metric_mid < best_metric:
                best_price, best_excess, best_metric, sol_best = mid, ex_mid, metric_mid, sol_mid
            info["price_evaluations"] = int(eval_count)
            info["price_evaluation_time_sec"] = float(eval_time)
            if metric_mid < tol:
                info["iterations"] = int(k)
                break
            if f_left * ex_mid <= 0:
                right, f_right = mid, ex_mid
                if method == "illinois" and side == -1:
                    f_left *= 0.5
                side = -1
            else:
                left, f_left = mid, ex_mid
                if method == "illinois" and side == 1:
                    f_right *= 0.5
                side = 1
            info["iterations"] = int(k)

    info["best_price"] = float(best_price)
    info["best_metric"] = float(best_metric)
    info["best_excess"] = float(best_excess)
    info["price_evaluations"] = int(eval_count)
    info["price_evaluation_time_sec"] = float(eval_time)
    if verbose:
        print(f"  Scalar refine: residual={best_metric:.3e} p={best_price:.4f}")
    return sol_best, np.array([best_price]), best_metric, info


def solve_markov_income_at_prices(
    p_eq: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = False,
    fast_stats: bool = False,
    SD: SimpleNamespace | None = None,
) -> SimpleNamespace:
    p = np.asarray(p_eq, dtype=float).reshape(-1).copy()
    r = P.user_cost_rate * p
    if SD is None:
        SD = precompute_shared(P, b_grid)
    t0 = time.perf_counter()
    V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, btime = solve_bellman_full_markov_income(
        r, p, P, b_grid, SD
    )
    t_bellman = time.perf_counter() - t0
    t0 = time.perf_counter()
    g, stats = forward_distribution_markov_income(
        bp_pol, hR_pol, tc, lp_j, fp, r, p, P, b_grid, SD, fast_stats=fast_stats, tenure_probs=tp
    )
    t_dist = time.perf_counter() - t0
    if fast_stats:
        sol = pack_fast_solution_markov_income(stats, p, P)
    else:
        sol = pack_solution_markov_income(V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, g, stats, P.w_hat, p, P)
    sol.b_grid = np.asarray(b_grid, dtype=float).copy()
    sol.timings = {
        "bellman_full": float(btime.get("bellman", t_bellman)),
        "distribution": float(t_dist),
        "n_full": 1,
        "n_eval": 0,
        "n_dist": 1,
        "income_process": "markov",
        "bellman_mode": "full_only",
    }
    if verbose:
        print(
            f"  Markov income fixed-price solve: own={100 * sol.own_rate:.1f}% "
            f"TFR={2 * sol.mean_parity:.2f}"
        )
    return sol


def make_income_type_params(P: SimpleNamespace, z_value: float, z_weight: float) -> SimpleNamespace:
    Pz = copy.deepcopy(P)
    Pz.use_income_types = False
    Pz.Nz = 1
    Pz.z_grid = np.array([float(z_value)])
    Pz.z_weights = np.array([1.0])
    Pz.income_type_z = float(z_value)
    Pz.income_type_weight = float(z_weight)
    Pz.N_target = float(getattr(P, "N_target", 1.0)) * float(z_weight)
    Pz.w_hat = np.asarray(P.w_hat, dtype=float).reshape(-1) * float(z_value)
    Pz.entry_shares = np.asarray(P.entry_shares, dtype=float).reshape(-1).copy()
    Pz.entry_by_loc = float(getattr(Pz, "E_total", 1.0 / Pz.J)) * Pz.entry_shares
    Pz = apply_overrides(Pz, {"w_hat": Pz.w_hat, "entry_shares": Pz.entry_shares})
    Pz.N_target = float(getattr(P, "N_target", 1.0)) * float(z_weight)
    Pz.income_type_z = float(z_value)
    Pz.income_type_weight = float(z_weight)
    return Pz


def solve_income_type_partial_equilibrium(
    p_fixed: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    z_grid, z_weights = income_type_values(P)
    type_solutions = []
    for z_value, z_weight in zip(z_grid, z_weights):
        Pz = make_income_type_params(P, float(z_value), float(z_weight))
        sol_z, _, _ = solve_partial_equilibrium_dt(
            p_fixed,
            Pz.w_hat,
            Pz.entry_shares,
            Pz,
            b_grid,
            verbose=False,
        )
        type_solutions.append(sol_z)
    sol = aggregate_income_type_solutions(type_solutions, z_grid, z_weights, P, p_fixed)
    sol.timings = {
        "income_type_solves": len(type_solutions),
        "bellman_full": float(sum(getattr(s, "timings", {}).get("bellman_full", 0.0) for s in type_solutions)),
        "distribution": float(sum(getattr(s, "timings", {}).get("distribution", 0.0) for s in type_solutions)),
    }
    return sol, P, p_fixed


def solve_income_type_equilibrium(
    p_init: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    p = np.asarray(p_init, dtype=float).reshape(-1)
    z_grid, z_weights = income_type_values(P)
    lam = float(getattr(P, "lambda_eq", 0.30))
    best_err = np.inf
    best_p = p.copy()
    best_solutions: list[SimpleNamespace] = []
    t_solve = 0.0
    final_err = np.nan

    if verbose:
        z_text = ", ".join(f"{z:.2f}:{w:.2f}" for z, w in zip(z_grid, z_weights))
        print(f"  Income types: {z_text}")

    for it in range(1, int(P.max_iter_eq) + 1):
        type_solutions = []
        t0 = time.perf_counter()
        for z_value, z_weight in zip(z_grid, z_weights):
            Pz = make_income_type_params(P, float(z_value), float(z_weight))
            sol_z, _, _ = solve_partial_equilibrium_dt(
                p,
                Pz.w_hat,
                Pz.entry_shares,
                Pz,
                b_grid,
                verbose=False,
            )
            type_solutions.append(sol_z)
        t_solve += time.perf_counter() - t0
        sol_it = aggregate_income_type_solutions(type_solutions, z_grid, z_weights, P, p)
        Hd = np.asarray(sol_it.housing_demand, dtype=float).reshape(-1)
        p_target = np.zeros(P.I)
        for i in range(P.I):
            hd = max(float(Hd[i]), float(P.housing_demand_floor_for_supply))
            p_target[i] = P.r_bar[i] * (hd / P.H0[i]) ** (1.0 / P.xi_supply[i]) / P.user_cost_rate
        err = float(np.max(np.abs(p_target - p) / np.maximum(np.abs(p), 1e-6)))
        final_err = err
        if err < best_err:
            best_err = err
            best_p = p.copy()
            best_solutions = type_solutions
        if verbose:
            print(
                f"  Z{it:3d}: ep={err:.4f} own={100 * sol_it.own_rate:.1f}% "
                f"TFR={2 * sol_it.mean_parity:.2f} p={','.join(f'{x:.3f}' for x in p)}"
            )
        if err < P.tol_eq:
            break
        p = p + lam * (p_target - p)
        if bool(getattr(P, "enforce_price_bounds", True)):
            p = np.clip(p, P.p_min, P.p_max)

    sol = aggregate_income_type_solutions(best_solutions, z_grid, z_weights, P, best_p)
    sol.timings = {
        "income_type_solves": int(len(z_grid)),
        "iterations_completed": int(it),
        "best_eq_error": float(best_err),
        "final_eq_error": float(final_err),
        "accepted": bool(best_err < P.tol_eq),
        "strict_converged": bool(best_err < P.tol_eq),
        "convergence_reason": "strict_tol" if best_err < P.tol_eq else "max_iter",
        "income_type_solve_time": float(t_solve),
    }
    sol.converged = bool(best_err < P.tol_eq)
    return sol, P, best_p


def aggregate_income_type_solutions(
    type_solutions: list[SimpleNamespace],
    z_grid: np.ndarray,
    z_weights: np.ndarray,
    P: SimpleNamespace,
    p: np.ndarray,
) -> SimpleNamespace:
    if not type_solutions:
        raise ValueError("No income-type solutions to aggregate.")
    total_mass = sum(float(getattr(s, "total_mass", 0.0)) for s in type_solutions)
    total_mass = max(total_mass, 1e-12)
    weights = np.array([float(getattr(s, "total_mass", 0.0)) / total_mass for s in type_solutions])
    g_total = sum(s.g for s in type_solutions)

    sol = SimpleNamespace(
        type_solutions=type_solutions,
        type_values=np.asarray(z_grid, dtype=float).copy(),
        type_weights=np.asarray(z_weights, dtype=float).copy(),
        g=g_total,
        p_eq=np.asarray(p, dtype=float).reshape(-1),
        owner_asset_price=np.asarray(p, dtype=float).reshape(-1),
        owner_user_cost=P.user_cost_rate * np.asarray(p, dtype=float).reshape(-1),
    )
    sol.V = type_solutions[0].V
    sol.c_pol = type_solutions[0].c_pol
    sol.hR_pol = type_solutions[0].hR_pol
    sol.bp_pol = type_solutions[0].bp_pol
    sol.tenure_choice = type_solutions[0].tenure_choice
    sol.tenure_probs = getattr(type_solutions[0], "tenure_probs", None)
    sol.loc_probs = type_solutions[0].loc_probs
    sol.fert_probs = type_solutions[0].fert_probs
    sol.fert_value = type_solutions[0].fert_value

    for name in ("housing_demand", "pop_share", "own_by_age", "fert_by_age", "parity_dist", "own_by_parity", "mean_housing_by_parity"):
        vals = [np.asarray(getattr(s, name), dtype=float) for s in type_solutions if hasattr(s, name)]
        if vals:
            setattr(sol, name, sum(w * v for w, v in zip(weights, vals)))

    scalar_names = [
        "own_rate",
        "mean_parity",
        "young_own_rate",
        "old_age_own_rate_6575",
        "old_age_own_rate_parents_6575",
        "old_age_own_rate_childless_6575",
        "old_age_parent_childless_gap_6575",
        "mean_age_first_birth",
        "wealth_to_income",
        "liquid_wealth_to_income",
        "young_liquid_wealth_to_income",
        "own_rate_parents",
        "own_rate_childless",
        "own_family_gap",
        "housing_increment_1to2",
    ]
    for name in scalar_names:
        vals = np.array([float(getattr(s, name, np.nan)) for s in type_solutions])
        finite = np.isfinite(vals)
        setattr(sol, name, float(np.sum(weights[finite] * vals[finite]) / max(np.sum(weights[finite]), 1e-12)) if np.any(finite) else np.nan)

    sol.total_mass = float(total_mass)
    sol.owner_demand_by_size = sum(weights[k] * np.asarray(type_solutions[k].owner_demand_by_size) for k in range(len(type_solutions)))
    sol.rental_demand_by_market = sum(weights[k] * np.asarray(type_solutions[k].rental_demand_by_market) for k in range(len(type_solutions)))
    sol.owner_demand_by_market = sum(weights[k] * np.asarray(type_solutions[k].owner_demand_by_market) for k in range(len(type_solutions)))
    sol.rental_demand_by_size = sol.rental_demand_by_market.copy()
    sol.housing_supply = P.H0 * (sol.owner_user_cost / P.r_bar) ** P.xi_supply
    sol.aggregate_rental_demand = float(np.sum(sol.rental_demand_by_market))
    sol.aggregate_owner_demand = float(np.sum(sol.owner_demand_by_market))
    sol.aggregate_housing_demand = float(sol.aggregate_rental_demand + sol.aggregate_owner_demand)
    sol.aggregate_housing_supply = float(np.sum(sol.housing_supply))
    sol.aggregate_housing_excess = float(sol.aggregate_housing_demand - sol.aggregate_housing_supply)
    sol.best_max_abs_rel_excess = float(
        np.max(np.abs((sol.rental_demand_by_market + sol.owner_demand_by_market - sol.housing_supply) / np.maximum(sol.housing_supply, 1e-12)))
    )
    sol.best_market_metric = sol.best_max_abs_rel_excess
    sol.converged = bool(sol.best_max_abs_rel_excess <= getattr(P, "tol_eq", 1e-4))
    sol.young_owner_rate = float(getattr(sol, "young_own_rate", np.nan))
    sol.old_owner_rate = float(getattr(sol, "old_age_own_rate_6575", np.nan))
    sol.mean_completed_fertility = float(getattr(sol, "mean_parity", np.nan))
    sol.childless_rate = float(sol.parity_dist[0]) if hasattr(sol, "parity_dist") and len(sol.parity_dist) else np.nan
    sol.own_rate_by_income_type = np.array([float(getattr(s, "own_rate", np.nan)) for s in type_solutions])
    sol.mean_fertility_by_income_type = np.array([float(getattr(s, "mean_parity", np.nan)) for s in type_solutions])
    sol.housing_demand_by_income_type = np.array([np.asarray(getattr(s, "housing_demand", np.zeros(P.I)), dtype=float) for s in type_solutions])
    return sol


def solve_partial_equilibrium_dt(
    p_fixed: np.ndarray,
    w_fixed: np.ndarray,
    entry_shares_fixed: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    p_eq = p_fixed.copy()
    w = w_fixed.copy()
    entry_shares = np.maximum(entry_shares_fixed.copy(), 0)
    entry_shares = entry_shares / entry_shares.sum()
    P.w_hat = w
    P = apply_overrides(P, {"w_hat": w, "entry_shares": entry_shares})
    P.entry_by_loc = P.E_total * entry_shares
    P.eq_iter = 1
    r = P.user_cost_rate * p_eq
    SD = precompute_shared(P, b_grid)
    V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, timings = solve_bellman_full(r, p_eq, P, b_grid, SD)
    g, stats = forward_distribution(
        bp_pol, hR_pol, tc, lp_j, fp, r, p_eq, P, b_grid, SD, fast_stats=False, tenure_probs=tp
    )
    sol = pack_solution(V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, g, stats, P.w_hat, p_eq, P)
    sol.timings = {"bellman_full": timings["bellman"], "distribution": timings.get("distribution", 0.0)}
    return sol, P, p_eq


def solve_equilibrium(
    p_init: np.ndarray, P: SimpleNamespace, b_grid: np.ndarray, verbose: bool = True
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    p = np.asarray(p_init, dtype=float).reshape(-1)
    entry_closure = uses_outside_option_closure(P)
    benchmark_norm_closure = uses_benchmark_normalized_outside_closure(P)
    renewal_closure = uses_renewal_valve_closure(P)
    scale_price_closure = uses_accounting_scale_price_closure(P) or renewal_closure
    if scale_price_closure and not benchmark_norm_closure and not (
        renewal_closure
        or
        bool(getattr(P, "outside_value_is_calibrated", False))
        or bool(getattr(P, "allow_uncalibrated_outside_value", False))
    ):
        raise ValueError(
            "population_closure='accounting_scale_prices' requires a calibrated outside value. "
            "Set P.outside_value from accounting_population_scale(..., calibrate_outside_value=True) "
            "and set P.outside_value_is_calibrated=True, or set P.allow_uncalibrated_outside_value=True "
            "for a deliberate experiment."
        )
    if entry_closure:
        entry_by_loc = np.maximum(np.asarray(P.entry_by_loc, dtype=float).reshape(-1), 0.0)
        if entry_by_loc.size != P.I or np.sum(entry_by_loc) <= 0:
            entry_shares0 = np.maximum(np.asarray(P.entry_shares, dtype=float).reshape(-1), 0.0)
            entry_shares0 = entry_shares0 / entry_shares0.sum() if entry_shares0.sum() > 0 else np.ones(P.I) / P.I
            entry_by_loc = float(getattr(P, "E_total", 1 / P.J)) * entry_shares0
        entry_shares, E_total = set_entry_masses(P, entry_by_loc)
    else:
        entry_shares = np.maximum(np.asarray(P.entry_shares, dtype=float).reshape(-1), 0)
        entry_shares = entry_shares / entry_shares.sum() if entry_shares.sum() > 0 else np.ones(P.I) / P.I
        E_total = P.E_total
        P.entry_shares = entry_shares
        P.entry_by_loc = E_total * entry_shares
        entry_by_loc = P.entry_by_loc.copy()

    lam = P.lambda_eq
    lp_damp = lam * np.ones(P.I)
    le_damp = lam * np.ones(P.I)
    dp_prev = np.zeros(P.I)
    de_prev = np.zeros(P.I)
    pts = p.copy()
    ets = entry_by_loc.copy() if entry_closure else entry_shares.copy()
    best_err = np.inf
    best_p = p.copy()
    best_entry = entry_by_loc.copy() if entry_closure else entry_shares.copy()
    best_snap: dict[str, Any] = {}

    t_full = t_eval = t_dist = 0.0
    n_full = n_eval = n_dist = 0
    SD = precompute_shared(P, b_grid)
    stored_bp = None
    howard_freq = max(1, int(round(getattr(P, "howard_freq", 5))))
    force_full = bool(getattr(P, "force_full_bellman", False))
    stall_count = 0
    prev_err = np.inf
    accepted_reason = "max_iter"
    accepted = False
    strict_converged = False
    iterations_completed = 0
    final_eq_error = np.nan
    best_iter = 0

    collect_trace = bool(getattr(P, "collect_ge_trace", False))
    ge_trace: list[dict] | None = [] if collect_trace else None

    if verbose:
        extra = " [full-only]" if force_full else ""
        print(f"  Howard eq: lam={lam:.2f}, tol={P.tol_eq:.1e}, freq={howard_freq}{extra}")

    for it in range(1, P.max_iter_eq + 1):
        r = P.user_cost_rate * p
        P.eq_iter = it
        # Howard policy iteration: every `howard_freq` GE iter (and the
        # first 3 / when stalled) we resolve the full Bellman optimization;
        # in between we re-evaluate at the stored savings policy `bp_pol`,
        # which is much cheaper. `force_full` disables the Howard cycle
        # for parity tests.
        do_full = force_full or it <= 3 or stored_bp is None or (it % howard_freq == 0) or stall_count >= 2

        if do_full:
            t0 = time.perf_counter()
            V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, _ = solve_bellman_full(r, p, P, b_grid, SD)
            t_full += time.perf_counter() - t0
            n_full += 1
            stored_bp = bp_pol.copy()
            stall_count = 0
            mode = "F"
        else:
            t0 = time.perf_counter()
            V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, _ = solve_bellman_eval(stored_bp, r, p, P, b_grid, SD)
            t_eval += time.perf_counter() - t0
            n_eval += 1
            mode = "E"

        t0 = time.perf_counter()
        g, stats = forward_distribution(
            bp_pol, hR_pol, tc, lp_j, fp, r, p, P, b_grid, SD, fast_stats=True, tenure_probs=tp
        )
        t_dist += time.perf_counter() - t0
        n_dist += 1

        scale_info = None
        if scale_price_closure:
            scale_state = SimpleNamespace(
                V=V,
                entry_rate=stats.entry_rate,
                entry_by_loc=stats.entry_by_loc,
                entrants_mature_by_loc=stats.entrants_mature_by_loc,
                entrants_mature_total=stats.entrants_mature_total,
                total_mass=stats.total_mass,
                housing_demand=stats.housing_demand,
            )
            if benchmark_norm_closure:
                scale_info = benchmark_normalized_outside_population_scale(scale_state, P, b_grid)
            elif renewal_closure:
                scale_info = renewal_population_scale(scale_state, P, b_grid)
            else:
                scale_info = accounting_population_scale(scale_state, P, b_grid)
            Hd = np.asarray(scale_info.implied_housing_demand, dtype=float)
            if not scale_info.finite_stationary_scale or not np.all(np.isfinite(Hd)):
                Hd = stats.housing_demand
        else:
            Hd = stats.housing_demand
        p_target = np.zeros(P.I)
        Hs = np.zeros(P.I)
        for i in range(P.I):
            hd = max(Hd[i], P.housing_demand_floor_for_supply)
            p_target[i] = P.r_bar[i] * (hd / P.H0[i]) ** (1 / P.xi_supply[i]) / P.user_cost_rate
            Hs[i] = P.H0[i] * (r[i] / P.r_bar[i]) ** P.xi_supply[i]
        if entry_closure:
            entry_target, entry_info = outside_option_entry_target(V, b_grid, stats, P)
        elif scale_price_closure and scale_info is not None and scale_info.finite_stationary_scale:
            if renewal_closure:
                entry_target = np.asarray(scale_info.conditional_entry_shares, dtype=float)
            else:
                entry_level_target = np.asarray(scale_info.implied_entry_by_loc, dtype=float)
                et_sum = float(np.sum(entry_level_target))
                entry_target = entry_level_target / et_sum if et_sum > 1e-14 else stats.mature_entry_shares.copy()
            entry_info = {
                "entry_values": getattr(scale_info, "entry_values", np.full(P.I, np.nan)),
                "city_entry_prob": getattr(scale_info, "city_entry_prob", np.full(P.I, np.nan)),
                "outside_entry_prob": getattr(scale_info, "outside_entry_prob", np.nan),
                "potential_entrant_mass": getattr(scale_info, "implied_potential_entrant_mass", np.nan),
                "outside_entry_mass": getattr(scale_info, "implied_outside_entry_mass", np.nan),
                "mature_cityborn_flow": scale_info.reference_mature_cityborn_flow,
                "scale_factor": scale_info.scale_factor,
                "implied_total_population": scale_info.implied_total_population,
            }
        else:
            entry_target = stats.mature_entry_shares.copy()
            entry_info = {}

        if P.adaptive_price_damping:
            fw = P.target_filter_weight
            pts = fw * p_target + (1 - fw) * pts
            if entry_closure and str(getattr(P, "entry_level_update", "log")).lower() == "log":
                floor = max(float(getattr(P, "entry_mass_floor", 1e-12)), 1e-14)
                ets = np.exp(fw * np.log(np.maximum(entry_target, floor)) + (1 - fw) * np.log(np.maximum(ets, floor)))
            else:
                ets = fw * entry_target + (1 - fw) * ets
            pt = pts
            et = ets
        else:
            pt = p_target
            et = entry_target

        err_p = float(np.max(np.abs(p_target - p) / np.maximum(np.abs(p), 1e-6)))
        if entry_closure:
            floor = max(float(getattr(P, "entry_mass_floor", 1e-12)), 1e-14)
            err_e = float(np.max(np.abs(np.log(np.maximum(entry_target, floor)) - np.log(np.maximum(entry_by_loc, floor)))))
        else:
            err_e = float(np.max(np.abs(entry_target - entry_shares)))
        ov = max(err_p, err_e)
        iterations_completed = it
        final_eq_error = ov

        if collect_trace:
            ge_trace.append({
                "iter": int(it),
                "p": p.tolist(),
                "entry_shares": entry_shares.tolist(),
                "entry_total": float(E_total),
                "err_p": float(err_p),
                "err_e": float(err_e),
                "err": float(ov),
                "own_rate": float(stats.own_rate),
                "TFR": float(2 * stats.mean_parity),
                "pop_share": stats.pop_share.tolist() if hasattr(stats.pop_share, "tolist") else list(stats.pop_share),
                "mode": mode,
            })
            if entry_closure:
                ge_trace[-1]["entry_target_total"] = float(np.sum(entry_target))
                ge_trace[-1]["potential_entrant_mass"] = float(entry_info.get("potential_entrant_mass", np.nan))
                ge_trace[-1]["outside_entry_prob"] = float(entry_info.get("outside_entry_prob", np.nan))
            if scale_price_closure:
                ge_trace[-1]["scale_factor"] = float(entry_info.get("scale_factor", np.nan))
                ge_trace[-1]["implied_total_population"] = float(entry_info.get("implied_total_population", np.nan))
                ge_trace[-1]["outside_entry_prob"] = float(entry_info.get("outside_entry_prob", np.nan))
                if scale_info is not None:
                    ge_trace[-1]["scale_denominator"] = float(scale_info.denominator)
                    ge_trace[-1]["scale_residual"] = float(scale_info.stationary_entry_residual)

        if ov < best_err:
            best_err = ov
            best_iter = it
            best_p = p.copy()
            best_entry = entry_by_loc.copy() if entry_closure else entry_shares.copy()
            best_snap = {
                "V": V,
                "c_pol": c_pol,
                "hR_pol": hR_pol,
                "bp_pol": bp_pol,
                "tc": tc,
                "tp": tp,
                "lp_j": lp_j,
                "fp": fp,
                "fv": fv,
                "g": g,
                "r": r.copy(),
                "p": p.copy(),
            }

        if verbose:
            pop_text = ",".join(f"{x:.2f}" for x in np.asarray(stats.pop_share).reshape(-1))
            print(
                f"  {mode}{it:3d}: ep={err_p:.4f} ee={err_e:.4f} "
                f"own={100 * stats.own_rate:.1f}% TFR={2 * stats.mean_parity:.2f} "
                f"pop=[{pop_text}]"
            )

        # Soft-accept: matches the MATLAB live solver's tiered exit. We
        # don't insist on `ov < tol_eq` once enough iterations have run
        # and `best_err` is close enough; otherwise oscillation around
        # tenure thresholds can prevent strict convergence.
        if ov < P.tol_eq:
            if verbose:
                print(f"  *** CONVERGED {it} iters ***")
            accepted_reason = "strict_tol"
            accepted = True
            strict_converged = True
            break
        if it > 10 and best_err <= P.tol_eq:
            if verbose:
                print("  *** NEAR-CONV ***")
            accepted_reason = "best_within_tol"
            accepted = True
            strict_converged = True
            break
        if it > 25 and best_err < 10 * P.tol_eq:
            if verbose:
                print(f"  *** ACCEPTED at {it} iters (best_err={best_err:.2e}) ***")
            accepted_reason = "soft_tol_10x"
            accepted = True
            break
        if it > 50 and best_err < 0.02:
            if verbose:
                print(f"  *** CAPPED at {it} iters (best_err={best_err:.2e}) ***")
            accepted_reason = "capped_best_err"
            accepted = True
            break

        stall_count = stall_count + 1 if ov >= prev_err * 0.999 else 0
        prev_err = ov

        if P.adaptive_price_damping:
            dp = pt - p
            entry_state = entry_by_loc if entry_closure else entry_shares
            if entry_closure and str(getattr(P, "entry_level_update", "log")).lower() == "log":
                floor = max(float(getattr(P, "entry_mass_floor", 1e-12)), 1e-14)
                max_step = max(float(getattr(P, "entry_log_step_max", 1.0)), 1e-6)
                de = np.log(np.maximum(et, floor)) - np.log(np.maximum(entry_state, floor))
                de = np.clip(de, -max_step, max_step)
            else:
                de = et - entry_state
            fp_mask = (np.abs(dp_prev) > 1e-12) & (np.sign(dp) != np.sign(dp_prev))
            fe_mask = (np.abs(de_prev) > 1e-12) & (np.sign(de) != np.sign(de_prev))
            lp_damp[fp_mask] = np.maximum(P.lambda_price_min, P.adaptive_damping_decay * lp_damp[fp_mask])
            lp_damp[~fp_mask] = np.minimum(P.lambda_price_max, P.adaptive_damping_grow * lp_damp[~fp_mask])
            le_damp[fe_mask] = np.maximum(P.lambda_entry_min, P.adaptive_damping_decay * le_damp[fe_mask])
            le_damp[~fe_mask] = np.minimum(P.lambda_entry_max, P.adaptive_damping_grow * le_damp[~fe_mask])
            if fp_mask.sum() + fe_mask.sum() >= 2 * P.I:
                lp_damp = np.maximum(P.lambda_price_min, P.cycle_guard_factor * lp_damp)
                le_damp = np.maximum(P.lambda_entry_min, P.cycle_guard_factor * le_damp)
            if it > 15 and ov > 1.05 * best_err:
                lp_damp = np.maximum(P.lambda_price_min, 0.9 * lp_damp)
                le_damp = np.maximum(P.lambda_entry_min, 0.9 * le_damp)
            p = p + lp_damp * dp
            if entry_closure:
                if str(getattr(P, "entry_level_update", "log")).lower() == "log":
                    floor = max(float(getattr(P, "entry_mass_floor", 1e-12)), 1e-14)
                    entry_by_loc = np.exp(np.log(np.maximum(entry_by_loc, floor)) + le_damp * de)
                else:
                    entry_by_loc = entry_by_loc + le_damp * de
            else:
                entry_shares = entry_shares + le_damp * de
            dp_prev = dp
            de_prev = de
        else:
            p = p + lam * (pt - p)
            if entry_closure:
                if str(getattr(P, "entry_level_update", "log")).lower() == "log":
                    floor = max(float(getattr(P, "entry_mass_floor", 1e-12)), 1e-14)
                    max_step = max(float(getattr(P, "entry_log_step_max", 1.0)), 1e-6)
                    dlog = np.log(np.maximum(et, floor)) - np.log(np.maximum(entry_by_loc, floor))
                    dlog = np.clip(dlog, -max_step, max_step)
                    entry_by_loc = np.exp(np.log(np.maximum(entry_by_loc, floor)) + lam * dlog)
                else:
                    entry_by_loc = entry_by_loc + lam * (et - entry_by_loc)
            else:
                entry_shares = entry_shares + lam * (et - entry_shares)

        if P.enforce_price_bounds:
            p = np.clip(p, P.p_min, P.p_max)
        if entry_closure:
            entry_by_loc = np.maximum(entry_by_loc, 0.0)
            if P.enforce_entry_share_floor:
                entry_floor = max(float(getattr(P, "entry_mass_floor", 0.0)), 0.0)
                if entry_floor > 0:
                    entry_by_loc = np.maximum(entry_by_loc, entry_floor)
            entry_shares, E_total = set_entry_masses(P, entry_by_loc)
        else:
            if P.enforce_entry_share_floor:
                entry_shares = np.maximum(entry_shares, P.entry_share_floor)
            entry_shares = np.maximum(entry_shares, 0)
            entry_shares = entry_shares / entry_shares.sum() if entry_shares.sum() > 0 else np.ones(P.I) / P.I
            P.entry_shares = entry_shares
            P.entry_by_loc = E_total * entry_shares
            entry_by_loc = P.entry_by_loc.copy()

    if entry_closure:
        entry_shares, E_total = set_entry_masses(P, best_entry)
    else:
        P.entry_shares = best_entry
        P.entry_by_loc = E_total * best_entry
    p_eq = best_p
    if verbose:
        print(f"\n  [TIME] Full Bellman: {n_full} calls, {t_full:.3f}s ({1000 * t_full / max(n_full, 1):.0f} ms/call)")
        print(f"  [TIME] Eval Bellman: {n_eval} calls, {t_eval:.3f}s ({1000 * t_eval / max(n_eval, 1):.0f} ms/call)")
        print(f"  [TIME] Distribution: {n_dist} calls, {t_dist:.3f}s ({1000 * t_dist / max(n_dist, 1):.0f} ms/call)")

    S = best_snap
    _, full_stats = forward_distribution(
        S["bp_pol"],
        S["hR_pol"],
        S["tc"],
        S["lp_j"],
        S["fp"],
        S["r"],
        S["p"],
        P,
        b_grid,
        SD,
        fast_stats=False,
        tenure_probs=S.get("tp"),
    )
    sol = pack_solution(
        S["V"],
        S["c_pol"],
        S["hR_pol"],
        S["bp_pol"],
        S["tc"],
        S.get("tp"),
        S["lp_j"],
        S["fp"],
        S["fv"],
        S["g"],
        full_stats,
        P.w_hat,
        S["p"],
        P,
    )
    sol.timings = {
        "bellman_full": t_full,
        "bellman_eval": t_eval,
        "distribution": t_dist,
        "n_full": n_full,
        "n_eval": n_eval,
        "n_dist": n_dist,
        "best_eq_error": best_err,
        "best_eq_iter": best_iter,
        "final_eq_error": final_eq_error,
        "iterations_completed": iterations_completed,
        "accepted": accepted,
        "strict_converged": strict_converged,
        "convergence_reason": accepted_reason,
        "population_closure": str(getattr(P, "population_closure", "normalized")),
    }
    if scale_price_closure:
        if benchmark_norm_closure:
            sol.accounting_scale = benchmark_normalized_outside_population_scale(sol, P, b_grid)
        elif renewal_closure:
            sol.accounting_scale = renewal_population_scale(sol, P, b_grid)
        else:
            sol.accounting_scale = accounting_population_scale(sol, P, b_grid)
    if collect_trace:
        sol.ge_trace = ge_trace
    return sol, P, p_eq


def precompute_shared(P: SimpleNamespace, b_grid: np.ndarray) -> SimpleNamespace:
    Nb = len(b_grid)
    nc = P.n_parity * P.n_child_states
    K = P.n_child_stages
    csm1 = K + 1
    c_bar = np.zeros((P.n_parity, P.n_child_states))
    h_bar = np.zeros((P.n_parity, P.n_child_states))
    psi_v = np.zeros((P.n_parity, P.n_child_states))
    for nn in range(P.n_parity):
        nk = nn
        for cs in range(P.n_child_states):
            kp = (cs >= 1) and (cs < csm1)
            if kp:
                c_bar[nn, cs] = P.c_bar_0 + P.c_bar_n * nk
                if str(P.child_housing_spec).lower() == "linear_only":
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_n * nk
                else:
                    h_bar[nn, cs] = P.h_bar_0 + P.h_bar_jump + P.h_bar_n * nk
                psi_v[nn, cs] = P.psi_child * nk
            else:
                c_bar[nn, cs] = P.c_bar_0
                h_bar[nn, cs] = P.h_bar_0

    triples = np.column_stack(
        [
            c_bar.reshape(-1, order="F"),
            h_bar.reshape(-1, order="F"),
            psi_v.reshape(-1, order="F"),
        ]
    )
    unique_triples, type_map = np.unique(triples, axis=0, return_inverse=True)
    birth_dp = np.zeros((P.n_parity, P.n_child_states, 1 + P.n_house, 1 + P.n_house), dtype=bool)
    for nn in range(P.n_parity):
        for cs in range(P.n_child_states):
            for to in range(1 + P.n_house):
                for tn in range(1 + P.n_house):
                    birth_dp[nn, cs, to, tn] = has_birth_dp_grant(P, nn, cs, to, tn)

    return SimpleNamespace(
        c_bar=c_bar,
        h_bar=h_bar,
        psi_v=psi_v,
        cb_flat=c_bar.reshape(1, nc, order="F"),
        hb_flat=h_bar.reshape(1, nc, order="F"),
        psi_flat=psi_v.reshape(1, nc, order="F"),
        nc=nc,
        b=b_grid.reshape(-1, 1),
        bp=b_grid.reshape(1, -1),
        phi_state=get_phi_state_matrix(P),
        phi_choice=get_phi_choice_tensor(P),
        n_types=unique_triples.shape[0],
        type_map=type_map,
        type_cb=unique_triples[:, 0],
        type_hb=unique_triples[:, 1],
        type_psi=unique_triples[:, 2],
        birth_dp=birth_dp,
        birth_entry_grant=get_birth_entry_grant_tensor(P),
    )


def solve_bellman_full(r_hat: np.ndarray, p_hat: np.ndarray, P: SimpleNamespace, b_grid: np.ndarray, SD: SimpleNamespace):
    return solve_bellman_core(r_hat, p_hat, P, b_grid, SD, stored_bp=None, eval_mode=False)


def solve_bellman_eval(
    stored_bp: np.ndarray, r_hat: np.ndarray, p_hat: np.ndarray, P: SimpleNamespace, b_grid: np.ndarray, SD: SimpleNamespace
):
    return solve_bellman_core(r_hat, p_hat, P, b_grid, SD, stored_bp=stored_bp, eval_mode=True)


def solve_bellman_full_markov_income(
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
):
    t0 = time.perf_counter()
    J = P.J
    I = P.I
    Nb = len(b_grid)
    nh = P.n_house
    nt = 1 + nh
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    z_grid, z_weights, Pi_z = income_transition_values(P)
    Nz = len(z_grid)
    beta = P.beta
    Rg = P.R_gross
    sigma = P.sigma
    alpha = P.alpha_cons
    oms = 1.0 - sigma
    owner_h_bar_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
    owner_size_cost = float(getattr(P, "owner_size_cost", 0.0))
    owner_size_cost_ref = float(getattr(P, "owner_size_cost_ref", 6.0))
    owner_size_cost_power = float(getattr(P, "owner_size_cost_power", 2.0))
    tenure_choice_kappa = max(float(getattr(P, "tenure_choice_kappa", 0.0)), 0.0)
    use_tenure_logit = tenure_choice_kappa > 0.0
    b = b_grid.reshape(-1, 1)
    use_full_kernel = NUMBA_AVAILABLE and bool(getattr(P, "use_full_kernel", True))
    cb_v = np.ascontiguousarray(SD.cb_flat.reshape(-1))
    hb_v = np.ascontiguousarray(SD.hb_flat.reshape(-1))
    psi_v_flat = np.ascontiguousarray(SD.psi_flat.reshape(-1))

    V = np.zeros((Nb, nt, I, J, Nz, npar, ncs))
    c_pol = np.zeros_like(V)
    hR_pol = np.zeros_like(V)
    bp_pol = np.ones_like(V)
    tenure_choice = np.zeros((Nb, nt, I, J, Nz, npar, ncs), dtype=np.int16)
    tenure_probs = (
        np.zeros((Nb, nt, I, J, Nz, npar, ncs, nt), dtype=np.float32)
        if use_tenure_logit
        else None
    )
    loc_probs = np.zeros((Nb, nt, I, I, J, Nz, npar, ncs))
    fert_probs = np.zeros((Nb, nt, I, J, Nz, npar))
    fert_value = np.zeros((Nb, nt, I, J, Nz))

    phi_choice = SD.phi_choice
    birth_entry_grant = SD.birth_entry_grant
    hcost = np.zeros((I, nt))
    heq = np.zeros((I, nt))
    dp_arr = np.zeros((I, nt, npar, ncs))
    bmo = np.zeros((I, nt, npar, ncs))
    hsrv = np.zeros((I, nt))
    ocst = np.zeros((I, nt))
    for i in range(I):
        for ten in range(1, nt):
            hs = P.H_own[ten - 1]
            hcost[i, ten] = p_hat[i] * hs
            heq[i, ten] = (1 - P.psi) * p_hat[i] * hs
            hsrv[i, ten] = P.chi * hs
            extra_size_cost = owner_size_cost * p_hat[i] * max(hs - owner_size_cost_ref, 0.0) ** owner_size_cost_power
            ocst[i, ten] = (P.delta + P.tau_H) * p_hat[i] * hs + extra_size_cost
            for nn in range(npar):
                for cs in range(ncs):
                    phi_ncs = phi_choice[i, ten, nn, cs]
                    dp_arr[i, ten, nn, cs] = (1 - phi_ncs) * hcost[i, ten]
                    bmo[i, ten, nn, cs] = -phi_ncs * hcost[i, ten]

    Vbq = np.zeros((Nb, nt, I, npar, ncs))
    for i in range(I):
        for ten in range(nt):
            hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    nk = get_completed_fertility(nn, cs, P)
                    Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)

    loc_shift = np.zeros((I, I))
    for io in range(I):
        for id_ in range(I):
            move_cost = P.mu_stay if id_ == io else P.mu_move
            loc_shift[io, id_] = P.E_loc[id_] - move_cost

    iidx = np.zeros((Nb, I, nt), dtype=np.int64)
    iwt = np.zeros((Nb, I, nt))
    for io in range(I):
        for to in range(nt):
            ba = np.clip(b_grid + heq[io, to], b_grid[0], b_grid[-1])
            idx, wt = interp_indices(b_grid, ba)
            iidx[:, io, to] = idx
            iwt[:, io, to] = wt

    gs_tol = 1e-3
    gs_alpha1 = (3 - math.sqrt(5)) / 2
    gs_alpha2 = (math.sqrt(5) - 1) / 2

    for j in range(J - 1, -1, -1):
        in_fert = (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end)
        for zz, z_value in enumerate(z_grid):
            if j == J - 1:
                Vnr = Vbq
            else:
                Vnr = np.zeros((Nb, nt, I, npar, ncs))
                for znext in range(Nz):
                    Vnr += Pi_z[zz, znext] * V[:, :, :, j + 1, znext, :, :]
            Vc = apply_child_aging(Vnr, P, Nb, nt, I, npar, ncs)
            Vd = np.zeros((Nb, nt, I, npar, ncs))
            cd = np.zeros_like(Vd)
            hd = np.zeros_like(Vd)
            bd = np.zeros_like(Vd)
            Vo_nc = np.zeros((Nb, nc))
            bp_nc = np.zeros((Nb, nc))

            for i in range(I):
                yj = income_at_state(P, i, j, float(z_value))
                ri = r_hat[i]
                Rv = Rg * b + yj
                hRmax = P.hR_max
                Vcr = flat_nc(Vc[:, 0, i, :, :], Nb, nc)
                Rv1d_full = np.ascontiguousarray(Rv[:, 0])
                if use_full_kernel:
                    bp_prev_r = np.zeros((Nb, nc))
                    has_prev_r = 0
                    Vo_nc, bp_nc, co_nc, ho_nc = full_renter_block_kernel(
                        Rv1d_full, Vcr, bp_prev_r, has_prev_r, b_grid,
                        cb_v, hb_v, psi_v_flat,
                        ri, hRmax, P.c_min, P.c_bar_0, P.h_bar_0,
                        alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol,
                    )
                else:
                    Kr = (alpha**alpha * ((1 - alpha) / ri) ** (1 - alpha)) ** oms
                    d_nc = SD.cb_flat + ri * SD.hb_flat
                    cap_nc = ri * (hRmax - SD.hb_flat) / (1 - alpha)
                    for c in range(nc):
                        Vbar = Vcr[:, c]
                        dc = d_nc[0, c]
                        pc = SD.psi_flat[0, c]
                        cc = cap_nc[0, c]
                        cb_c = SD.cb_flat[0, c]
                        hb_c = SD.hb_flat[0, c]
                        ht_cap_c = max(hRmax - hb_c, 1e-10)
                        lo = np.maximum(np.zeros(Nb), b_grid[0])
                        hi = np.maximum((Rv[:, 0] - dc - 1e-6), 0.0)
                        bp, val = golden_renter(
                            lo, hi, Rv[:, 0], Vbar, b_grid, dc, pc, cc, cb_c, hb_c,
                            ri, hRmax, ht_cap_c, Kr, alpha, oms, beta,
                            gs_alpha1, gs_alpha2, gs_tol,
                        )
                        bp_nc[:, c] = bp
                        Vo_nc[:, c] = val
                    surplus_nc = Rv - d_nc - bp_nc
                    ct_nc = alpha * np.maximum(surplus_nc, 1e-10)
                    ht_nc = (1 - alpha) / ri * np.maximum(surplus_nc, 1e-10)
                    cm = (SD.hb_flat + ht_nc) > hRmax
                    if np.any(cm):
                        ct_cap = np.maximum(Rv - SD.cb_flat - ri * hRmax - bp_nc, 1e-10)
                        hcap = np.tile(np.maximum(hRmax - SD.hb_flat, 1e-10), (Nb, 1))
                        ct_nc[cm] = ct_cap[cm]
                        ht_nc[cm] = hcap[cm]
                    co_nc = SD.cb_flat + np.maximum(ct_nc, P.c_min)
                    ho_nc = SD.hb_flat + np.maximum(ht_nc, 0.01)
                    bad = surplus_nc <= 1e-10
                    co_nc[bad] = P.c_bar_0 + P.c_min
                    ho_nc[bad] = P.h_bar_0 + 0.01
                Vd[:, 0, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                bd[:, 0, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
                cd[:, 0, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)
                hd[:, 0, i, :, :] = unflat_nc(ho_nc, Nb, npar, ncs)

                for ten in range(1, nt):
                    oc = ocst[i, ten]
                    hsv = hsrv[i, ten]
                    Vco = flat_nc(Vc[:, ten, i, :, :], Nb, nc)
                    if use_full_kernel:
                        bf_v = np.ascontiguousarray(bmo[i, ten, :, :].reshape(-1, order="F"))
                        bp_prev_o = np.zeros((Nb, nc))
                        has_prev_o = 0
                        Vo_nc, bp_nc, co_nc = full_owner_block_kernel(
                            Rv1d_full, Vco, bp_prev_o, has_prev_o, b_grid,
                            cb_v, hb_v, psi_v_flat, bf_v,
                            oc, hsv, owner_h_bar_scale, P.c_min,
                            alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol,
                        )
                    else:
                        for c in range(nc):
                            Vbar = Vco[:, c]
                            cb_c = SD.cb_flat[0, c]
                            pc = SD.psi_flat[0, c]
                            ht_c = max(hsv - owner_h_bar_scale * SD.hb_flat[0, c], 1e-10)
                            Ko_c = ht_c ** ((1 - alpha) * oms)
                            nn_c_1 = math.ceil((c + 1) / ncs)
                            cs_c_1 = (c + 1) - (nn_c_1 - 1) * ncs
                            bf_c = bmo[i, ten, nn_c_1 - 1, cs_c_1 - 1]
                            lo = np.maximum(bf_c, b_grid[0]) * np.ones(Nb)
                            hi = np.maximum(Rv[:, 0] - oc - cb_c - 1e-6, lo)
                            bp, val = golden_owner(
                                lo, hi, Rv[:, 0], Vbar, b_grid, oc, cb_c, pc,
                                Ko_c, alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol,
                            )
                            bp_nc[:, c] = bp
                            Vo_nc[:, c] = val
                        co_nc = SD.cb_flat + np.maximum(Rv - oc - SD.cb_flat - bp_nc, P.c_min)
                    Vd[:, ten, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                    bd[:, ten, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
                    cd[:, ten, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)

            c_pol[:, :, :, j, zz, :, :] = cd
            hR_pol[:, :, :, j, zz, :, :] = hd
            bp_pol[:, :, :, j, zz, :, :] = bd

            dp_choice = dp_arr
            if bool(getattr(P, "use_pti_constraint", False)):
                income_j = np.array([income_at_state(P, i, j, float(z_value)) for i in range(I)], dtype=float)
                dp_choice = pti_adjusted_downpayment(dp_arr, hcost, income_j, P, b_grid)

            if use_tenure_logit and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
                VH, tcj, prj = tenure_logit_kernel(
                    Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant, tenure_choice_kappa
                )
                tenure_probs[:, :, :, j, zz, :, :, :] = prj
            elif (not use_tenure_logit) and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
                VH, tcj = tenure_choice_kernel(
                    Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant
                )
            else:
                VH = np.zeros((Nb, nt, I, npar, ncs))
                tcj = np.zeros((Nb, nt, I, npar, ncs), dtype=np.int16)
                for id_ in range(I):
                    for to in range(nt):
                        sp = heq[id_, to] if to > 0 else 0.0
                        Vopt = np.zeros((Nb, npar, ncs, nt))
                        if to == 0:
                            Vopt[:, :, :, 0] = Vd[:, 0, id_, :, :]
                        else:
                            ba = np.maximum(b_grid + sp, 0.0)
                            Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)
                        for tn in range(1, nt):
                            hc = hcost[id_, tn]
                            Vow = Vd[:, tn, id_, :, :]
                            if to == tn:
                                Vopt[:, :, :, tn] = Vow
                            elif to == 0:
                                bab = b_grid - hc
                                Vb = interp_on_grid(b_grid, Vow, bab)
                                for nn in range(npar):
                                    for cs in range(ncs):
                                        dpn = dp_choice[id_, tn, nn, cs]
                                        bmn = bmo[id_, tn, nn, cs]
                                        if SD.birth_dp[nn, cs, to, tn]:
                                            bag = np.maximum(bab, bmn)
                                            Vb[:, nn, cs] = interp_vector(b_grid, Vow[:, nn, cs], bag)
                                        elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                            gfix = birth_entry_grant[id_, tn, nn, cs]
                                            babg = bab + gfix
                                            Vg = interp_vector(b_grid, Vow[:, nn, cs], babg)
                                            inf_m = ((b_grid + gfix) < dpn) | (babg < bmn)
                                            Vg[inf_m] = -1e10
                                            Vb[:, nn, cs] = Vg
                                        else:
                                            inf_m = (b_grid < dpn) | (bab < bmn)
                                            Vb[inf_m, nn, cs] = -1e10
                                Vopt[:, :, :, tn] = Vb
                            else:
                                bar = b_grid + sp - hc
                                Vrs = interp_on_grid(b_grid, Vow, bar)
                                for nn in range(npar):
                                    for cs in range(ncs):
                                        dpn = dp_choice[id_, tn, nn, cs]
                                        bmn = bmo[id_, tn, nn, cs]
                                        dpc = dpn - sp
                                        inf_m = (b_grid < dpc) | (bar < bmn)
                                        Vrs[inf_m, nn, cs] = -1e10
                                Vopt[:, :, :, tn] = Vrs
                        tc = np.argmax(Vopt, axis=3)
                        if use_tenure_logit:
                            ls, pr = logsumexp(Vopt / tenure_choice_kappa, axis=3)
                            VH[:, to, id_, :, :] = tenure_choice_kappa * ls
                            tenure_probs[:, to, id_, j, zz, :, :, :] = pr.astype(np.float32)
                        else:
                            VH[:, to, id_, :, :] = np.max(Vopt, axis=3)
                        tcj[:, to, id_, :, :] = tc
            tenure_choice[:, :, :, j, zz, :, :] = tcj

            kl = P.kappa_loc
            if NUMBA_AVAILABLE and bool(getattr(P, "use_loc_kernel", True)):
                VI, lpj = location_logit_kernel(VH, iidx, iwt, loc_shift, kl)
            else:
                VI = np.zeros((Nb, nt, I, npar, ncs))
                lpj = np.zeros((Nb, nt, I, I, npar, ncs))
                for io in range(I):
                    for to in range(nt):
                        Va = np.zeros((Nb, I, npar, ncs))
                        Va[:, io, :, :] = VH[:, to, io, :, :]
                        idx = iidx[:, io, to]
                        wt = iwt[:, io, to]
                        for id_ in range(I):
                            if id_ == io:
                                continue
                            Vdst = VH[:, 0, id_, :, :]
                            Va[:, id_, :, :] = (1 - wt)[:, None, None] * Vdst[idx, :, :] + wt[:, None, None] * Vdst[idx + 1, :, :]
                        la = Va.copy()
                        for id_ in range(I):
                            la[:, id_, :, :] += loc_shift[io, id_]
                        la = la / kl
                        ls, pr = logsumexp(la, axis=1)
                        VI[:, to, io, :, :] = kl * ls
                        lpj[:, to, io, :, :, :] = pr
            loc_probs[:, :, :, :, j, zz, :, :] = lpj

            if in_fert:
                Vfa = np.zeros((Nb, nt, I, npar))
                Vfa[:, :, :, 0] = VI[:, :, :, 0, 0]
                for nn in range(1, npar):
                    Vfa[:, :, :, nn] = VI[:, :, :, nn, 1]
                lf = Vfa / P.kappa_fert
                ls, pr = logsumexp(lf, axis=3)
                fert_probs[:, :, :, j, zz, :] = pr
                fert_value[:, :, :, j, zz] = P.kappa_fert * ls
                V[:, :, :, j, zz, 0, 0] = fert_value[:, :, :, j, zz]
                V[:, :, :, j, zz, 1:, :] = VI[:, :, :, 1:, :]
                V[:, :, :, j, zz, 0, 1:] = VI[:, :, :, 0, 1:]
            else:
                V[:, :, :, j, zz, :, :] = VI

    return (
        V,
        c_pol,
        hR_pol,
        bp_pol,
        tenure_choice,
        tenure_probs,
        loc_probs,
        fert_probs,
        fert_value,
        {"bellman": time.perf_counter() - t0},
    )


def solve_bellman_core(
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    stored_bp: np.ndarray | None,
    eval_mode: bool,
):
    # Backward induction over age `j`. At each `j` we solve (per i, ten):
    # savings choice via golden-section (full mode) or plug-in at
    # stored_bp (eval mode); then tenure choice; then location logit;
    # then fertility logit at fertile ages. State arrays are indexed
    # (b, ten, i, j, parity n, child_state cs); the `(n, cs)` pair is
    # often flattened to a single column index in F (column-major) order
    # for compatibility with the MATLAB reference.
    t0 = time.perf_counter()
    J = P.J
    I = P.I
    Nb = len(b_grid)
    nh = P.n_house
    nt = 1 + nh
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    beta = P.beta
    Rg = P.R_gross
    sigma = P.sigma
    alpha = P.alpha_cons
    oms = 1.0 - sigma
    owner_h_bar_scale = float(getattr(P, "owner_h_bar_scale", 1.0))
    owner_size_cost = float(getattr(P, "owner_size_cost", 0.0))
    owner_size_cost_ref = float(getattr(P, "owner_size_cost_ref", 6.0))
    owner_size_cost_power = float(getattr(P, "owner_size_cost_power", 2.0))
    tenure_choice_kappa = max(float(getattr(P, "tenure_choice_kappa", 0.0)), 0.0)
    use_tenure_logit = tenure_choice_kappa > 0.0
    b = b_grid.reshape(-1, 1)
    col_offset = Nb * np.arange(nc)
    b_lo = b_grid[0]
    b_hi = b_grid[-1]
    stored_idx = None
    stored_wt = None
    if eval_mode and stored_bp is not None:
        stored_idx, stored_wt = interp_indices(b_grid, np.clip(stored_bp, b_lo, b_hi))
    use_eval_kernel = (
        eval_mode
        and NUMBA_AVAILABLE
        and stored_idx is not None
        and bool(getattr(P, "use_eval_kernel", True))
    )
    use_full_kernel = (
        (not eval_mode) and NUMBA_AVAILABLE and bool(getattr(P, "use_full_kernel", True))
    )
    cb_v = np.ascontiguousarray(SD.cb_flat.reshape(-1))
    hb_v = np.ascontiguousarray(SD.hb_flat.reshape(-1))
    psi_v_flat = np.ascontiguousarray(SD.psi_flat.reshape(-1))

    V = np.zeros((Nb, nt, I, J, npar, ncs))
    c_pol = np.zeros_like(V)
    hR_pol = np.zeros_like(V)
    bp_pol = np.ones_like(V)
    tenure_choice = np.zeros((Nb, nt, I, J, npar, ncs), dtype=np.int16)
    tenure_probs = (
        np.zeros((Nb, nt, I, J, npar, ncs, nt), dtype=np.float32)
        if use_tenure_logit
        else None
    )
    loc_probs = np.zeros((Nb, nt, I, I, J, npar, ncs))
    fert_probs = np.zeros((Nb, nt, I, J, npar))
    fert_value = np.zeros((Nb, nt, I, J))

    phi_choice = SD.phi_choice
    birth_entry_grant = SD.birth_entry_grant
    hcost = np.zeros((I, nt))
    heq = np.zeros((I, nt))
    dp_arr = np.zeros((I, nt, npar, ncs))
    bmo = np.zeros((I, nt, npar, ncs))
    hsrv = np.zeros((I, nt))
    ocst = np.zeros((I, nt))
    for i in range(I):
        for ten in range(1, nt):
            hs = P.H_own[ten - 1]
            hcost[i, ten] = p_hat[i] * hs
            heq[i, ten] = (1 - P.psi) * p_hat[i] * hs
            hsrv[i, ten] = P.chi * hs
            extra_size_cost = owner_size_cost * p_hat[i] * max(hs - owner_size_cost_ref, 0.0) ** owner_size_cost_power
            ocst[i, ten] = (P.delta + P.tau_H) * p_hat[i] * hs + extra_size_cost
            for nn in range(npar):
                for cs in range(ncs):
                    # phi is the financed share, so the down-payment
                    # threshold is (1 - phi) * hcost and the borrowing
                    # limit is -phi * hcost.
                    phi_ncs = phi_choice[i, ten, nn, cs]
                    dp_arr[i, ten, nn, cs] = (1 - phi_ncs) * hcost[i, ten]
                    bmo[i, ten, nn, cs] = -phi_ncs * hcost[i, ten]

    Vbq = np.zeros((Nb, nt, I, npar, ncs))
    for i in range(I):
        for ten in range(nt):
            hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0
            for nn in range(npar):
                for cs in range(ncs):
                    nk = get_completed_fertility(nn, cs, P)
                    Vbq[:, ten, i, nn, cs] = bequest_utility_vec(b_grid + hv, nk, P)

    loc_shift = np.zeros((I, I))
    for io in range(I):
        for id_ in range(I):
            move_cost = P.mu_stay if id_ == io else P.mu_move
            loc_shift[io, id_] = P.E_loc[id_] - move_cost

    iidx = np.zeros((Nb, I, nt), dtype=np.int64)
    iwt = np.zeros((Nb, I, nt))
    for io in range(I):
        for to in range(nt):
            ba = np.clip(b_grid + heq[io, to], b_grid[0], b_grid[-1])
            idx, wt = interp_indices(b_grid, ba)
            iidx[:, io, to] = idx
            iwt[:, io, to] = wt

    Vd = np.zeros((Nb, nt, I, npar, ncs))
    cd = np.zeros_like(Vd)
    hd = np.zeros_like(Vd)
    bd = np.zeros_like(Vd)
    Vo_nc = np.zeros((Nb, nc))
    bp_nc = np.zeros((Nb, nc))

    gs_tol = 1e-3
    gs_alpha1 = (3 - math.sqrt(5)) / 2
    gs_alpha2 = (math.sqrt(5) - 1) / 2

    for j in range(J - 1, -1, -1):
        in_fert = (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end)
        if j == J - 1:
            Vnr = Vbq
        else:
            Vnr = V[:, :, :, j + 1, :, :]
        Vc = apply_child_aging(Vnr, P, Nb, nt, I, npar, ncs)

        for i in range(I):
            yj = P.income[i, j]
            ri = r_hat[i]
            Rv = Rg * b + yj
            if not eval_mode:
                hRmax = P.hR_max
                Vcr = flat_nc(Vc[:, 0, i, :, :], Nb, nc)
                Rv1d_full = np.ascontiguousarray(Rv[:, 0])
                if use_full_kernel:
                    if j < J - 1:
                        bp_prev_r = flat_nc(bd[:, 0, i, :, :], Nb, nc)
                        has_prev_r = 1
                    else:
                        bp_prev_r = np.zeros((Nb, nc))
                        has_prev_r = 0
                    Vo_nc, bp_nc, co_nc, ho_nc = full_renter_block_kernel(
                        Rv1d_full, Vcr, bp_prev_r, has_prev_r, b_grid,
                        cb_v, hb_v, psi_v_flat,
                        ri, hRmax, P.c_min, P.c_bar_0, P.h_bar_0,
                        alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol,
                    )
                else:
                    Kr = (alpha**alpha * ((1 - alpha) / ri) ** (1 - alpha)) ** oms
                    aoms = alpha * oms
                    d_nc = SD.cb_flat + ri * SD.hb_flat
                    cap_nc = ri * (hRmax - SD.hb_flat) / (1 - alpha)
                    bp_prev_r = flat_nc(bd[:, 0, i, :, :], Nb, nc) if j < J - 1 else None
                    for c in range(nc):
                        Vbar = Vcr[:, c]
                        dc = d_nc[0, c]
                        pc = SD.psi_flat[0, c]
                        cc = cap_nc[0, c]
                        cb_c = SD.cb_flat[0, c]
                        hb_c = SD.hb_flat[0, c]
                        ht_cap_c = max(hRmax - hb_c, 1e-10)
                        lo = np.maximum(np.zeros(Nb), b_grid[0])
                        hi = np.maximum((Rv[:, 0] - dc - 1e-6), 0.0)
                        if j < J - 1:
                            lo = np.maximum(lo, bp_prev_r[:, c] - 2.0)
                            hi = np.minimum(hi, bp_prev_r[:, c] + 2.0)
                            lo = np.maximum(lo, 0.0)
                            hi = np.maximum(hi, lo)
                        if NUMBA_AVAILABLE:
                            bp, val = golden_renter_kernel(
                                lo, hi, Rv[:, 0], Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol
                            )
                        else:
                            bp, val = golden_renter(
                                lo, hi, Rv[:, 0], Vbar, b_grid, dc, pc, cc, cb_c, hb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol
                            )
                        bp_nc[:, c] = bp
                        Vo_nc[:, c] = val

                    surplus_nc = Rv - d_nc - bp_nc
                    ct_nc = alpha * np.maximum(surplus_nc, 1e-10)
                    ht_nc = (1 - alpha) / ri * np.maximum(surplus_nc, 1e-10)
                    cm = (SD.hb_flat + ht_nc) > hRmax
                    if np.any(cm):
                        ct_cap = np.maximum(Rv - SD.cb_flat - ri * hRmax - bp_nc, 1e-10)
                        hcap = np.tile(np.maximum(hRmax - SD.hb_flat, 1e-10), (Nb, 1))
                        ct_nc[cm] = ct_cap[cm]
                        ht_nc[cm] = hcap[cm]
                    co_nc = SD.cb_flat + np.maximum(ct_nc, P.c_min)
                    ho_nc = SD.hb_flat + np.maximum(ht_nc, 0.01)
                    bad = surplus_nc <= 1e-10
                    co_nc[bad] = P.c_bar_0 + P.c_min
                    ho_nc[bad] = P.h_bar_0 + 0.01
                Vd[:, 0, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                bd[:, 0, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
                cd[:, 0, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)
                hd[:, 0, i, :, :] = unflat_nc(ho_nc, Nb, npar, ncs)

                for ten in range(1, nt):
                    oc = ocst[i, ten]
                    hsv = hsrv[i, ten]
                    Vco = flat_nc(Vc[:, ten, i, :, :], Nb, nc)
                    if use_full_kernel:
                        bf_v = np.ascontiguousarray(bmo[i, ten, :, :].reshape(-1, order="F"))
                        if j < J - 1:
                            bp_prev_o = flat_nc(bd[:, ten, i, :, :], Nb, nc)
                            has_prev_o = 1
                        else:
                            bp_prev_o = np.zeros((Nb, nc))
                            has_prev_o = 0
                        Vo_nc, bp_nc, co_nc = full_owner_block_kernel(
                            Rv1d_full, Vco, bp_prev_o, has_prev_o, b_grid,
                            cb_v, hb_v, psi_v_flat, bf_v,
                            oc, hsv, owner_h_bar_scale, P.c_min,
                            alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol,
                        )
                    else:
                        bp_prev_o = flat_nc(bd[:, ten, i, :, :], Nb, nc) if j < J - 1 else None
                        for c in range(nc):
                            Vbar = Vco[:, c]
                            cb_c = SD.cb_flat[0, c]
                            pc = SD.psi_flat[0, c]
                            ht_c = max(hsv - owner_h_bar_scale * SD.hb_flat[0, c], 1e-10)
                            Ko_c = ht_c ** ((1 - alpha) * oms)
                            nn_c_1 = math.ceil((c + 1) / ncs)
                            cs_c_1 = (c + 1) - (nn_c_1 - 1) * ncs
                            bf_c = bmo[i, ten, nn_c_1 - 1, cs_c_1 - 1]
                            lo = np.maximum(bf_c, b_grid[0]) * np.ones(Nb)
                            hi = np.maximum(Rv[:, 0] - oc - cb_c - 1e-6, lo)
                            if j < J - 1:
                                lo = np.maximum(lo, bp_prev_o[:, c] - 2.0)
                                hi = np.minimum(hi, bp_prev_o[:, c] + 2.0)
                                lo = np.maximum(lo, bf_c)
                                hi = np.maximum(hi, lo)
                            if NUMBA_AVAILABLE:
                                bp, val = golden_owner_kernel(
                                    lo, hi, Rv[:, 0], Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol
                                )
                            else:
                                bp, val = golden_owner(
                                    lo, hi, Rv[:, 0], Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol
                                )
                            bp_nc[:, c] = bp
                            Vo_nc[:, c] = val
                        co_nc = SD.cb_flat + np.maximum(Rv - oc - SD.cb_flat - bp_nc, P.c_min)
                    Vd[:, ten, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                    bd[:, ten, i, :, :] = unflat_nc(bp_nc, Nb, npar, ncs)
                    cd[:, ten, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)
            else:
                bpv = flat_nc(stored_bp[:, 0, i, j, :, :], Nb, nc)
                Vcr_nc = flat_nc(Vc[:, 0, i, :, :], Nb, nc)
                Rv1d = np.ascontiguousarray(Rv[:, 0]) if use_eval_kernel else None
                if use_eval_kernel:
                    idx_nc = flat_nc(stored_idx[:, 0, i, j, :, :], Nb, nc).astype(np.int64, copy=False)
                    wt_nc = flat_nc(stored_wt[:, 0, i, j, :, :], Nb, nc)
                    Vo_nc, co_nc, ho_nc = eval_renter_block_kernel(
                        Rv1d, bpv, Vcr_nc, idx_nc, wt_nc, cb_v, hb_v, psi_v_flat,
                        ri, P.hR_max, P.c_min, P.c_bar_0, P.h_bar_0, alpha, oms, beta,
                    )
                else:
                    if NUMBA_AVAILABLE and stored_idx is not None and stored_wt is not None:
                        idx_nc = flat_nc(stored_idx[:, 0, i, j, :, :], Nb, nc)
                        wt_nc = flat_nc(stored_wt[:, 0, i, j, :, :], Nb, nc)
                        Vcbp = interp_cols_preidx_kernel(Vcr_nc, idx_nc, wt_nc)
                    elif NUMBA_AVAILABLE:
                        Vcbp = interp_cols_kernel(b_grid, Vcr_nc, np.clip(bpv, b_lo, b_hi))
                    else:
                        Vcbp = interp_cols(b_grid, Vcr_nc, np.clip(bpv, b_lo, b_hi))
                    Kr_ev = (alpha**alpha * ((1 - alpha) / ri) ** (1 - alpha)) ** oms
                    d_nc = SD.cb_flat + ri * SD.hb_flat
                    cap_nc_ev = ri * (P.hR_max - SD.hb_flat) / (1 - alpha)
                    surplus_nc = Rv - d_nc - bpv
                    ss = np.maximum(surplus_nc, 1e-10)
                    Vo_nc = Kr_ev * ss**oms / oms + SD.psi_flat + beta * Vcbp
                    cm = surplus_nc > cap_nc_ev
                    if np.any(cm):
                        ct_cap = np.maximum(Rv - SD.cb_flat - ri * P.hR_max - bpv, 1e-10)
                        ht_cap_v = np.tile(np.maximum(P.hR_max - SD.hb_flat, 1e-10), (Nb, 1))
                        comp_cap = ct_cap**alpha * ht_cap_v ** (1 - alpha)
                        u_cap = comp_cap**oms / oms + SD.psi_flat
                        Vo_nc[cm] = u_cap[cm] + beta * Vcbp[cm]
                    Vo_nc[surplus_nc <= 1e-10] = -1e10
                    ct_nc = alpha * np.maximum(surplus_nc, 1e-10)
                    ht_nc = (1 - alpha) / ri * np.maximum(surplus_nc, 1e-10)
                    if np.any(cm):
                        ct_nc[cm] = ct_cap[cm]
                        ht_nc[cm] = ht_cap_v[cm]
                    co_nc = SD.cb_flat + np.maximum(ct_nc, P.c_min)
                    ho_nc = SD.hb_flat + np.maximum(ht_nc, 0.01)
                    bad = surplus_nc <= 1e-10
                    co_nc[bad] = P.c_bar_0 + P.c_min
                    ho_nc[bad] = P.h_bar_0 + 0.01
                Vd[:, 0, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                bd[:, 0, i, :, :] = unflat_nc(bpv, Nb, npar, ncs)
                cd[:, 0, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)
                hd[:, 0, i, :, :] = unflat_nc(ho_nc, Nb, npar, ncs)

                for ten in range(1, nt):
                    oc = ocst[i, ten]
                    hsv = hsrv[i, ten]
                    bpv_o = flat_nc(stored_bp[:, ten, i, j, :, :], Nb, nc)
                    Vco_nc = flat_nc(Vc[:, ten, i, :, :], Nb, nc)
                    if use_eval_kernel:
                        idx_nc_o = flat_nc(stored_idx[:, ten, i, j, :, :], Nb, nc).astype(np.int64, copy=False)
                        wt_nc_o = flat_nc(stored_wt[:, ten, i, j, :, :], Nb, nc)
                        Vo_nc, co_nc = eval_owner_block_kernel(
                            Rv1d, bpv_o, Vco_nc, idx_nc_o, wt_nc_o, cb_v, hb_v, psi_v_flat,
                            oc, hsv, owner_h_bar_scale, P.c_min, alpha, oms, beta,
                        )
                    else:
                        if NUMBA_AVAILABLE and stored_idx is not None and stored_wt is not None:
                            idx_nc = flat_nc(stored_idx[:, ten, i, j, :, :], Nb, nc)
                            wt_nc = flat_nc(stored_wt[:, ten, i, j, :, :], Nb, nc)
                            Vcbpo = interp_cols_preidx_kernel(Vco_nc, idx_nc, wt_nc)
                        elif NUMBA_AVAILABLE:
                            Vcbpo = interp_cols_kernel(b_grid, Vco_nc, np.clip(bpv_o, b_lo, b_hi))
                        else:
                            Vcbpo = interp_cols(b_grid, Vco_nc, np.clip(bpv_o, b_lo, b_hi))
                        ct_o = np.maximum(Rv - oc - SD.cb_flat - bpv_o, 1e-10)
                        ht_o = np.maximum(hsv - owner_h_bar_scale * SD.hb_flat, 1e-10)
                        Ko_ev = ht_o ** ((1 - alpha) * oms)
                        Vo_nc = Ko_ev * ct_o ** (alpha * oms) / oms + SD.psi_flat + beta * Vcbpo
                        Vo_nc[ct_o <= 1e-10] = -1e10
                        co_nc = SD.cb_flat + np.maximum(ct_o, P.c_min)
                    Vd[:, ten, i, :, :] = unflat_nc(Vo_nc, Nb, npar, ncs)
                    bd[:, ten, i, :, :] = unflat_nc(bpv_o, Nb, npar, ncs)
                    cd[:, ten, i, :, :] = unflat_nc(co_nc, Nb, npar, ncs)

        c_pol[:, :, :, j, :, :] = cd
        hR_pol[:, :, :, j, :, :] = hd
        bp_pol[:, :, :, j, :, :] = bd

        dp_choice = dp_arr
        if bool(getattr(P, "use_pti_constraint", False)):
            dp_choice = pti_adjusted_downpayment(dp_arr, hcost, P.income[:, j], P, b_grid)

        if use_tenure_logit and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
            VH, tcj, prj = tenure_logit_kernel(
                Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant, tenure_choice_kappa
            )
            tenure_probs[:, :, :, j, :, :, :] = prj
        elif (not use_tenure_logit) and NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
            VH, tcj = tenure_choice_kernel(
                Vd, b_grid, heq, hcost, dp_choice, bmo, SD.birth_dp, birth_entry_grant
            )
        else:
            VH = np.zeros((Nb, nt, I, npar, ncs))
            tcj = np.zeros((Nb, nt, I, npar, ncs), dtype=np.int16)
            for id_ in range(I):
                for to in range(nt):
                    sp = heq[id_, to] if to > 0 else 0.0
                    Vopt = np.zeros((Nb, npar, ncs, nt))
                    if to == 0:
                        Vopt[:, :, :, 0] = Vd[:, 0, id_, :, :]
                    else:
                        ba = np.maximum(b_grid + sp, 0.0)
                        Vopt[:, :, :, 0] = interp_on_grid(b_grid, Vd[:, 0, id_, :, :], ba)
                    for tn in range(1, nt):
                        hc = hcost[id_, tn]
                        Vow = Vd[:, tn, id_, :, :]
                        if to == tn:
                            Vopt[:, :, :, tn] = Vow
                        elif to == 0:
                            bab = b_grid - hc
                            Vb = interp_on_grid(b_grid, Vow, bab)
                            for nn in range(npar):
                                for cs in range(ncs):
                                    dpn = dp_choice[id_, tn, nn, cs]
                                    bmn = bmo[id_, tn, nn, cs]
                                    if SD.birth_dp[nn, cs, to, tn]:
                                        bag = np.maximum(bab, bmn)
                                        Vb[:, nn, cs] = interp_vector(b_grid, Vow[:, nn, cs], bag)
                                    elif birth_entry_grant[id_, tn, nn, cs] > 0:
                                        gfix = birth_entry_grant[id_, tn, nn, cs]
                                        babg = bab + gfix
                                        Vg = interp_vector(b_grid, Vow[:, nn, cs], babg)
                                        inf_m = ((b_grid + gfix) < dpn) | (babg < bmn)
                                        Vg[inf_m] = -1e10
                                        Vb[:, nn, cs] = Vg
                                    else:
                                        inf_m = (b_grid < dpn) | (bab < bmn)
                                        Vb[inf_m, nn, cs] = -1e10
                            Vopt[:, :, :, tn] = Vb
                        else:
                            bar = b_grid + sp - hc
                            Vrs = interp_on_grid(b_grid, Vow, bar)
                            for nn in range(npar):
                                for cs in range(ncs):
                                    dpn = dp_choice[id_, tn, nn, cs]
                                    bmn = bmo[id_, tn, nn, cs]
                                    dpc = dpn - sp
                                    inf_m = (b_grid < dpc) | (bar < bmn)
                                    Vrs[inf_m, nn, cs] = -1e10
                            Vopt[:, :, :, tn] = Vrs
                    tc = np.argmax(Vopt, axis=3)
                    if use_tenure_logit:
                        ls, pr = logsumexp(Vopt / tenure_choice_kappa, axis=3)
                        VH[:, to, id_, :, :] = tenure_choice_kappa * ls
                        tenure_probs[:, to, id_, j, :, :, :] = pr.astype(np.float32)
                    else:
                        VH[:, to, id_, :, :] = np.max(Vopt, axis=3)
                    tcj[:, to, id_, :, :] = tc
        tenure_choice[:, :, :, j, :, :] = tcj

        kl = P.kappa_loc
        if NUMBA_AVAILABLE and bool(getattr(P, "use_loc_kernel", True)):
            VI, lpj = location_logit_kernel(VH, iidx, iwt, loc_shift, kl)
        else:
            VI = np.zeros((Nb, nt, I, npar, ncs))
            lpj = np.zeros((Nb, nt, I, I, npar, ncs))
            for io in range(I):
                for to in range(nt):
                    Va = np.zeros((Nb, I, npar, ncs))
                    Va[:, io, :, :] = VH[:, to, io, :, :]
                    idx = iidx[:, io, to]
                    wt = iwt[:, io, to]
                    for id_ in range(I):
                        if id_ == io:
                            continue
                        Vdst = VH[:, 0, id_, :, :]
                        Va[:, id_, :, :] = (1 - wt)[:, None, None] * Vdst[idx, :, :] + wt[:, None, None] * Vdst[idx + 1, :, :]
                    la = Va.copy()
                    for id_ in range(I):
                        la[:, id_, :, :] += loc_shift[io, id_]
                    la = la / kl
                    ls, pr = logsumexp(la, axis=1)
                    VI[:, to, io, :, :] = kl * ls
                    lpj[:, to, io, :, :, :] = pr
        loc_probs[:, :, :, :, j, :, :] = lpj

        if in_fert:
            Vfa = np.zeros((Nb, nt, I, npar))
            Vfa[:, :, :, 0] = VI[:, :, :, 0, 0]
            for nn in range(1, npar):
                Vfa[:, :, :, nn] = VI[:, :, :, nn, 1]
            lf = Vfa / P.kappa_fert
            ls, pr = logsumexp(lf, axis=3)
            fert_probs[:, :, :, j, :] = pr
            fert_value[:, :, :, j] = P.kappa_fert * ls
            V[:, :, :, j, 0, 0] = fert_value[:, :, :, j]
            V[:, :, :, j, 1:, :] = VI[:, :, :, 1:, :]
            V[:, :, :, j, 0, 1:] = VI[:, :, :, 0, 1:]
        else:
            V[:, :, :, j, :, :] = VI

    return (
        V,
        c_pol,
        hR_pol,
        bp_pol,
        tenure_choice,
        tenure_probs,
        loc_probs,
        fert_probs,
        fert_value,
        {"bellman": time.perf_counter() - t0},
    )


def golden_renter(lo, hi, Rv, Vbar, b_grid, dc, pc, cc, cb_c, hb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, a1, a2, tol):
    d = hi - lo
    x1 = lo + a1 * d
    x2 = lo + a2 * d
    f1 = eval_renter(x1, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
    f2 = eval_renter(x2, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
    d = a1 * a2 * d
    while np.any(d > tol):
        bt = f2 >= f1
        xe = np.clip(np.where(bt, x2 + d, x1 - d), lo, hi)
        fe = eval_renter(xe, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
        x1n = np.where(bt, x2, xe)
        f1n = np.where(bt, f2, fe)
        x2n = np.where(bt, xe, x1)
        f2n = np.where(bt, fe, f1)
        d = d * a2
        x1, x2, f1, f2 = x1n, x2n, f1n, f2n
    bt = f2 >= f1
    return np.where(bt, x2, x1), np.maximum(f1, f2)


def eval_renter(bp, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta):
    surplus = Rv - dc - bp
    f = Kr * np.maximum(surplus, 1e-10) ** oms / oms + pc + beta * interp_vector(b_grid, Vbar, bp)
    cm = surplus > cc
    if np.any(cm):
        ct = np.maximum(Rv[cm] - cb_c - ri * hRmax - bp[cm], 1e-10)
        f[cm] = (ct**alpha * ht_cap_c ** (1 - alpha)) ** oms / oms + pc + beta * interp_vector(b_grid, Vbar, bp[cm])
    f[surplus <= 1e-10] = -1e10
    return f


def golden_owner(lo, hi, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta, a1, a2, tol):
    d = hi - lo
    x1 = lo + a1 * d
    x2 = lo + a2 * d
    f1 = eval_owner(x1, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta)
    f2 = eval_owner(x2, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta)
    d = a1 * a2 * d
    while np.any(d > tol):
        bt = f2 >= f1
        xe = np.clip(np.where(bt, x2 + d, x1 - d), lo, hi)
        fe = eval_owner(xe, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta)
        x1n = np.where(bt, x2, xe)
        f1n = np.where(bt, f2, fe)
        x2n = np.where(bt, xe, x1)
        f2n = np.where(bt, fe, f1)
        d = d * a2
        x1, x2, f1, f2 = x1n, x2n, f1n, f2n
    bt = f2 >= f1
    return np.where(bt, x2, x1), np.maximum(f1, f2)


def eval_owner(bp, Rv, Vbar, b_grid, oc, cb_c, pc, Ko_c, alpha, oms, beta):
    ct = np.maximum(Rv - oc - cb_c - bp, 1e-10)
    f = Ko_c * ct ** (alpha * oms) / oms + pc + beta * interp_vector(b_grid, Vbar, bp)
    f[ct <= 1e-10] = -1e10
    return f


def forward_distribution(
    bp_pol: np.ndarray,
    hR_pol: np.ndarray,
    tenure_choice: np.ndarray,
    loc_probs: np.ndarray,
    fert_probs: np.ndarray,
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    fast_stats: bool = False,
    tenure_probs: np.ndarray | None = None,
) -> tuple[np.ndarray, SimpleNamespace]:
    J = P.J
    I = P.I
    Nb = len(b_grid)
    nt = 1 + P.n_house
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    bmin = b_grid[0]
    bmax = b_grid[-1]
    g = np.zeros((Nb, nt, I, J, npar, ncs))
    ie = int(np.argmax(b_grid >= P.b_entry_fixed)) if np.any(b_grid >= P.b_entry_fixed) else 0
    for i in range(I):
        g[ie, 0, i, 0, 0, 0] = P.entry_by_loc[i]
    total_births = 0.0
    births_by_loc = np.zeros(I)
    entrants_mature_by_loc = np.zeros(I)
    entrants_mature_total = 0.0
    event_horizon = 3
    birth_es3_pre_sum = birth_es3_post_sum = birth_es3_mass = 0.0
    addchild_es3_one_sum = addchild_es3_one_mass = 0.0
    addchild_es3_two_plus_sum = addchild_es3_two_plus_mass = 0.0
    onechild_es3_pre_sum = onechild_es3_post_sum = onechild_es3_mass = 0.0
    twoplus_es3_pre_sum = twoplus_es3_post_sum = twoplus_es3_mass = 0.0

    hc = np.zeros((I, nt))
    he = np.zeros((I, nt))
    for i in range(I):
        for ten in range(1, nt):
            hs = P.H_own[ten - 1]
            hc[i, ten] = p_hat[i] * hs
            he[i, ten] = (1 - P.psi) * p_hat[i] * hs

    phi_choice = SD.phi_choice
    use_lin = str(P.kfe_wealth_interp).lower() == "linear"
    if not use_lin:
        raise NotImplementedError("Only linear distribution interpolation is implemented in the Python port.")

    lmm_idx = np.zeros((I, nt, Nb), dtype=np.int64)
    lmm_wt = np.zeros((I, nt, Nb))
    for io in range(I):
        for to in range(nt):
            ba = np.clip(b_grid + he[io, to], bmin, bmax)
            lmm_idx[io, to, :], lmm_wt[io, to, :] = interp_indices(b_grid, ba)

    tmx_idx = np.zeros((I, nt, nt, npar, ncs, Nb), dtype=np.int64)
    tmx_wt = np.zeros((I, nt, nt, npar, ncs, Nb))
    for nn in range(npar):
        for cs in range(ncs):
            for id_ in range(I):
                for to in range(nt):
                    sp = he[id_, to]
                    for tn in range(nt):
                        pn = phi_choice[id_, tn, nn, cs] if tn > 0 else 1.0
                        if tn == to:
                            bf = b_grid.copy()
                        elif to == 0 and tn > 0:
                            bf = np.clip(np.maximum(b_grid - hc[id_, tn], -pn * hc[id_, tn]), bmin, bmax)
                        elif to > 0 and tn == 0:
                            bf = np.clip(np.maximum(b_grid + sp, 0.0), bmin, bmax)
                        else:
                            bf = np.clip(np.maximum(b_grid + sp - hc[id_, tn], -pn * hc[id_, tn]), bmin, bmax)
                        tmx_idx[id_, to, tn, nn, cs, :], tmx_wt[id_, to, tn, nn, cs, :] = interp_indices(b_grid, bf)

    K = P.n_child_stages
    csm1 = K + 1
    csm2 = K + 2
    ust = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    Pia = P.Pi_child if ust else None

    if (
        tenure_probs is None
        and fast_stats
        and NUMBA_AVAILABLE
        and bool(getattr(P, "use_compiled_forward_distribution", True))
    ):
        bp_idx, bp_wt = interp_indices(b_grid, np.clip(bp_pol, bmin, bmax))
        pia_arg = np.asarray(Pia if Pia is not None else np.zeros((ncs, ncs, npar)), dtype=float)
        g, total_births, births_by_loc, entrants_mature_by_loc, entrants_mature_total = forward_distribution_fast_kernel(
            np.asarray(P.entry_by_loc, dtype=float),
            fert_probs,
            loc_probs,
            tenure_choice.astype(np.int64, copy=False),
            bp_idx,
            bp_wt,
            lmm_idx,
            lmm_wt,
            tmx_idx,
            tmx_wt,
            pia_arg,
            Nb,
            nt,
            I,
            J,
            npar,
            ncs,
            ie,
            int(P.A_f_start),
            int(P.A_f_end),
            int(K),
            bool(ust),
        )
        tm = float(np.sum(g))
        if tm > 1e-12 and normalize_population_mass(P):
            sc = P.N_target / tm
            g *= sc
            total_births *= sc
            births_by_loc *= sc
            entrants_mature_by_loc *= sc
            entrants_mature_total *= sc

        stats = compute_eq_stats(g, P, b_grid, p_hat, hR_pol)
        stats.total_births_kfe = total_births
        stats.births_by_loc = births_by_loc
        stats.entry_by_loc = np.sum(g[:, :, :, 0, :, :], axis=(0, 1, 3, 4))
        stats.entry_rate = float(np.sum(g[:, :, :, 0, :, :]))
        stats.total_mass = float(np.sum(g))
        stats.entrants_mature_by_loc = entrants_mature_by_loc
        stats.entrants_mature_total = entrants_mature_total
        stats.mature_entry_shares = entrants_mature_by_loc / max(entrants_mature_total, 1e-12)
        stats.housing_increment_0to1_eventstudy_t3 = 0.0
        stats.housing_increment_1to2_proxy_t3 = 0.0
        stats.housing_increment_0to1_onechild_eventstudy_t3 = 0.0
        stats.housing_increment_0to2plus_eventstudy_t3 = 0.0
        stats.housing_event_horizon = 3
        return g, stats

    event_horizon = 3
    compute_event_stats = not fast_stats
    birth_es3_pre_sum = birth_es3_post_sum = birth_es3_mass = 0.0
    addchild_es3_one_sum = addchild_es3_one_mass = 0.0
    addchild_es3_two_plus_sum = addchild_es3_two_plus_mass = 0.0
    onechild_es3_pre_sum = onechild_es3_post_sum = onechild_es3_mass = 0.0
    twoplus_es3_pre_sum = twoplus_es3_post_sum = twoplus_es3_mass = 0.0

    for j in range(J - 1):
        if (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end):
            gc = g[:, :, :, j, 0, 0]
            pa = fert_probs[:, :, :, j, :]
            pv = np.arange(npar).reshape(1, 1, 1, npar)
            ba_j = gc * np.sum(pa * pv, axis=3)
            total_births += float(np.sum(ba_j))
            for i in range(I):
                births_by_loc[i] += float(np.sum(ba_j[:, :, i]))
            mbp = gc[:, :, :, None] * pa
            gpf = np.zeros((Nb, nt, I, npar, ncs))
            gpf[:, :, :, 0, 0] = mbp[:, :, :, 0]
            gpf[:, :, :, 1:, 1] = mbp[:, :, :, 1:]
            g[:, :, :, j, 0, 0] = 0.0
            g[:, :, :, j, :, :] += gpf

            birth_mass = np.sum(mbp[:, :, :, 1:], axis=3)
            birth_mass_total = float(np.sum(birth_mass))
            if compute_event_stats and birth_mass_total > 1e-12 and (j + event_horizon) < J:
                pre_h = mean_housing_childless_weighted(birth_mass, j, hR_pol, P)
                birth_cohort = np.zeros((Nb, nt, I, npar, ncs))
                birth_cohort[:, :, :, 1:, 1] = mbp[:, :, :, 1:]
                birth_cohort = advance_cohort_horizon(
                    birth_cohort,
                    j,
                    event_horizon,
                    loc_probs,
                    tenure_choice,
                    tenure_probs,
                    bp_pol,
                    P,
                    b_grid,
                    SD,
                    lmm_idx,
                    lmm_wt,
                    tmx_idx,
                    tmx_wt,
                    ust,
                    Pia,
                )
                post_h = mean_housing_distribution(birth_cohort, j + event_horizon, hR_pol, P)
                birth_es3_pre_sum += birth_mass_total * pre_h
                birth_es3_post_sum += birth_mass_total * post_h
                birth_es3_mass += birth_mass_total

                one_child = np.zeros_like(birth_cohort)
                one_child[:, :, :, 1, :] = birth_cohort[:, :, :, 1, :]
                one_mass = float(np.sum(one_child))
                if one_mass > 1e-12:
                    one_child_birth_mass = mbp[:, :, :, 1]
                    one_child_birth_mass_total = float(np.sum(one_child_birth_mass))
                    addchild_es3_one_sum += one_mass * mean_housing_distribution(one_child, j + event_horizon, hR_pol, P)
                    addchild_es3_one_mass += one_mass
                    if one_child_birth_mass_total > 1e-12:
                        one_pre = mean_housing_childless_weighted(one_child_birth_mass, j, hR_pol, P)
                        onechild_es3_pre_sum += one_child_birth_mass_total * one_pre
                        onechild_es3_post_sum += one_mass * mean_housing_distribution(one_child, j + event_horizon, hR_pol, P)
                        onechild_es3_mass += one_child_birth_mass_total

                if npar >= 3:
                    two_plus_birth_mass = np.sum(mbp[:, :, :, 2:], axis=3)
                    two_plus_birth_mass_total = float(np.sum(two_plus_birth_mass))
                    two_plus = np.zeros_like(birth_cohort)
                    two_plus[:, :, :, 2:, :] = birth_cohort[:, :, :, 2:, :]
                    two_plus_mass = float(np.sum(two_plus))
                    if two_plus_mass > 1e-12:
                        addchild_es3_two_plus_sum += two_plus_mass * mean_housing_distribution(two_plus, j + event_horizon, hR_pol, P)
                        addchild_es3_two_plus_mass += two_plus_mass
                        if two_plus_birth_mass_total > 1e-12:
                            two_pre = mean_housing_childless_weighted(two_plus_birth_mass, j, hR_pol, P)
                            twoplus_es3_pre_sum += two_plus_birth_mass_total * two_pre
                            twoplus_es3_post_sum += two_plus_mass * mean_housing_distribution(two_plus, j + event_horizon, hR_pol, P)
                            twoplus_es3_mass += two_plus_birth_mass_total

        gj = g[:, :, :, j, :, :]
        gpl = np.zeros((Nb, nt, I, npar, ncs))
        for io in range(I):
            for to in range(nt):
                go = flat_nc(gj[:, to, io, :, :], Nb, nc)
                if np.sum(go) < 1e-15:
                    continue
                po = np.reshape(loc_probs[:, to, io, :, j, :, :], (Nb, I, nc), order="F")
                sp = po[:, io, :]
                gpl[:, to, io, :, :] += unflat_nc(go * sp, Nb, npar, ncs)
                idx = lmm_idx[io, to, :]
                wt = lmm_wt[io, to, :]
                for id_ in range(I):
                    if id_ == io:
                        continue
                    mp = go * po[:, id_, :]
                    if use_compiled_scatter:
                        moved = scatter_cols_sameidx_kernel(idx, wt, mp, Nb)
                    else:
                        moved = scatter_redistribute_cols_sameidx(idx, wt, mp, Nb)
                    gpl[:, 0, id_, :, :] += unflat_nc(moved, Nb, npar, ncs)

        gpt = np.zeros((Nb, nt, I, npar, ncs))
        for nn in range(npar):
            for id_ in range(I):
                for to in range(nt):
                    gs = gpl[:, to, id_, nn, :]
                    if np.sum(gs) < 1e-15:
                        continue
                    for tn in range(nt):
                        if tenure_probs is None:
                            tcs = tenure_choice[:, to, id_, j, nn, :]
                            mk = tcs == tn
                            if not np.any(mk):
                                continue
                            mt = gs * mk
                        else:
                            pr = tenure_probs[:, to, id_, j, nn, :, tn]
                            mt = gs * pr
                        if np.sum(mt) < 1e-15:
                            continue
                        rd = np.zeros((Nb, ncs))
                        for cs in range(ncs):
                            idx = tmx_idx[id_, to, tn, nn, cs, :]
                            wt = tmx_wt[id_, to, tn, nn, cs, :]
                            if use_compiled_scatter:
                                rd[:, cs] = scatter_vec_kernel(idx, wt, mt[:, cs], Nb)
                            else:
                                rd[:, cs] = scatter_redistribute(idx, wt, mt[:, cs], Nb)
                        gpt[:, tn, id_, nn, :] += rd

        gps = np.zeros((Nb, nt, I, npar, ncs))
        for i in range(I):
            for ten in range(nt):
                gf = flat_nc(gpt[:, ten, i, :, :], Nb, nc)
                bpv = flat_nc(bp_pol[:, ten, i, j, :, :], Nb, nc)
                bpc = np.clip(bpv, bmin, bmax)
                idx, wt = interp_indices(b_grid, bpc)
                if use_compiled_scatter:
                    g_new = scatter_cols_kernel(idx, wt, gf, Nb)
                else:
                    g_new = scatter_redistribute_cols(idx, wt, gf, Nb)
                gps[:, ten, i, :, :] = unflat_nc(g_new, Nb, npar, ncs)

        for nn in range(npar):
            for cs in range(ncs):
                gp = gps[:, :, :, nn, cs]
                if ust:
                    Pi = Pia[:, :, nn]
                    if cs == K and nn >= 1:
                        pm = Pi[cs, csm1] if nn == 1 else Pi[cs, csm2]
                        if pm > 0:
                            nk = nn
                            for im in range(I):
                                fi = nk * pm * float(np.sum(gp[:, :, im]))
                                entrants_mature_by_loc[im] += fi
                                entrants_mature_total += fi
                    for csn in range(ncs):
                        wt = Pi[cs, csn]
                        if wt > 0:
                            g[:, :, :, j + 1, nn, csn] += wt * gp
                else:
                    if cs == 0:
                        csn = 0
                    elif cs >= csm1:
                        csn = cs
                    elif cs < K:
                        csn = cs + 1
                    else:
                        csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                    if cs == K and csn >= csm1 and nn >= 1:
                        nk = nn
                        for im in range(I):
                            fi = nk * float(np.sum(gp[:, :, im]))
                            entrants_mature_by_loc[im] += fi
                            entrants_mature_total += fi
                    g[:, :, :, j + 1, nn, csn] += gp

    tm = float(np.sum(g))
    if tm > 1e-12 and normalize_population_mass(P):
        sc = P.N_target / tm
        g *= sc
        total_births *= sc
        births_by_loc *= sc
        entrants_mature_by_loc *= sc
        entrants_mature_total *= sc

    stats = compute_eq_stats(g, P, b_grid, p_hat, hR_pol) if fast_stats else compute_statistics(g, fert_probs, loc_probs, P, b_grid, p_hat, hR_pol)
    stats.total_births_kfe = total_births
    stats.births_by_loc = births_by_loc
    stats.entry_by_loc = np.sum(g[:, :, :, 0, :, :], axis=(0, 1, 3, 4))
    stats.entry_rate = float(np.sum(g[:, :, :, 0, :, :]))
    stats.total_mass = float(np.sum(g))
    stats.entrants_mature_by_loc = entrants_mature_by_loc
    stats.entrants_mature_total = entrants_mature_total
    stats.mature_entry_shares = entrants_mature_by_loc / max(entrants_mature_total, 1e-12)
    stats.housing_increment_0to1_eventstudy_t3 = (
        birth_es3_post_sum / birth_es3_mass - birth_es3_pre_sum / birth_es3_mass if birth_es3_mass > 1e-12 else 0.0
    )
    stats.housing_increment_1to2_proxy_t3 = (
        addchild_es3_two_plus_sum / addchild_es3_two_plus_mass - addchild_es3_one_sum / addchild_es3_one_mass
        if addchild_es3_one_mass > 1e-12 and addchild_es3_two_plus_mass > 1e-12
        else 0.0
    )
    stats.housing_increment_0to1_onechild_eventstudy_t3 = (
        onechild_es3_post_sum / onechild_es3_mass - onechild_es3_pre_sum / onechild_es3_mass if onechild_es3_mass > 1e-12 else 0.0
    )
    stats.housing_increment_0to2plus_eventstudy_t3 = (
        twoplus_es3_post_sum / twoplus_es3_mass - twoplus_es3_pre_sum / twoplus_es3_mass if twoplus_es3_mass > 1e-12 else 0.0
    )
    stats.housing_event_horizon = event_horizon
    return g, stats


def forward_distribution_markov_income(
    bp_pol: np.ndarray,
    hR_pol: np.ndarray,
    tenure_choice: np.ndarray,
    loc_probs: np.ndarray,
    fert_probs: np.ndarray,
    r_hat: np.ndarray,
    p_hat: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    SD: SimpleNamespace,
    fast_stats: bool = False,
    tenure_probs: np.ndarray | None = None,
) -> tuple[np.ndarray, SimpleNamespace]:
    J = P.J
    I = P.I
    Nb = len(b_grid)
    nt = 1 + P.n_house
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    z_grid, z_weights, Pi_z = income_transition_values(P)
    Nz = len(z_grid)
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    bmin = b_grid[0]
    bmax = b_grid[-1]
    g = np.zeros((Nb, nt, I, J, Nz, npar, ncs))
    ie = int(np.argmax(b_grid >= P.b_entry_fixed)) if np.any(b_grid >= P.b_entry_fixed) else 0
    for i in range(I):
        for zz in range(Nz):
            g[ie, 0, i, 0, zz, 0, 0] = P.entry_by_loc[i] * z_weights[zz]
    total_births = 0.0
    births_by_loc = np.zeros(I)
    entrants_mature_by_loc = np.zeros(I)
    entrants_mature_total = 0.0
    event_horizon = 3
    birth_es3_pre_sum = birth_es3_post_sum = birth_es3_mass = 0.0
    addchild_es3_one_sum = addchild_es3_one_mass = 0.0
    addchild_es3_two_plus_sum = addchild_es3_two_plus_mass = 0.0
    onechild_es3_pre_sum = onechild_es3_post_sum = onechild_es3_mass = 0.0
    twoplus_es3_pre_sum = twoplus_es3_post_sum = twoplus_es3_mass = 0.0

    hc = np.zeros((I, nt))
    he = np.zeros((I, nt))
    for i in range(I):
        for ten in range(1, nt):
            hs = P.H_own[ten - 1]
            hc[i, ten] = p_hat[i] * hs
            he[i, ten] = (1 - P.psi) * p_hat[i] * hs

    phi_choice = SD.phi_choice
    lmm_idx = np.zeros((I, nt, Nb), dtype=np.int64)
    lmm_wt = np.zeros((I, nt, Nb))
    for io in range(I):
        for to in range(nt):
            ba = np.clip(b_grid + he[io, to], bmin, bmax)
            lmm_idx[io, to, :], lmm_wt[io, to, :] = interp_indices(b_grid, ba)

    tmx_idx = np.zeros((I, nt, nt, npar, ncs, Nb), dtype=np.int64)
    tmx_wt = np.zeros((I, nt, nt, npar, ncs, Nb))
    for nn in range(npar):
        for cs in range(ncs):
            for id_ in range(I):
                for to in range(nt):
                    sp = he[id_, to]
                    for tn in range(nt):
                        pn = phi_choice[id_, tn, nn, cs] if tn > 0 else 1.0
                        if tn == to:
                            bf = b_grid.copy()
                        elif to == 0 and tn > 0:
                            bf = np.clip(np.maximum(b_grid - hc[id_, tn], -pn * hc[id_, tn]), bmin, bmax)
                        elif to > 0 and tn == 0:
                            bf = np.clip(np.maximum(b_grid + sp, 0.0), bmin, bmax)
                        else:
                            bf = np.clip(np.maximum(b_grid + sp - hc[id_, tn], -pn * hc[id_, tn]), bmin, bmax)
                        tmx_idx[id_, to, tn, nn, cs, :], tmx_wt[id_, to, tn, nn, cs, :] = interp_indices(b_grid, bf)

    K = P.n_child_stages
    csm1 = K + 1
    csm2 = K + 2
    ust = bool(P.use_stochastic_aging and hasattr(P, "Pi_child"))
    Pia = P.Pi_child if ust else None

    for j in range(J - 1):
        if (j + 1 >= P.A_f_start) and (j + 1 <= P.A_f_end):
            for zz in range(Nz):
                gc = g[:, :, :, j, zz, 0, 0]
                pa = fert_probs[:, :, :, j, zz, :]
                pv = np.arange(npar).reshape(1, 1, 1, npar)
                ba_j = gc * np.sum(pa * pv, axis=3)
                total_births += float(np.sum(ba_j))
                for i in range(I):
                    births_by_loc[i] += float(np.sum(ba_j[:, :, i]))
                mbp = gc[:, :, :, None] * pa
                gpf = np.zeros((Nb, nt, I, npar, ncs))
                gpf[:, :, :, 0, 0] = mbp[:, :, :, 0]
                gpf[:, :, :, 1:, 1] = mbp[:, :, :, 1:]
                g[:, :, :, j, zz, 0, 0] = 0.0
                g[:, :, :, j, zz, :, :] += gpf

                birth_mass = np.sum(mbp[:, :, :, 1:], axis=3)
                birth_mass_total = float(np.sum(birth_mass))
                if not fast_stats and birth_mass_total > 1e-12 and (j + event_horizon) < J:
                    birth_weight = np.zeros((Nb, nt, I, Nz))
                    birth_weight[:, :, :, zz] = birth_mass
                    pre_h = mean_housing_childless_weighted_markov(birth_weight, j, hR_pol, P)
                    birth_cohort = np.zeros((Nb, nt, I, Nz, npar, ncs))
                    birth_cohort[:, :, :, zz, 1:, 1] = mbp[:, :, :, 1:]
                    birth_cohort = advance_cohort_horizon_markov_income(
                        birth_cohort,
                        j,
                        event_horizon,
                        loc_probs,
                        tenure_choice,
                        tenure_probs,
                        bp_pol,
                        P,
                        b_grid,
                        SD,
                        lmm_idx,
                        lmm_wt,
                        tmx_idx,
                        tmx_wt,
                        ust,
                        Pia,
                        Pi_z,
                    )
                    post_h = mean_housing_distribution_markov(birth_cohort, j + event_horizon, hR_pol, P)
                    birth_es3_pre_sum += birth_mass_total * pre_h
                    birth_es3_post_sum += birth_mass_total * post_h
                    birth_es3_mass += birth_mass_total

                    one_child = np.zeros_like(birth_cohort)
                    one_child[:, :, :, :, 1, :] = birth_cohort[:, :, :, :, 1, :]
                    one_mass = float(np.sum(one_child))
                    if one_mass > 1e-12:
                        one_child_birth_mass = mbp[:, :, :, 1]
                        one_child_birth_mass_total = float(np.sum(one_child_birth_mass))
                        addchild_es3_one_sum += one_mass * mean_housing_distribution_markov(
                            one_child, j + event_horizon, hR_pol, P
                        )
                        addchild_es3_one_mass += one_mass
                        if one_child_birth_mass_total > 1e-12:
                            one_weight = np.zeros((Nb, nt, I, Nz))
                            one_weight[:, :, :, zz] = one_child_birth_mass
                            one_pre = mean_housing_childless_weighted_markov(one_weight, j, hR_pol, P)
                            onechild_es3_pre_sum += one_child_birth_mass_total * one_pre
                            onechild_es3_post_sum += one_mass * mean_housing_distribution_markov(
                                one_child, j + event_horizon, hR_pol, P
                            )
                            onechild_es3_mass += one_child_birth_mass_total

                    if npar >= 3:
                        two_plus_birth_mass = np.sum(mbp[:, :, :, 2:], axis=3)
                        two_plus_birth_mass_total = float(np.sum(two_plus_birth_mass))
                        two_plus = np.zeros_like(birth_cohort)
                        two_plus[:, :, :, :, 2:, :] = birth_cohort[:, :, :, :, 2:, :]
                        two_plus_mass = float(np.sum(two_plus))
                        if two_plus_mass > 1e-12:
                            addchild_es3_two_plus_sum += two_plus_mass * mean_housing_distribution_markov(
                                two_plus, j + event_horizon, hR_pol, P
                            )
                            addchild_es3_two_plus_mass += two_plus_mass
                            if two_plus_birth_mass_total > 1e-12:
                                two_weight = np.zeros((Nb, nt, I, Nz))
                                two_weight[:, :, :, zz] = two_plus_birth_mass
                                two_pre = mean_housing_childless_weighted_markov(two_weight, j, hR_pol, P)
                                twoplus_es3_pre_sum += two_plus_birth_mass_total * two_pre
                                twoplus_es3_post_sum += two_plus_mass * mean_housing_distribution_markov(
                                    two_plus, j + event_horizon, hR_pol, P
                                )
                                twoplus_es3_mass += two_plus_birth_mass_total

        gj = g[:, :, :, j, :, :, :]
        gpl = np.zeros((Nb, nt, I, Nz, npar, ncs))
        for zz in range(Nz):
            for io in range(I):
                for to in range(nt):
                    go = flat_nc(gj[:, to, io, zz, :, :], Nb, nc)
                    if np.sum(go) < 1e-15:
                        continue
                    po = np.reshape(loc_probs[:, to, io, :, j, zz, :, :], (Nb, I, nc), order="F")
                    sp = po[:, io, :]
                    gpl[:, to, io, zz, :, :] += unflat_nc(go * sp, Nb, npar, ncs)
                    idx = lmm_idx[io, to, :]
                    wt = lmm_wt[io, to, :]
                    for id_ in range(I):
                        if id_ == io:
                            continue
                        mp = go * po[:, id_, :]
                        if use_compiled_scatter:
                            moved = scatter_cols_sameidx_kernel(idx, wt, mp, Nb)
                        else:
                            moved = scatter_redistribute_cols_sameidx(idx, wt, mp, Nb)
                        gpl[:, 0, id_, zz, :, :] += unflat_nc(moved, Nb, npar, ncs)

        gpt = np.zeros((Nb, nt, I, Nz, npar, ncs))
        for zz in range(Nz):
            for nn in range(npar):
                for id_ in range(I):
                    for to in range(nt):
                        gs = gpl[:, to, id_, zz, nn, :]
                        if np.sum(gs) < 1e-15:
                            continue
                        for tn in range(nt):
                            if tenure_probs is None:
                                tcs = tenure_choice[:, to, id_, j, zz, nn, :]
                                mk = tcs == tn
                                if not np.any(mk):
                                    continue
                                mt = gs * mk
                            else:
                                pr = tenure_probs[:, to, id_, j, zz, nn, :, tn]
                                mt = gs * pr
                            if np.sum(mt) < 1e-15:
                                continue
                            rd = np.zeros((Nb, ncs))
                            for cs in range(ncs):
                                idx = tmx_idx[id_, to, tn, nn, cs, :]
                                wt = tmx_wt[id_, to, tn, nn, cs, :]
                                if use_compiled_scatter:
                                    rd[:, cs] = scatter_vec_kernel(idx, wt, mt[:, cs], Nb)
                                else:
                                    rd[:, cs] = scatter_redistribute(idx, wt, mt[:, cs], Nb)
                            gpt[:, tn, id_, zz, nn, :] += rd

        gps = np.zeros((Nb, nt, I, Nz, npar, ncs))
        for zz in range(Nz):
            for i in range(I):
                for ten in range(nt):
                    gf = flat_nc(gpt[:, ten, i, zz, :, :], Nb, nc)
                    bpv = flat_nc(bp_pol[:, ten, i, j, zz, :, :], Nb, nc)
                    bpc = np.clip(bpv, bmin, bmax)
                    idx, wt = interp_indices(b_grid, bpc)
                    if use_compiled_scatter:
                        g_new = scatter_cols_kernel(idx, wt, gf, Nb)
                    else:
                        g_new = scatter_redistribute_cols(idx, wt, gf, Nb)
                    gps[:, ten, i, zz, :, :] = unflat_nc(g_new, Nb, npar, ncs)

        for zz in range(Nz):
            for nn in range(npar):
                for cs in range(ncs):
                    gp = gps[:, :, :, zz, nn, cs]
                    if ust:
                        Pi = Pia[:, :, nn]
                        if cs == K and nn >= 1:
                            pm = Pi[cs, csm1] if nn == 1 else Pi[cs, csm2]
                            if pm > 0:
                                nk = nn
                                for im in range(I):
                                    fi = nk * pm * float(np.sum(gp[:, :, im]))
                                    entrants_mature_by_loc[im] += fi
                                    entrants_mature_total += fi
                        for csn in range(ncs):
                            wt_child = Pi[cs, csn]
                            if wt_child > 0:
                                for zn in range(Nz):
                                    g[:, :, :, j + 1, zn, nn, csn] += Pi_z[zz, zn] * wt_child * gp
                    else:
                        if cs == 0:
                            csn = 0
                        elif cs >= csm1:
                            csn = cs
                        elif cs < K:
                            csn = cs + 1
                        else:
                            csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                        if cs == K and csn >= csm1 and nn >= 1:
                            nk = nn
                            for im in range(I):
                                fi = nk * float(np.sum(gp[:, :, im]))
                                entrants_mature_by_loc[im] += fi
                                entrants_mature_total += fi
                        for zn in range(Nz):
                            g[:, :, :, j + 1, zn, nn, csn] += Pi_z[zz, zn] * gp

    tm = float(np.sum(g))
    if tm > 1e-12 and normalize_population_mass(P):
        sc = P.N_target / tm
        g *= sc
        total_births *= sc
        births_by_loc *= sc
        entrants_mature_by_loc *= sc
        entrants_mature_total *= sc

    if fast_stats:
        stats = compute_markov_eq_stats(g, P, b_grid, p_hat, hR_pol)
    else:
        stats = compute_markov_statistics(g, fert_probs, loc_probs, P, b_grid, p_hat, hR_pol)
    stats.total_births_kfe = total_births
    stats.births_by_loc = births_by_loc
    stats.entry_by_loc = np.sum(g[:, :, :, 0, :, :, :], axis=(0, 1, 3, 4, 5))
    stats.entry_rate = float(np.sum(g[:, :, :, 0, :, :, :]))
    stats.total_mass = float(np.sum(g))
    stats.entrants_mature_by_loc = entrants_mature_by_loc
    stats.entrants_mature_total = entrants_mature_total
    stats.mature_entry_shares = entrants_mature_by_loc / max(entrants_mature_total, 1e-12)
    stats.housing_increment_0to1_eventstudy_t3 = (
        birth_es3_post_sum / birth_es3_mass - birth_es3_pre_sum / birth_es3_mass if birth_es3_mass > 1e-12 else 0.0
    )
    stats.housing_increment_1to2_proxy_t3 = (
        addchild_es3_two_plus_sum / addchild_es3_two_plus_mass - addchild_es3_one_sum / addchild_es3_one_mass
        if addchild_es3_one_mass > 1e-12 and addchild_es3_two_plus_mass > 1e-12
        else 0.0
    )
    stats.housing_increment_0to1_onechild_eventstudy_t3 = (
        onechild_es3_post_sum / onechild_es3_mass - onechild_es3_pre_sum / onechild_es3_mass
        if onechild_es3_mass > 1e-12
        else 0.0
    )
    stats.housing_increment_0to2plus_eventstudy_t3 = (
        twoplus_es3_post_sum / twoplus_es3_mass - twoplus_es3_pre_sum / twoplus_es3_mass
        if twoplus_es3_mass > 1e-12
        else 0.0
    )
    stats.housing_event_horizon = event_horizon
    return g, stats


def collapse_markov_policy(policy: np.ndarray, g: np.ndarray, z_weights: np.ndarray) -> np.ndarray:
    fallback = np.tensordot(policy, z_weights, axes=([4], [0]))
    den = np.sum(g, axis=4)
    num = np.sum(policy * g, axis=4)
    out = fallback.copy()
    mask = den > 1e-15
    out[mask] = num[mask] / den[mask]
    return out


def collapse_markov_fertility_probs(fp: np.ndarray, g: np.ndarray, z_weights: np.ndarray) -> np.ndarray:
    fallback = np.tensordot(fp, z_weights, axes=([4], [0]))
    mass = g[:, :, :, :, :, 0, 0]
    den = np.sum(mass, axis=4)
    num = np.sum(fp * mass[:, :, :, :, :, None], axis=4)
    out = fallback.copy()
    mask = den > 1e-15
    out[mask, :] = num[mask, :] / den[mask, None]
    return out


def collapse_markov_location_probs(lp: np.ndarray, g: np.ndarray, z_weights: np.ndarray) -> np.ndarray:
    fallback = np.tensordot(lp, z_weights, axes=([5], [0]))
    mass = g[:, :, :, None, :, :, :, :]
    den = np.sum(mass, axis=5)
    num = np.sum(lp * mass, axis=5)
    return np.divide(num, den, out=fallback.copy(), where=den > 1e-15)


def advance_cohort_horizon_markov_income(
    g_in,
    start_age,
    horizon,
    loc_probs,
    tenure_choice,
    tenure_probs,
    bp_pol,
    P,
    b_grid,
    SD,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    ust,
    Pia,
    Pi_z,
):
    g_out = g_in
    for step in range(1, horizon + 1):
        age_idx = start_age + step - 1
        if age_idx >= P.J - 1:
            break
        g_out = advance_cohort_one_period_markov_income(
            g_out,
            age_idx,
            loc_probs,
            tenure_choice,
            tenure_probs,
            bp_pol,
            P,
            b_grid,
            SD,
            lmm_idx,
            lmm_wt,
            tmx_idx,
            tmx_wt,
            ust,
            Pia,
            Pi_z,
        )
    return g_out


def advance_cohort_one_period_markov_income(
    gj,
    j,
    loc_probs,
    tenure_choice,
    tenure_probs,
    bp_pol,
    P,
    b_grid,
    SD,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    ust,
    Pia,
    Pi_z,
):
    Nb = len(b_grid)
    nt = 1 + P.n_house
    I = P.I
    Nz = gj.shape[3]
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    K = P.n_child_stages
    csm1 = K + 1
    csm2 = K + 2

    gpl = np.zeros((Nb, nt, I, Nz, npar, ncs))
    for zz in range(Nz):
        for io in range(I):
            for to in range(nt):
                go = flat_nc(gj[:, to, io, zz, :, :], Nb, nc)
                if np.sum(go) < 1e-15:
                    continue
                po = np.reshape(loc_probs[:, to, io, :, j, zz, :, :], (Nb, I, nc), order="F")
                sp = po[:, io, :]
                gpl[:, to, io, zz, :, :] += unflat_nc(go * sp, Nb, npar, ncs)
                idx = lmm_idx[io, to, :]
                wt = lmm_wt[io, to, :]
                for id_ in range(I):
                    if id_ == io:
                        continue
                    mp = go * po[:, id_, :]
                    if use_compiled_scatter:
                        moved = scatter_cols_sameidx_kernel(idx, wt, mp, Nb)
                    else:
                        moved = scatter_redistribute_cols_sameidx(idx, wt, mp, Nb)
                    gpl[:, 0, id_, zz, :, :] += unflat_nc(moved, Nb, npar, ncs)

    gpt = np.zeros((Nb, nt, I, Nz, npar, ncs))
    for zz in range(Nz):
        for nn in range(npar):
            for id_ in range(I):
                for to in range(nt):
                    gs = gpl[:, to, id_, zz, nn, :]
                    if np.sum(gs) < 1e-15:
                        continue
                    for tn in range(nt):
                        if tenure_probs is None:
                            tcs = tenure_choice[:, to, id_, j, zz, nn, :]
                            mk = tcs == tn
                            if not np.any(mk):
                                continue
                            mt = gs * mk
                        else:
                            pr = tenure_probs[:, to, id_, j, zz, nn, :, tn]
                            mt = gs * pr
                        if np.sum(mt) < 1e-15:
                            continue
                        rd = np.zeros((Nb, ncs))
                        for cs in range(ncs):
                            idx = tmx_idx[id_, to, tn, nn, cs, :]
                            wt = tmx_wt[id_, to, tn, nn, cs, :]
                            if use_compiled_scatter:
                                rd[:, cs] = scatter_vec_kernel(idx, wt, mt[:, cs], Nb)
                            else:
                                rd[:, cs] = scatter_redistribute(idx, wt, mt[:, cs], Nb)
                        gpt[:, tn, id_, zz, nn, :] += rd

    gps = np.zeros((Nb, nt, I, Nz, npar, ncs))
    for zz in range(Nz):
        for i in range(I):
            for ten in range(nt):
                gf = flat_nc(gpt[:, ten, i, zz, :, :], Nb, nc)
                bpv = flat_nc(bp_pol[:, ten, i, j, zz, :, :], Nb, nc)
                idx, wt = interp_indices(b_grid, np.clip(bpv, b_grid[0], b_grid[-1]))
                if use_compiled_scatter:
                    g_new = scatter_cols_kernel(idx, wt, gf, Nb)
                else:
                    g_new = scatter_redistribute_cols(idx, wt, gf, Nb)
                gps[:, ten, i, zz, :, :] = unflat_nc(g_new, Nb, npar, ncs)

    g_next = np.zeros_like(gj)
    for zz in range(Nz):
        for nn in range(npar):
            for cs in range(ncs):
                gp = gps[:, :, :, zz, nn, cs]
                if ust:
                    Pi = Pia[:, :, nn]
                    for csn in range(ncs):
                        wt_child = Pi[cs, csn]
                        if wt_child > 0:
                            for zn in range(Nz):
                                g_next[:, :, :, zn, nn, csn] += Pi_z[zz, zn] * wt_child * gp
                else:
                    if cs == 0:
                        csn = 0
                    elif cs >= csm1:
                        csn = cs
                    elif cs < K:
                        csn = cs + 1
                    else:
                        csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                    for zn in range(Nz):
                        g_next[:, :, :, zn, nn, csn] += Pi_z[zz, zn] * gp
    return g_next


def mean_housing_childless_weighted_markov(weight_dist, j, hR_pol, P):
    Nb, nt, I, Nz = weight_dist.shape
    th = mn = 0.0
    for zz in range(Nz):
        for i in range(I):
            for ten in range(nt):
                gs = weight_dist[:, ten, i, zz]
                mh = float(np.sum(gs))
                if mh < 1e-15:
                    continue
                if ten == 0:
                    th += float(np.sum(gs * hR_pol[:, ten, i, j, zz, 0, 0]))
                else:
                    th += mh * P.H_own[ten - 1]
                mn += mh
    return th / max(mn, 1e-12)


def mean_housing_distribution_markov(g_dist, j, hR_pol, P):
    Nb, nt, I, Nz, npar, ncs = g_dist.shape
    th = mn = 0.0
    for zz in range(Nz):
        for nn in range(npar):
            for cs in range(ncs):
                for i in range(I):
                    for ten in range(nt):
                        gs = g_dist[:, ten, i, zz, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        if ten == 0:
                            th += float(np.sum(gs * hR_pol[:, ten, i, j, zz, nn, cs]))
                        else:
                            th += mh * P.H_own[ten - 1]
                        mn += mh
    return th / max(mn, 1e-12)


def compute_markov_statistics(
    g: np.ndarray,
    fp: np.ndarray,
    lp: np.ndarray,
    P: SimpleNamespace,
    bg: np.ndarray,
    ph: np.ndarray,
    hR: np.ndarray,
) -> SimpleNamespace:
    z_grid, z_weights, Pi_z = income_transition_values(P)
    g_total = np.sum(g, axis=4)
    hR_total = collapse_markov_policy(hR, g, z_weights)
    fp_total = collapse_markov_fertility_probs(fp, g, z_weights)
    lp_total = collapse_markov_location_probs(lp, g, z_weights)
    stats = compute_statistics(g_total, fp_total, lp_total, P, bg, ph, hR_total)
    Nz = len(z_grid)
    nt = 1 + P.n_house
    stats.income_state_mass = np.zeros(Nz)
    stats.own_rate_by_income_type = np.zeros(Nz)
    stats.mean_fertility_by_income_type = np.zeros(Nz)
    stats.housing_demand_by_income_type = np.zeros((Nz, P.I))
    for zz in range(Nz):
        gz = g[:, :, :, :, zz, :, :]
        mz = float(np.sum(gz))
        stats.income_state_mass[zz] = mz
        stats.own_rate_by_income_type[zz] = float(np.sum(gz[:, 1:, :, :, :, :]) / max(mz, 1e-12))
        mp = float(np.sum(g[:, :, :, P.A_f_end :, zz, :, :]))
        if mp > 1e-12:
            mean_n = 0.0
            for nn in range(P.n_parity):
                mean_n += nn * float(np.sum(g[:, :, :, P.A_f_end :, zz, nn, :])) / mp
            stats.mean_fertility_by_income_type[zz] = mean_n
        for i in range(P.I):
            Hd = 0.0
            for j in range(P.J):
                for nn in range(P.n_parity):
                    for cs in range(P.n_child_states):
                        Hd += float(np.sum(g[:, 0, i, j, zz, nn, cs] * hR[:, 0, i, j, zz, nn, cs]))
                        for ten in range(1, nt):
                            Hd += float(np.sum(g[:, ten, i, j, zz, nn, cs]) * P.H_own[ten - 1])
            stats.housing_demand_by_income_type[zz, i] = Hd / max(housing_demand_normalizer(P), 1e-12)

    worker_income = worker_mass = 0.0
    young_income = young_mass = young_liquid = 0.0
    payroll_tax_revenue = 0.0
    period_scale = float(getattr(P, "period_years", getattr(P, "da", 1.0))) if bool(getattr(P, "scale_flows_to_period", False)) else 1.0
    a25s = age_to_index(P, 25)
    aye = age_to_index(P, 35)
    for j in range(P.J):
        for i in range(P.I):
            for zz, z_value in enumerate(z_grid):
                yj = income_at_state(P, i, j, float(z_value))
                mass = float(np.sum(g[:, :, i, j, zz, :, :]))
                if j < P.J_R:
                    worker_income += yj * mass
                    worker_mass += mass
                    payroll_tax_revenue += period_scale * P.tau_pay * P.w_hat[i] * P.income_age_profile[j] * float(z_value) * mass
                if a25s <= j <= aye:
                    gm = g[:, 0, i, j, zz, 0, 0]
                    mh = float(np.sum(gm))
                    if mh > 1e-15:
                        young_income += yj * mh
                        young_mass += mh
                        young_liquid += float(np.sum(gm * bg))
    stats.mean_income = worker_income / max(worker_mass, 1e-12)
    stats.young_childless_renter_income = young_income / max(young_mass, 1e-12)
    stats.young_liquid_wealth = young_liquid / max(young_mass, 1e-12)
    stats.young_liquid_wealth_to_income = young_liquid / max(young_income, 1e-12)
    stats.wealth_to_income = getattr(stats, "mean_wealth_4555", 0.0) / max(stats.mean_income, 1e-12)
    stats.liquid_wealth_to_income = getattr(stats, "liquid_wealth_4555", 0.0) / max(stats.mean_income, 1e-12)
    stats.payroll_tax_revenue = payroll_tax_revenue
    stats.pension_outlays = P.pension * stats.retiree_mass_total
    stats.pension_budget_residual = stats.payroll_tax_revenue - stats.pension_outlays
    stats.implied_balanced_pension = stats.payroll_tax_revenue / max(stats.retiree_mass_total, 1e-12)
    stats.income_transition = Pi_z.copy()
    return stats


def compute_markov_eq_stats(g: np.ndarray, P: SimpleNamespace, bg: np.ndarray, ph: np.ndarray, hR: np.ndarray) -> SimpleNamespace:
    """Minimal Markov-income statistics needed for price clearing."""

    _ = bg, ph
    nt = 1 + P.n_house
    npar = P.n_parity
    tm = float(np.sum(g))
    stats = SimpleNamespace()
    stats.own_rate = float(np.sum(g[:, 1:, :, :, :, :, :]) / max(tm, 1e-12))
    stats.pop_share = np.zeros(P.I)
    stats.housing_demand = np.zeros(P.I)
    norm = housing_demand_normalizer(P)
    for i in range(P.I):
        gi = g[:, :, i, :, :, :, :]
        stats.pop_share[i] = float(np.sum(gi)) / max(tm, 1e-12)
        renter_demand = float(np.sum(gi[:, 0, :, :, :, :] * hR[:, 0, i, :, :, :, :]))
        owner_demand = 0.0
        for ten in range(1, nt):
            owner_demand += float(np.sum(gi[:, ten, :, :, :, :])) * float(P.H_own[ten - 1])
        stats.housing_demand[i] = (renter_demand + owner_demand) / max(norm, 1e-12)
    mp = float(np.sum(g[:, :, :, P.A_f_end :, :, :, :]))
    stats.parity_dist = np.zeros(npar)
    for nn in range(npar):
        stats.parity_dist[nn] = np.sum(g[:, :, :, P.A_f_end :, :, nn, :]) / max(mp, 1e-12)
    stats.mean_parity = float(np.sum(np.arange(npar) * stats.parity_dist))
    return stats


def compute_eq_stats(g: np.ndarray, P: SimpleNamespace, bg: np.ndarray, ph: np.ndarray, hR: np.ndarray) -> SimpleNamespace:
    J = P.J
    I = P.I
    nt = 1 + P.n_house
    npar = P.n_parity
    ncs = P.n_child_states
    tm = float(np.sum(g))
    stats = SimpleNamespace()
    stats.own_rate = float(np.sum(g[:, 1:, :, :, :, :]) / max(tm, 1e-12))
    stats.pop_share = np.zeros(I)
    stats.housing_demand = np.zeros(I)
    for i in range(I):
        pi = float(np.sum(g[:, :, i, :, :, :]))
        stats.pop_share[i] = pi / max(tm, 1e-12)
        Hd = 0.0
        for j in range(J):
            for nn in range(npar):
                for cs in range(ncs):
                    Hd += float(np.sum(g[:, 0, i, j, nn, cs] * hR[:, 0, i, j, nn, cs]))
                    for ten in range(1, nt):
                        Hd += float(np.sum(g[:, ten, i, j, nn, cs]) * P.H_own[ten - 1])
        stats.housing_demand[i] = Hd / housing_demand_normalizer(P)
    mp = float(np.sum(g[:, :, :, P.A_f_end :, :, :]))
    stats.parity_dist = np.zeros(npar)
    for nn in range(npar):
        stats.parity_dist[nn] = np.sum(g[:, :, :, P.A_f_end :, nn, :]) / max(mp, 1e-12)
    stats.mean_parity = float(np.sum(np.arange(npar) * stats.parity_dist))
    return stats


def age_to_index(P: SimpleNamespace, age: float) -> int:
    idx = int(round((float(age) - float(P.age_start)) / max(float(P.da), 1e-12)))
    return int(np.clip(idx, 0, P.J - 1))


def compute_statistics(g: np.ndarray, fp: np.ndarray, lp: np.ndarray, P: SimpleNamespace, bg: np.ndarray, ph: np.ndarray, hR: np.ndarray) -> SimpleNamespace:
    J = P.J
    I = P.I
    Nb = len(bg)
    nt = 1 + P.n_house
    npar = P.n_parity
    ncs = P.n_child_states
    tm = float(np.sum(g))
    stats = SimpleNamespace()
    stats.own_rate = float(np.sum(g[:, 1:, :, :, :, :]) / max(tm, 1e-12))
    stats.pop_share = np.zeros(I)
    stats.own_by_loc = np.zeros(I)
    stats.housing_demand = np.zeros(I)
    for i in range(I):
        pi = float(np.sum(g[:, :, i, :, :, :]))
        stats.pop_share[i] = pi / max(tm, 1e-12)
        stats.own_by_loc[i] = np.sum(g[:, 1:, i, :, :, :]) / max(pi, 1e-12)
        Hd = 0.0
        for j in range(J):
            for nn in range(npar):
                for cs in range(ncs):
                    Hd += float(np.sum(g[:, 0, i, j, nn, cs] * hR[:, 0, i, j, nn, cs]))
                    for ten in range(1, nt):
                        Hd += float(np.sum(g[:, ten, i, j, nn, cs]) * P.H_own[ten - 1])
        stats.housing_demand[i] = Hd / housing_demand_normalizer(P)

    stats.worker_mass_by_loc = np.array([np.sum(g[:, :, i, : P.J_R, :, :]) for i in range(I)])
    stats.retiree_mass_by_loc = np.array([np.sum(g[:, :, i, P.J_R :, :, :]) for i in range(I)])
    stats.worker_mass_total = float(np.sum(stats.worker_mass_by_loc))
    stats.retiree_mass_total = float(np.sum(stats.retiree_mass_by_loc))
    stats.parity_dist = np.zeros(npar)
    mp = float(np.sum(g[:, :, :, P.A_f_end :, :, :]))
    for nn in range(npar):
        stats.parity_dist[nn] = np.sum(g[:, :, :, P.A_f_end :, nn, :]) / max(mp, 1e-12)
    stats.mean_parity = float(np.sum(np.arange(npar) * stats.parity_dist))
    stats.own_by_parity = np.zeros(npar)
    for nn in range(npar):
        mn = float(np.sum(g[:, :, :, :, nn, :]))
        stats.own_by_parity[nn] = np.sum(g[:, 1:, :, :, nn, :]) / max(mn, 1e-12)
    stats.mean_parity_by_loc = np.zeros(I)
    stats.frac_childless_by_loc = np.zeros(I)
    for i in range(I):
        mip = float(np.sum(g[:, :, i, P.A_f_end :, :, :]))
        if mip > 1e-12:
            stats.mean_parity_by_loc[i] = sum(
                nn * np.sum(g[:, :, i, P.A_f_end :, nn, :]) / mip for nn in range(npar)
            )
            stats.frac_childless_by_loc[i] = np.sum(g[:, :, i, P.A_f_end :, 0, :]) / mip
    stats.own_by_age = np.zeros(J)
    for jj in range(J):
        gj = g[:, :, :, jj, :, :]
        mj = float(np.sum(gj))
        if mj > 1e-12:
            stats.own_by_age[jj] = np.sum(gj[:, 1:, :, :, :]) / mj
    stats.child_state_dist = np.zeros((J, ncs))
    for jj in range(J):
        mj = float(np.sum(g[:, :, :, jj, :, :]))
        if mj > 1e-12:
            for cs in range(ncs):
                stats.child_state_dist[jj, cs] = np.sum(g[:, :, :, jj, :, cs]) / mj
    stats.fert_by_age = np.zeros(J)
    for j in range(P.A_f_start - 1, P.A_f_end):
        mj = float(np.sum(g[:, :, :, j, 0, 0]))
        if mj > 1e-12:
            En = 0.0
            for i in range(I):
                for ten in range(nt):
                    gs = g[:, ten, i, j, 0, 0]
                    nz = gs > 1e-15
                    if not np.any(nz):
                        continue
                    pr = fp[:, ten, i, j, :]
                    En += float(np.sum(gs[nz] * (pr[nz, :] @ np.arange(npar))))
            stats.fert_by_age[j] = En / mj

    a22s = age_to_index(P, 22)
    a25s = age_to_index(P, 25)
    a45e = age_to_index(P, 45)
    a30s = age_to_index(P, 30)
    a55e = age_to_index(P, 55)
    a65s = age_to_index(P, 65)
    a75e = age_to_index(P, 75)
    asw = age_to_index(P, 45)
    aew = age_to_index(P, 55)
    newparent_cs = list(range(1, min(P.n_child_stages + 1, 3)))

    ti = sum(P.income[i, 0] * stats.worker_mass_by_loc[i] for i in range(I))
    tmw = float(np.sum(stats.worker_mass_by_loc))
    mean_income = ti / max(tmw, 1e-12)
    tw = tm4 = 0.0
    for jj in range(asw, aew + 1):
        for i in range(I):
            for ten in range(nt):
                heq = (1 - P.psi) * ph[i] * P.H_own[ten - 1] if ten > 0 else 0.0
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, i, jj, nn, cs]
                        tw += float(np.sum(gs * (bg + heq)))
                        tm4 += float(np.sum(gs))
    stats.mean_wealth_4555 = tw / max(tm4, 1e-12)
    stats.mean_income = mean_income
    stats.wealth_to_income = stats.mean_wealth_4555 / max(mean_income, 1e-12)
    stats.entry_mass_by_loc = np.array([np.sum(g[:, :, i, 0, :, :]) for i in range(I)])
    append_pension_budget_stats(stats, g, P)

    tmv = tas = 0.0
    for j in range(J - 1):
        for io in range(I):
            for ten in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, io, j, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        sp = lp[:, ten, io, io, j, nn, cs]
                        tmv += float(np.sum(gs * (1 - sp)))
                        tas += mh
    stats.migration_rate = tmv / max(tas, 1e-12)
    tmv2 = tas2 = 0.0
    for j in range(a22s, min(a45e, J - 1) + 1):
        for io in range(I):
            for ten in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, io, j, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        sp = lp[:, ten, io, io, j, nn, cs]
                        tmv2 += float(np.sum(gs * (1 - sp)))
                        tas2 += mh
    stats.migration_rate_2245 = tmv2 / max(tas2, 1e-12)

    tl = tml = 0.0
    for jj in range(asw, aew + 1):
        for i in range(I):
            for ten in range(nt):
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, i, jj, nn, cs]
                        tl += float(np.sum(gs * bg))
                        tml += float(np.sum(gs))
    stats.liquid_wealth_4555 = tl / max(tml, 1e-12)
    stats.liquid_wealth_to_income = stats.liquid_wealth_4555 / max(mean_income, 1e-12)

    aye = age_to_index(P, 35)
    young_block = g[:, :, :, a25s : aye + 1, :, :]
    yt = float(np.sum(young_block))
    yo = float(np.sum(young_block[:, 1:, :, :, :, :]))
    stats.young_own_rate = yo / max(yt, 1e-12)
    ylw = yinc = ycm = 0.0
    for jj in range(a25s, aye + 1):
        for i in range(I):
            gs = g[:, 0, i, jj, 0, 0]
            mh = float(np.sum(gs))
            if mh < 1e-15:
                continue
            ylw += float(np.sum(gs * bg))
            yinc += P.income[i, jj] * mh
            ycm += mh
    stats.young_liquid_wealth = ylw / max(ycm, 1e-12)
    stats.young_childless_renter_income = yinc / max(ycm, 1e-12)
    stats.young_liquid_wealth_to_income = ylw / max(yinc, 1e-12)

    tba = tfb = 0.0
    for j in range(P.A_f_start - 1, P.A_f_end):
        ra = P.age_start + j * P.da
        for i in range(I):
            for ten in range(nt):
                gs = g[:, ten, i, j, 0, 0]
                nz = gs > 1e-15
                if not np.any(nz):
                    continue
                pb = 1 - fp[nz, ten, i, j, 0]
                wb = float(np.sum(gs[nz] * pb))
                tba += ra * wb
                tfb += wb
    stats.mean_age_first_birth = tba / max(tfb, 1e-12)
    mg1 = float(np.sum(g[:, :, :, P.A_f_end :, 1:, :]))
    stats.parity_progression_1to2 = float(np.sum(g[:, :, :, P.A_f_end :, 2:, :]) / max(mg1, 1e-12))

    ic = 1 if I > 1 else 0
    mcy = float(np.sum(g[:, :, ic, a25s : aye + 1, 0, :]))
    mct = float(np.sum(g[:, :, :, a25s : aye + 1, 0, :]))
    stats.center_share_childless_young = mcy / max(mct, 1e-12)
    stats.center_share_parents = float(np.sum(g[:, :, ic, :, 1:, :]) / max(np.sum(g[:, :, :, :, 1:, :]), 1e-12))
    mpo = float(np.sum(g[:, 1:, :, :, 1:, :]))
    mpa = float(np.sum(g[:, :, :, :, 1:, :]))
    stats.own_rate_parents = mpo / max(mpa, 1e-12)
    mco = float(np.sum(g[:, 1:, :, :, 0, :]))
    mca = float(np.sum(g[:, :, :, :, 0, :]))
    stats.own_rate_childless = mco / max(mca, 1e-12)
    stats.own_family_gap = stats.own_rate_parents - stats.own_rate_childless

    a34e = age_to_index(P, 34)
    a35s = age_to_index(P, 35)
    a44e = age_to_index(P, 44)
    prime_mass = float(np.sum(g[:, :, :, a30s : a55e + 1, :, :]))
    prime_owner = float(np.sum(g[:, 1:, :, a30s : a55e + 1, :, :]))
    stats.own_rate_3055 = prime_owner / max(prime_mass, 1e-12)
    early_mass = float(np.sum(g[:, :, :, a25s : a34e + 1, :, :]))
    early_owner = float(np.sum(g[:, 1:, :, a25s : a34e + 1, :, :]))
    stats.own_rate_2534 = early_owner / max(early_mass, 1e-12)
    mid_mass = float(np.sum(g[:, :, :, a35s : a44e + 1, :, :]))
    mid_owner = float(np.sum(g[:, 1:, :, a35s : a44e + 1, :, :]))
    stats.own_rate_3544 = mid_owner / max(mid_mass, 1e-12)
    if I > 1:
        prime_mass_p = float(np.sum(g[:, :, 0, a30s : a55e + 1, :, :]))
        prime_mass_c = float(np.sum(g[:, :, 1, a30s : a55e + 1, :, :]))
        prime_owner_p = float(np.sum(g[:, 1:, 0, a30s : a55e + 1, :, :]))
        prime_owner_c = float(np.sum(g[:, 1:, 1, a30s : a55e + 1, :, :]))
        stats.own_gradient_3055 = prime_owner_p / max(prime_mass_p, 1e-12) - prime_owner_c / max(prime_mass_c, 1e-12)
    else:
        stats.own_gradient_3055 = 0.0
    nonparent_mass_2245 = float(np.sum(g[:, :, :, a22s : a45e + 1, 0, 0]))
    stats.center_share_nonparents_2245 = (
        float(np.sum(g[:, :, 1, a22s : a45e + 1, 0, 0]) / max(nonparent_mass_2245, 1e-12)) if I > 1 else 1.0
    )
    nonparent_mass_3055 = float(np.sum(g[:, :, :, a30s : a55e + 1, 0, 0]))
    nonparent_owner_3055 = float(np.sum(g[:, 1:, :, a30s : a55e + 1, 0, 0]))
    stats.own_rate_nonparents_3055 = nonparent_owner_3055 / max(nonparent_mass_3055, 1e-12)
    if not newparent_cs:
        stats.center_share_newparents_2245 = 0.0
        stats.own_rate_newparents_3055 = 0.0
    else:
        newparent_mass_2245 = float(np.sum(g[:, :, :, a22s : a45e + 1, 1:, newparent_cs]))
        stats.center_share_newparents_2245 = (
            float(np.sum(g[:, :, 1, a22s : a45e + 1, 1:, newparent_cs]) / max(newparent_mass_2245, 1e-12))
            if I > 1
            else 1.0
        )
        newparent_mass_3055 = float(np.sum(g[:, :, :, a30s : a55e + 1, 1:, newparent_cs]))
        newparent_owner_3055 = float(np.sum(g[:, 1:, :, a30s : a55e + 1, 1:, newparent_cs]))
        stats.own_rate_newparents_3055 = newparent_owner_3055 / max(newparent_mass_3055, 1e-12)
    stats.own_gap_newparent_nonparent_3055 = stats.own_rate_newparents_3055 - stats.own_rate_nonparents_3055

    old_mass = float(np.sum(g[:, :, :, a65s : a75e + 1, :, :]))
    old_owner = float(np.sum(g[:, 1:, :, a65s : a75e + 1, :, :]))
    stats.old_age_own_rate_6575 = old_owner / max(old_mass, 1e-12)
    old_parent_mass = float(np.sum(g[:, :, :, a65s : a75e + 1, 1:, :]))
    old_parent_owner = float(np.sum(g[:, 1:, :, a65s : a75e + 1, 1:, :]))
    old_childless_mass = float(np.sum(g[:, :, :, a65s : a75e + 1, 0, 0]))
    old_childless_owner = float(np.sum(g[:, 1:, :, a65s : a75e + 1, 0, 0]))
    stats.old_age_own_rate_parents_6575 = old_parent_owner / max(old_parent_mass, 1e-12)
    stats.old_age_own_rate_childless_6575 = old_childless_owner / max(old_childless_mass, 1e-12)
    stats.old_age_parent_childless_gap_6575 = stats.old_age_own_rate_parents_6575 - stats.old_age_own_rate_childless_6575
    old_nonhousing_wealth = old_income = 0.0
    old_parent_nonhousing_wealth = old_parent_income = 0.0
    old_childless_nonhousing_wealth = old_childless_income = 0.0
    old_total_wealth = old_parent_total_wealth = old_childless_total_wealth = 0.0
    old_nonhousing_ratio_vals: list[np.ndarray] = []
    old_nonhousing_ratio_wts: list[np.ndarray] = []
    old_total_ratio_vals: list[np.ndarray] = []
    old_total_ratio_wts: list[np.ndarray] = []
    old_parent_nonhousing_ratio_vals: list[np.ndarray] = []
    old_parent_nonhousing_ratio_wts: list[np.ndarray] = []
    old_childless_nonhousing_ratio_vals: list[np.ndarray] = []
    old_childless_nonhousing_ratio_wts: list[np.ndarray] = []
    old_parent_total_ratio_vals: list[np.ndarray] = []
    old_parent_total_ratio_wts: list[np.ndarray] = []
    old_childless_total_ratio_vals: list[np.ndarray] = []
    old_childless_total_ratio_wts: list[np.ndarray] = []
    for jj in range(a65s, a75e + 1):
        for i in range(I):
            yj = float(P.income[i, jj])
            for ten in range(nt):
                home_equity = (1 - P.psi) * ph[i] * P.H_own[ten - 1] if ten > 0 else 0.0
                for nn in range(npar):
                    for cs in range(ncs):
                        gs = g[:, ten, i, jj, nn, cs]
                        mass = float(np.sum(gs))
                        if mass < 1e-15:
                            continue
                        fin = float(np.sum(gs * bg))
                        total = fin + mass * home_equity
                        income = yj * mass
                        old_nonhousing_wealth += fin
                        old_total_wealth += total
                        old_income += income
                        positive = gs > 0
                        if np.any(positive) and yj > 0:
                            wts = gs[positive]
                            nonhousing_ratio = bg[positive] / yj
                            total_ratio = (bg[positive] + home_equity) / yj
                            old_nonhousing_ratio_vals.append(nonhousing_ratio)
                            old_nonhousing_ratio_wts.append(wts)
                            old_total_ratio_vals.append(total_ratio)
                            old_total_ratio_wts.append(wts)
                        if nn > 0:
                            old_parent_nonhousing_wealth += fin
                            old_parent_total_wealth += total
                            old_parent_income += income
                            if np.any(positive) and yj > 0:
                                old_parent_nonhousing_ratio_vals.append(nonhousing_ratio)
                                old_parent_nonhousing_ratio_wts.append(wts)
                                old_parent_total_ratio_vals.append(total_ratio)
                                old_parent_total_ratio_wts.append(wts)
                        elif cs == 0:
                            old_childless_nonhousing_wealth += fin
                            old_childless_total_wealth += total
                            old_childless_income += income
                            if np.any(positive) and yj > 0:
                                old_childless_nonhousing_ratio_vals.append(nonhousing_ratio)
                                old_childless_nonhousing_ratio_wts.append(wts)
                                old_childless_total_ratio_vals.append(total_ratio)
                                old_childless_total_ratio_wts.append(wts)
    stats.old_nonhousing_wealth_to_income_6575 = old_nonhousing_wealth / max(old_income, 1e-12)
    stats.old_total_wealth_to_income_6575 = old_total_wealth / max(old_income, 1e-12)
    parent_nonhousing_ratio = old_parent_nonhousing_wealth / max(old_parent_income, 1e-12)
    childless_nonhousing_ratio = old_childless_nonhousing_wealth / max(old_childless_income, 1e-12)
    parent_total_ratio = old_parent_total_wealth / max(old_parent_income, 1e-12)
    childless_total_ratio = old_childless_total_wealth / max(old_childless_income, 1e-12)
    stats.old_parent_nonhousing_wealth_to_income_6575 = parent_nonhousing_ratio
    stats.old_childless_nonhousing_wealth_to_income_6575 = childless_nonhousing_ratio
    stats.old_parent_childless_nonhousing_wealth_to_income_gap_6575 = parent_nonhousing_ratio - childless_nonhousing_ratio
    stats.old_parent_total_wealth_to_income_6575 = parent_total_ratio
    stats.old_childless_total_wealth_to_income_6575 = childless_total_ratio
    stats.old_parent_childless_total_wealth_to_income_gap_6575 = parent_total_ratio - childless_total_ratio
    stats.old_nonhousing_wealth_to_income_median_6575 = weighted_median_from_cells(
        old_nonhousing_ratio_vals, old_nonhousing_ratio_wts
    )
    stats.old_total_wealth_to_income_median_6575 = weighted_median_from_cells(old_total_ratio_vals, old_total_ratio_wts)
    old_parent_nonhousing_median = weighted_median_from_cells(
        old_parent_nonhousing_ratio_vals, old_parent_nonhousing_ratio_wts
    )
    old_childless_nonhousing_median = weighted_median_from_cells(
        old_childless_nonhousing_ratio_vals, old_childless_nonhousing_ratio_wts
    )
    old_parent_total_median = weighted_median_from_cells(old_parent_total_ratio_vals, old_parent_total_ratio_wts)
    old_childless_total_median = weighted_median_from_cells(old_childless_total_ratio_vals, old_childless_total_ratio_wts)
    stats.old_parent_nonhousing_wealth_to_income_median_6575 = old_parent_nonhousing_median
    stats.old_childless_nonhousing_wealth_to_income_median_6575 = old_childless_nonhousing_median
    stats.old_parent_childless_nonhousing_wealth_to_income_median_gap_6575 = (
        old_parent_nonhousing_median - old_childless_nonhousing_median
    )
    stats.old_parent_total_wealth_to_income_median_6575 = old_parent_total_median
    stats.old_childless_total_wealth_to_income_median_6575 = old_childless_total_median
    stats.old_parent_childless_total_wealth_to_income_median_gap_6575 = (
        old_parent_total_median - old_childless_total_median
    )
    owner_mass_2545 = float(np.sum(g[:, 1:, :, a25s : a45e + 1, :, :]))
    stats.owner_neg_liquid_share_2545 = float(np.sum(g[bg < 0, 1:, :, a25s : a45e + 1, :, :]) / max(owner_mass_2545, 1e-12))
    owner_mass_2534 = float(np.sum(g[:, 1:, :, a25s : a34e + 1, :, :]))
    stats.owner_neg_liquid_share_2534 = float(np.sum(g[bg < 0, 1:, :, a25s : a34e + 1, :, :]) / max(owner_mass_2534, 1e-12))

    dep_last = P.n_child_stages
    renter_vals: list[np.ndarray] = []
    renter_wts: list[np.ndarray] = []
    owner_vals: list[np.ndarray] = []
    owner_wts: list[np.ndarray] = []
    for j in range(a25s, a45e + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    if current_child_bin_dt(nn, cs, dep_last) != 2:
                        continue
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr):
                        renter_vals.append(hr[kr])
                        renter_wts.append(gr[kr])
                    for ten in range(1, nt):
                        go = g[:, ten, i, j, nn, cs]
                        ko = go > 0
                        if np.any(ko):
                            owner_vals.append(P.H_own[ten - 1] * np.ones(np.count_nonzero(ko)))
                            owner_wts.append(go[ko])
    stats.prime_childless_renter_median_rooms = weighted_median_from_cells(renter_vals, renter_wts)
    stats.prime_childless_owner_median_rooms = weighted_median_from_cells(owner_vals, owner_wts)
    renter_mass_3055_childless = 0.0
    renter_rooms_3055_childless = 0.0
    renter_rooms_ge6_3055_childless = 0.0
    owner_mass_3055_childless = 0.0
    owner_rooms_3055_childless = 0.0
    owner_rooms_ge6_3055_childless = 0.0
    owner_mass_3055_parent = 0.0
    owner_rooms_3055_parent = 0.0
    renter_mass_3055_parent = 0.0
    renter_rooms_3055_parent = 0.0
    for j in range(a30s, a55e + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    child_bin = current_child_bin_dt(nn, cs, dep_last)
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr):
                        wr = gr[kr]
                        rr = hr[kr]
                        mass = float(np.sum(wr))
                        if child_bin == 2:
                            renter_mass_3055_childless += mass
                            renter_rooms_3055_childless += float(np.sum(wr * rr))
                            renter_rooms_ge6_3055_childless += float(np.sum(wr[rr >= 6.0 - 1e-8]))
                        elif child_bin > 2:
                            renter_mass_3055_parent += mass
                            renter_rooms_3055_parent += float(np.sum(wr * rr))
                    for ten in range(1, nt):
                        go = g[:, ten, i, j, nn, cs]
                        mo = float(np.sum(go))
                        if mo <= 1e-15:
                            continue
                        rooms = float(P.H_own[ten - 1])
                        if child_bin == 2:
                            owner_mass_3055_childless += mo
                            owner_rooms_3055_childless += mo * rooms
                            if rooms >= 6.0 - 1e-8:
                                owner_rooms_ge6_3055_childless += mo
                        elif child_bin > 2:
                            owner_mass_3055_parent += mo
                            owner_rooms_3055_parent += mo * rooms
    stats.prime30_55_childless_renter_mean_rooms = (
        renter_rooms_3055_childless / max(renter_mass_3055_childless, 1e-12)
    )
    stats.prime30_55_childless_owner_mean_rooms = (
        owner_rooms_3055_childless / max(owner_mass_3055_childless, 1e-12)
    )
    stats.prime30_55_childless_owner_minus_renter_mean_rooms = (
        stats.prime30_55_childless_owner_mean_rooms - stats.prime30_55_childless_renter_mean_rooms
    )
    stats.prime30_55_childless_renter_share_rooms_ge6 = (
        renter_rooms_ge6_3055_childless / max(renter_mass_3055_childless, 1e-12)
    )
    stats.prime30_55_childless_owner_share_rooms_ge6 = (
        owner_rooms_ge6_3055_childless / max(owner_mass_3055_childless, 1e-12)
    )
    parent_owner_mean = owner_rooms_3055_parent / max(owner_mass_3055_parent, 1e-12)
    parent_renter_mean = renter_rooms_3055_parent / max(renter_mass_3055_parent, 1e-12)
    stats.prime30_55_parent_owner_minus_renter_mean_rooms = parent_owner_mean - parent_renter_mean
    owner_all_mass_2545 = 0.0
    owner_le6_mass_2545 = 0.0
    owner_7to8_mass_2545 = 0.0
    owner_ge9_mass_2545 = 0.0
    renter_cap_mass_2545_all = 0.0
    renter_cap_mass_2545_child0 = 0.0
    renter_cap_mass_2545_child1 = 0.0
    renter_mass_2545_all = 0.0
    renter_mass_2545_child0 = 0.0
    renter_mass_2545_child1 = 0.0
    renter_rooms_2545_child0 = 0.0
    renter_rooms_2545_child1 = 0.0
    owner_mass_2545_child0 = 0.0
    owner_mass_2545_child1 = 0.0
    owner_rooms_2545_child0 = 0.0
    owner_rooms_2545_child1 = 0.0
    owner_rung_mass_2545 = np.zeros(P.n_house)
    for j in range(a25s, a45e + 1):
        for i in range(I):
            for nn in range(npar):
                for cs in range(ncs):
                    child_bin = current_child_bin_dt(nn, cs, dep_last)
                    gr = g[:, 0, i, j, nn, cs]
                    hr = hR[:, 0, i, j, nn, cs]
                    kr = (gr > 0) & np.isfinite(hr) & (hr > 0)
                    if np.any(kr):
                        wr = gr[kr]
                        rr = hr[kr]
                        mass = float(np.sum(wr))
                        cap_mass = float(np.sum(wr[rr >= float(P.hR_max) - 1e-8]))
                        renter_mass_2545_all += mass
                        renter_cap_mass_2545_all += cap_mass
                        if child_bin == 2:
                            renter_mass_2545_child0 += mass
                            renter_cap_mass_2545_child0 += cap_mass
                            renter_rooms_2545_child0 += float(np.sum(wr * rr))
                        elif child_bin == 3:
                            renter_mass_2545_child1 += mass
                            renter_cap_mass_2545_child1 += cap_mass
                            renter_rooms_2545_child1 += float(np.sum(wr * rr))
                    for ten in range(1, nt):
                        go = g[:, ten, i, j, nn, cs]
                        mo = float(np.sum(go))
                        if mo <= 1e-15:
                            continue
                        rooms = float(P.H_own[ten - 1])
                        owner_all_mass_2545 += mo
                        owner_rung_mass_2545[ten - 1] += mo
                        if rooms <= 6.0 + 1e-8:
                            owner_le6_mass_2545 += mo
                        elif rooms <= 8.0 + 1e-8:
                            owner_7to8_mass_2545 += mo
                        else:
                            owner_ge9_mass_2545 += mo
                        if child_bin == 2:
                            owner_mass_2545_child0 += mo
                            owner_rooms_2545_child0 += mo * rooms
                        elif child_bin == 3:
                            owner_mass_2545_child1 += mo
                            owner_rooms_2545_child1 += mo * rooms
    for rung in range(P.n_house):
        setattr(
            stats,
            f"owner25_45_rung{rung + 1}_share",
            float(owner_rung_mass_2545[rung] / max(owner_all_mass_2545, 1e-12)),
        )
    stats.owner25_45_rooms_le6_share = float(owner_le6_mass_2545 / max(owner_all_mass_2545, 1e-12))
    stats.owner25_45_rooms_7to8_share = float(owner_7to8_mass_2545 / max(owner_all_mass_2545, 1e-12))
    stats.owner25_45_rooms_ge9_share = float(owner_ge9_mass_2545 / max(owner_all_mass_2545, 1e-12))
    stats.renter25_45_all_cap_share = float(renter_cap_mass_2545_all / max(renter_mass_2545_all, 1e-12))
    stats.renter25_45_current0_cap_share = float(renter_cap_mass_2545_child0 / max(renter_mass_2545_child0, 1e-12))
    stats.renter25_45_current1_cap_share = float(renter_cap_mass_2545_child1 / max(renter_mass_2545_child1, 1e-12))
    stats.renter25_45_current0_mean_rooms = float(renter_rooms_2545_child0 / max(renter_mass_2545_child0, 1e-12))
    stats.renter25_45_current1_mean_rooms = float(renter_rooms_2545_child1 / max(renter_mass_2545_child1, 1e-12))
    stats.owner25_45_current0_mean_rooms = float(owner_rooms_2545_child0 / max(owner_mass_2545_child0, 1e-12))
    stats.owner25_45_current1_mean_rooms = float(owner_rooms_2545_child1 / max(owner_mass_2545_child1, 1e-12))
    stats.mean_housing_by_parity = np.zeros(npar)
    for nn in range(npar):
        th = mn = 0.0
        for i in range(I):
            for ten in range(nt):
                for j in range(J):
                    for cs in range(ncs):
                        gs = g[:, ten, i, j, nn, cs]
                        mh = float(np.sum(gs))
                        if mh < 1e-15:
                            continue
                        if ten == 0:
                            th += float(np.sum(gs * hR[:, ten, i, j, nn, cs]))
                        else:
                            th += mh * P.H_own[ten - 1]
                        mn += mh
        stats.mean_housing_by_parity[nn] = th / max(mn, 1e-12)
    stats.housing_increment_1to2 = (
        stats.mean_housing_by_parity[2] - stats.mean_housing_by_parity[1] if npar >= 3 else 0.0
    )
    return stats


def append_pension_budget_stats(stats: SimpleNamespace, g: np.ndarray, P: SimpleNamespace) -> None:
    income_profile = P.income_age_profile
    period_scale = float(getattr(P, "period_years", getattr(P, "da", 1.0))) if bool(getattr(P, "scale_flows_to_period", False)) else 1.0
    payroll_tax_revenue = 0.0
    for i in range(P.I):
        for j in range(P.J_R):
            mass_ij = float(np.sum(g[:, :, i, j, :, :]))
            payroll_tax_revenue += period_scale * P.tau_pay * P.w_hat[i] * income_profile[j] * mass_ij
    pension_outlays = P.pension * stats.retiree_mass_total
    stats.payroll_tax_revenue = payroll_tax_revenue
    stats.pension_outlays = pension_outlays
    stats.pension_budget_residual = payroll_tax_revenue - pension_outlays
    stats.implied_balanced_pension = payroll_tax_revenue / max(stats.retiree_mass_total, 1e-12)
    stats.pension = P.pension


def advance_cohort_horizon(
    g_in,
    start_age,
    horizon,
    loc_probs,
    tenure_choice,
    tenure_probs,
    bp_pol,
    P,
    b_grid,
    SD,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    ust,
    Pia,
):
    g_out = g_in
    for step in range(1, horizon + 1):
        age_idx = start_age + step - 1
        if age_idx >= P.J - 1:
            break
        g_out = advance_cohort_one_period(
            g_out,
            age_idx,
            loc_probs,
            tenure_choice,
            tenure_probs,
            bp_pol,
            P,
            b_grid,
            SD,
            lmm_idx,
            lmm_wt,
            tmx_idx,
            tmx_wt,
            ust,
            Pia,
        )
    return g_out


def advance_cohort_one_period(
    gj,
    j,
    loc_probs,
    tenure_choice,
    tenure_probs,
    bp_pol,
    P,
    b_grid,
    SD,
    lmm_idx,
    lmm_wt,
    tmx_idx,
    tmx_wt,
    ust,
    Pia,
):
    Nb = len(b_grid)
    nt = 1 + P.n_house
    I = P.I
    npar = P.n_parity
    ncs = P.n_child_states
    nc = SD.nc
    use_compiled_scatter = NUMBA_AVAILABLE and bool(getattr(P, "use_numba_scatter", False))
    K = P.n_child_stages
    csm1 = K + 1
    csm2 = K + 2
    gpl = np.zeros((Nb, nt, I, npar, ncs))
    for io in range(I):
        for to in range(nt):
            go = flat_nc(gj[:, to, io, :, :], Nb, nc)
            if np.sum(go) < 1e-15:
                continue
            po = np.reshape(loc_probs[:, to, io, :, j, :, :], (Nb, I, nc), order="F")
            sp = po[:, io, :]
            gpl[:, to, io, :, :] += unflat_nc(go * sp, Nb, npar, ncs)
            idx = lmm_idx[io, to, :]
            wt = lmm_wt[io, to, :]
            for id_ in range(I):
                if id_ == io:
                    continue
                mp = go * po[:, id_, :]
                if use_compiled_scatter:
                    moved = scatter_cols_sameidx_kernel(idx, wt, mp, Nb)
                else:
                    moved = scatter_redistribute_cols_sameidx(idx, wt, mp, Nb)
                gpl[:, 0, id_, :, :] += unflat_nc(moved, Nb, npar, ncs)

    gpt = np.zeros((Nb, nt, I, npar, ncs))
    for nn in range(npar):
        for id_ in range(I):
            for to in range(nt):
                gs = gpl[:, to, id_, nn, :]
                if np.sum(gs) < 1e-15:
                    continue
                for tn in range(nt):
                    if tenure_probs is None:
                        tcs = tenure_choice[:, to, id_, j, nn, :]
                        mk = tcs == tn
                        if not np.any(mk):
                            continue
                        mt = gs * mk
                    else:
                        pr = tenure_probs[:, to, id_, j, nn, :, tn]
                        mt = gs * pr
                    if np.sum(mt) < 1e-15:
                        continue
                    rd = np.zeros((Nb, ncs))
                    for cs in range(ncs):
                        idx = tmx_idx[id_, to, tn, nn, cs, :]
                        wt = tmx_wt[id_, to, tn, nn, cs, :]
                        if use_compiled_scatter:
                            rd[:, cs] = scatter_vec_kernel(idx, wt, mt[:, cs], Nb)
                        else:
                            rd[:, cs] = scatter_redistribute(idx, wt, mt[:, cs], Nb)
                    gpt[:, tn, id_, nn, :] += rd

    gps = np.zeros((Nb, nt, I, npar, ncs))
    for i in range(I):
        for ten in range(nt):
            gf = flat_nc(gpt[:, ten, i, :, :], Nb, nc)
            bpv = flat_nc(bp_pol[:, ten, i, j, :, :], Nb, nc)
            idx, wt = interp_indices(b_grid, np.clip(bpv, b_grid[0], b_grid[-1]))
            if use_compiled_scatter:
                g_new = scatter_cols_kernel(idx, wt, gf, Nb)
            else:
                g_new = scatter_redistribute_cols(idx, wt, gf, Nb)
            gps[:, ten, i, :, :] = unflat_nc(g_new, Nb, npar, ncs)

    g_next = np.zeros((Nb, nt, I, npar, ncs))
    if ust:
        for nn in range(npar):
            Pi = Pia[:, :, nn]
            for i in range(I):
                for ten in range(nt):
                    gin = gps[:, ten, i, nn, :]
                    gout = gin @ Pi
                    g_next[:, ten, i, nn, :] = gout
    else:
        for nn in range(npar):
            for cs in range(ncs):
                if cs == 0:
                    csn = 0
                elif cs >= csm1:
                    csn = cs
                elif cs < K:
                    csn = cs + 1
                else:
                    csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                g_next[:, :, :, nn, csn] += gps[:, :, :, nn, cs]
    return g_next


def mean_housing_childless_weighted(weight_dist, j, hR_pol, P):
    Nb, nt, I = weight_dist.shape
    th = mn = 0.0
    for i in range(I):
        for ten in range(nt):
            gs = weight_dist[:, ten, i]
            mh = float(np.sum(gs))
            if mh < 1e-15:
                continue
            if ten == 0:
                th += float(np.sum(gs * hR_pol[:, ten, i, j, 0, 0]))
            else:
                th += mh * P.H_own[ten - 1]
            mn += mh
    return th / max(mn, 1e-12)


def mean_housing_distribution(g_dist, j, hR_pol, P):
    Nb, nt, I, npar, ncs = g_dist.shape
    th = mn = 0.0
    for nn in range(npar):
        for cs in range(ncs):
            for i in range(I):
                for ten in range(nt):
                    gs = g_dist[:, ten, i, nn, cs]
                    mh = float(np.sum(gs))
                    if mh < 1e-15:
                        continue
                    if ten == 0:
                        th += float(np.sum(gs * hR_pol[:, ten, i, j, nn, cs]))
                    else:
                        th += mh * P.H_own[ten - 1]
                    mn += mh
    return th / max(mn, 1e-12)


def apply_child_aging(Vn, P, Nb, nt, I, npar, ncs):
    Vc = np.zeros((Nb, nt, I, npar, ncs))
    K = P.n_child_stages
    if P.use_stochastic_aging and hasattr(P, "Pi_child"):
        Pa = P.Pi_child
        for nn in range(npar):
            Pi = Pa[:, :, nn]
            Vnn = np.reshape(Vn[:, :, :, nn, :], (-1, ncs), order="F")
            Vc[:, :, :, nn, :] = np.reshape(Vnn @ Pi.T, (Nb, nt, I, ncs), order="F")
    else:
        csm1 = K + 1
        csm2 = K + 2
        for nn in range(npar):
            for cs in range(ncs):
                if cs == 0:
                    csn = 0
                elif cs >= csm1:
                    csn = cs
                elif cs < K:
                    csn = cs + 1
                else:
                    csn = 0 if nn == 0 else csm1 if nn == 1 else csm2
                Vc[:, :, :, nn, cs] = Vn[:, :, :, nn, csn]
    return Vc


def pack_solution_markov_income(
    V,
    c,
    h,
    bp,
    tc,
    tp,
    lp,
    fp,
    fv,
    g,
    st: SimpleNamespace,
    w,
    p,
    P: SimpleNamespace,
) -> SimpleNamespace:
    z_grid, z_weights, Pi_z = income_transition_values(P)
    sol = SimpleNamespace(
        V=V,
        c_pol=c,
        hR_pol=h,
        bp_pol=bp,
        tenure_choice=tc,
        tenure_probs=tp,
        loc_probs=lp,
        fert_probs=fp,
        fert_value=fv,
        g=g,
        g_collapsed=np.sum(g, axis=4),
        w_hat=w,
        p_eq=p,
        type_values=z_grid.copy(),
        type_weights=z_weights.copy(),
        income_transition=Pi_z.copy(),
    )
    for key, value in vars(st).items():
        setattr(sol, key, value)
    norm = housing_demand_normalizer(P)
    owner_by_size = np.zeros(P.n_house)
    renter_by_market = np.zeros(P.I)
    owner_by_market = np.zeros(P.I)
    for i in range(P.I):
        renter_by_market[i] = float(np.sum(g[:, 0, i, :, :, :, :] * h[:, 0, i, :, :, :, :])) / max(norm, 1e-12)
        for ten in range(1, 1 + P.n_house):
            demand = float(np.sum(g[:, ten, i, :, :, :, :])) * float(P.H_own[ten - 1]) / max(norm, 1e-12)
            owner_by_size[ten - 1] += demand
            owner_by_market[i] += demand
    user_cost = P.user_cost_rate * np.asarray(p, dtype=float).reshape(-1)
    supply = P.H0 * (user_cost / P.r_bar) ** P.xi_supply
    sol.owner_user_cost = user_cost
    sol.owner_asset_price = np.asarray(p, dtype=float).reshape(-1)
    sol.rental_demand_by_market = renter_by_market
    sol.owner_demand_by_market = owner_by_market
    sol.owner_demand_by_size = owner_by_size
    sol.rental_demand_by_size = renter_by_market.copy()
    sol.housing_supply = supply
    sol.aggregate_rental_demand = float(np.sum(renter_by_market))
    sol.aggregate_owner_demand = float(np.sum(owner_by_market))
    sol.aggregate_housing_demand = float(sol.aggregate_rental_demand + sol.aggregate_owner_demand)
    sol.aggregate_housing_supply = float(np.sum(supply))
    sol.aggregate_housing_excess = float(sol.aggregate_housing_demand - sol.aggregate_housing_supply)
    sol.best_max_abs_rel_excess = float(
        np.max(np.abs((renter_by_market + owner_by_market - supply) / np.maximum(supply, 1e-12)))
    )
    sol.best_market_metric = sol.best_max_abs_rel_excess
    sol.converged = bool(sol.best_max_abs_rel_excess <= getattr(P, "tol_eq", 1e-4))
    sol.young_owner_rate = float(getattr(sol, "young_own_rate", np.nan))
    sol.old_owner_rate = float(getattr(sol, "old_age_own_rate_6575", np.nan))
    sol.mean_completed_fertility = float(getattr(sol, "mean_parity", np.nan))
    sol.childless_rate = float(getattr(sol, "parity_dist", np.array([np.nan]))[0])
    if not hasattr(sol, "own_rate_by_income_type"):
        sol.own_rate_by_income_type = np.full(len(z_grid), np.nan)
    if not hasattr(sol, "mean_fertility_by_income_type"):
        sol.mean_fertility_by_income_type = np.full(len(z_grid), np.nan)
    if not hasattr(sol, "housing_demand_by_income_type"):
        sol.housing_demand_by_income_type = np.full((len(z_grid), P.I), np.nan)
    return sol


def pack_fast_solution_markov_income(st: SimpleNamespace, p: np.ndarray, P: SimpleNamespace) -> SimpleNamespace:
    p_vec = np.asarray(p, dtype=float).reshape(-1)
    user_cost = P.user_cost_rate * p_vec
    supply = P.H0 * (user_cost / P.r_bar) ** P.xi_supply
    demand = np.asarray(st.housing_demand, dtype=float).reshape(-1)
    sol = SimpleNamespace()
    for key, value in vars(st).items():
        setattr(sol, key, value)
    sol.p_eq = p_vec
    sol.owner_user_cost = user_cost
    sol.owner_asset_price = p_vec
    sol.housing_supply = supply
    sol.aggregate_housing_demand = float(np.sum(demand))
    sol.aggregate_housing_supply = float(np.sum(supply))
    sol.aggregate_housing_excess = float(sol.aggregate_housing_demand - sol.aggregate_housing_supply)
    sol.best_max_abs_rel_excess = float(
        np.max(np.abs((demand - supply) / np.maximum(supply, 1e-12)))
    )
    sol.best_market_metric = sol.best_max_abs_rel_excess
    sol.converged = bool(sol.best_max_abs_rel_excess <= getattr(P, "tol_eq", 1e-4))
    sol.mean_completed_fertility = float(getattr(sol, "mean_parity", np.nan))
    sol.childless_rate = float(getattr(sol, "parity_dist", np.array([np.nan]))[0])
    return sol


def pack_solution(V, c, h, bp, tc, tp, lp, fp, fv, g, st: SimpleNamespace, w, p, P: SimpleNamespace) -> SimpleNamespace:
    sol = SimpleNamespace(
        V=V,
        c_pol=c,
        hR_pol=h,
        bp_pol=bp,
        tenure_choice=tc,
        tenure_probs=tp,
        loc_probs=lp,
        fert_probs=fp,
        fert_value=fv,
        g=g,
        w_hat=w,
        p_eq=p,
    )
    for key, value in vars(st).items():
        setattr(sol, key, value)
    norm = housing_demand_normalizer(P)
    owner_by_size = np.zeros(P.n_house)
    renter_by_market = np.zeros(P.I)
    owner_by_market = np.zeros(P.I)
    for i in range(P.I):
        renter_by_market[i] = float(np.sum(g[:, 0, i, :, :, :] * h[:, 0, i, :, :, :])) / max(norm, 1e-12)
        for ten in range(1, 1 + P.n_house):
            demand = float(np.sum(g[:, ten, i, :, :, :])) * float(P.H_own[ten - 1]) / max(norm, 1e-12)
            owner_by_size[ten - 1] += demand
            owner_by_market[i] += demand
    user_cost = P.user_cost_rate * np.asarray(p, dtype=float).reshape(-1)
    supply = P.H0 * (user_cost / P.r_bar) ** P.xi_supply
    sol.owner_user_cost = user_cost
    sol.owner_asset_price = np.asarray(p, dtype=float).reshape(-1)
    sol.rental_demand_by_market = renter_by_market
    sol.owner_demand_by_market = owner_by_market
    sol.owner_demand_by_size = owner_by_size
    sol.rental_demand_by_size = renter_by_market.copy()
    sol.housing_supply = supply
    sol.aggregate_rental_demand = float(np.sum(renter_by_market))
    sol.aggregate_owner_demand = float(np.sum(owner_by_market))
    sol.aggregate_housing_demand = float(sol.aggregate_rental_demand + sol.aggregate_owner_demand)
    sol.aggregate_housing_supply = float(np.sum(supply))
    sol.aggregate_housing_excess = float(sol.aggregate_housing_demand - sol.aggregate_housing_supply)
    sol.best_max_abs_rel_excess = float(
        np.max(np.abs((renter_by_market + owner_by_market - supply) / np.maximum(supply, 1e-12)))
    )
    sol.best_market_metric = sol.best_max_abs_rel_excess
    sol.converged = bool(sol.best_max_abs_rel_excess <= getattr(P, "tol_eq", 1e-4))
    sol.young_owner_rate = float(getattr(sol, "young_own_rate", np.nan))
    sol.old_owner_rate = float(getattr(sol, "old_age_own_rate_6575", np.nan))
    sol.mean_completed_fertility = float(getattr(sol, "mean_parity", np.nan))
    sol.childless_rate = float(getattr(sol, "parity_dist", np.array([np.nan]))[0])
    return sol


def get_completed_fertility(nn: int, cs: int, P: SimpleNamespace) -> int:
    K = P.n_child_stages
    if cs == 0:
        return 0
    if cs == K + 1:
        return 1
    if cs == K + 2:
        return 2
    return nn


def bequest_utility_vec(b, nk, P):
    b_gross = np.maximum(b, 0.0)
    estate_tax_rate = min(max(float(getattr(P, "estate_tax_rate", 0.0)), 0.0), 0.999)
    estate_tax_exemption = max(float(getattr(P, "estate_tax_exemption", 0.0)), 0.0)
    taxable = np.maximum(b_gross - estate_tax_exemption, 0.0)
    b = np.maximum(b_gross - estate_tax_rate * taxable, 0.0)
    scale = P.theta0 * max(1 + P.theta_n * nk, 0)
    if abs(P.sigma - 1) < 1e-6:
        return scale * np.log(P.theta1 + b)
    return scale * (P.theta1 + b) ** (1 - P.sigma) / (1 - P.sigma)


def pti_adjusted_downpayment(dp_arr, hcost, income, P, b_grid):
    """Optional underwriting screen based on actual transaction debt.

    The core collateral constraint is `dp_arr`: cash available before purchase
    must cover `(1 - phi) * pH`, and branch liquid wealth may not fall below
    `-phi * pH`. If PTI is enabled, use actual transaction debt
    `D=max(pH-cash, 0)`, not the maximum allowed LTV debt, so extra cash can
    relax the payment screen.
    """
    out = np.array(dp_arr, dtype=float, copy=True)
    pti_limit = max(float(getattr(P, "pti_limit", 1.0)), 0.0)
    q = max(float(getattr(P, "q", 0.0)), 0.0)
    tau_h = max(float(getattr(P, "tau_H", 0.0)), 0.0)
    incomes = np.asarray(income, dtype=float).reshape(-1)
    big_dp = max(float(b_grid[-1]) + 10.0 * np.max(np.maximum(hcost, 1.0)), 1e8)
    for i in range(out.shape[0]):
        y = float(incomes[i]) if i < incomes.size and np.isfinite(incomes[i]) else 0.0
        for ten in range(1, out.shape[1]):
            house_cost = float(hcost[i, ten])
            tax_payment = tau_h * house_cost
            allowed_debt_payment = pti_limit * y - tax_payment
            if allowed_debt_payment < 0.0:
                out[i, ten, :, :] = big_dp
            elif q > 1e-12:
                max_debt_pti = allowed_debt_payment / q
                min_cash_pti = house_cost - max_debt_pti
                out[i, ten, :, :] = np.maximum(out[i, ten, :, :], min_cash_pti)
    return out


def has_birth_dp_grant(P, nn, cs, to, tn):
    if not bool(getattr(P, "birth_dp_grant", False)):
        return False
    if to != 0 or tn <= 0:
        return False
    return (nn >= 1) and (cs == 1)


def get_phi_state_matrix(P):
    ps = np.tile(P.phi.reshape(-1, 1), (1, P.n_child_states))
    if not bool(getattr(P, "parent_dp_waiver", False)):
        return ps
    kp = get_parent_target_child_states(P)
    po = getattr(P, "parent_dp_waiver_phi", 1.0)
    if P.n_parity >= 2:
        ps[1:, kp] = np.maximum(ps[1:, kp], po)
    return ps


def get_phi_choice_tensor(P):
    nt = 1 + P.n_house
    pc = np.tile(P.phi.reshape(1, 1, P.n_parity, 1), (P.I, nt, 1, P.n_child_states))
    if not bool(getattr(P, "parent_dp_waiver", False)) or P.n_parity < 2:
        return pc
    po = getattr(P, "parent_dp_waiver_phi", 1.0)
    loc_idx = get_parent_target_locations(P)
    ten_idx = get_parent_target_owner_tenures(P)
    cs_idx = get_parent_target_child_states(P)
    if len(loc_idx) == 0 or len(ten_idx) == 0 or not np.any(cs_idx):
        return pc
    pc[np.ix_(loc_idx, ten_idx, np.arange(1, P.n_parity), np.where(cs_idx)[0])] = np.maximum(
        pc[np.ix_(loc_idx, ten_idx, np.arange(1, P.n_parity), np.where(cs_idx)[0])], po
    )
    return pc


def get_parent_target_locations(P):
    loc_idx = np.arange(P.I)
    vals = getattr(P, "parent_dp_waiver_locations", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        loc_idx = np.asarray(vals, dtype=int) - 1
    return loc_idx[(loc_idx >= 0) & (loc_idx < P.I)]


def get_parent_target_owner_tenures(P):
    ten_idx = np.arange(1, 1 + P.n_house)
    vals = getattr(P, "parent_dp_waiver_owner_rungs", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        ten_idx = np.asarray(vals, dtype=int)
    return ten_idx[(ten_idx >= 1) & (ten_idx <= P.n_house)]


def get_parent_target_child_states(P):
    kp = np.zeros(P.n_child_states, dtype=bool)
    if bool(getattr(P, "parent_dp_waiver_birth_state_only", False)):
        if P.n_child_states >= 2:
            kp[1] = True
    else:
        kp[1 : P.n_child_stages + 1] = True
    return kp


def get_birth_entry_grant_tensor(P):
    nt = 1 + P.n_house
    bg = np.zeros((P.I, nt, P.n_parity, P.n_child_states))
    if not bool(getattr(P, "birth_entry_grant", False)):
        return bg
    grant = getattr(P, "birth_entry_grant_amount", 0.0)
    if not (np.isscalar(grant) and np.isfinite(grant) and grant > 0) or P.n_parity < 2 or P.n_child_states < 2:
        return bg
    loc_idx = np.arange(P.I)
    vals = getattr(P, "birth_entry_grant_locations", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        loc_idx = np.asarray(vals, dtype=int) - 1
    loc_idx = loc_idx[(loc_idx >= 0) & (loc_idx < P.I)]
    ten_idx = np.arange(1, nt)
    vals = getattr(P, "birth_entry_grant_owner_rungs", None)
    if vals is not None and len(np.atleast_1d(vals)) > 0:
        ten_idx = np.asarray(vals, dtype=int)
    ten_idx = ten_idx[(ten_idx >= 1) & (ten_idx < nt)]
    if len(loc_idx) == 0 or len(ten_idx) == 0:
        return bg
    bg[np.ix_(loc_idx, ten_idx, np.arange(1, P.n_parity), np.array([1]))] = grant
    return bg


def current_child_bin_dt(nn, cs, dep_last):
    if cs == 0 or cs > dep_last:
        return 2
    current_n = max(nn, 0)
    if current_n <= 0:
        return 2
    if current_n == 1:
        return 3
    return 4
