"""Solver for the Python port of `run_model_cp_dt.m`."""

from __future__ import annotations

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
    return str(getattr(P, "population_closure", "normalized")).lower() in ACCOUNTING_SCALE_PRICE_CLOSURES


def uses_renewal_valve_closure(P: SimpleNamespace) -> bool:
    return str(getattr(P, "population_closure", "normalized")).lower() in RENEWAL_VALVE_CLOSURES


def uses_calibrated_renewal_valve_closure(P: SimpleNamespace) -> bool:
    mode = str(getattr(P, "population_closure", "normalized")).lower()
    return mode in RENEWAL_VALVE_CALIBRATED_CLOSURES or bool(getattr(P, "renewal_calibrate_outside_flow", False))


def normalize_population_mass(P: SimpleNamespace) -> bool:
    if uses_outside_option_closure(P):
        return False
    return bool(getattr(P, "normalize_population_mass", True))


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
        print("  DT LIFECYCLE MODEL PYTHON PORT")
        print(f"  J={P.J}, beta={P.beta:.4f}, R={P.R_gross:.4f}, Nb={P.Nb}")
        print(f"  alpha={P.alpha_cons:.2f}, sigma={P.sigma:.2f}, chi={P.chi:.2f}, hR_max={P.hR_max:.1f}")
        print(f"  kappa_loc={P.kappa_loc:.2f}, kappa_fert={P.kappa_fert:.2f}")
        print(f"  location_choice_form={P.location_choice_form}")
        print("=" * 60)

    b_grid = make_grid(P)
    p_init = P.r_bar / P.user_cost_rate
    if hasattr(P, "p_init_override") and P.p_init_override is not None:
        p_init = np.asarray(P.p_init_override, dtype=float).reshape(-1)

    solve_mode = str(getattr(P, "solve_mode", "ge")).lower()
    do_pe = solve_mode in ("pe", "partial", "partial_equilibrium", "fixed")
    if do_pe:
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
        print("\n" + "=" * 20 + " RESULTS " + "=" * 20)
        print(f"TFR={2 * sol.mean_parity:.2f}, Own={100 * sol.own_rate:.1f}%, Prices=[{p_eq[0]:.3f},{p_eq[1]:.3f}]")
        print(
            f"Pop=[{sol.pop_share[0]:.2f},{sol.pop_share[1]:.2f}], "
            f"Gradient={sol.mean_parity_by_loc[0] - sol.mean_parity_by_loc[1]:.3f}"
        )
        print(f"Total: {elapsed:.1f} sec")
    return sol, P, p_eq


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
    V, c_pol, hR_pol, bp_pol, tc, lp_j, fp, fv, timings = solve_bellman_full(r, p_eq, P, b_grid, SD)
    g, stats = forward_distribution(bp_pol, hR_pol, tc, lp_j, fp, r, p_eq, P, b_grid, SD, fast_stats=False)
    sol = pack_solution(V, c_pol, hR_pol, bp_pol, tc, lp_j, fp, fv, g, stats, P.w_hat, p_eq)
    sol.timings = {"bellman_full": timings["bellman"], "distribution": timings.get("distribution", 0.0)}
    return sol, P, p_eq


def solve_equilibrium(
    p_init: np.ndarray, P: SimpleNamespace, b_grid: np.ndarray, verbose: bool = True
) -> tuple[SimpleNamespace, SimpleNamespace, np.ndarray]:
    p = np.asarray(p_init, dtype=float).reshape(-1)
    entry_closure = uses_outside_option_closure(P)
    renewal_closure = uses_renewal_valve_closure(P)
    scale_price_closure = uses_accounting_scale_price_closure(P) or renewal_closure
    if scale_price_closure and not (
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
            V, c_pol, hR_pol, bp_pol, tc, lp_j, fp, fv, _ = solve_bellman_full(r, p, P, b_grid, SD)
            t_full += time.perf_counter() - t0
            n_full += 1
            stored_bp = bp_pol.copy()
            stall_count = 0
            mode = "F"
        else:
            t0 = time.perf_counter()
            V, c_pol, hR_pol, bp_pol, tc, lp_j, fp, fv, _ = solve_bellman_eval(stored_bp, r, p, P, b_grid, SD)
            t_eval += time.perf_counter() - t0
            n_eval += 1
            mode = "E"

        t0 = time.perf_counter()
        g, stats = forward_distribution(bp_pol, hR_pol, tc, lp_j, fp, r, p, P, b_grid, SD, fast_stats=True)
        t_dist += time.perf_counter() - t0
        n_dist += 1

        scale_info = None
        if scale_price_closure:
            scale_state = SimpleNamespace(
                V=V,
                entry_rate=stats.entry_rate,
                entrants_mature_by_loc=stats.entrants_mature_by_loc,
                entrants_mature_total=stats.entrants_mature_total,
                total_mass=stats.total_mass,
                housing_demand=stats.housing_demand,
            )
            if renewal_closure:
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
                "lp_j": lp_j,
                "fp": fp,
                "fv": fv,
                "g": g,
                "r": r.copy(),
                "p": p.copy(),
            }

        if verbose:
            print(
                f"  {mode}{it:3d}: ep={err_p:.4f} ee={err_e:.4f} "
                f"own={100 * stats.own_rate:.1f}% TFR={2 * stats.mean_parity:.2f} "
                f"pop=[{stats.pop_share[0]:.2f},{stats.pop_share[1]:.2f}]"
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
        S["bp_pol"], S["hR_pol"], S["tc"], S["lp_j"], S["fp"], S["r"], S["p"], P, b_grid, SD, fast_stats=False
    )
    sol = pack_solution(S["V"], S["c_pol"], S["hR_pol"], S["bp_pol"], S["tc"], S["lp_j"], S["fp"], S["fv"], S["g"], full_stats, P.w_hat, S["p"])
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
        sol.accounting_scale = renewal_population_scale(sol, P, b_grid) if renewal_closure else accounting_population_scale(sol, P, b_grid)
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
            ocst[i, ten] = (P.delta + P.tau_H) * p_hat[i] * hs
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
                            oc, hsv, P.c_min,
                            alpha, oms, beta, gs_alpha1, gs_alpha2, gs_tol,
                        )
                    else:
                        bp_prev_o = flat_nc(bd[:, ten, i, :, :], Nb, nc) if j < J - 1 else None
                        for c in range(nc):
                            Vbar = Vco[:, c]
                            cb_c = SD.cb_flat[0, c]
                            pc = SD.psi_flat[0, c]
                            ht_c = max(hsv - SD.hb_flat[0, c], 1e-10)
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
                            oc, hsv, P.c_min, alpha, oms, beta,
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
                        ht_o = np.maximum(hsv - SD.hb_flat, 1e-10)
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

        if NUMBA_AVAILABLE and bool(getattr(P, "use_tenure_kernel", True)):
            VH, tcj = tenure_choice_kernel(
                Vd, b_grid, heq, hcost, dp_arr, bmo, SD.birth_dp, birth_entry_grant
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
                                    dpn = dp_arr[id_, tn, nn, cs]
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
                                    dpn = dp_arr[id_, tn, nn, cs]
                                    bmn = bmo[id_, tn, nn, cs]
                                    dpc = dpn - sp
                                    inf_m = (b_grid < dpc) | (bar < bmn)
                                    Vrs[inf_m, nn, cs] = -1e10
                            Vopt[:, :, :, tn] = Vrs
                    tc = np.argmax(Vopt, axis=3)
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

    return V, c_pol, hR_pol, bp_pol, tenure_choice, loc_probs, fert_probs, fert_value, {"bellman": time.perf_counter() - t0}


def golden_renter(lo, hi, Rv, Vbar, b_grid, dc, pc, cc, cb_c, hb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta, a1, a2, tol):
    d = hi - lo
    x1 = lo + a1 * d
    x2 = lo + a2 * d
    f1 = eval_renter(x1, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
    f2 = eval_renter(x2, Rv, Vbar, b_grid, dc, pc, cc, cb_c, ri, hRmax, ht_cap_c, Kr, alpha, oms, beta)
    d = a1 * a2 * d
    while np.any(d > tol):
        bt = f2 >= f1
        xe = np.where(bt, x2 + d, x1 - d)
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
        xe = np.where(bt, x2 + d, x1 - d)
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

    if fast_stats and NUMBA_AVAILABLE and bool(getattr(P, "use_compiled_forward_distribution", True)):
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
                    tcs = tenure_choice[:, to, id_, j, nn, :]
                    for tn in range(nt):
                        mk = tcs == tn
                        if not np.any(mk):
                            continue
                        mt = gs * mk
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

    a22s = max(0, round(22 - P.age_start))
    a25s = max(0, round(25 - P.age_start))
    a45e = min(J - 1, round(45 - P.age_start))
    a30s = max(0, round(30 - P.age_start))
    a55e = min(J - 1, round(55 - P.age_start))
    a65s = max(0, round(65 - P.age_start))
    a75e = min(J - 1, round(75 - P.age_start))
    asw = max(0, round(45 - P.age_start))
    aew = min(J - 1, round(55 - P.age_start))
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

    aye = min(J - 1, round(35 - P.age_start))
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

    ic = 1
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

    a34e = min(J - 1, round(34 - P.age_start))
    a35s = max(0, round(35 - P.age_start))
    a44e = min(J - 1, round(44 - P.age_start))
    prime_mass = float(np.sum(g[:, :, :, a30s : a55e + 1, :, :]))
    prime_owner = float(np.sum(g[:, 1:, :, a30s : a55e + 1, :, :]))
    stats.own_rate_3055 = prime_owner / max(prime_mass, 1e-12)
    early_mass = float(np.sum(g[:, :, :, a25s : a34e + 1, :, :]))
    early_owner = float(np.sum(g[:, 1:, :, a25s : a34e + 1, :, :]))
    stats.own_rate_2534 = early_owner / max(early_mass, 1e-12)
    mid_mass = float(np.sum(g[:, :, :, a35s : a44e + 1, :, :]))
    mid_owner = float(np.sum(g[:, 1:, :, a35s : a44e + 1, :, :]))
    stats.own_rate_3544 = mid_owner / max(mid_mass, 1e-12)
    prime_mass_p = float(np.sum(g[:, :, 0, a30s : a55e + 1, :, :]))
    prime_mass_c = float(np.sum(g[:, :, 1, a30s : a55e + 1, :, :]))
    prime_owner_p = float(np.sum(g[:, 1:, 0, a30s : a55e + 1, :, :]))
    prime_owner_c = float(np.sum(g[:, 1:, 1, a30s : a55e + 1, :, :]))
    stats.own_gradient_3055 = prime_owner_p / max(prime_mass_p, 1e-12) - prime_owner_c / max(prime_mass_c, 1e-12)
    nonparent_mass_2245 = float(np.sum(g[:, :, :, a22s : a45e + 1, 0, 0]))
    stats.center_share_nonparents_2245 = float(
        np.sum(g[:, :, 1, a22s : a45e + 1, 0, 0]) / max(nonparent_mass_2245, 1e-12)
    )
    nonparent_mass_3055 = float(np.sum(g[:, :, :, a30s : a55e + 1, 0, 0]))
    nonparent_owner_3055 = float(np.sum(g[:, 1:, :, a30s : a55e + 1, 0, 0]))
    stats.own_rate_nonparents_3055 = nonparent_owner_3055 / max(nonparent_mass_3055, 1e-12)
    if not newparent_cs:
        stats.center_share_newparents_2245 = 0.0
        stats.own_rate_newparents_3055 = 0.0
    else:
        newparent_mass_2245 = float(np.sum(g[:, :, :, a22s : a45e + 1, 1:, newparent_cs]))
        stats.center_share_newparents_2245 = float(
            np.sum(g[:, :, 1, a22s : a45e + 1, 1:, newparent_cs]) / max(newparent_mass_2245, 1e-12)
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
    payroll_tax_revenue = 0.0
    for i in range(P.I):
        for j in range(P.J_R):
            mass_ij = float(np.sum(g[:, :, i, j, :, :]))
            payroll_tax_revenue += P.tau_pay * P.w_hat[i] * income_profile[j] * mass_ij
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
            g_out, age_idx, loc_probs, tenure_choice, bp_pol, P, b_grid, SD, lmm_idx, lmm_wt, tmx_idx, tmx_wt, ust, Pia
        )
    return g_out


def advance_cohort_one_period(gj, j, loc_probs, tenure_choice, bp_pol, P, b_grid, SD, lmm_idx, lmm_wt, tmx_idx, tmx_wt, ust, Pia):
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
                tcs = tenure_choice[:, to, id_, j, nn, :]
                for tn in range(nt):
                    mk = tcs == tn
                    if not np.any(mk):
                        continue
                    mt = gs * mk
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


def pack_solution(V, c, h, bp, tc, lp, fp, fv, g, st: SimpleNamespace, w, p) -> SimpleNamespace:
    sol = SimpleNamespace(
        V=V,
        c_pol=c,
        hR_pol=h,
        bp_pol=bp,
        tenure_choice=tc,
        loc_probs=lp,
        fert_probs=fp,
        fert_value=fv,
        g=g,
        w_hat=w,
        p_eq=p,
    )
    for key, value in vars(st).items():
        setattr(sol, key, value)
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
    b = np.maximum(b, 0.0)
    scale = P.theta0 * max(1 + P.theta_n * nk, 0)
    if abs(P.sigma - 1) < 1e-6:
        return scale * np.log(P.theta1 + b)
    return scale * (P.theta1 + b) ** (1 - P.sigma) / (1 - P.sigma)


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
