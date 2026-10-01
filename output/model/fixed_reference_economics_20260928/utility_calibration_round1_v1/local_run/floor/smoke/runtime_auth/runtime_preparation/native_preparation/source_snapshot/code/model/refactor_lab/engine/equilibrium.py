"""Equilibrium stage of the extracted solver (bodies byte-identical; see split_receipt.json)."""
from __future__ import annotations
import copy
import math
import time
from types import SimpleNamespace
from typing import Any
import numpy as np
from . import joint_nested
from .adult_entry import adjusted_births, potential_entry_households
from .warm_price import search_warm_price
from .child_preferences import apply_child_preferences
from .parameters import (bequest_utility_net_active, child_earnings_multiplier, child_earnings_penalty_active, children_at_home_count, estate_housing_value, estate_receiver_active, estate_transfer_at_age, get_fecundity_by_age, independent_child_maturation_active, mortgage_stay_floor_active, parent_age_maturation_active, readiness_childless_states, readiness_cumulative_probability, readiness_gate_active, readiness_settled_state, readiness_transition_hazard, rental_wedge_active, unsecured_debt_floor)
from .kernels import (NUMBA_AVAILABLE, full_owner_block_kernel, full_renter_block_kernel, location_logit_kernel, scatter_cols_kernel, scatter_cols_sameidx_kernel, scatter_vec_kernel, tenure_choice_kernel, tenure_logit_kernel)
from .utils import (decode_flat_family_state, flat_nc, interp_indices, interp_on_grid, interp_vector, logsumexp, make_grid, make_value_interp, scatter_redistribute, scatter_redistribute_cols, scatter_redistribute_cols_sameidx, unflat_nc, weighted_median_from_cells, weighted_quantile)
from .shared import (ACCOUNTING_SCALE_PRICE_CLOSURES, BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES, RENEWAL_VALVE_CALIBRATED_CLOSURES, RENEWAL_VALVE_CLOSURES, housing_demand_normalizer, income_transition_values, precompute_shared)
from .household import (solve_bellman_full_markov_income)
from .distribution import (forward_distribution_markov_income, pack_fast_solution_markov_income, pack_solution_markov_income, upgrade_fast_markov_solution)


def uses_accounting_scale_price_closure(P: SimpleNamespace) -> bool:
    mode = str(getattr(P, "population_closure", "normalized")).lower()
    return mode in ACCOUNTING_SCALE_PRICE_CLOSURES or mode in BENCHMARK_NORMALIZED_OUTSIDE_CLOSURES


def uses_renewal_valve_closure(P: SimpleNamespace) -> bool:
    return str(getattr(P, "population_closure", "normalized")).lower() in RENEWAL_VALVE_CLOSURES


def uses_calibrated_renewal_valve_closure(P: SimpleNamespace) -> bool:
    mode = str(getattr(P, "population_closure", "normalized")).lower()
    return mode in RENEWAL_VALVE_CALIBRATED_CLOSURES or bool(getattr(P, "renewal_calibrate_outside_flow", False))


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


def markov_market_housing_demand(
    sol: SimpleNamespace,
    P: SimpleNamespace,
    b_grid: np.ndarray,
) -> tuple[np.ndarray, SimpleNamespace | None]:
    """Housing demand used to clear the active Markov-income equilibrium.

    The Markov-income branch solves a normalized stationary composition.  A
    renewal closure must therefore scale that composition's demand before the
    price is updated.  Previously this branch ignored the already-implemented
    renewal accounting and always cleared on normalized demand.

    The value returned in ``sol.housing_demand`` remains demand per normalized
    adult household.  The first return value is aggregate market demand after
    applying the stationary scale factor.
    """

    raw_demand = np.asarray(getattr(sol, "housing_demand", np.zeros(P.I)), dtype=float).reshape(-1)
    if uses_renewal_valve_closure(P):
        scale_info = renewal_population_scale(sol, P, b_grid)
        return np.asarray(scale_info.implied_housing_demand, dtype=float).reshape(-1), scale_info
    if uses_accounting_scale_price_closure(P):
        raise NotImplementedError(
            "The Markov-income equilibrium currently supports the transparent renewal-valve "
            "scale closure only; its value-based outside-option closures require entrant values "
            "to be carried through the fast price loop."
        )
    return raw_demand, None


def attach_markov_market_accounting(
    sol: SimpleNamespace,
    P: SimpleNamespace,
    b_grid: np.ndarray,
) -> SimpleNamespace:
    """Attach scale and actual market residuals to an accepted Markov solve."""

    market_demand, scale_info = markov_market_housing_demand(sol, P, b_grid)
    if scale_info is None:
        return sol
    supply = np.asarray(sol.housing_supply, dtype=float).reshape(-1)
    excess = market_demand - supply
    metric = float(np.max(np.abs(excess) / np.maximum(np.abs(supply), 1e-12)))
    sol.accounting_scale = scale_info
    sol.market_housing_demand = market_demand
    sol.market_aggregate_housing_demand = float(np.sum(market_demand))
    sol.market_housing_excess = excess
    sol.market_aggregate_housing_excess = float(np.sum(excess))
    sol.best_max_abs_rel_excess = metric
    sol.best_market_metric = metric
    sol.converged = bool(metric <= getattr(P, "tol_eq", 1e-4))
    return sol


def solve_markov_income_equilibrium(
    p_init: np.ndarray,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
    *,
    warm_price_state: dict[str, float | None] | None = None,
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
    price_cache: dict[tuple[float, ...], SimpleNamespace] = {}
    cache_hits = 0

    if verbose:
        z_text = ", ".join(f"{z:.2f}:{w:.2f}" for z, w in zip(z_grid, z_weights))
        print(f"  Markov income states: {z_text}")
        print(f"  Markov transition rows: {np.array2string(Pi_z, precision=3)}")

    SD_shared = precompute_shared(P, b_grid)
    equilibrium_method = str(getattr(P, "markov_equilibrium_method", "legacy_damped")).lower()
    use_direct_scalar = P.I == 1 and equilibrium_method in {"direct", "direct_brent", "brent"}
    warm_info: dict[str, Any] | None = None
    if warm_price_state is not None:
        if not use_direct_scalar:
            raise ValueError("Warm prices require the one-market direct equilibrium solver")
        if set(warm_price_state) - {"price", "slope"}:
            raise ValueError("Warm price state accepts only price and slope")
        warm_info = {"used": False, "initial_price": warm_price_state.get("price"),
                     "initial_slope": warm_price_state.get("slope"), "fallback_required": False}

    def fast_solve(price: np.ndarray) -> SimpleNamespace:
        nonlocal cache_hits, t_solve
        price_arr = np.asarray(price, dtype=float).reshape(-1)
        key = tuple(float(value) for value in price_arr)
        cached = price_cache.get(key)
        if cached is not None:
            cache_hits += 1
            return cached
        started = time.perf_counter()
        solution = solve_markov_income_at_prices(
            price_arr, P, b_grid, verbose=False, fast_stats=True, SD=SD_shared,
            retain_payload=True,
        )
        t_solve += time.perf_counter() - started
        summary = copy.copy(solution)
        if hasattr(summary, "_model_payload"):
            delattr(summary, "_model_payload")
        price_cache[key] = summary
        return solution

    use_warm_scalar = warm_price_state is not None and warm_price_state.get("price") is not None
    if use_warm_scalar:
        def warm_evaluate(price: float) -> tuple[float, float, SimpleNamespace]:
            P.eq_iter = len(price_cache) + 1
            sol = fast_solve(np.array([price]))
            demand, _ = markov_market_housing_demand(sol, P, b_grid)
            supply = float(np.asarray(sol.housing_supply).reshape(-1)[0])
            excess = float(demand[0]) - supply
            return excess, abs(excess) / max(abs(supply), 1e-12), sol

        price, best_sol, best_err, warm_info = search_warm_price(
            warm_evaluate,
            initial_price=float(warm_price_state["price"]),
            initial_slope=warm_price_state.get("slope"),
            lower_bound=float(getattr(P, "p_min", 1e-4)),
            upper_bound=float(getattr(P, "p_max", 100.0)),
            tolerance=float(P.tol_eq),
        )
        best_p, final_err = np.array([price]), best_err
        iterations_completed = int(warm_info["price_evaluations"])
    elif use_direct_scalar:
        P.eq_iter = 1
        sol_it = fast_solve(p)
        market_demand, _ = markov_market_housing_demand(sol_it, P, b_grid)
        demand = float(market_demand[0])
        supply = float(np.asarray(sol_it.housing_supply).reshape(-1)[0])
        best_err = abs(demand - supply) / max(abs(supply), 1e-12)
        best_p, best_sol, final_err, iterations_completed = p.copy(), sol_it, best_err, 1
    for it in range(1, 0 if use_direct_scalar else int(P.max_iter_eq) + 1):
        P.eq_iter = it
        t0 = time.perf_counter()
        sol_it = fast_solve(p)
        t_solve += 0.0
        Hd, _ = markov_market_housing_demand(sol_it, P, b_grid)
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
        best_sol = fast_solve(best_p)
    scalar_refine_info: dict[str, Any] = {"used": False}
    if (P.I == 1 and bool(getattr(P, "scalar_market_refine", True))
            and not (use_warm_scalar and best_err < P.tol_eq)):
        best_sol, best_p, best_err, scalar_refine_info = refine_one_market_markov_income(
            best_p,
            best_sol,
            best_err,
            P,
            b_grid,
            verbose=verbose,
            SD=SD_shared,
            price_cache=price_cache,
        )
    if hasattr(best_sol, "_model_payload"):
        best_sol = upgrade_fast_markov_solution(best_sol, P, b_grid, SD_shared)
    else:
        best_sol = solve_markov_income_at_prices(best_p, P, b_grid, verbose=False, fast_stats=False, SD=SD_shared)
    best_sol = attach_markov_market_accounting(best_sol, P, b_grid)
    if str(getattr(P, "adult_entry_clock", "child_departure")) == "split_birth_vintage":
        best_sol.adult_entry_stationary_residual = (
            float(best_sol.entry_rate) - float(best_sol.adult_entry_potential_total)
        )
        best_sol.adult_entry_stationary_relative_gap = abs(
            best_sol.adult_entry_stationary_residual
        ) / max(
            float(best_sol.entry_rate), float(best_sol.adult_entry_potential_total), 1e-12
        )
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
        "equilibrium_method": equilibrium_method,
        "unique_fast_price_evaluations": int(len(price_cache)),
        "price_cache_hits": int(cache_hits + scalar_refine_info.get("cache_hits", 0)),
        "scalar_market_refine": scalar_refine_info,
    }
    best_sol.converged = bool(best_err < P.tol_eq)
    if warm_price_state is not None:
        # Publish only certified starts. Retain lightweight cache summaries,
        # never whole household payloads, between normalization trials.
        if best_sol.converged:
            points = []
            for key, summary in price_cache.items():
                demand, _ = markov_market_housing_demand(summary, P, b_grid)
                excess = float(demand[0]) - float(np.asarray(summary.housing_supply).reshape(-1)[0])
                if math.isfinite(excess):
                    points.append((key[0], excess))
            slope = warm_price_state.get("slope")
            if len(points) >= 2 and points[-1][0] != points[-2][0]:
                slope = (points[-1][1] - points[-2][1]) / (points[-1][0] - points[-2][0])
                if not math.isfinite(slope):
                    slope = None
            warm_price_state.update(price=float(best_p[0]), slope=slope)
        best_sol.timings["warm_price_search"] = {
            **warm_info, "accepted_state": dict(warm_price_state) if best_sol.converged else None,
        }
    return best_sol, P, best_p


def refine_one_market_markov_income(
    best_p: np.ndarray,
    best_sol: SimpleNamespace,
    best_err: float,
    P: SimpleNamespace,
    b_grid: np.ndarray,
    verbose: bool = True,
    SD: SimpleNamespace | None = None,
    price_cache: dict[tuple[float, ...], SimpleNamespace] | None = None,
) -> tuple[SimpleNamespace, np.ndarray, float, dict[str, Any]]:
    p0 = float(np.asarray(best_p, dtype=float).reshape(-1)[0])
    p_min = float(getattr(P, "p_min", 1e-4))
    p_max = float(getattr(P, "p_max", 100.0))
    expand = max(float(getattr(P, "scalar_market_refine_expand", 1.5)), 1.05)
    max_expand = max(0, int(getattr(P, "scalar_market_refine_max_expand", 8)))
    direct_search = str(getattr(P, "markov_equilibrium_method", "legacy_damped")).lower() in {"direct", "direct_brent", "brent"}
    max_iter = max(1, int(getattr(P, "scalar_market_refine_iter", 24)))
    if direct_search:
        max_iter = max(max_iter, int(getattr(P, "scalar_market_direct_iter", 32)))
    tol = float(getattr(P, "tol_eq", 1e-4))
    eval_count = 0
    eval_time = 0.0
    cache_hits = 0
    cache = price_cache if price_cache is not None else {}

    def eval_price(price: float) -> tuple[float, float, SimpleNamespace]:
        nonlocal eval_count, eval_time, cache_hits
        key = (float(price),)
        sol = cache.get(key)
        if sol is None:
            t_eval = time.perf_counter()
            sol = solve_markov_income_at_prices(np.array([price]), P, b_grid, verbose=False, fast_stats=True, SD=SD, retain_payload=True)
            eval_time += time.perf_counter() - t_eval
            eval_count += 1
            summary = copy.copy(sol)
            if hasattr(summary, "_model_payload"):
                delattr(summary, "_model_payload")
            cache[key] = summary
        else:
            cache_hits += 1
        market_demand, _ = markov_market_housing_demand(sol, P, b_grid)
        demand = float(market_demand[0])
        supply = float(np.asarray(sol.housing_supply).reshape(-1)[0])
        excess = demand - supply
        metric = abs(excess) / max(abs(supply), 1e-12)
        return excess, metric, sol

    initial_market_demand, _ = markov_market_housing_demand(best_sol, P, b_grid)
    best_excess = float(initial_market_demand[0] - np.asarray(best_sol.housing_supply).reshape(-1)[0])
    best_metric = float(best_err)
    best_price = p0
    sol_best = best_sol

    expansions = 0
    directional_fallback_used = False
    if direct_search:
        lo = hi = p0
        ex_lo = ex_hi = best_excess
        metric_lo = metric_hi = best_metric
        sol_lo = sol_hi = best_sol
        direct_expand = max(float(getattr(P, "scalar_market_direct_expand", 1.10)), 1.01)
        search_lower = best_excess < 0.0
        while ex_lo * ex_hi > 0 and expansions < max_expand:
            expansions += 1
            candidate = max(p_min, lo / direct_expand) if search_lower else min(p_max, hi * direct_expand)
            excess, metric, solution = eval_price(candidate)
            if search_lower:
                lo, ex_lo = candidate, excess
            else:
                hi, ex_hi = candidate, excess
            if metric < best_metric:
                best_price, best_excess, best_metric, sol_best = candidate, excess, metric, solution
            if candidate <= p_min or candidate >= p_max:
                break
        if ex_lo * ex_hi > 0:
            directional_fallback_used = True
            lo, hi = max(p_min, p0 / expand), min(p_max, p0 * expand)
            ex_lo, metric_lo, sol_lo = eval_price(lo)
            ex_hi, metric_hi, sol_hi = eval_price(hi)
    else:
        lo, hi = max(p_min, p0 / expand), min(p_max, p0 * expand)
        ex_lo, metric_lo, sol_lo = eval_price(lo)
        ex_hi, metric_hi, sol_hi = eval_price(hi)
    for price, excess, metric, sol in ((lo, ex_lo, metric_lo, sol_lo), (hi, ex_hi, metric_hi, sol_hi)):
        if metric < best_metric:
            best_price, best_excess, best_metric, sol_best = price, excess, metric, sol
    while ex_lo * ex_hi > 0 and expansions < max_expand:
        expansions += 1
        lo, hi = max(p_min, lo / expand), min(p_max, hi * expand)
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
        "cache_hits": int(cache_hits),
        "directional_fallback_used": bool(directional_fallback_used),
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
            # Degenerate secant/IQI denominators are possible on a flat or
            # repeated excess-demand bracket.  Preserve the safeguarded
            # invariant: do not divide, force the bisection candidate below.
            degenerate_secant = fb == fa
            if fa != fc and fb != fc and not degenerate_secant:
                mid = (
                    a * fb * fc / ((fa - fb) * (fa - fc))
                    + b_pt * fa * fc / ((fb - fa) * (fb - fc))
                    + c_pt * fa * fb / ((fc - fa) * (fc - fb))
                )
            elif not degenerate_secant:
                mid = b_pt - fb * (b_pt - a) / (fb - fa)
            else:
                mid = 0.5 * (a + b_pt)
            lo_gate = min((3.0 * a + b_pt) / 4.0, b_pt)
            hi_gate = max((3.0 * a + b_pt) / 4.0, b_pt)
            use_bisect = (
                degenerate_secant
                or not (lo_gate < mid < hi_gate)
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
    info["cache_hits"] = int(cache_hits)
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
    retain_payload: bool = False,
) -> SimpleNamespace:
    if str(getattr(P, "adult_entry_clock", "child_departure")) == "split_birth_vintage":
        if int(P.I) != 1 or str(getattr(P, "population_closure", "normalized")) != "normalized":
            raise ValueError("Split birth-vintage stationary entry requires the closed one-market normalized closure")
        if float(P.period_years) != 4.0:
            raise ValueError("Split birth-vintage entry requires four-year model periods")
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
        bp_pol, hR_pol, tc, lp_j, fp, V, r, p, P, b_grid, SD, fast_stats=fast_stats, tenure_probs=tp,
        bp_pol_stay=getattr(P, "_bp_pol_stay", None),
    )
    t_dist = time.perf_counter() - t0
    if fast_stats:
        sol = pack_fast_solution_markov_income(stats, p, P)
    else:
        sol = pack_solution_markov_income(V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, g, stats, P.w_hat, p, P)
    if str(getattr(P, "adult_entry_clock", "child_departure")) == "split_birth_vintage":
        sol.adult_entry_stationary_residual = (
            float(sol.entry_rate) - float(sol.adult_entry_potential_total)
        )
        sol.adult_entry_stationary_relative_gap = abs(sol.adult_entry_stationary_residual) / max(
            float(sol.entry_rate), float(sol.adult_entry_potential_total), 1e-12
        )
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
    # Sequential births use this Bellman-to-KFE side channel.  It is carried
    # with any retained payload so an accepted cached price cannot reuse the
    # probabilities from a later rejected price evaluation.
    if fast_stats and retain_payload:
        sol._model_payload = (V, c_pol, hR_pol, bp_pol, tc, tp, lp_j, fp, fv, r, p, P._fert2_probs.copy())
        sol.joint_choice = getattr(P, "_joint_choice", None)
        sol._bp_pol_stay = getattr(P, "_bp_pol_stay", None)
        sol._c_pol_stay = getattr(P, "_c_pol_stay", None)
    if verbose:
        print(
            f"  Markov income fixed-price solve: own={100 * sol.own_rate:.1f}% "
            f"TFR={2 * sol.mean_parity:.2f}"
        )
    return sol
