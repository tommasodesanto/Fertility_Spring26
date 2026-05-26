"""SMM objective for the Python DT center-periphery port."""

from __future__ import annotations

import copy
from types import SimpleNamespace
from typing import Any

import numpy as np

from .parameters import asdict
from .solver import run_model_cp_dt
from .theta import PARAM_NAMES, apply_theta, theta_bounds


def smm_objective_dt(
    theta: np.ndarray,
    targets: dict[str, float],
    weights: dict[str, float],
    P_base: Any,
    inversion_targets: dict[str, float],
    verbose: bool = False,
) -> tuple[float, SimpleNamespace, SimpleNamespace | None, bool, SimpleNamespace | None, np.ndarray | None]:
    theta = np.asarray(theta, dtype=float).reshape(-1)
    try:
        lb, ub = theta_bounds(theta)
    except ValueError:
        return 1e6, SimpleNamespace(), None, False, None, None

    if np.any(theta < lb) or np.any(theta > ub) or np.any(~np.isfinite(theta)):
        return 1e6, SimpleNamespace(), None, False, None, None

    P = SimpleNamespace(**copy.deepcopy(asdict(P_base)))
    P = apply_theta(P, theta, PARAM_NAMES[: len(theta)])

    max_inv_iter = int(inversion_targets.get("max_iter", 3))
    inv_tol = float(inversion_targets.get("tol", 0.03))
    inv_speed = float(inversion_targets.get("speed", 0.7))
    inv_loss_weight = float(inversion_targets.get("loss_weight", 0.0))
    freeze_geometry = bool(inversion_targets.get("freeze", False))
    if freeze_geometry:
        E_C = float(inversion_targets["fixed_E_C"])
        r_bar_C = float(inversion_targets["fixed_r_bar_C"])
        max_inv_iter = 1
    else:
        E_C = float(P.E_loc[1])
        r_bar_C = float(P.r_bar[1])

    best_inv_err = np.inf
    best_err_pop = np.inf
    best_err_rent = np.inf
    best_sol = None
    best_P = None
    best_p_eq = None

    for inv_iter in range(1, max_inv_iter + 1):
        P_iter = SimpleNamespace(**copy.deepcopy(asdict(P)))
        P_iter.E_loc = np.array([P.E_loc[0], E_C], dtype=float)
        P_iter.r_bar = np.array([P.r_bar[0], r_bar_C], dtype=float)
        try:
            sol_i, P_i, p_eq_i = run_model_cp_dt(P_iter, verbose=verbose)
            solve_ok = True
        except Exception as exc:  # pragma: no cover - exercised by benchmark scripts
            if verbose:
                print(f"  [INV iter {inv_iter}] Solve failed: {exc}")
            solve_ok = False
        if not solve_ok:
            if inv_iter == 1:
                return 1e6, SimpleNamespace(), None, False, None, None
            break

        model_pop_C = sol_i.pop_share[1]
        model_rent_ratio = (P_i.user_cost_rate * p_eq_i[1]) / (P_i.user_cost_rate * p_eq_i[0])
        err_pop = abs(model_pop_C - inversion_targets["pop_share_C"]) / max(abs(inversion_targets["pop_share_C"]), 1e-6)
        err_rent = abs(model_rent_ratio - inversion_targets["rent_ratio"]) / max(abs(inversion_targets["rent_ratio"]), 1e-6)
        inv_err = max(err_pop, err_rent)
        if verbose:
            print(
                f"  [INV {inv_iter}] pop_C={model_pop_C:.3f} ({inversion_targets['pop_share_C']:.3f}, {100 * err_pop:.1f}%), "
                f"rr={model_rent_ratio:.2f} ({inversion_targets['rent_ratio']:.2f}, {100 * err_rent:.1f}%)"
            )
        if inv_err < best_inv_err:
            best_inv_err = inv_err
            best_err_pop = err_pop
            best_err_rent = err_rent
            best_sol = sol_i
            best_P = P_i
            best_p_eq = p_eq_i
        if inv_err < inv_tol or freeze_geometry or inv_iter == max_inv_iter:
            break
        if model_pop_C > 1e-6 and model_pop_C < 1 - 1e-6:
            target_odds = inversion_targets["pop_share_C"] / max(1 - inversion_targets["pop_share_C"], 1e-6)
            model_odds = model_pop_C / max(1 - model_pop_C, 1e-6)
            E_C = E_C + P.kappa_loc * inv_speed * np.log(target_odds / max(model_odds, 1e-6))
        if model_rent_ratio > 1e-6:
            r_bar_C = r_bar_C * (inversion_targets["rent_ratio"] / model_rent_ratio) ** inv_speed
        E_C = float(np.clip(E_C, -5.0, 5.0))
        r_bar_C = float(np.clip(r_bar_C, 0.01, 0.30))

    if best_sol is None:
        return 1e6, SimpleNamespace(), None, False, None, None

    sol = best_sol
    P_out = best_P
    p_eq = best_p_eq
    converged = bool(best_inv_err < inv_tol)
    moments = extract_moments(sol, P_out, p_eq, best_inv_err, best_err_pop, best_err_rent, converged)

    loss = 0.0
    for fn, target_val in targets.items():
        if hasattr(moments, fn) and fn in weights:
            model_val = getattr(moments, fn)
            dev = (model_val - target_val) / abs(target_val) if abs(target_val) > 1e-6 else model_val - target_val
            loss += weights[fn] * dev**2
    if not np.isfinite(loss):
        return 1e6, moments, sol, False, P_out, p_eq
    inv_loss = inv_loss_weight * (best_err_pop**2 + best_err_rent**2)
    loss = loss + inv_loss if np.isfinite(inv_loss) else 1e6
    moments.inv_loss = inv_loss
    moments.inv_converged = float(best_inv_err < inv_tol)
    return float(loss), moments, sol, converged, P_out, p_eq


def extract_moments(
    sol: SimpleNamespace,
    P_out: SimpleNamespace,
    p_eq: np.ndarray,
    inv_err: float,
    err_pop: float,
    err_rent: float,
    inv_converged: bool,
) -> SimpleNamespace:
    moments = SimpleNamespace()
    moments.location_inv_err = inv_err
    moments.location_inv_err_pop = err_pop
    moments.location_inv_err_rent = err_rent
    moments.location_inv_converged = inv_converged
    moments.tfr = 2 * sol.mean_parity
    moments.childless_rate = sol.parity_dist[0]
    moments.mean_age_first_birth = sol.mean_age_first_birth
    moments.tfr_gradient = 2 * (sol.mean_parity_by_loc[0] - sol.mean_parity_by_loc[1])
    moments.parity_progression_1to2 = sol.parity_progression_1to2
    moments.own_rate = sol.own_rate_3055
    moments.own_rate_2534 = getattr(sol, "own_rate_2534", np.nan)
    moments.own_rate_3544 = getattr(sol, "own_rate_3544", np.nan)
    moments.own_gradient = sol.own_gradient_3055
    moments.own_family_gap = sol.own_gap_newparent_nonparent_3055
    moments.own_lifecycle_slope = sol.old_age_own_rate_6575 - moments.own_rate_2534
    moments.prime_childless_renter_median_rooms = sol.prime_childless_renter_median_rooms
    moments.prime_childless_owner_median_rooms = sol.prime_childless_owner_median_rooms
    moments.housing_increment_0to1 = getattr(
        sol,
        "housing_increment_0to1_eventstudy_t3",
        sol.mean_housing_by_parity[1] - sol.mean_housing_by_parity[0] if len(sol.mean_housing_by_parity) >= 2 else 0.0,
    )
    moments.housing_increment_0to1_onechild = getattr(sol, "housing_increment_0to1_onechild_eventstudy_t3", np.nan)
    moments.housing_increment_1to2 = getattr(sol, "housing_increment_1to2_proxy_t3", sol.housing_increment_1to2)
    moments.housing_increment_0to2plus = getattr(sol, "housing_increment_0to2plus_eventstudy_t3", np.nan)
    moments.young_liquid_wealth_to_income = sol.young_liquid_wealth_to_income
    moments.center_share_nonparents = sol.center_share_nonparents_2245
    moments.center_share_newparents = sol.center_share_newparents_2245
    moments.migration_rate = sol.migration_rate_2245
    moments.old_age_own_rate = sol.old_age_own_rate_6575
    moments.old_age_parent_childless_gap = sol.old_age_parent_childless_gap_6575
    moments.owner_neg_liquid_share_2545 = getattr(sol, "owner_neg_liquid_share_2545", np.nan)
    moments.inv_E_C = P_out.E_loc[1]
    moments.inv_r_bar_C = P_out.r_bar[1]
    moments.inv_pop_C = sol.pop_share[1]
    moments.inv_rent_ratio = (P_out.user_cost_rate * p_eq[1]) / (P_out.user_cost_rate * p_eq[0])
    moments.inv_err_pop = err_pop
    moments.inv_err_rent = err_rent
    moments.inv_max_error = inv_err
    return moments
