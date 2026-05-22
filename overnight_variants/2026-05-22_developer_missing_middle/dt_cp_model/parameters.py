"""Parameter and calibration setup for the Python DT center-periphery model."""

from __future__ import annotations

import copy
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np


def asdict(obj: Any) -> dict[str, Any]:
    if obj is None:
        return {}
    if isinstance(obj, Mapping):
        return dict(obj)
    if isinstance(obj, SimpleNamespace):
        return vars(obj).copy()
    return dict(vars(obj))


def clone_ns(obj: SimpleNamespace) -> SimpleNamespace:
    return copy.deepcopy(obj)


def setup_parameters() -> SimpleNamespace:
    P = SimpleNamespace()
    P.J = 60
    P.da = 1.0
    P.age_start = 18
    P.J_R = 47
    P.A_f_start = 1
    P.A_f_end = 25
    P.A_m = 18
    P.n_parity = 3
    P.use_stochastic_aging = True
    P.stage_durations = np.array([2.0, 4.0, 6.0, 6.0])
    P.n_child_stages = len(P.stage_durations)
    P.n_child_states = P.n_child_stages + 3
    P.Pi_child = make_child_transition_matrix_with_matured(P.stage_durations, P.n_parity)

    P.kappa_fert = 4.5
    P.kappa_loc = 2.0
    P.eps_fert = P.kappa_fert
    P.eps_loc = P.kappa_loc
    P.alpha_cons = 0.70
    P.c_bar_0 = 0.10
    P.c_bar_n = 0.12
    P.h_bar_0 = 4.0
    P.h_bar_jump = 0.75
    P.h_bar_n = 0.60
    P.child_housing_spec = "jump_plus_linear"
    P.kappa_h_base = 0.40
    P.kappa_h_slope = 0.0
    P.psi_child = 0.07
    P.theta0 = 0.53
    P.theta_n = 0.25
    P.theta1 = 0.01
    P.u_bar = 0.0
    P.b_entry_fixed = 0.0
    P.beta = 0.96
    P.rho = 1 / P.beta - 1
    P.gamma = 0.0
    P.sigma = 2.0
    P.q = 0.04
    P.delta = 0.02
    P.tau_H = 0.01
    P.R_gross = 1 + P.q
    P.chi = 1.10
    P.psi = 0.06
    P.phi = 0.80 * np.ones(P.n_parity)
    P.parent_dp_waiver = False
    P.parent_dp_waiver_phi = 1.0
    P.parent_dp_waiver_locations = np.array([], dtype=int)
    P.parent_dp_waiver_owner_rungs = np.array([], dtype=int)
    P.parent_dp_waiver_birth_state_only = False
    P.birth_dp_grant = False
    P.c_min = 0.01
    P.birth_entry_grant = False
    P.birth_entry_grant_amount = 0.0
    P.birth_entry_grant_locations = np.array([], dtype=int)
    P.birth_entry_grant_owner_rungs = np.array([], dtype=int)

    P.n_house = 6
    P.h_own_min = 4.0
    P.h_own_max = 11.0
    P.H_own = np.linspace(P.h_own_min, P.h_own_max, P.n_house)
    P.hR_max = 8.0
    P.I = 2
    P.w_hat = np.array([1.0, 1.0])
    P.E_loc = np.array([0.0, 0.30])
    P.income_age_breaks = np.array([18.0, 25.0, 35.0, 45.0, 55.0])
    P.income_age_values = np.array([0.565, 0.838, 1.000, 0.985, 0.935])
    P.mu_stay = 0.0
    P.mu_move = 5.0
    P.mu_move_parent = 5.0
    P.mu_age_decay = 0.0
    P.location_choice_form = "additive_due"
    P.N_0 = np.array([0.5, 0.5])
    P.entry_shares = P.N_0 / P.N_0.sum()
    P.N_target = 1.0
    P.E_total = 1 / P.J
    P.entry_by_loc = P.E_total * P.entry_shares
    P.population_closure = "normalized"
    P.normalize_population_mass = True
    P.outside_entry_flow = P.E_total
    P.outside_entry_shares = P.entry_shares.copy()
    P.local_birth_entry_weight = 1.0
    P.renewal_retention = 1.0
    P.renewal_calibrate_outside_flow = False
    P.renewal_target_total_population = P.N_target
    P.outside_value = 0.0
    P.outside_value_is_calibrated = False
    P.allow_uncalibrated_outside_value = False
    P.kappa_entry = P.kappa_loc
    P.entry_level_update = "log"
    P.entry_log_step_max = 1.0
    P.entry_mass_floor = 1e-12
    P.r_bar = np.array([0.04, 0.08])
    P.H0 = np.array([6.2, 5.3])
    P.eta_supply = np.array([2.0, 0.8])
    P.xi_supply = P.eta_supply.copy()
    P.p_min = 0.01
    P.p_max = 30.0
    P.alpha_price = 0.35
    P.rho_hat = P.rho
    P.user_cost_rate = P.q + P.delta + P.tau_H
    P.tau_pay = 0.179
    P.pension_mode = "balanced_stationary"
    P.pension = np.nan
    P.pension_by_loc = np.full(P.I, np.nan)
    P = set_income_given_w_and_pension(P)

    P.Nb = 80
    P.b_min = -35.0
    P.b_max = 100.0
    P.b_grid_power = 1.5
    P.n_sub = 1
    P.max_iter_eq = 200
    P.tol_eq = 1e-4
    P.lambda_eq = 0.30
    P.adaptive_price_damping = True
    P.lambda_price_min = 0.005
    P.lambda_price_max = 0.35
    P.lambda_entry_min = 0.005
    P.lambda_entry_max = 0.35
    P.adaptive_damping_decay = 0.70
    P.adaptive_damping_grow = 1.03
    P.cycle_guard_factor = 0.85
    P.target_filter_weight = 0.35
    P.enforce_price_bounds = True
    P.enforce_entry_share_floor = True
    P.entry_share_floor = 1e-4
    P.housing_demand_floor_for_supply = 1e-6
    P.report_clamp_hits = True
    P.kfe_wealth_interp = "linear"
    P.interp_method = "linear"
    P.has_lookup_mex = False
    return P


def build_calibration_setup(setup_mode: str = "benchmark") -> SimpleNamespace:
    setup_mode = (setup_mode or "benchmark").lower()

    targets = {
        "tfr": 1.70,
        "childless_rate": 0.15,
        "mean_age_first_birth": 26.0,
        "tfr_gradient": 0.133,
        "own_rate": 0.627,
        "own_gradient": 0.170,
        "own_family_gap": 0.110,
        "prime_childless_renter_median_rooms": 4.0,
        "prime_childless_owner_median_rooms": 6.0,
        "housing_increment_0to1": 0.66443467,
        "housing_increment_1to2": 0.56581378,
        "young_liquid_wealth_to_income": 0.60,
        "center_share_nonparents": 0.494,
        "center_share_newparents": 0.416,
        "migration_rate": 0.032,
        "old_age_own_rate": 0.863,
        "old_age_parent_childless_gap": 0.070,
    }
    weights = {
        "tfr": 12.0,
        "childless_rate": 12.0,
        "mean_age_first_birth": 12.0,
        "tfr_gradient": 10.0,
        "own_rate": 12.0,
        "own_gradient": 12.0,
        "own_family_gap": 10.0,
        "prime_childless_renter_median_rooms": 10.0,
        "prime_childless_owner_median_rooms": 10.0,
        "housing_increment_0to1": 8.0,
        "housing_increment_1to2": 2.0,
        "young_liquid_wealth_to_income": 12.0,
        "center_share_nonparents": 8.0,
        "center_share_newparents": 8.0,
        "migration_rate": 8.0,
        "old_age_own_rate": 10.0,
        "old_age_parent_childless_gap": 10.0,
    }
    inversion_targets = {
        "pop_share_C": 0.45,
        "rent_ratio": 1.14,
        "max_iter": 6,
        "tol": 0.01,
        "speed": 0.5,
        "loss_weight": 100.0,
    }

    P_base = SimpleNamespace()
    P_base.J = 60
    P_base.sigma = 2.0
    P_base.alpha_cons = 0.70
    P_base.c_bar_0 = 0.10
    P_base.h_bar_0 = 4.0
    P_base.h_bar_jump = 0.75
    P_base.q = 0.04
    P_base.delta = 0.02
    P_base.tau_H = 0.01
    P_base.psi = 0.06
    P_base.tau_pay = 0.179
    P_base.theta0 = 0.53
    P_base.theta_n = 0.25
    P_base.theta1 = 0.01
    P_base.w_hat = np.array([1.0, 1.0])
    P_base.income_age_breaks = np.array([18.0, 25.0, 35.0, 45.0, 55.0])
    P_base.income_age_values = np.array([0.565, 0.838, 1.000, 0.985, 0.935])
    P_base.eta_supply = np.array([2.0, 0.8])
    P_base.H0 = np.array([5.27, 4.505])
    P_base.E_loc = np.array([0.0, 0.30])
    P_base.I = len(P_base.E_loc)
    P_base.r_bar = np.array([0.04, 0.09])
    P_base.mu_stay = 0.0
    P_base.location_choice_form = "additive_due"
    P_base.n_parity = 4
    P_base.phi = 0.80 * np.ones(P_base.n_parity)
    P_base.n_house = 6
    P_base.H_own = np.linspace(4.0, 11.0, P_base.n_house)
    P_base.h_own_min = 4.0
    P_base.h_own_max = 11.0
    P_base.hR_max = 8.0
    P_base.Nb = 80
    P_base.b_min = -35.0
    P_base.b_max = 100.0
    P_base.b_grid_power = 1.5
    P_base.max_iter_eq = 200
    P_base.tol_eq = 5e-4
    P_base.setup_mode = setup_mode

    if setup_mode in ("benchmark", "full"):
        pass
    elif setup_mode in ("fast", "fast_explore", "explore"):
        P_base.n_house = 4
        P_base.H_own = np.array([4.0, 6.0, 8.5, 11.0])
        P_base.h_own_min = P_base.H_own[0]
        P_base.h_own_max = P_base.H_own[-1]
        P_base.Nb = 60
        P_base.max_iter_eq = 120
        P_base.tol_eq = 1e-3
        inversion_targets["max_iter"] = 4
        inversion_targets["tol"] = 0.02
    else:
        raise ValueError(f"Unknown setup_mode: {setup_mode}")

    lb = np.array([0.90, 0.00, 0.02, 0.20, 0.10, 0.05, 1.25, 1.01, 0.50, 0.00, 0.10, 0.00, 1.50])
    ub = np.array([0.99, 1.00, 0.12, 2.50, 1.00, 0.25, 8.00, 1.10, 4.00, 10.00, 1.50, 0.75, 4.50])
    x0 = np.array([0.955, 0.15, 0.060, 1.20, 0.70, 0.12, 3.00, 1.09, 1.50, 5.00, 0.53, 0.25, 4.00])
    names = [
        "beta",
        "b_entry_fixed",
        "psi_child",
        "h_bar_jump",
        "h_bar_n",
        "c_bar_n",
        "kappa_fert",
        "chi",
        "kappa_loc",
        "mu_move",
        "theta0",
        "theta_n",
        "h_bar_0",
    ]
    return SimpleNamespace(
        targets=targets,
        weights=weights,
        inversion_targets=inversion_targets,
        P_base=P_base,
        names=names,
        lb=lb,
        ub=ub,
        x0=x0,
    )


def apply_overrides(P: SimpleNamespace, overrides: Any | None) -> SimpleNamespace:
    explicit_mapping = isinstance(overrides, Mapping)
    od = asdict(overrides)
    for key, value in od.items():
        setattr(P, key, _coerce_value(value))

    if "eps_loc" in od:
        P.kappa_loc = P.eps_loc
        P.eps_loc = P.kappa_loc
    if "eps_fert" in od:
        P.kappa_fert = P.eps_fert
        P.eps_fert = P.kappa_fert
    if "kappa_loc" in od:
        P.eps_loc = P.kappa_loc
    if "kappa_fert" in od:
        P.eps_fert = P.kappa_fert
    if "H_bar" in od:
        P.H0 = np.asarray(P.H_bar, dtype=float)
    if "n_parity" in od:
        P.n_parity = int(P.n_parity)
        if not hasattr(P, "phi") or len(np.atleast_1d(P.phi)) != P.n_parity:
            P.phi = 0.80 * np.ones(P.n_parity)
        P.Pi_child = make_child_transition_matrix_with_matured(P.stage_durations, P.n_parity)
    if "H_own" in od:
        P.H_own = np.asarray(P.H_own, dtype=float)
        P.n_house = len(P.H_own)
        P.h_own_min = P.H_own[0]
        P.h_own_max = P.H_own[-1]
    if "J" in od:
        P.J = int(P.J)
        P.E_total = 1 / P.J
    if "I" in od:
        P.I = int(P.I)
    if "Nb" in od:
        P.Nb = int(P.Nb)
    if "phi" in od:
        P.phi = np.asarray(P.phi, dtype=float).reshape(-1)
        if P.phi.size == 1:
            P.phi = P.phi.item() * np.ones(P.n_parity)
    if "eta_supply" in od:
        P.eta_supply = np.asarray(P.eta_supply, dtype=float)
        P.xi_supply = P.eta_supply.copy()
    if "w_hat" in od:
        P.w_hat = np.asarray(P.w_hat, dtype=float)
    if "entry_shares" in od:
        e0 = np.asarray(P.entry_shares, dtype=float).reshape(-1)
        e0 = np.maximum(e0, 0)
        P.entry_shares = e0 / e0.sum() if e0.sum() > 0 else np.ones(P.I) / P.I
    if "N_0" in od:
        P.N_0 = np.asarray(P.N_0, dtype=float)
    if "E_loc" in od:
        P.E_loc = np.asarray(P.E_loc, dtype=float)
        P.I = len(P.E_loc)
    if "r_bar" in od:
        P.r_bar = np.asarray(P.r_bar, dtype=float)
    if "H0" in od:
        P.H0 = np.asarray(P.H0, dtype=float)

    P.user_cost_rate = P.q + P.delta + P.tau_H
    P.R_gross = 1 + P.q
    P.rho = 1 / P.beta - 1
    P.rho_hat = P.rho
    P = set_income_given_w_and_pension(P)
    if not hasattr(P, "population_closure"):
        P.population_closure = "normalized"
    if not hasattr(P, "normalize_population_mass"):
        P.normalize_population_mass = True
    if "population_closure" in od:
        mode = str(P.population_closure).lower()
        if mode in ("outside_option", "outside_option_local_births", "local_births_outside", "open_city"):
            P.normalize_population_mass = False
    if not hasattr(P, "outside_entry_flow") or ("J" in od and "outside_entry_flow" not in od):
        P.outside_entry_flow = P.E_total
    if not hasattr(P, "outside_entry_shares"):
        P.outside_entry_shares = P.entry_shares.copy()
    if "outside_entry_shares" in od:
        oes = np.maximum(np.asarray(P.outside_entry_shares, dtype=float).reshape(-1), 0.0)
        P.outside_entry_shares = oes / oes.sum() if oes.sum() > 0 else np.ones(P.I) / P.I
    if not hasattr(P, "local_birth_entry_weight"):
        P.local_birth_entry_weight = 1.0
    if not hasattr(P, "renewal_retention"):
        P.renewal_retention = 1.0
    if not hasattr(P, "renewal_calibrate_outside_flow"):
        P.renewal_calibrate_outside_flow = False
    if not hasattr(P, "renewal_target_total_population"):
        P.renewal_target_total_population = getattr(P, "N_target", 1.0)
    if not hasattr(P, "outside_value"):
        P.outside_value = 0.0
    if explicit_mapping and "outside_value" in od:
        P.outside_value_is_calibrated = True
    if not hasattr(P, "outside_value_is_calibrated"):
        P.outside_value_is_calibrated = False
    if not hasattr(P, "allow_uncalibrated_outside_value"):
        P.allow_uncalibrated_outside_value = False
    if not hasattr(P, "kappa_entry") or ("kappa_loc" in od and "kappa_entry" not in od):
        P.kappa_entry = P.kappa_loc
    if "entry_by_loc" in od:
        entry = np.maximum(np.asarray(P.entry_by_loc, dtype=float).reshape(-1), 0.0)
        if entry.size != P.I:
            raise ValueError("entry_by_loc must have length P.I")
        P.entry_by_loc = entry
        P.E_total = float(np.sum(entry))
        P.entry_shares = entry / P.E_total if P.E_total > 0 else np.ones(P.I) / P.I
    else:
        P.entry_by_loc = P.E_total * P.entry_shares
    return P


def finalize_location_choice_spec(P: SimpleNamespace) -> SimpleNamespace:
    form = getattr(P, "location_choice_form", "legacy_multiplicative")
    form = str(form).lower()
    if form in ("additive_due", "additive", "due"):
        P.location_choice_form = "additive_due"
        P.location_choice_legacy_converted = False
    elif form in ("legacy_multiplicative", "multiplicative"):
        if np.any(np.asarray(P.E_loc) <= 0):
            raise ValueError("Legacy multiplicative location amenities must be strictly positive.")
        if getattr(P, "mu_stay", 1.0) <= 0 or P.mu_move <= 0:
            raise ValueError("Legacy multiplicative moving wedges must be strictly positive.")
        P.E_loc = P.kappa_loc * np.log(P.E_loc)
        P.mu_stay = -P.kappa_loc * np.log(P.mu_stay)
        P.mu_move = -P.kappa_loc * np.log(P.mu_move)
        if hasattr(P, "mu_move_parent") and P.mu_move_parent > 0:
            P.mu_move_parent = -P.kappa_loc * np.log(P.mu_move_parent)
        P.location_choice_form = "additive_due"
        P.location_choice_legacy_converted = True
    else:
        raise ValueError(f"Unknown location_choice_form: {form}")

    if not hasattr(P, "mu_stay"):
        P.mu_stay = 0.0
    return P


def set_income_given_w_and_pension(P: SimpleNamespace) -> SimpleNamespace:
    income_profile = get_income_age_profile(P)
    P.income_age_profile = income_profile
    P.pension = resolve_pension_value(P, income_profile)
    P.pension_by_loc = P.pension * np.ones(P.I)
    P.income = np.zeros((P.I, P.J))
    for i in range(P.I):
        for j in range(P.J):
            if j < P.J_R:
                P.income[i, j] = (1 - P.tau_pay) * P.w_hat[i] * income_profile[j]
            else:
                P.income[i, j] = P.pension
    return P


def resolve_pension_value(P: SimpleNamespace, income_profile: np.ndarray) -> float:
    mode = str(getattr(P, "pension_mode", "balanced_stationary")).lower()
    if mode in ("balanced_stationary", "balanced", "paygo", "fixed"):
        avg_worker_income = compute_stationary_worker_income(P, income_profile)
        retiree_ratio = P.J_R / max(P.J - P.J_R, 1)
        return P.tau_pay * avg_worker_income * retiree_ratio
    if mode in ("legacy_fixed", "legacy"):
        return P.tau_pay * 1.0 * (P.J_R / max(P.J - P.J_R, 1))
    if mode in ("manual", "custom"):
        pension = getattr(P, "pension", np.nan)
        if not np.isfinite(pension):
            raise ValueError("Manual pension_mode requires finite P.pension.")
        return float(pension)
    raise ValueError(f"Unknown pension_mode: {mode}")


def compute_stationary_worker_income(P: SimpleNamespace, income_profile: np.ndarray) -> float:
    worker_profile = income_profile[: P.J_R]
    if worker_profile.size == 0:
        return 0.0
    loc_weights = np.asarray(getattr(P, "entry_shares", np.ones_like(P.w_hat)), dtype=float)
    if loc_weights.size != len(P.w_hat) or loc_weights.sum() <= 0:
        loc_weights = np.ones(len(P.w_hat))
    loc_weights = loc_weights / loc_weights.sum()
    avg_wage = np.sum(loc_weights * P.w_hat)
    return float(avg_wage * np.mean(worker_profile))


def get_income_age_profile(P: SimpleNamespace) -> np.ndarray:
    income_profile = np.ones(P.J)
    age_breaks = np.asarray(getattr(P, "income_age_breaks", [18, 25, 35, 45, 55]), dtype=float)
    age_values = np.asarray(getattr(P, "income_age_values", [0.565, 0.838, 1.0, 0.985, 0.935]), dtype=float)
    age_vec = P.age_start + np.arange(P.J) * P.da
    for k, val in enumerate(age_values):
        age_hi = age_breaks[k + 1] if k < len(age_values) - 1 else P.age_start + P.J_R * P.da
        mask = (age_vec >= age_breaks[k]) & (age_vec < age_hi)
        income_profile[mask] = val
    return income_profile


def make_child_transition_matrix_with_matured(stage_durations: np.ndarray, n_parity: int) -> np.ndarray:
    sd = np.asarray(stage_durations, dtype=float)
    K = len(sd)
    ncs = K + 3
    Pa = np.zeros((ncs, ncs, n_parity))
    for nn in range(n_parity):
        Pi = np.zeros((ncs, ncs))
        Pi[0, 0] = 1.0
        for k in range(K):
            cs = k + 1
            pa = min(max(1.0 / sd[k], 0.0), 1.0)
            Pi[cs, cs] = 1.0 - pa
            if k < K - 1:
                Pi[cs, cs + 1] = pa
            else:
                if nn == 0:
                    Pi[cs, 0] = pa
                elif nn == 1:
                    Pi[cs, K + 1] = pa
                else:
                    Pi[cs, K + 2] = pa
        Pi[K + 1, K + 1] = 1.0
        Pi[K + 2, K + 2] = 1.0
        Pa[:, :, nn] = Pi
    return Pa


def _coerce_value(value: Any) -> Any:
    if isinstance(value, list | tuple):
        return np.asarray(value, dtype=float)
    return copy.deepcopy(value)
