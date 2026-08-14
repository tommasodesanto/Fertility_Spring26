"""Regression tests for renewal scaling in the active Markov-income solver."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from intergen_eqscale_seq_optimized.solver import (
    attach_markov_market_accounting,
    markov_market_housing_demand,
    run_model_cp_dt,
)


def _renewal_parameters() -> SimpleNamespace:
    return SimpleNamespace(
        I=1,
        N_target=1.0,
        population_closure="renewal_valve",
        normalize_population_mass=True,
        renewal_retention=1.0,
        renewal_calibrate_outside_flow=False,
        renewal_target_total_population=1.0,
        outside_entry_flow=0.2,
        outside_entry_shares=np.array([1.0]),
        entry_shares=np.array([1.0]),
        tol_eq=1e-8,
    )


def test_markov_market_demand_applies_renewal_scale() -> None:
    parameters = _renewal_parameters()
    baseline = SimpleNamespace(
        total_mass=1.0,
        entry_rate=1.0,
        entrants_mature_total=0.8,
        entrants_mature_by_loc=np.array([0.8]),
        housing_demand=np.array([4.0]),
    )
    policy = SimpleNamespace(**vars(baseline))
    policy.entrants_mature_total = 0.9
    policy.entrants_mature_by_loc = np.array([0.9])

    baseline_demand, baseline_scale = markov_market_housing_demand(
        baseline, parameters, np.array([0.0, 1.0])
    )
    policy_demand, policy_scale = markov_market_housing_demand(
        policy, parameters, np.array([0.0, 1.0])
    )

    np.testing.assert_allclose(baseline_scale.scale_factor, 1.0)
    np.testing.assert_allclose(baseline_demand, np.array([4.0]))
    np.testing.assert_allclose(policy_scale.scale_factor, 2.0)
    np.testing.assert_allclose(policy_demand, np.array([8.0]))


def test_attached_market_residual_uses_scaled_not_per_unit_demand() -> None:
    parameters = _renewal_parameters()
    solution = SimpleNamespace(
        total_mass=1.0,
        entry_rate=1.0,
        entrants_mature_total=0.9,
        entrants_mature_by_loc=np.array([0.9]),
        housing_demand=np.array([4.0]),
        housing_supply=np.array([8.0]),
    )

    attached = attach_markov_market_accounting(
        solution, parameters, np.array([0.0, 1.0])
    )

    np.testing.assert_allclose(attached.housing_demand, np.array([4.0]))
    np.testing.assert_allclose(attached.market_housing_demand, np.array([8.0]))
    assert attached.best_max_abs_rel_excess < 1e-12
    assert attached.converged


def test_tiny_markov_renewal_baseline_nests_normalized_equilibrium() -> None:
    common = {
        "J": 6,
        "J_R": 5,
        "A_f_start": 1,
        "A_f_end": 4,
        "Nb": 20,
        "b_min": 0.0,
        "b_core_lo": 0.0,
        "b_core_hi": 3.0,
        "b_mid_hi": 6.0,
        "b_max": 10.0,
        "n_house": 2,
        "H_own": np.array([2.0, 4.0]),
        "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]),
        "solve_mode": "ge",
        "use_income_types": True,
        "income_type_transition": "persistent",
        "z_grid": np.array([0.8, 1.2]),
        "z_weights": np.array([0.5, 0.5]),
        "Pi_z": np.array([[0.9, 0.1], [0.1, 0.9]]),
        "entry_wealth_mode": "scalar",
        "b_entry_fixed": 0.5,
        "entry_wealth_spread_nodes": 1,
        "c_bar_0": 0.05,
        "c_bar_n": 0.02,
        "h_bar_0": 0.25,
        "h_bar_jump": 0.10,
        "h_bar_n": 0.05,
        "lambda_d": 0.0,
        "use_full_kernel": True,
        "use_tenure_kernel": True,
        "use_loc_kernel": True,
        "tenure_choice_kappa": 0.0,
        "max_iter_eq": 20,
        "tol_eq": 2e-4,
    }
    normalized, _, normalized_price = run_model_cp_dt(common, verbose=False)
    rho = min(0.5, 0.5 * float(normalized.entry_rate) / max(float(normalized.entrants_mature_total), 1e-12))
    outside_flow = float(normalized.entry_rate) - rho * float(normalized.entrants_mature_total)
    renewal, _, renewal_price = run_model_cp_dt(
        {
            **common,
            "population_closure": "renewal_valve",
            "renewal_retention": rho,
            "outside_entry_flow": outside_flow,
            "renewal_calibrate_outside_flow": False,
        },
        verbose=False,
    )

    assert hasattr(renewal, "accounting_scale")
    np.testing.assert_allclose(renewal.accounting_scale.scale_factor, 1.0, atol=5e-4, rtol=0.0)
    np.testing.assert_allclose(renewal_price, normalized_price, atol=5e-4, rtol=0.0)
    assert renewal.best_max_abs_rel_excess <= renewal.timings["best_eq_error"] + 1e-10
