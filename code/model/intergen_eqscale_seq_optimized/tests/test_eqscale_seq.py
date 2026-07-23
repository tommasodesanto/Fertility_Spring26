from __future__ import annotations

import numpy as np

from intergen_housing_fertility.kernels import eval_renter_scalar as production_eval_renter_scalar
from intergen_eqscale_seq_optimized.kernels import eval_renter_scalar
from intergen_housing_fertility.solver import run_model_cp_dt as run_production
from intergen_eqscale_seq_optimized.local_panel import income_process_overrides
from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters
from intergen_eqscale_seq_optimized.solver import (
    InfeasibleThetaError,
    precompute_shared,
    run_model_cp_dt as run_fork,
)


def _tiny_markov() -> dict[str, object]:
    return {
        "J": 6, "J_R": 5, "A_f_start": 1, "A_f_end": 4,
        "Nb": 20, "b_min": 0.0, "b_core_lo": 0.0, "b_core_hi": 3.0,
        "b_mid_hi": 6.0, "b_max": 10.0, "n_house": 2,
        "H_own": np.array([2.0, 4.0]), "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]), "solve_mode": "pe",
        "p_fixed": np.array([1.0]), "w_fixed": np.array([1.0]),
        "entry_shares_fixed": np.array([1.0]), "use_income_types": True,
        "income_type_transition": "persistent", "z_grid": np.array([0.8, 1.2]),
        "z_weights": np.array([0.5, 0.5]), "Pi_z": np.array([[0.9, 0.1], [0.1, 0.9]]),
        "entry_wealth_mode": "scalar", "b_entry_fixed": 0.5,
        "entry_wealth_spread_nodes": 1, "c_bar_0": 0.05, "c_bar_n": 0.02,
        "h_bar_0": 0.25, "h_bar_jump": 0.10, "h_bar_n": 0.05,
        "lambda_d": 0.0, "use_full_kernel": True, "use_tenure_kernel": True,
        "use_loc_kernel": True, "tenure_choice_kappa": 0.0,
    }


def test_default_nesting_chain_is_bitwise_production() -> None:
    fork, _, _ = run_fork(_tiny_markov(), verbose=False)
    prod, _, _ = run_production(_tiny_markov(), verbose=False)
    assert np.array_equal(fork.V, prod.V)
    assert np.array_equal(fork.g, prod.g)


def test_sequential_accounting_pi_one() -> None:
    sol, P, _ = run_fork({**_tiny_markov(), "sequential_births": True}, verbose=False)
    np.testing.assert_allclose(P._second_births_by_age, P._second_attempts_by_age, atol=1e-12, rtol=0.0)
    age_mass = np.sum(sol.g, axis=(0, 1, 2, 4, 5, 6))
    np.testing.assert_allclose(age_mass, np.full(P.J, 1.0 / P.J), atol=1e-10, rtol=0.0)
    assert 0.0 <= sol.parity_progression_1to2_flow <= 1.0


def test_sequential_second_birth_hazard_respects_fecundity() -> None:
    sol, _, _ = run_fork(
        {**_tiny_markov(), "sequential_births": True, "fecundity_omega1": 0.5, "fecundity_omega2": 0.0},
        verbose=False,
    )
    np.testing.assert_allclose(
        sol.second_birth_hazard_by_age, 0.5 * sol.second_attempt_hazard_by_age, atol=1e-12, rtol=0.0
    )


def test_eqscale_renter_allocation_and_childless_invariance() -> None:
    base = {**_tiny_markov(), "preference_spec": "eqscale", "delta_alpha": 0.1, "gamma_e": 0.5}
    sol, P, _ = run_fork(base, verbose=False)
    # Choose an interior renter/parent node.  The F-order column is n + cs*npar.
    b, ten, i, j, z, n, cs = 10, 0, 0, 2, 0, 1, 1
    c = sol.c_pol[b, ten, i, j, z, n, cs]
    h = sol.hR_pol[b, ten, i, j, z, n, cs]
    al = np.clip(P.alpha_cons - (P.delta_alpha_jump + P.delta_alpha * n), 0.05, 0.95)
    assert c > P.c_min and h < P.hR_max
    np.testing.assert_allclose(h / c, ((1.0 - al) / al) / P.user_cost_rate, rtol=1e-7, atol=1e-7)
    no_shift, _, _ = run_fork({**base, "delta_alpha": 0.0}, verbose=False)
    # At the terminal period no continuation/fertility channel remains, so
    # the childless column isolates the shifter's direct incidence exactly.
    np.testing.assert_array_equal(sol.c_pol[:, :, :, -1, :, 0, :], no_shift.c_pol[:, :, :, -1, :, 0, :])
    np.testing.assert_array_equal(sol.hR_pol[:, :, :, -1, :, 0, :], no_shift.hR_pol[:, :, :, -1, :, 0, :])


def test_eqscale_form_default_is_bitwise_inert() -> None:
    base = {**_tiny_markov(), "preference_spec": "eqscale", "delta_alpha": 0.1, "gamma_e": 0.5}
    implicit, _, _ = run_fork(base, verbose=False)
    explicit, _, _ = run_fork({**base, "eqscale_form": "linear"}, verbose=False)
    assert np.array_equal(implicit.V, explicit.V)
    assert np.array_equal(implicit.g, explicit.g)


def test_eqscale_power_form_scale_values() -> None:
    power = apply_overrides(setup_parameters(), {"preference_spec": "eqscale", "eqscale_form": "power"})
    power_shared = precompute_shared(power, np.array([0.0, 1.0]))
    parity_one_parent = 1 + power.n_parity
    np.testing.assert_allclose(power_shared.escale_flat[0, parity_one_parent], 1.35**0.7, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(power_shared.escale_flat[0, 0], 1.0, atol=1e-12, rtol=0.0)

    sqrt = apply_overrides(setup_parameters(), {"preference_spec": "eqscale", "eqscale_form": "sqrt"})
    sqrt_shared = precompute_shared(sqrt, np.array([0.0, 1.0]))
    np.testing.assert_allclose(sqrt_shared.escale_flat[0, parity_one_parent], 1.5**0.5, atol=1e-12, rtol=0.0)
    np.testing.assert_allclose(sqrt_shared.escale_flat[0, 0], 1.0, atol=1e-12, rtol=0.0)


def test_eqscale_power_form_changes_parent_cells_only_at_terminal_age() -> None:
    base = {
        **_tiny_markov(),
        "preference_spec": "eqscale",
        "delta_alpha": 0.0,
        "delta_alpha_jump": 0.0,
        "gamma_e": 0.0,
    }
    linear, _, _ = run_fork({**base, "eqscale_form": "linear"}, verbose=False)
    power, _, _ = run_fork({**base, "eqscale_form": "power"}, verbose=False)
    assert not np.array_equal(linear.V, power.V)
    np.testing.assert_array_equal(linear.c_pol[:, :, :, -1, :, 0, :], power.c_pol[:, :, :, -1, :, 0, :])


def test_eqscale_sigma_point_twenty_is_feasible() -> None:
    try:
        income = income_process_overrides(3, process="rouwenhorst", annual_innovation_sd=0.20, annual_rho=0.90)
        run_fork({**_tiny_markov(), **income, "preference_spec": "eqscale"}, verbose=False)
    except InfeasibleThetaError as exc:  # pragma: no cover - assertion message is the contract
        raise AssertionError("eqscale sigma=0.20 must remain feasible") from exc


def test_sequential_preserves_childless_mass() -> None:
    sol, P, _ = run_fork(
        {**_tiny_markov(), "sequential_births": True, "kappa_fert": 8.0, "fecundity_omega1": 0.0},
        verbose=False,
    )
    for j in range(P.A_f_start - 1, P.A_f_end):
        assert float(np.sum(sol.g[:, :, :, j, :, 0, 0])) > 0.0
    age_mass = np.sum(sol.g, axis=(0, 1, 2, 4, 5, 6))
    np.testing.assert_allclose(age_mass, np.full(P.J, 1.0 / P.J), atol=1e-12, rtol=0.0)


def test_no_same_period_second_attempt() -> None:
    _, P, _ = run_fork({**_tiny_markov(), "sequential_births": True}, verbose=False)
    assert P._second_births_by_age[P.A_f_start - 1] == 0.0


def test_second_birth_hazard_bounded() -> None:
    sol, _, _ = run_fork(
        {**_tiny_markov(), "sequential_births": True, "fecundity_omega1": 0.5, "fecundity_omega2": 0.0},
        verbose=False,
    )
    assert np.all((0.0 <= sol.second_attempt_hazard_by_age) & (sol.second_attempt_hazard_by_age <= 1.0))
    assert np.all((0.0 <= sol.second_birth_hazard_by_age) & (sol.second_birth_hazard_by_age <= 1.0))
    np.testing.assert_allclose(
        sol.second_birth_hazard_by_age, 0.5 * sol.second_attempt_hazard_by_age, atol=1e-12, rtol=0.0
    )


def test_sequential_first_birth_event_study_is_measured() -> None:
    sol, _, _ = run_fork({**_tiny_markov(), "sequential_births": True}, verbose=False)
    assert sol.housing_increment_0to1_eventstudy_t3 > 0.0


def test_eval_renter_scalar_bitwise_vs_production_across_rents() -> None:
    production = production_eval_renter_scalar.py_func
    fork = eval_renter_scalar.py_func
    bg = np.array([0.0, 1.0, 2.0])
    Vbar = np.array([0.2, 0.3, 0.5])
    oms = -1.0
    for ri in np.geomspace(1e-3, 2.0, 200):
        for alpha, Rv, cc in ((0.3, 1.1, 2.0), (0.7, 4.0, 0.5)):
            Kr = (alpha ** alpha * ((1.0 - alpha) / ri) ** (1.0 - alpha)) ** oms
            args = (0.0, Rv, Vbar, bg, 0.1, 0.05, cc, 0.1, ri, 1.0, 0.9, Kr, alpha, oms, 0.95)
            assert fork(*args, 1.0) == production(*args)
