from __future__ import annotations

import numpy as np
import pytest

from intergen_housing_fertility.solver import run_model_cp_dt as run_production
from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    get_fecundity_by_age,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import run_model_cp_dt as run_fork


def _tiny_markov_overrides() -> dict[str, object]:
    return {
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
        "solve_mode": "pe",
        "p_fixed": np.array([1.0]),
        "w_fixed": np.array([1.0]),
        "entry_shares_fixed": np.array([1.0]),
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
    }


def test_fecundity_schedule_nesting_and_terminal_age() -> None:
    default = apply_overrides(setup_parameters(), {"J": 10, "A_f_end": 10})
    default.fecundity_terminal_age = 22.0
    np.testing.assert_array_equal(get_fecundity_by_age(default), np.ones(default.J))

    active = apply_overrides(
        setup_parameters(),
        {
            "J": 10,
            "fecundity_omega1": 0.1,
            "fecundity_omega2": 0.1,
            "fecundity_terminal_age": 42.0,
        },
    )
    schedule = get_fecundity_by_age(active)
    ages = active.age_start + np.arange(active.J) * active.da
    assert np.all((0.0 <= schedule) & (schedule <= 1.0))
    assert np.all(np.diff(schedule[ages < 42.0]) <= 0.0)
    assert np.all(schedule[ages >= 42.0] == 0.0)


def test_default_fork_is_bitwise_nested_in_production() -> None:
    overrides = _tiny_markov_overrides()
    fork, _, _ = run_fork(overrides, verbose=False)
    production, _, _ = run_production(overrides, verbose=False)
    assert np.array_equal(fork.V, production.V)
    assert np.array_equal(fork.g, production.g)
    assert fork.total_births_kfe == production.total_births_kfe
    assert fork.mean_parity == production.mean_parity


def test_attempt_gate_realizes_half_of_attempt_hazard_and_conserves_mass() -> None:
    baseline, _, _ = run_fork(_tiny_markov_overrides(), verbose=False)
    overrides = {
        **_tiny_markov_overrides(),
        "fecundity_omega1": 0.5,
        "fecundity_omega2": 0.0,
    }
    gated, P, _ = run_fork(overrides, verbose=False)
    np.testing.assert_allclose(
        gated.first_birth_hazard_by_age,
        0.5 * gated.attempt_hazard_by_age,
        atol=0.0,
        rtol=0.0,
    )
    assert 0.0 < gated.total_births_kfe < baseline.total_births_kfe
    age_mass = np.sum(gated.g, axis=(0, 1, 2, 4, 5, 6))
    np.testing.assert_allclose(age_mass, np.full(P.J, 1.0 / P.J), atol=1e-10, rtol=0.0)


def test_active_gate_rejects_nonmarkov_path() -> None:
    overrides = {
        **_tiny_markov_overrides(),
        "use_income_types": False,
        "fecundity_omega1": 0.5,
    }
    with pytest.raises(NotImplementedError, match="fecundity/sequential births: markov-income path only"):
        run_fork(overrides, verbose=False)
