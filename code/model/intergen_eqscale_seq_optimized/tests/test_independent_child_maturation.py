"""Tests for independent stochastic maturation without child-age states."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    make_independent_child_count_transition_matrix,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import (
    current_child_bin_dt,
    get_birth_entry_grant_tensor,
    get_completed_fertility,
    get_parent_target_child_states,
    income_at_state,
    markov_grant_outlays,
    precompute_shared,
    property_tax_revenue_from_distribution,
    run_model_cp_dt,
)
from intergen_eqscale_seq_optimized.tests.test_literal_parity import L4, _seq_eqscale


def repaired_parameters(**extra: object) -> SimpleNamespace:
    P = setup_parameters()
    overrides: dict[str, object] = {
        "n_parity": 4,
        "child_state_mode": "independent_count",
        "A_m": 18.0,
        "use_stochastic_aging": True,
    }
    overrides.update(extra)
    return apply_overrides(P, overrides)


def test_two_child_transition_matches_independent_binomial_formula() -> None:
    matrix = make_independent_child_count_transition_matrix(18.0 / 4.0, 4)
    np.testing.assert_allclose(
        matrix[2, :, 2],
        np.array([4.0 / 81.0, 28.0 / 81.0, 49.0 / 81.0, 0.0]),
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_allclose(matrix.sum(axis=1), 1.0, rtol=0.0, atol=1e-15)
    expected_survivors = sum(nxt * matrix[2, nxt, 2] for nxt in range(4))
    assert np.isclose(expected_survivors, 2.0 * 7.0 / 9.0)
    assert np.isclose(4.0 / (1.0 - matrix[1, 1, 1]), 18.0)


def test_repair_reuses_existing_dimension_and_preserves_parity_after_exit() -> None:
    P = repaired_parameters()
    assert P.n_child_states == P.n_parity == 4
    assert get_completed_fertility(3, 0, P) == 3
    assert current_child_bin_dt(3, 2, 1, 3, P.child_state_mode) == 3
    assert current_child_bin_dt(3, 3, 1, 3, P.child_state_mode) == 4
    assert current_child_bin_dt(3, 0, 1, 3, P.child_state_mode) == 2


def test_preferences_depend_on_children_at_home_not_ever_born() -> None:
    P = repaired_parameters(
        preference_spec="eqscale",
        eqscale_form="linear",
        gamma_e=0.4,
        delta_alpha=0.05,
        delta_alpha_jump=0.02,
    )
    shared = precompute_shared(P, np.linspace(-1.0, 2.0, 5))
    assert shared.c_bar[3, 0] == 0.0
    assert shared.psi_v[3, 0] == 0.0
    assert shared.psi_v[3, 2] == 2.0 * P.psi_child
    alpha = shared.alpha_flat.reshape(P.n_parity, P.n_child_states, order="F")
    assert np.isclose(alpha[3, 2], P.alpha_cons - P.delta_alpha_jump - 2.0 * P.delta_alpha)


def test_policy_eligibility_covers_every_valid_parent_count() -> None:
    P = repaired_parameters(
        I=1,
        n_house=1,
        H_own=np.array([6.0]),
        parent_dp_waiver=True,
        birth_entry_grant=True,
        birth_entry_grant_amount=0.4,
    )
    np.testing.assert_array_equal(get_parent_target_child_states(P), [False, True, True, True])
    grant = get_birth_entry_grant_tensor(P)
    assert grant[0, 1, 1, 1] == 0.4
    assert grant[0, 1, 2, 1] == 0.4
    assert grant[0, 1, 2, 2] == 0.4
    assert grant[0, 1, 1, 2] == 0.0


def test_parity_progression_remains_available_after_older_child_leaves() -> None:
    solution, parameters, _ = run_model_cp_dt(
        _seq_eqscale(**L4, child_state_mode="independent_count"),
        verbose=False,
    )
    assert parameters._fert2_probs.ndim == 8
    # A parity-one household with no child left at home can still attempt the
    # second birth; the newborn then enters at-home count one.
    assert np.any(parameters._fert2_probs[..., 1, 0, 0] > 0.0)
    assert float(np.sum(solution.g[..., 1, 0])) > 0.0


def test_transfer_and_budget_accounting_are_available_in_repaired_package() -> None:
    P = SimpleNamespace(
        income=np.array([[2.0, 1.0]]),
        J=2,
        J_R=1,
        retirement_income_z_scale=0.0,
        property_tax_lump_sum_transfer=0.3,
    )
    assert income_at_state(P, 0, 0, 1.5) == 3.3
    assert income_at_state(P, 0, 1, 1.5) == 1.3

    P = SimpleNamespace(I=1, n_house=1, H_own=np.array([6.0]), tau_H=0.08)
    distribution = np.zeros((1, 2, 1, 1, 1, 1, 1))
    distribution[0, 0, 0, 0, 0, 0, 0] = 0.4
    distribution[0, 1, 0, 0, 0, 0, 0] = 0.6
    rental_policy = np.zeros_like(distribution)
    rental_policy[0, 0, 0, 0, 0, 0, 0] = 4.0
    revenue = property_tax_revenue_from_distribution(
        distribution, rental_policy, np.array([0.75]), P
    )
    assert np.isclose(revenue, 0.08 * 0.75 * (0.4 * 4.0 + 0.6 * 6.0))


def test_grant_outlays_sum_all_valid_at_home_counts() -> None:
    P = SimpleNamespace(
        I=1,
        J=1,
        n_parity=4,
        n_child_states=4,
        child_state_mode="independent_count",
    )
    distribution = np.zeros((1, 2, 1, 1, 1, 4, 4))
    distribution[0, 0, 0, 0, 0, 2, 1] = 0.25
    distribution[0, 0, 0, 0, 0, 2, 2] = 0.50
    tenure_choice = np.zeros((1, 2, 1, 1, 1, 4, 4), dtype=int)
    tenure_probs = np.zeros((1, 2, 1, 1, 1, 4, 4, 2))
    tenure_probs[0, 0, 0, 0, 0, 2, 1, :] = [0.2, 0.8]
    tenure_probs[0, 0, 0, 0, 0, 2, 2, :] = [0.5, 0.5]
    birth_entry_grant = np.zeros((1, 2, 4, 4))
    birth_entry_grant[0, 1, 2, 1:3] = 0.4
    birth_dp = np.zeros((4, 4, 2, 2), dtype=bool)
    shared = SimpleNamespace(birth_entry_grant=birth_entry_grant, birth_dp=birth_dp)
    recipient_mass, outlays = markov_grant_outlays(
        distribution, tenure_choice, tenure_probs, P, shared
    )
    assert np.isclose(recipient_mass, 0.25 * 0.8 + 0.50 * 0.5)
    assert np.isclose(outlays, 0.4 * recipient_mass)
