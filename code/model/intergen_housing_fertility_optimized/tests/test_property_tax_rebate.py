from types import SimpleNamespace

import numpy as np

from intergen_housing_fertility_optimized.solver import (
    income_at_state,
    markov_grant_outlays,
    property_tax_revenue_from_distribution,
)


def test_lump_sum_transfer_enters_worker_and_retiree_budgets() -> None:
    parameters = SimpleNamespace(
        income=np.array([[2.0, 1.0]]),
        J=2,
        J_R=1,
        retirement_income_z_scale=0.0,
        property_tax_lump_sum_transfer=0.3,
    )
    assert income_at_state(parameters, 0, 0, 1.5) == 3.3
    assert income_at_state(parameters, 0, 1, 1.5) == 1.3


def test_property_tax_revenue_uses_rental_and_owner_housing() -> None:
    parameters = SimpleNamespace(I=1, n_house=1, H_own=np.array([6.0]), tau_H=0.08)
    distribution = np.zeros((1, 2, 1, 1, 1, 1, 1))
    distribution[0, 0, 0, 0, 0, 0, 0] = 0.4
    distribution[0, 1, 0, 0, 0, 0, 0] = 0.6
    rental_policy = np.zeros_like(distribution)
    rental_policy[0, 0, 0, 0, 0, 0, 0] = 4.0
    revenue = property_tax_revenue_from_distribution(
        distribution,
        rental_policy,
        np.array([0.75]),
        parameters,
    )
    assert np.isclose(revenue, 0.08 * 0.75 * (0.4 * 4.0 + 0.6 * 6.0))


def test_grant_outlays_count_only_eligible_renter_purchases() -> None:
    parameters = SimpleNamespace(I=1, J=1, n_parity=2)
    distribution = np.zeros((1, 2, 1, 1, 1, 2, 2))
    distribution[0, 0, 0, 0, 0, 1, 1] = 0.25
    tenure_choice = np.zeros((1, 2, 1, 1, 1, 2, 2), dtype=int)
    tenure_probs = np.zeros((1, 2, 1, 1, 1, 2, 2, 2))
    tenure_probs[0, 0, 0, 0, 0, 1, 1, :] = np.array([0.2, 0.8])
    birth_entry_grant = np.zeros((1, 2, 2, 2))
    birth_entry_grant[0, 1, 1, 1] = 0.4
    birth_dp = np.zeros((2, 2, 2, 2), dtype=bool)
    shared = SimpleNamespace(birth_entry_grant=birth_entry_grant, birth_dp=birth_dp)
    recipient_mass, outlays = markov_grant_outlays(
        distribution,
        tenure_choice,
        tenure_probs,
        parameters,
        shared,
    )
    assert np.isclose(recipient_mass, 0.20)
    assert np.isclose(outlays, 0.08)
