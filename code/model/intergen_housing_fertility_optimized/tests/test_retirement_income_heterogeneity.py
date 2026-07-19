from types import SimpleNamespace

import numpy as np
import pytest

from intergen_housing_fertility_optimized.parameters import apply_overrides, setup_parameters
from intergen_housing_fertility_optimized.solver import (
    add_old_wealth_income_moments,
    annual_gross_income_at_state,
    income_at_state,
)


def test_default_retirement_income_is_common_across_income_states() -> None:
    P = SimpleNamespace(income=np.array([[2.0, 3.0]]), J=2, J_R=1)
    assert income_at_state(P, 0, 1, 0.6) == pytest.approx(3.0)
    assert income_at_state(P, 0, 1, 1.4) == pytest.approx(3.0)


def test_retirement_income_scale_loads_persistent_state_and_preserves_mean() -> None:
    P = SimpleNamespace(
        income=np.array([[2.0, 3.0]]),
        J=2,
        J_R=1,
        retirement_income_z_scale=1.5,
    )
    z = np.array([0.6, 1.0, 1.4])
    weights = np.array([0.25, 0.5, 0.25])
    incomes = np.array([income_at_state(P, 0, 1, value) for value in z])
    assert incomes == pytest.approx(3.0 * (1.0 + 1.5 * (z - 1.0)))
    assert float(weights @ incomes) == pytest.approx(3.0)


def test_annual_gross_retirement_income_uses_same_multiplier_without_payroll_grossup() -> None:
    P = SimpleNamespace(
        income=np.array([[2.0, 8.0]]),
        J=2,
        J_R=1,
        period_years=4.0,
        tau_pay=0.2,
        retirement_income_z_scale=1.0,
    )
    assert annual_gross_income_at_state(P, 0, 1, 1.25) == pytest.approx(2.5)


def test_invalid_retirement_income_scale_is_rejected() -> None:
    with pytest.raises(ValueError, match="nonpositive retirement income"):
        apply_overrides(setup_parameters(), {"retirement_income_z_scale": 10.0})


def test_old_age_ratio_family_preserves_markov_income_states() -> None:
    P = SimpleNamespace(
        J=1,
        I=1,
        age_start=65.0,
        da=4.0,
        J_R=0,
        income=np.array([[10.0]]),
        period_years=1.0,
        tau_pay=0.0,
        retirement_income_z_scale=1.0,
        z_grid=np.array([0.5, 1.5]),
        n_parity=2,
        n_child_states=1,
        psi=0.0,
        H_own=np.array([], dtype=float),
    )
    # Equal mass at retirement incomes 5 and 15, for childless and parents.
    distribution = np.ones((1, 1, 1, 1, 2, 2, 1))
    stats = SimpleNamespace()
    add_old_wealth_income_moments(
        stats,
        distribution,
        P,
        np.array([10.0]),
        np.array([1.0]),
    )
    expected_mean = (2.0 + 2.0 / 3.0) / 2.0
    expected_median = 2.0 / 3.0
    for name in (
        "old_nonhousing_wealth_to_income_6575",
        "old_total_wealth_to_income_6575",
        "old_parent_nonhousing_wealth_to_income_6575",
        "old_childless_nonhousing_wealth_to_income_6575",
        "old_parent_total_wealth_to_income_6575",
        "old_childless_total_wealth_to_income_6575",
    ):
        assert getattr(stats, name) == pytest.approx(expected_mean)
    for name in (
        "old_nonhousing_wealth_to_income_median_6575",
        "old_total_wealth_to_income_median_6575",
        "old_parent_nonhousing_wealth_to_income_median_6575",
        "old_childless_nonhousing_wealth_to_income_median_6575",
        "old_parent_total_wealth_to_income_median_6575",
        "old_childless_total_wealth_to_income_median_6575",
    ):
        assert getattr(stats, name) == pytest.approx(expected_median)
