from types import SimpleNamespace

import numpy as np
import pytest

from intergen_housing_fertility.parameters import apply_overrides, setup_parameters
from intergen_housing_fertility.solver import annual_gross_income_at_state, income_at_state


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
