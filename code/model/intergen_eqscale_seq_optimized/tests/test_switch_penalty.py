"""Tests for the default-off children-at-home earnings penalty switch.

Economics: a per-period time cost of children at home, proportional to
earnings. With entries [t0, t1, t2, t3] for at-home counts m = 0, 1, 2, 3+,
working-age after-tax earnings become y^d(a,z)(1 - t(m)) in the household
budget. Pensions, the payroll-tax base, and the PAYGO balance use
unpenalized earnings by construction (see also the comment on
``child_earnings_multiplier`` in parameters.py).
"""

from __future__ import annotations

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    child_earnings_multiplier,
    children_at_home_count,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import (
    income_at_state,
    penalized_income_at_state,
    run_model_cp_dt,
)
from intergen_eqscale_seq_optimized.tests.test_eqscale_seq import _tiny_markov


PENALTY = [0.0, 0.2, 0.2, 0.2]


def test_off_bitwise_identical() -> None:
    """Zeros leave every solved array bit for bit unchanged."""
    base, _, _ = run_model_cp_dt(_tiny_markov(), verbose=False)
    explicit, _, _ = run_model_cp_dt(
        {**_tiny_markov(), "child_earnings_penalty": [0.0, 0.0, 0.0, 0.0]},
        verbose=False,
    )
    assert np.array_equal(base.V, explicit.V)
    assert np.array_equal(base.g, explicit.g)
    assert np.array_equal(base.c_pol, explicit.c_pol)
    assert np.array_equal(base.hR_pol, explicit.hR_pol)
    assert np.array_equal(base.bp_pol, explicit.bp_pol)
    assert float(base.mean_income) == float(explicit.mean_income)


def test_working_age_resources_scaled_retirement_untouched() -> None:
    """With [0, 0.2, 0.2, 0.2], resources at m >= 1 are 0.8x while working."""
    P = apply_overrides(setup_parameters(), {**_tiny_markov(), "child_earnings_penalty": PENALTY})
    j_work, j_ret = 0, int(P.J_R)
    assert j_ret < int(P.J)
    for m in (1, 2, 3):
        assert penalized_income_at_state(P, 0, j_work, 1.0, m) == pytest.approx(
            0.8 * income_at_state(P, 0, j_work, 1.0), rel=0.0, abs=0.0
        )
        assert penalized_income_at_state(P, 0, j_ret, 1.0, m) == pytest.approx(
            income_at_state(P, 0, j_ret, 1.0), rel=0.0, abs=0.0
        )
    assert penalized_income_at_state(P, 0, j_work, 1.0, 0) == pytest.approx(
        income_at_state(P, 0, j_work, 1.0), rel=0.0, abs=0.0
    )
    # The penalty binds in equilibrium: values move when it is turned on.
    base, _, _ = run_model_cp_dt(_tiny_markov(), verbose=False)
    penalized, _, _ = run_model_cp_dt(
        {**_tiny_markov(), "child_earnings_penalty": PENALTY}, verbose=False
    )
    assert not np.array_equal(base.V, penalized.V)
    assert np.all(np.isfinite(penalized.V))


def test_at_home_count_mapping() -> None:
    """At-home counts: shared clock flags dependent stages; counts cap at 3."""
    P = apply_overrides(setup_parameters(), _tiny_markov())
    assert children_at_home_count(0, 0, P) == 0
    assert children_at_home_count(1, 1, P) == 1
    assert children_at_home_count(0, 1, P) == 0
    assert child_earnings_multiplier(P, 0, 5) == pytest.approx(
        child_earnings_multiplier(P, 0, 3), rel=0.0, abs=0.0
    )
    assert child_earnings_multiplier(P, int(P.J_R), 2) == pytest.approx(1.0)


def test_forward_step_conserves_mass() -> None:
    """The penalty changes budgets, not transitions: mass is conserved."""
    sol, P, _ = run_model_cp_dt(
        {**_tiny_markov(), "child_earnings_penalty": PENALTY}, verbose=False
    )
    assert float(np.sum(sol.g)) == pytest.approx(1.0, abs=1e-9)
    # Fiscal objects built on unpenalized earnings do not move with children.
    assert float(P.pension) == pytest.approx(
        float(apply_overrides(setup_parameters(), _tiny_markov()).pension)
    )


def test_penalty_rejects_bad_values() -> None:
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"child_earnings_penalty": [0.0, 0.2]})
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"child_earnings_penalty": [0.0, 0.2, 0.2, 1.0]})
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"child_earnings_penalty": [0.0, -0.1, 0.0, 0.0]})


def test_scalar_penalty_broadcasts() -> None:
    P = apply_overrides(setup_parameters(), {"child_earnings_penalty": 0.1})
    assert np.array_equal(np.asarray(P.child_earnings_penalty), np.full(4, 0.1))
