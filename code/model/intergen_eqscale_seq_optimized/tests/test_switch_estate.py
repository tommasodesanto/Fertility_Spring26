"""Tests for the default-off estates-paid-to-households switch.

Economics: a dying household's net estate (liquid holdings plus housing net
of the selling cost) is paid as an equal per-household lump sum to everyone
aged 45-65. The transfer enters income exactly like the property-tax rebate
(through ``income_at_state``), so it is in the budget but not in the
down-payment test (which reads ``dp_arr``/``bmo``, never income). Off,
estates enter bequest utility gross and leave the economy.
"""

from __future__ import annotations

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.parameters import (
    apply_overrides,
    estate_recipient_age_indices,
    setup_parameters,
)
from intergen_eqscale_seq_optimized.solver import (
    income_at_state,
    penalized_income_at_state,
    run_model_cp_dt,
)
from intergen_eqscale_seq_optimized.tests.test_eqscale_seq import _tiny_markov


def _tiny_estate() -> dict[str, object]:
    cfg = _tiny_markov()
    cfg.update(
        {
            "J": 12,
            "J_R": 9,
            "Nb": 12,
            "b_core_hi": 3.0,
            "b_mid_hi": 6.0,
            "b_max": 10.0,
            "n_house": 1,
            "H_own": np.array([4.0]),
            "H0": np.array([4.0]),
        }
    )
    return cfg


def test_off_bitwise_identical() -> None:
    """Off (and explicit off) leaves every solved array bit for bit unchanged."""
    base, _, _ = run_model_cp_dt(_tiny_markov(), verbose=False)
    explicit, _, _ = run_model_cp_dt(
        {**_tiny_markov(), "estate_receiver": "none", "bequest_net_of_selling_cost": False},
        verbose=False,
    )
    assert np.array_equal(base.V, explicit.V)
    assert np.array_equal(base.g, explicit.g)
    assert np.array_equal(base.c_pol, explicit.c_pol)
    assert np.array_equal(base.hR_pol, explicit.hR_pol)
    assert np.array_equal(base.bp_pol, explicit.bp_pol)


def test_bequest_net_zero_selling_cost_is_inert() -> None:
    """With psi=0 net equals gross, so the utility switch is inert."""
    cfg = {**_tiny_markov(), "psi": 0.0}
    gross, _, _ = run_model_cp_dt(cfg, verbose=False)
    net, _, _ = run_model_cp_dt({**cfg, "bequest_net_of_selling_cost": True}, verbose=False)
    assert np.array_equal(gross.V, net.V)
    assert np.array_equal(gross.g, net.g)


def test_bequest_net_moves_values_with_positive_selling_cost() -> None:
    """With psi>0 the utility-side switch revalues bequests and moves V."""
    base, _, _ = run_model_cp_dt(_tiny_markov(), verbose=False)
    net, _, _ = run_model_cp_dt(
        {**_tiny_markov(), "bequest_net_of_selling_cost": True}, verbose=False
    )
    assert not np.array_equal(base.V, net.V)
    assert np.all(np.isfinite(net.V))


def test_transfer_enters_income_like_rebate() -> None:
    """The lump sum adds to income at ages 45-65 only, unscaled by penalties."""
    P_base = apply_overrides(setup_parameters(), _tiny_estate())
    P_on = apply_overrides(
        setup_parameters(),
        {
            **_tiny_estate(),
            "estate_receiver": "ages_45_65",
            "estate_lump_sum_transfer": 0.5,
            "child_earnings_penalty": [0.0, 0.2, 0.2, 0.2],
        },
    )
    recipients = set(int(k) for k in estate_recipient_age_indices(P_on))
    assert len(recipients) > 0
    for j in range(int(P_on.J)):
        assert income_at_state(P_on, 0, j, 1.0) == (
            income_at_state(P_base, 0, j, 1.0) + (0.5 if j in recipients else 0.0)
        )
    j_rec = min(recipients)
    if j_rec < int(P_on.J_R):
        penalized_on = penalized_income_at_state(P_on, 0, j_rec, 1.0, 1)
        P_nopen = apply_overrides(
            setup_parameters(),
            {**_tiny_estate(), "child_earnings_penalty": [0.0, 0.2, 0.2, 0.2]},
        )
        assert penalized_on == penalized_income_at_state(P_nopen, 0, j_rec, 1.0, 1) + 0.5


def test_rejects_bad_values() -> None:
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"estate_receiver": "everyone"})
    with pytest.raises(ValueError):
        apply_overrides(setup_parameters(), {"estate_lump_sum_transfer": -0.1})


def test_on_flow_clears_and_mass_conserved() -> None:
    """On: paid flow equals generated flow to 1e-10; mass conserved; smoke."""
    sol, P, _ = run_model_cp_dt(
        {**_tiny_estate(), "estate_receiver": "ages_45_65", "bequest_net_of_selling_cost": True},
        verbose=False,
    )
    assert bool(getattr(sol, "estate_converged", False))
    assert abs(float(sol.estate_flow_paid) - float(sol.estate_flow_generated)) <= 1e-10
    assert float(sol.estate_transfer) >= 0.0
    assert float(sol.estate_recipient_mass) > 0.0
    assert float(np.sum(sol.g)) == pytest.approx(1.0, abs=1e-9)
    assert np.all(np.isfinite(sol.V))
    assert float(P.estate_lump_sum_transfer) == float(sol.estate_transfer)


def test_on_without_utility_net_also_solves() -> None:
    """The transfer and the utility-side revaluation separate cleanly."""
    sol, _, _ = run_model_cp_dt(
        {**_tiny_estate(), "estate_receiver": "ages_45_65"},
        verbose=False,
    )
    assert bool(getattr(sol, "estate_converged", False))
    assert abs(float(sol.estate_flow_paid) - float(sol.estate_flow_generated)) <= 1e-10
    assert float(np.sum(sol.g)) == pytest.approx(1.0, abs=1e-9)
