"""Focused accounting tests; run on Torch, never on the author's Mac."""
import inspect
from types import SimpleNamespace

import numpy as np
import pytest

import e5f_estate_receiver_adapter as adapter


def parameters():
    return SimpleNamespace(I=1, J=17, J_R=12, age_start=18., da=4.,
        period_years=4., psi=.06, psi_child=.135, tau_pay=.08751017424959717,
        H_own=np.array([2.]), survival_probs=np.ones(16), use_age_survival=True,
        use_postdecision_current_distribution=True)


def fixture(P, *, future_transfer=0.):
    g = np.zeros((1, 2, 1, 17, 1, 1, 1))
    g[0, 0, 0, 7] = .5  # recipient cell starts at age46
    g[0, 1, 0, 16] = .5  # terminal owner age82
    bp = np.full_like(g, 1. + .05 * future_transfer)
    gross = .5 * (1. + .05 * future_transfer + 2.)
    return SimpleNamespace(g=g, bp_pol=bp, annual_bequest_flow=gross/4.)


def test_estate_uses_post_saving_wealth_net_housing_and_period_units():
    P = adapter.configure(parameters(), "net_valuation_transfer", 2.88)
    receipt = adapter.estate_accounts(fixture(P), P, np.array([1.]))
    assert receipt["generated_gross_period"] == pytest.approx(1.5)
    assert receipt["generated_net_period"] == pytest.approx(1.44)
    assert receipt["paid_period"] == pytest.approx(1.44)
    assert receipt["residual"] == pytest.approx(0.)
    assert receipt["recipient_ages"] == [46., 50., 54., 58., 62.]


def test_earlier_mortality_and_negative_estates_are_handled_once():
    P = adapter.configure(parameters(), "control")
    P.survival_probs[7] = .8
    sol = fixture(P)
    sol.bp_pol[0, 0, 0, 7] = 5.
    sol.bp_pol[0, 1, 0, 16] = -3.  # insolvent estate does not fund recipients
    sol.annual_bequest_flow = .5 / 4.
    receipt = adapter.estate_accounts(sol, P, [1.])
    assert receipt["generated_gross_period"] == pytest.approx(.5)
    assert receipt["generated_net_period"] == pytest.approx(.5)
    assert receipt["paid_period"] == 0.


def test_rejects_mismatched_death_observer_and_double_receivers():
    P = parameters()
    sol = fixture(P)
    sol.annual_bequest_flow *= 2
    with pytest.raises(RuntimeError, match="native gross-flow"):
        adapter.estate_accounts(sol, P, [1.])
    P.estate_receiver = "ages_45_65"
    with pytest.raises(ValueError, match="Reference already"):
        adapter.configure(P, "control")


def solve_bellman_full_markov_income(r_hat, p_hat, P, b_grid, SD, continuation_V=None):
    # Minimal source fixture retains the exact production estate-table anchor.
    income_for_purchase = 0.
    bmo_purchase = income_for_purchase
    out = []
    for i in range(P.I):
        for ten in range(2):
            hv = p_hat[i] * P.H_own[ten - 1] if ten > 0 else 0.0
            out.append(hv)
    return out, bmo_purchase


def test_hooks_leave_control_unchanged_and_exclude_receipts_from_earnings(tmp_path):
    # A module-like shared namespace matches the frozen solver's global lookup.
    model = SimpleNamespace()
    model.__dict__["np"] = np
    exec(inspect.getsource(solve_bellman_full_markov_income), model.__dict__)
    model.income_at_state = lambda P, i, j, z: 4. * z
    exec("def annual_gross_income_at_state(P,i,j,z):\n"
         "    return income_at_state(P,i,j,z)/4./(1.-P.tau_pay)\n", model.__dict__)
    P = parameters()
    original = model.solve_bellman_full_markov_income(None, [1.], P, None, None)
    # inspect needs a file-backed function, as in the real generated runtime.
    model.solve_bellman_full_markov_income = solve_bellman_full_markov_income
    adapter.install(model, tmp_path)
    control = adapter.configure(P, "control")
    receiving = adapter.configure(P, "net_valuation_transfer", .5)
    assert model.solve_bellman_full_markov_income(None, [1.], control, None, None) == original
    assert model.solve_bellman_full_markov_income(None, [1.], receiving, None, None)[0] == [0., 1.88]
    assert model.income_at_state(receiving, 0, 7, 1.) == 4.5
    assert model.income_at_state(receiving, 0, 0, 1.) == 4.
    assert model.annual_gross_income_at_state(receiving, 0, 7, 1.) == pytest.approx(1./(1.-P.tau_pay))
    assert model.annual_gross_income_at_state(control, 0, 7, 1.) == pytest.approx(1./(1.-P.tau_pay))
    with pytest.raises(RuntimeError, match="exactly once"):
        adapter.install(model, tmp_path)


def test_exact_case_loop_balances_receipts_and_preserves_preferences(tmp_path):
    import time
    model = SimpleNamespace(_estate_probe_installed=True)
    base = parameters()
    calls = []

    def native(**kwargs):
        P = kwargs["parameters"]
        calls.append((P.estate_probe_case, P.psi_child, kwargs["payroll_tax"]))
        sol = fixture(P, future_transfer=P.estate_probe_transfer)
        return sol, P, np.array([1.]), {"fiscal_gate": True}

    for case in adapter.CASES:
        _, result_P, _, _, receipt = adapter.solve_case(
            model=model, native_solver=native, parameters=base, b_grid=np.array([0.]),
            initial_prices=np.array([1.]), case=case, deadline_epoch=time.time()+30.,
            output_dir=tmp_path/case)
        assert result_P.psi_child == base.psi_child
        assert receipt["converged"]
        if case == adapter.TRANSFER_CASE:
            assert abs(receipt["residual"]) <= 1e-10
            assert receipt["transfer"] == pytest.approx(2.88/.95, abs=1e-9)
        else:
            assert receipt["native_solves"] == 1
    assert all(psi == base.psi_child and tax == base.tau_pay for _, psi, tax in calls)
    assert base.__dict__.get("estate_probe_case") is None


def test_budget_and_no_recipient_fail_closed(tmp_path):
    import time
    model = SimpleNamespace(_estate_probe_installed=True)
    def must_not_run(**kwargs):
        raise AssertionError("Expired budget called native solver")
    with pytest.raises(TimeoutError):
        adapter.solve_case(model=model, native_solver=must_not_run,
            parameters=parameters(), b_grid=[0.], initial_prices=[1.], case="control",
            deadline_epoch=time.time()-1., output_dir=tmp_path/"expired")
    with pytest.raises(ValueError, match="forty"):
        adapter.solve_case(model=model, native_solver=must_not_run,
            parameters=parameters(), b_grid=[0.], initial_prices=[1.], case="control",
            deadline_epoch=time.time()+1., output_dir=tmp_path/"oversized", max_solves=41)
