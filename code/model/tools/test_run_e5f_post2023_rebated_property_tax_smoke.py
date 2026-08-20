#!/usr/bin/env python3
"""Focused contract tests for the exact rebated property-tax smoke."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

import run_e5f_post2023_rebated_property_tax_smoke as tax


def test_case_menu_excludes_purchase_grant() -> None:
    assert tuple(tax.CASES) == (
        "status-quo-tax1-unrebated",
        "tax1-equal-rebate",
        "tax2-equal-rebate",
    )
    assert not tax.CASES["status-quo-tax1-unrebated"].rebate_revenue
    assert tax.CASES["tax1-equal-rebate"].rebate_revenue
    assert tax.CASES["tax2-equal-rebate"].annual_tax_rate == 0.02


def test_tax_setter_disables_grants_and_updates_user_cost() -> None:
    parameters = SimpleNamespace(period_years=4, q=0.03, delta=0.02)
    tax.set_tax_policy(parameters, tax.CASES["tax2-equal-rebate"], 0.25)
    assert parameters.tau_H == 0.08
    assert parameters.user_cost_rate == 0.13
    assert parameters.property_tax_lump_sum_transfer == 0.25
    assert parameters.birth_entry_grant is False
    assert parameters.birth_entry_grant_amount == 0.0
    assert parameters.birth_entry_grant_locations.size == 0
    assert parameters.birth_entry_grant_owner_rungs.size == 0


def test_fiscal_ledger_equal_rebate_identity(monkeypatch) -> None:
    evaluation = SimpleNamespace(
        g_current=np.ones((2,)),
        policy=SimpleNamespace(hR_pol=np.array([1.0]), price=np.array([2.0])),
    )
    parameters = SimpleNamespace(property_tax_lump_sum_transfer=0.5)
    monkeypatch.setattr(
        tax.calendar.model,
        "property_tax_revenue_from_distribution",
        lambda *args: 1.0,
    )
    ledger = tax.fiscal_ledger(evaluation, parameters)
    assert ledger["property_tax_revenue"] == 1.0
    assert ledger["equal_transfer_outlays"] == 1.0
    assert ledger["government_budget_residual"] == 0.0


def test_renewal_lag_contract() -> None:
    result = tax.renewal_lag_gate()
    assert result["passed"]
    assert result["waiting_slots"] == 4
    assert result["birth_to_population_effect_years"] == 20
    assert np.isclose(result["conversion"], 1.0 / 2.1)


if __name__ == "__main__":
    test_case_menu_excludes_purchase_grant()
    test_tax_setter_disables_grants_and_updates_user_cost()
    test_renewal_lag_contract()
    print("POST2023_REBATED_PROPERTY_TAX_TESTS_PASS tests=3")
