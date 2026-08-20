#!/usr/bin/env python3
"""Focused tests for the Coven-style property-tax impact smoke."""

from __future__ import annotations

import numpy as np

import run_e5f_post2023_coven_property_tax_smoke as coven


def test_coven_policy_contract_is_two_rebated_tax_cases() -> None:
    assert tuple(coven.CASES) == (
        "rebated-tax1-baseline",
        "rebated-tax2-reform",
    )
    assert all(case.rebate_revenue for case in coven.CASES.values())
    assert coven.CASES["rebated-tax1-baseline"].annual_tax_rate == 0.01
    assert coven.CASES["rebated-tax2-reform"].annual_tax_rate == 0.02


def test_supply_reanchor_holds_inherited_price_and_stock() -> None:
    for elasticity in (0.232, 1.75):
        rule, gate = coven.reanchor_supply_rule(
            inherited_price=0.6,
            baseline_demand=7.5,
            elasticity=elasticity,
        )
        assert gate["passed"]
        assert gate["relative_identity_gap"] == 0.0
        assert rule.initial_price == 0.6
        assert rule.initial_stock == 7.5
        assert rule.elasticity == elasticity
        assert gate["supply_elasticity"] == elasticity
        assert np.isclose(rule.quantity(np.array([0.6]))[0], 7.5)


def test_supply_reanchor_rejects_invalid_elasticity() -> None:
    for elasticity in (0.0, -0.1, float("nan")):
        try:
            coven.reanchor_supply_rule(
                inherited_price=0.6,
                baseline_demand=7.5,
                elasticity=elasticity,
            )
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid elasticity accepted: {elasticity}")


if __name__ == "__main__":
    test_coven_policy_contract_is_two_rebated_tax_cases()
    test_supply_reanchor_holds_inherited_price_and_stock()
    test_supply_reanchor_rejects_invalid_elasticity()
    print("POST2023_COVEN_PROPERTY_TAX_TESTS_PASS tests=3")
