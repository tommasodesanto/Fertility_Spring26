#!/usr/bin/env python3
"""Focused contract tests for the post-2023 policy mechanism driver."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

import run_dynamic_population_transition as calendar
import run_e5f_post2023_policy_mechanisms as policy


def test_policy_menu_is_small_and_predeclared() -> None:
    assert tuple(policy.POLICIES) == (
        "baseline",
        "supply-plus-20",
        "dependent-child-ltv95",
        "combined",
        "property-tax-2pct-no-rebate",
    )
    assert policy.POLICIES["baseline"].supply_multiplier == 1.0
    assert policy.POLICIES["supply-plus-20"].supply_multiplier == 1.2
    assert policy.POLICIES["dependent-child-ltv95"].dependent_child_ltv95
    assert policy.POLICIES["combined"].dependent_child_ltv95
    assert policy.POLICIES["baseline"].annual_property_tax_rate == 0.01
    assert (
        policy.POLICIES["property-tax-2pct-no-rebate"].annual_property_tax_rate
        == 0.02
    )


def test_supply_policy_shifts_only_the_intercept() -> None:
    source = calendar.HousingSupplyRule(
        mode="static-elastic",
        initial_price=2.0,
        initial_stock=5.0,
        elasticity=1.75,
    )
    shifted = policy.policy_supply_rule(source, policy.POLICIES["supply-plus-20"])
    assert shifted.initial_price == source.initial_price
    assert shifted.elasticity == source.elasticity
    assert shifted.initial_stock == 6.0
    assert np.isclose(shifted.quantity(np.array([2.0]))[0], 6.0)
    assert source.initial_stock == 5.0


def test_credit_policy_is_dependent_child_not_birth_only() -> None:
    parameters = SimpleNamespace(period_years=4.0, q=0.02, delta=0.03)
    policy.apply_policy(parameters, policy.POLICIES["dependent-child-ltv95"])
    assert parameters.parent_dp_waiver is True
    assert parameters.parent_dp_waiver_phi == 0.95
    assert parameters.parent_dp_waiver_birth_state_only is False
    assert parameters.parent_dp_waiver_locations.size == 0
    assert parameters.parent_dp_waiver_owner_rungs.size == 0
    policy.apply_policy(parameters, policy.POLICIES["baseline"])
    assert parameters.parent_dp_waiver is False
    assert parameters.parent_dp_waiver_phi == 1.0


def test_property_tax_policy_is_unrebated_and_has_no_grant() -> None:
    parameters = SimpleNamespace(period_years=4.0, q=0.02, delta=0.03)
    policy.apply_policy(
        parameters, policy.POLICIES["property-tax-2pct-no-rebate"]
    )
    assert parameters.tau_H == 0.08
    assert parameters.user_cost_rate == 0.13
    assert parameters.property_tax_lump_sum_transfer == 0.0
    assert parameters.birth_entry_grant is False
    assert parameters.birth_entry_grant_amount == 0.0
    assert parameters.parent_dp_waiver is False


if __name__ == "__main__":
    test_policy_menu_is_small_and_predeclared()
    test_supply_policy_shifts_only_the_intercept()
    test_credit_policy_is_dependent_child_not_birth_only()
    test_property_tax_policy_is_unrebated_and_has_no_grant()
    print("POST2023_POLICY_MECHANISM_TESTS_PASS tests=4")
