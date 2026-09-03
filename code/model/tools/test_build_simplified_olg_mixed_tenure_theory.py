#!/usr/bin/env python3
"""Focused first-order-condition tests for the simplified OLG construction."""

from __future__ import annotations

import numpy as np

import build_simplified_olg_mixed_tenure_theory as theory


def solved_endpoints() -> tuple[dict, dict]:
    initial = theory.solve_steady_state(theory.P.theta_initial, 0.34)
    final = theory.solve_steady_state(theory.P.theta_final, 0.48)
    return initial, final


def test_end_of_period_rent_timing_matches_stationary_user_cost() -> None:
    for price in (0.25, 0.50, 0.75):
        rent = (
            (theory.P.gross_return + theory.P.tau_p) * price - price
        )
        zero_profit_residual = (
            theory.P.q * (rent + price)
            - price
            - theory.P.q * theory.P.tau_p * price
        )
        stationary_rent_rate = (
            theory.P.gross_return - 1.0 + theory.P.tau_p
        )
        assert abs(zero_profit_residual) < 1e-12
        assert np.isclose(rent / price, stationary_rent_rate, atol=1e-12)


def test_capped_old_renter_reduction_and_euler_equation() -> None:
    for steady_state in solved_endpoints():
        renter = steady_state["young"]["renter"]
        price = steady_state["price"]
        transfer = steady_state["transfer"]
        rent = (
            (theory.P.gross_return + theory.P.tau_p) * price - price
        )
        rental_user_cost = theory.P.q * rent
        saving_weight = theory.P.beta * (1.0 + theory.P.omega_b)
        lifetime_wealth = (
            theory.P.income
            + theory.P.liquid_wealth
            + transfer
            + theory.P.q * transfer
        )
        effective_wealth = (
            lifetime_wealth
            - theory.P.q**2 * rent * theory.P.renter_cap
        )
        lifetime_budget_left = (
            (1.0 + saving_weight) * renter["usable_consumption"]
            + theory.P.chi * renter["fertility"]
            + rental_user_cost * renter["housing"]
        )
        expected_saving = (
            saving_weight * renter["usable_consumption"]
            - theory.P.q * transfer
            + theory.P.q**2 * rent * theory.P.renter_cap
        )
        budget_residual = (
            renter["gross_consumption"]
            + renter["saving"]
            + rental_user_cost * renter["housing"]
            - (theory.P.income + theory.P.liquid_wealth + transfer)
        )
        euler_residual = (
            1.0 / renter["usable_consumption"]
            - theory.P.beta
            * theory.P.gross_return
            / renter["future_old"]["consumption"]
        )
        family_space = (
            renter["housing"] - theory.P.kappa * renter["fertility"]
        )
        fertility_residual = (
            steady_state["theta"] / renter["fertility"]
            - theory.P.chi / renter["usable_consumption"]
            - theory.P.alpha * theory.P.kappa / family_space
        )
        old = renter["future_old"]
        old_budget_residual = (
            old["consumption"]
            + theory.P.q * old["estate"]
            + rental_user_cost * old["housing"]
            - (theory.P.gross_return * renter["saving"] + transfer)
        )
        old_estate_foc_residual = (
            theory.P.q * old["estate"]
            - theory.P.omega_b * old["consumption"]
        )

        assert renter["old_branch"] == "renter_cap"
        assert renter["future_old"]["cap_binds"]
        assert renter["future_old"]["cap_shadow_wedge"] > 0.0
        assert np.isclose(lifetime_budget_left, effective_wealth, atol=1e-12)
        assert np.isclose(renter["saving"], expected_saving, atol=1e-12)
        assert abs(budget_residual) < 1e-12
        assert np.isclose(
            renter["gross_consumption"],
            renter["usable_consumption"]
            + theory.P.chi * renter["fertility"],
            atol=1e-12,
        )
        assert abs(euler_residual / (1.0 / renter["usable_consumption"])) < 1e-12
        assert abs(fertility_residual) < 1e-12
        assert abs(renter["child_price_identity_residual"]) < 1e-12
        assert abs(old_budget_residual) < 1e-12
        assert abs(old_estate_foc_residual) < 1e-12


def test_owner_reduction_and_euler_equation() -> None:
    for steady_state in solved_endpoints():
        owner = steady_state["young"]["owner"]
        price = steady_state["price"]
        transfer = steady_state["transfer"]
        budget_residual = (
            owner["gross_consumption"]
            + owner["saving"]
            + ((1.0 - theory.P.phi) + theory.P.q * theory.P.tau_p)
            * price
            * owner["housing"]
            - (theory.P.income + theory.P.liquid_wealth + transfer)
        )
        euler_residual = (
            1.0 / owner["usable_consumption"]
            - theory.P.beta
            * theory.P.gross_return
            / owner["future_old"]["consumption"]
        )
        family_space = owner["housing"] - theory.P.kappa * owner["fertility"]
        fertility_residual = (
            steady_state["theta"] / owner["fertility"]
            - theory.P.chi / owner["usable_consumption"]
            - theory.P.alpha * theory.P.kappa / family_space
        )
        old = owner["future_old"]
        old_budget_residual = (
            old["consumption"]
            + theory.P.q * old["estate"]
            + old["incumbent_cost"] * old["housing"]
            - (
                owner["next_financial_wealth"]
                + (price - old["gain_tax"]) * owner["housing"]
                + transfer
            )
        )
        old_housing_foc_residual = (
            theory.P.gamma * old["consumption"] / old["housing"]
            - old["incumbent_cost"]
        )
        old_estate_foc_residual = (
            theory.P.q * old["estate"]
            - theory.P.omega_b * old["consumption"]
        )

        assert owner["old_branch"] == "interior"
        assert owner["future_old"]["retention_slack"] > 0.0
        assert owner["future_old"]["estate_slack"] > 0.0
        assert owner["access_wedge"] > 0.0
        assert abs(budget_residual) < 1e-12
        assert np.isclose(
            owner["gross_consumption"],
            owner["usable_consumption"]
            + theory.P.chi * owner["fertility"],
            atol=1e-12,
        )
        assert abs(euler_residual / (1.0 / owner["usable_consumption"])) < 1e-12
        assert abs(fertility_residual) < 1e-12
        assert abs(owner["child_price_identity_residual"]) < 1e-12
        assert abs(old_budget_residual) < 1e-12
        assert abs(old_housing_foc_residual) < 1e-12
        assert abs(old_estate_foc_residual) < 1e-12


def test_endpoint_market_and_fiscal_conditions() -> None:
    for steady_state in solved_endpoints():
        cohort_mass = steady_state["cohort_mass"]
        housing_residual = (
            cohort_mass * steady_state["housing_per_cohort_pair"]
            - theory.P.housing_stock
        )
        fiscal_residual = (
            2.0 * cohort_mass * steady_state["transfer"]
            - theory.P.q
            * theory.P.tau_p
            * steady_state["price"]
            * theory.P.housing_stock
        )
        replacement_residual = (
            theory.P.nu * steady_state["young"]["average_fertility"] - 1.0
        )

        assert abs(housing_residual) < 1e-12
        assert abs(fiscal_residual) < 1e-10
        assert abs(replacement_residual) < 1e-10


if __name__ == "__main__":
    test_end_of_period_rent_timing_matches_stationary_user_cost()
    test_capped_old_renter_reduction_and_euler_equation()
    test_owner_reduction_and_euler_equation()
    test_endpoint_market_and_fiscal_conditions()
    print("SIMPLIFIED_OLG_TESTS_PASS tests=4")
