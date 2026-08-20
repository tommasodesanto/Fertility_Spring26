#!/usr/bin/env python3
"""Focused tests for the short-run rebated property-tax path."""

from __future__ import annotations

import numpy as np

import run_e5f_post2023_coven_property_tax_path as path
import run_e5f_post2023_no_policy_continuations as baseline


def state() -> baseline.DynamicState:
    return baseline.DynamicState(
        g_pre=np.array([[[1.0, 2.0]]]),
        scheduled_entries=[1.0, 2.0, 3.0, 4.0],
        scheduled_raw_entries=[0.5, 1.0, 1.5, 2.0],
        price_guess=0.6,
        initial_policy=None,
    )


def test_copy_state_is_independent() -> None:
    original = state()
    copied = path.copy_state(original)
    copied.g_pre[0, 0, 0] = 99.0
    copied.scheduled_entries[0] = 99.0
    assert original.g_pre[0, 0, 0] == 1.0
    assert original.scheduled_entries[0] == 1.0


def test_period_gate_accepts_clean_row_and_rejects_fiscal_failure() -> None:
    row = {
        "calendar_year": 2023,
        "relative_market_residual": 1.0e-6,
        "government_budget_residual": 1.0e-6,
        "mass_accounting_residual": 1.0e-12,
        "nonfinite_distribution_count": 0,
        "min_distribution_mass": 0.0,
        "feasibility_frontier_projection_mass": 0.0,
    }
    gates = path.gate_period_row(row, state())
    assert gates["calendar_year"] == 2023
    bad = dict(row, government_budget_residual=1.0)
    try:
        path.gate_period_row(bad, state())
    except RuntimeError as exc:
        assert "Fiscal gate failed" in str(exc)
    else:
        raise AssertionError("Fiscal failure was accepted")


def test_effect_rows_use_matched_dates_and_correct_units() -> None:
    base = [{
        "calendar_year": 2035,
        "years_from_2023": 12,
        "asset_price": 2.0,
        "housing_user_cost": 0.4,
        "equal_transfer_period_units": 1.0,
        "owner_rate": 0.50,
        "dependent_child_owner_rate": 0.40,
        "topcode_adjusted_births_per_adult": 0.08,
        "adult_population": 1.0,
    }]
    reform = [{
        "calendar_year": 2035,
        "years_from_2023": 12,
        "asset_price": 1.8,
        "housing_user_cost": 0.5,
        "equal_transfer_period_units": 1.5,
        "owner_rate": 0.52,
        "dependent_child_owner_rate": 0.43,
        "topcode_adjusted_births_per_adult": 0.084,
        "adult_population": 1.0,
    }]
    result = path.build_effect_rows(base, reform)[0]
    assert np.isclose(result["asset_price_percent_change"], -10.0)
    assert np.isclose(result["housing_user_cost_percent_change"], 25.0)
    assert np.isclose(result["equal_transfer_percent_change"], 50.0)
    assert np.isclose(result["owner_rate_pp_change"], 2.0)
    assert np.isclose(result["dependent_child_owner_rate_pp_change"], 3.0)
    assert np.isclose(result["births_per_adult_percent_change"], 5.0)
    assert result["adult_population_percent_change"] == 0.0


if __name__ == "__main__":
    test_copy_state_is_independent()
    test_period_gate_accepts_clean_row_and_rejects_fiscal_failure()
    test_effect_rows_use_matched_dates_and_correct_units()
    print("POST2023_REBATED_PROPERTY_TAX_PATH_TESTS_PASS tests=3")
