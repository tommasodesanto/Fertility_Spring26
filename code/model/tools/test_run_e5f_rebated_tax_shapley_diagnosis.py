#!/usr/bin/env python3
"""Focused tests for the rebated-tax Shapley decomposition."""

from __future__ import annotations

import itertools
import csv
from pathlib import Path

import numpy as np

import run_e5f_rebated_tax_shapley_diagnosis as diagnosis


def test_shapley_recovers_additive_components() -> None:
    cells = {
        cell: 4.0 + 2.0 * cell[0] - 3.0 * cell[1] + 5.0 * cell[2]
        for cell in itertools.product((0, 1), repeat=3)
    }
    result = diagnosis.shapley_decomposition(cells)
    assert np.isclose(result["tax_rate"], 2.0)
    assert np.isclose(result["asset_price"], -3.0)
    assert np.isclose(result["equal_rebate"], 5.0)


def test_shapley_splits_interaction_symmetrically() -> None:
    cells = {
        cell: float(cell[0] * cell[1])
        for cell in itertools.product((0, 1), repeat=3)
    }
    result = diagnosis.shapley_decomposition(cells)
    assert np.isclose(result["tax_rate"], 0.5)
    assert np.isclose(result["asset_price"], 0.5)
    assert np.isclose(result["equal_rebate"], 0.0)


def test_policy_endpoints_accept_exact_smoke_case_names(tmp_path: Path) -> None:
    path = tmp_path / "impact_results.csv"
    fields = (
        "policy_case",
        "calendar_year",
        "asset_price",
        "equal_transfer_period_units",
        "annual_property_tax_rate",
        "topcode_adjusted_births_per_adult",
        "owner_rate",
        "dependent_child_owner_rate",
        "housing_demand_per_adult",
    )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for index, name in enumerate(("tax1-equal-rebate", "tax2-equal-rebate")):
            writer.writerow(
                {
                    "policy_case": name,
                    "calendar_year": 2023,
                    "asset_price": 1.0 - 0.1 * index,
                    "equal_transfer_period_units": 0.2 + 0.1 * index,
                    "annual_property_tax_rate": 0.01 + 0.01 * index,
                    "topcode_adjusted_births_per_adult": 0.08 + 0.001 * index,
                    "owner_rate": 0.6 - 0.01 * index,
                    "dependent_child_owner_rate": 0.55 - 0.01 * index,
                    "housing_demand_per_adult": 6.0 - 0.1 * index,
                }
            )
    endpoints = diagnosis.read_policy_endpoints(path)
    assert endpoints["rebated-tax1-baseline"]["annual_tax_rate"] == 0.01
    assert endpoints["rebated-tax2-reform"]["annual_tax_rate"] == 0.02


if __name__ == "__main__":
    test_shapley_recovers_additive_components()
    test_shapley_splits_interaction_symmetrically()
    print("REBATED_TAX_SHAPLEY_TESTS_PASS tests=2_direct_plus_1_pytest_fixture")
