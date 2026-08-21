#!/usr/bin/env python3
"""Focused tests for the rebated-tax Shapley decomposition."""

from __future__ import annotations

import itertools

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


if __name__ == "__main__":
    test_shapley_recovers_additive_components()
    test_shapley_splits_interaction_symmetrically()
    print("REBATED_TAX_SHAPLEY_TESTS_PASS tests=2")
