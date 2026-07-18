from __future__ import annotations

import math

import numpy as np

from intergen_seq_fertility.kernels import (
    full_owner_block_kernel,
    full_renter_block_kernel,
)
from intergen_seq_fertility.parameters import apply_overrides, setup_parameters
from intergen_seq_fertility.solver import precompute_shared


DEAD_VALUE_CUTOFF = -1e9
GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
DEBT_GRID = np.array([-1.0, -0.1, 0.0, 0.5, 1.0])
CBAR = 0.4
HBAR = 0.5
RI = 0.2
ALPHA = 0.7
OMS = -1.0
BETA = 0.9
GS_ALPHA1 = (3.0 - math.sqrt(5.0)) / 2.0
GS_ALPHA2 = (math.sqrt(5.0) - 1.0) / 2.0


def _python_kernel(function):
    return getattr(function, "py_func", function)


def _renter_call(
    resources: np.ndarray,
    guarantee: float,
    *,
    test_resources: np.ndarray | None = None,
    grid: np.ndarray = GRID,
):
    kernel = _python_kernel(full_renter_block_kernel)
    resources = np.asarray(resources, dtype=float)
    if test_resources is None:
        test_resources = resources
    return kernel(
        resources,
        np.asarray(test_resources, dtype=float),
        np.zeros((grid.size, 1)),
        np.zeros((grid.size, 1)),
        0,
        grid,
        np.array([CBAR]),
        np.array([HBAR]),
        np.array([0.0]),
        np.array([guarantee]),
        RI,
        3.0,
        0.04,
        CBAR,
        HBAR,
        ALPHA,
        OMS,
        BETA,
        0.0,
        0.0,
        GS_ALPHA1,
        GS_ALPHA2,
        1e-8,
    )


def _owner_call(
    resources: np.ndarray,
    guarantee: float,
    *,
    test_resources: np.ndarray | None = None,
    grid: np.ndarray = GRID,
    collateral_floor: float = 0.0,
):
    kernel = _python_kernel(full_owner_block_kernel)
    resources = np.asarray(resources, dtype=float)
    if test_resources is None:
        test_resources = resources
    return kernel(
        resources,
        np.asarray(test_resources, dtype=float),
        np.zeros((grid.size, 1)),
        np.zeros((grid.size, 1)),
        0,
        grid,
        np.array([CBAR]),
        np.array([HBAR]),
        np.array([0.0]),
        np.array([guarantee]),
        np.array([collateral_floor]),
        0.1,
        2.0,
        1.0,
        1.0,
        0.04,
        ALPHA,
        OMS,
        BETA,
        0.0,
        0.0,
        GS_ALPHA1,
        GS_ALPHA2,
        1e-8,
    )


def test_precompute_shared_builds_transfer_guarantees_in_f_order() -> None:
    params = apply_overrides(
        setup_parameters(),
        {"transfer_floor_G0": 0.52, "transfer_floor_Gn": 0.40},
    )
    shared = precompute_shared(params, GRID)

    assert shared.g_bar[2, 1] == 0.52 + 0.40 * 2
    assert shared.g_bar[2, 0] == 0.52
    np.testing.assert_array_equal(
        shared.gb_flat,
        shared.g_bar.reshape(1, shared.nc, order="F"),
    )


def test_renter_kernel_is_identical_off_and_matches_no_floor_formula() -> None:
    resources = np.array([0.25, 0.75, 1.25, 1.75, 2.25])
    first = _renter_call(resources, 0.0)
    duplicate = _renter_call(resources, 0.0)
    for actual, repeated in zip(first, duplicate):
        assert np.array_equal(actual, repeated)

    values, bp, consumption, housing = first
    cell = 2
    dc = CBAR + RI * HBAR
    surplus = resources[cell] - dc - bp[cell, 0]
    cap = RI * (3.0 - HBAR) / (1.0 - ALPHA)
    assert 0.0 < surplus < cap
    kr = (ALPHA**ALPHA * ((1.0 - ALPHA) / RI) ** (1.0 - ALPHA)) ** OMS
    expected_value = kr * surplus**OMS / OMS
    expected_consumption = CBAR + max(ALPHA * surplus, 0.04)
    expected_housing = HBAR + max((1.0 - ALPHA) * surplus / RI, 0.01)
    np.testing.assert_allclose(values[cell, 0], expected_value, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(consumption[cell, 0], expected_consumption, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(housing[cell, 0], expected_housing, rtol=0.0, atol=1e-12)


def test_income_poor_negative_asset_node_is_revived_but_cap_still_binds() -> None:
    income = 0.35
    resources = DEBT_GRID + income
    test_resources = np.maximum(DEBT_GRID, 0.0) + income
    values_off, _, _, _ = _renter_call(
        resources,
        0.0,
        test_resources=test_resources,
        grid=DEBT_GRID,
    )
    values_on, _, consumption_on, _ = _renter_call(
        resources,
        0.7,
        test_resources=test_resources,
        grid=DEBT_GRID,
    )

    revived = 1
    assert DEBT_GRID[revived] < 0.0
    assert test_resources[revived] < 0.7
    assert 0.7 > CBAR + RI * HBAR
    assert values_off[revived, 0] <= DEAD_VALUE_CUTOFF
    assert values_on[revived, 0] > DEAD_VALUE_CUTOFF
    assert consumption_on[revived, 0] >= CBAR

    capped = 0
    assert resources[capped] + 0.7 < CBAR + RI * HBAR
    assert values_on[capped, 0] <= DEAD_VALUE_CUTOFF


def test_mortgage_leverage_does_not_trigger_transfer() -> None:
    income = 1.0
    resources = DEBT_GRID + income
    test_resources = np.maximum(DEBT_GRID, 0.0) + income
    without_floor = _owner_call(
        resources,
        0.0,
        test_resources=test_resources,
        grid=DEBT_GRID,
        collateral_floor=-1.0,
    )
    with_floor = _owner_call(
        resources,
        0.7,
        test_resources=test_resources,
        grid=DEBT_GRID,
        collateral_floor=-1.0,
    )

    assert resources[0] < 0.1 + CBAR
    assert test_resources[0] > 0.7
    for actual, no_floor in zip(with_floor, without_floor):
        assert np.array_equal(actual, no_floor)


def test_owner_kernel_is_identical_off_and_floor_revives_node() -> None:
    resources = np.array([0.25, 0.75, 1.25, 1.75, 2.25])
    first = _owner_call(resources, 0.0)
    duplicate = _owner_call(resources, 0.0)
    for actual, repeated in zip(first, duplicate):
        assert np.array_equal(actual, repeated)

    values, bp, consumption = first
    cell = 2
    ct = resources[cell] - 0.1 - CBAR - bp[cell, 0]
    ko = (2.0 - HBAR) ** ((1.0 - ALPHA) * OMS)
    expected_value = ko * ct ** (ALPHA * OMS) / OMS
    np.testing.assert_allclose(values[cell, 0], expected_value, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(consumption[cell, 0], CBAR + max(ct, 0.04), rtol=0.0, atol=1e-12)

    values_on, _, consumption_on = _owner_call(resources, 0.7)
    assert values[0, 0] <= DEAD_VALUE_CUTOFF
    assert values_on[0, 0] > DEAD_VALUE_CUTOFF
    assert consumption_on[0, 0] >= CBAR
