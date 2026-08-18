"""Tests for the default-off E5F child-room-floor experiment."""

from __future__ import annotations

import numpy as np
import pytest

from intergen_eqscale_seq_optimized.e5f_floor_profile import (
    E5F_DOMAIN,
    E5F_FIXED,
    e5f_metadata,
    e5f_overrides,
)
from intergen_eqscale_seq_optimized.calibration import extract_moments
from intergen_eqscale_seq_optimized.e5_profile import E5_TARGETS
from intergen_eqscale_seq_optimized.kernels import eval_owner_block_kernel, eval_renter_scalar
from intergen_eqscale_seq_optimized.parameters import apply_overrides, setup_parameters
from intergen_eqscale_seq_optimized.solver import precompute_shared, run_model_cp_dt
from intergen_eqscale_seq_optimized.tests.test_literal_parity import L4, _seq_eqscale


def _repaired_tiny(**extra: object) -> dict[str, object]:
    overrides = {
        **_seq_eqscale(**L4),
        "child_state_mode": "independent_count",
        "A_m": 18.0,
        "use_stochastic_aging": True,
        "alpha_cons": 0.733,
    }
    overrides.update(extra)
    return overrides


def _parameters(**extra: object):
    P = setup_parameters()
    overrides = {
        **L4,
        "preference_spec": "eqscale",
        "eqscale_form": "power",
        "child_state_mode": "independent_count",
        "A_m": 18.0,
        "use_stochastic_aging": True,
        "alpha_cons": 0.733,
    }
    overrides.update(extra)
    return apply_overrides(P, overrides)


def test_e5f_contract_is_nine_parameters_against_same_twelve_targets() -> None:
    assert len(E5F_DOMAIN) == 9
    assert {name for name, *_ in E5F_DOMAIN}.isdisjoint({"delta_alpha", "delta_alpha_jump"})
    assert E5F_DOMAIN[-1] == ("hbar_child_rooms", 0.10, 1.80, "log")
    assert E5F_FIXED["alpha_cons"] == 0.733
    assert E5F_FIXED["delta_alpha"] == E5F_FIXED["delta_alpha_jump"] == 0.0
    assert e5f_metadata()["target_count"] == 12


def test_switch_off_and_zero_floor_are_bitwise_nested() -> None:
    base = _repaired_tiny(delta_alpha=0.05, delta_alpha_jump=0.03)
    repaired, repaired_parameters, _ = run_model_cp_dt(base, verbose=False)
    switch_off, switch_off_parameters, _ = run_model_cp_dt(
        {**base, "child_room_floor": False, "hbar_child_rooms": 1.2}, verbose=False
    )
    zero_floor, zero_floor_parameters, _ = run_model_cp_dt(
        {**base, "child_room_floor": True, "hbar_child_rooms": 0.0}, verbose=False
    )
    for candidate, candidate_parameters in (
        (switch_off, switch_off_parameters),
        (zero_floor, zero_floor_parameters),
    ):
        assert np.array_equal(repaired.V, candidate.V)
        assert np.array_equal(repaired.g, candidate.g)
        assert np.array_equal(repaired.c_pol, candidate.c_pol)
        assert np.array_equal(repaired.hR_pol, candidate.hR_pol)
        assert np.array_equal(repaired.bp_pol, candidate.bp_pol)
        assert np.array_equal(repaired.tenure_choice, candidate.tenure_choice)
        repaired_moments = extract_moments(repaired, repaired_parameters)
        candidate_moments = extract_moments(candidate, candidate_parameters)
        for name in E5_TARGETS:
            assert repaired_moments[name] == candidate_moments[name] or (
                np.isnan(repaired_moments[name]) and np.isnan(candidate_moments[name])
            )


def test_floor_is_keyed_to_children_currently_at_home() -> None:
    P = _parameters(
        **{
            **e5f_overrides(),
            "hbar_child_rooms": 0.6,
            "delta_alpha": 0.0,
            "delta_alpha_jump": 0.0,
        }
    )
    shared = precompute_shared(P, np.linspace(-1.0, 2.0, 5))
    np.testing.assert_allclose(
        shared.h_bar[3], np.array([0.0, 0.6, 1.2, 1.8]), rtol=0.0, atol=1e-15
    )
    assert shared.h_bar[1, 0] == shared.h_bar[3, 0] == 0.0
    assert shared.h_bar[1, 1] == shared.h_bar[3, 1] == 0.6
    assert np.all(shared.c_bar == 0.0)
    assert np.all(shared.alpha_flat == 0.733)


def test_first_child_jump_separates_level_from_additional_child_slope() -> None:
    P = _parameters(
        **{
            **e5f_overrides(),
            "hbar_child_rooms": 0.3,
            "hbar_first_child_jump": 0.5,
            "delta_alpha": 0.0,
            "delta_alpha_jump": 0.0,
        }
    )
    shared = precompute_shared(P, np.linspace(-1.0, 2.0, 5))
    np.testing.assert_allclose(
        shared.h_bar[3], np.array([0.0, 0.8, 1.1, 1.4]), rtol=0.0, atol=1e-15
    )
    # Once parenthood begins, each additional child changes the floor by the
    # original slope while the first-child jump is paid only once.
    np.testing.assert_allclose(np.diff(shared.h_bar[3, 1:]), 0.3, rtol=0.0, atol=1e-15)


def test_renter_utility_respects_floor_and_diverges_at_boundary() -> None:
    alpha, sigma, rent = 0.733, 2.0, 1.0
    oms = 1.0 - sigma
    floor = 1.2
    Kr = (alpha**alpha * ((1.0 - alpha) / rent) ** (1.0 - alpha)) ** oms
    grid = np.array([0.0, 1.0])
    continuation = np.zeros(2)

    def value(residual_expenditure: float) -> float:
        return float(
            eval_renter_scalar(
                0.0, floor + residual_expenditure, continuation, grid,
                floor, 0.0, 100.0, 0.0, rent, 6.0, 4.8,
                Kr, alpha, oms, 0.0, 1.0,
            )
        )

    assert value(1e-4) < value(1e-2) < value(1.0)
    assert value(0.0) == -1e10


def test_owner_rungs_below_floor_are_dead_and_receive_no_mass() -> None:
    floor = 1.1
    owner_values, owner_consumption = eval_owner_block_kernel(
        np.array([3.0, 4.0]),
        np.zeros((2, 1)),
        np.zeros((2, 1)),
        np.zeros((2, 1), dtype=np.int64),
        np.zeros((2, 1)),
        np.zeros(1),
        np.array([2.2]),
        np.zeros(1),
        0.0,
        2.0,
        1.0,
        1.0,
        0.01,
        0.733,
        -1.0,
        0.0,
        1,
    )
    assert np.all(owner_values == -1e10)
    assert np.all(np.isfinite(owner_consumption))
    solution, _, _ = run_model_cp_dt(
        _repaired_tiny(
            **{
                **e5f_overrides(),
                "hbar_child_rooms": floor,
                "delta_alpha": 0.0,
                "delta_alpha_jump": 0.0,
            }
        ),
        verbose=False,
    )
    assert np.all(np.isfinite(solution.V))
    assert np.all(np.isfinite(solution.g))
    # The first owner rung supplies two service units.  It is infeasible for
    # states with two or three children at home because m*1.1 >= 2.2.
    assert float(np.sum(solution.g[:, 1, ..., 2])) == 0.0
    assert float(np.sum(solution.g[:, 1, ..., 3])) == 0.0
    for m in (1, 2, 3):
        renter_mass = solution.g[:, 0, ..., m]
        realized_h = solution.hR_pol[:, 0, ..., m]
        if float(np.sum(renter_mass)) > 0.0:
            assert np.all(realized_h[renter_mass > 0.0] > m * floor)


def test_floor_cap_guard_accepts_domain_upper_bound_and_refuses_collision() -> None:
    accepted = _parameters(child_room_floor=True, hbar_child_rooms=1.8)
    assert accepted.hbar_child_rooms == 1.8
    with pytest.raises(ValueError, match="rental-service cap"):
        _parameters(child_room_floor=True, hbar_child_rooms=2.0)
    with pytest.raises(ValueError, match="rental-service cap"):
        _parameters(
            child_room_floor=True,
            hbar_child_rooms=1.7,
            hbar_first_child_jump=1.0,
        )
