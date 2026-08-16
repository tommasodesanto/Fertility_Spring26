from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import run_e5f_transition_calibration as transition_calibration


SOURCE_THETA = {
    "beta": 0.9787070083905574,
    "psi_child": 0.0,
    "kappa_fert": 2.9868403225209876,
    "kappa_fert_continuation": 1.7133240302362214,
    "chi": 1.0369353042009588,
    "H0": 9.873032233988866,
    "theta0": 0.15935051857132962,
    "theta1": 0.09373031820629316,
    "hbar_child_rooms": 0.6136054836829904,
}


def test_transition_unit_round_trip() -> None:
    unit = transition_calibration.transition_unit_from_candidate(
        SOURCE_THETA, -0.2419
    )
    theta, terminal = transition_calibration.candidate_from_transition_unit(
        SOURCE_THETA, unit
    )
    for name in SOURCE_THETA:
        assert np.isclose(theta[name], SOURCE_THETA[name])
    assert np.isclose(terminal, -0.2419)


def test_latin_hypercube_is_deterministic_and_stratified() -> None:
    first = transition_calibration.latin_hypercube(8, 3, 17)
    second = transition_calibration.latin_hypercube(8, 3, 17)
    assert np.array_equal(first, second)
    assert np.all((first >= 0.0) & (first <= 1.0))
    for column in range(first.shape[1]):
        strata = np.floor(first[:, column] * 8).astype(int)
        assert sorted(strata.tolist()) == list(range(8))


def test_panel_anchor_is_the_declared_working_point() -> None:
    args = SimpleNamespace(
        panel_task_id=1,
        panel_size=64,
        panel_seed=2026081601,
        panel_design="mixed",
        panel_local_radius=0.18,
        panel_center_json=None,
    )
    theta, terminal, metadata = transition_calibration.panel_candidate(
        SOURCE_THETA, args
    )
    assert metadata["design"] == "anchor"
    assert metadata["task_id"] == 1
    assert np.isclose(terminal, -0.2419)
    for name in SOURCE_THETA:
        assert np.isclose(theta[name], SOURCE_THETA[name])


def test_coordinate_panel_moves_one_dimension_at_a_time() -> None:
    args = SimpleNamespace(
        panel_task_id=2,
        panel_size=19,
        panel_seed=2026081602,
        panel_design="coordinate",
        panel_local_radius=0.04,
        panel_center_json=None,
    )
    _, terminal, metadata = transition_calibration.panel_candidate(
        SOURCE_THETA, args
    )
    center = transition_calibration.transition_unit_from_candidate(
        SOURCE_THETA, -0.2419
    )
    candidate = np.asarray(metadata["unit_vector"])
    changed = np.flatnonzero(~np.isclose(candidate, center))
    assert changed.tolist() == [0]
    assert np.isclose(candidate[0], center[0] - 0.04)
    assert np.isclose(terminal, -0.2419)


def test_repaired_profile_adds_cost_and_measured_income() -> None:
    theta = dict(SOURCE_THETA)
    domain, overrides, metadata = transition_calibration.activate_model_profile(
        transition_calibration.REPAIRED_MODEL_PROFILE,
        theta,
    )
    assert metadata["name"] == transition_calibration.REPAIRED_MODEL_PROFILE
    assert theta["first_birth_fixed_cost"] == 0.0
    assert domain[-1][0] == "first_birth_fixed_cost"
    assert overrides["permanent_income_levels_enabled"] is True
    assert len(np.asarray(overrides["z_grid"])) == 15
