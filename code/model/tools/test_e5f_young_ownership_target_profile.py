from __future__ import annotations

import pytest

from intergen_eqscale_seq_optimized.e5_profile import e5_target_system
from intergen_eqscale_seq_optimized.e5_young_ownership_profile import (
    E5_YOUNG_OWNERSHIP_CALIBRATION_SCALE,
    E5_YOUNG_OWNERSHIP_TARGET,
    E5_YOUNG_OWNERSHIP_WEIGHT,
    e5_target_system_for_profile,
    e5_young_ownership_target_system,
)


def test_diagnostic_profile_appends_one_row_without_changing_e5() -> None:
    baseline = e5_target_system()
    diagnostic = e5_young_ownership_target_system()
    assert diagnostic.moment_names[:-1] == baseline.moment_names
    assert diagnostic.target_values[:-1] == baseline.target_values
    assert diagnostic.weights[:-1] == baseline.weights
    assert diagnostic.moment_names[-1] == "own_rate_2534"
    assert diagnostic.target_values[-1] == pytest.approx(E5_YOUNG_OWNERSHIP_TARGET)
    assert diagnostic.weights[-1] == pytest.approx(E5_YOUNG_OWNERSHIP_WEIGHT)
    assert E5_YOUNG_OWNERSHIP_CALIBRATION_SCALE == pytest.approx(
        0.05 * E5_YOUNG_OWNERSHIP_TARGET
    )


def test_target_selector_preserves_production_default() -> None:
    baseline = e5_target_system_for_profile("baseline")
    assert baseline.count == 12
    assert baseline.fingerprint == (
        "3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497"
    )
    diagnostic = e5_target_system_for_profile("young-ownership-overid-v1")
    assert diagnostic.count == 13
    assert diagnostic.fingerprint != baseline.fingerprint


def test_target_selector_fails_closed() -> None:
    with pytest.raises(ValueError, match="Unknown target profile"):
        e5_target_system_for_profile("typo")
