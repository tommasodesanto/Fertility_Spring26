#!/usr/bin/env python3
"""Focused algebra tests for the full-model demographic-closure smoke."""

from __future__ import annotations

import numpy as np

import run_e5f_historical_demographic_closure_smoke as smoke


def synthetic_distribution(age_masses: np.ndarray) -> np.ndarray:
    values = np.zeros((1, 1, 1, age_masses.size, 1, 1, 1), dtype=float)
    values[0, 0, 0, :, 0, 0, 0] = age_masses
    return values


def test_age_mass_target_preserves_within_age_composition() -> None:
    values = np.zeros((2, 1, 1, 3, 1, 1, 1), dtype=float)
    values[:, 0, 0, 0, 0, 0, 0] = (1.0, 3.0)
    values[:, 0, 0, 1, 0, 0, 0] = (2.0, 2.0)
    values[:, 0, 0, 2, 0, 0, 0] = (3.0, 1.0)
    target = np.array([2.0, 8.0, 1.0])
    out, audit = smoke.apply_age_mass_target(values, target)
    assert np.allclose(smoke.age_masses(out), target)
    for age in range(3):
        before = values[:, 0, 0, age, 0, 0, 0]
        after = out[:, 0, 0, age, 0, 0, 0]
        assert np.allclose(before / before.sum(), after / after.sum())
    assert abs(audit["net_flow"] - (target.sum() - values.sum())) < 1e-14


def test_age_mass_target_rejects_negative_mass() -> None:
    values = synthetic_distribution(np.ones(3))
    try:
        smoke.apply_age_mass_target(values, np.array([1.0, -0.1, 1.0]))
    except RuntimeError:
        pass
    else:
        raise AssertionError("Negative target age mass was not rejected")


def test_age_mass_target_rejects_empty_source_cell() -> None:
    values = synthetic_distribution(np.array([1.0, 0.0, 1.0]))
    try:
        smoke.apply_age_mass_target(values, np.ones(3))
    except RuntimeError:
        pass
    else:
        raise AssertionError("Empty source age cell was not rejected")


if __name__ == "__main__":
    test_age_mass_target_preserves_within_age_composition()
    test_age_mass_target_rejects_negative_mass()
    test_age_mass_target_rejects_empty_source_cell()
    print("HISTORICAL_DEMOGRAPHIC_CLOSURE_TESTS_PASS tests=3")
