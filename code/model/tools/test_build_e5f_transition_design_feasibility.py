#!/usr/bin/env python3
"""Focused tests for the transition-design demographic accounting."""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

import build_e5f_transition_design_feasibility as design


def test_aging_operator_moves_survivors_one_cell() -> None:
    operator = design.aging_operator()
    assert operator.shape == (17, 17)
    assert np.allclose(operator[0], 0.0)
    assert np.allclose(operator[1:, :-1], np.diag(np.diag(operator[1:, :-1])))
    assert np.all((operator >= 0.0) & (operator <= 1.0))


def test_fit_group_profile_recovers_known_law() -> None:
    operator = design.aging_operator()
    true_profile = np.zeros(17)
    expected = {}
    for value, (label, indices) in zip(
        (0.01, 0.02, 0.003, -0.0002, -0.0004),
        design.FIVE_GROUPS,
        strict=True,
    ):
        true_profile[list(indices)] = value
        expected[label] = value
    targets = [np.full(17, 1.0 / 17.0)]
    for _ in range(4):
        current = targets[-1]
        targets.append(operator @ current + true_profile * np.sum(current))
    fitted, rows = design.fit_group_profile(targets, operator, design.FIVE_GROUPS)
    assert np.allclose(fitted, true_profile, atol=1e-14, rtol=0.0)
    recovered = {
        row["age_group"]: row["net_rate_per_age_cell_per_four_years"]
        for row in rows
    }
    assert set(recovered) == set(expected)
    for label in expected:
        assert np.isclose(recovered[label], expected[label], atol=1e-14, rtol=0.0)


def test_no_bridge_uses_only_inherited_youngest_entry() -> None:
    operator = design.aging_operator()
    initial = np.arange(1.0, 18.0)
    inherited = 0.123
    states = design.simulate_no_bridge(initial, operator, inherited)
    assert len(states) == len(design.YEARS)
    for current, following in zip(states[:-1], states[1:], strict=True):
        expected = operator @ current
        expected[0] = inherited
        assert np.allclose(following, expected)


def test_output_writer_refuses_empty_csv() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "empty.csv"
        try:
            design.write_csv(path, [])
        except ValueError:
            pass
        else:
            raise AssertionError("Empty CSV was not rejected")


if __name__ == "__main__":
    test_aging_operator_moves_survivors_one_cell()
    test_fit_group_profile_recovers_known_law()
    test_no_bridge_uses_only_inherited_youngest_entry()
    test_output_writer_refuses_empty_csv()
    print("TRANSITION_DESIGN_FEASIBILITY_TESTS_PASS tests=4")
