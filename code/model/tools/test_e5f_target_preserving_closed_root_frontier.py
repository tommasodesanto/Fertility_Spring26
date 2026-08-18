#!/usr/bin/env python3
"""Focused algebra tests for the target-preserving slope profiler."""

from __future__ import annotations

import run_e5f_target_preserving_closed_root_frontier as profile


def _trial(delta: float, gap: float) -> dict:
    target = 1.918
    return {
        "row": {
            "preference_change": delta,
            "target_gap": gap,
            "terminal_completed_fertility": target + gap,
        }
    }


def test_factor_contract() -> None:
    assert profile.factor_for_index(1) == 0.25
    assert profile.factor_for_index(4) == 1.0
    assert profile.factor_for_index(7) == 2.0


def test_bracket_and_secant() -> None:
    trials = [_trial(-0.5, -0.2), _trial(-0.2, 0.1)]
    bracket = profile.target_bracket(trials)
    assert bracket is not None
    delta = profile.safeguarded_secant(*bracket)
    assert abs(delta + 0.3) < 1e-12


def test_monotonicity_gate() -> None:
    assert profile.monotone_trials([_trial(-0.5, -0.1), _trial(-0.2, 0.1)])
    assert not profile.monotone_trials([_trial(-0.5, 0.1), _trial(-0.2, -0.1)])


def test_choose_best() -> None:
    best = profile.choose_best([_trial(-0.5, -0.1), _trial(-0.2, 0.01)])
    assert best["row"]["preference_change"] == -0.2


if __name__ == "__main__":
    test_factor_contract()
    test_bracket_and_secant()
    test_monotonicity_gate()
    test_choose_best()
    print("TARGET_PRESERVING_PROFILE_TESTS_PASS tests=4")
