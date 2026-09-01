#!/usr/bin/env python3
"""Focused tests for the positive-kappa H128 mechanism diagnostics."""

from __future__ import annotations

import csv
import gzip
import math
import tempfile
import unittest
from pathlib import Path

import collect_e5f_post2023_kappa_mechanism_decomposition as collector
import run_e5f_person_policy_fertility_mechanism_diagnostic as diagnostic


def _write_states(path, masses, probabilities) -> None:
    fields = [
        *collector.STATE_KEY,
        "age_group",
        "pre_fertility_mass",
        "attempt_probability",
        "fecundity_probability",
        "realized_birth_probability",
        "birth_child_units_if_realized",
        "expected_birth_child_units",
        "birth_value_gap_attempt_minus_wait",
    ]
    rows = []
    for index, (mass, probability) in enumerate(zip(masses, probabilities)):
        rows.append(
            {
                "calendar_year": 2023,
                "wealth_index": index,
                "origin_tenure": index,
                "market": 0,
                "age_index": index,
                "income_state": 0,
                "parity": index,
                "child_state": 0,
                "age_group": "25_34" if index == 0 else "35_44",
                "pre_fertility_mass": mass,
                "attempt_probability": probability,
                "fecundity_probability": 1.0,
                "realized_birth_probability": probability,
                "birth_child_units_if_realized": 1.0,
                "expected_birth_child_units": probability,
                "birth_value_gap_attempt_minus_wait": 0.01,
            }
        )
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class MechanismDecompositionTests(unittest.TestCase):
    def test_binary_value_gap_inverts_binary_logit(self) -> None:
        scale = 1.7
        gap = 0.23
        probability = 1.0 / (1.0 + math.exp(-gap / scale))
        self.assertTrue(
            math.isclose(
                diagnostic.binary_value_gap(probability, scale),
                gap,
                rel_tol=0.0,
                abs_tol=1.0e-14,
            )
        )
        self.assertEqual(diagnostic.binary_value_gap(0.0, scale), -math.inf)
        self.assertEqual(diagnostic.binary_value_gap(1.0, scale), math.inf)

    def test_plain_language_state_groups(self) -> None:
        self.assertEqual(diagnostic.age_group(24.0), "under_25")
        self.assertEqual(diagnostic.age_group(30.0), "25_34")
        self.assertEqual(diagnostic.age_group(40.0), "35_44")
        self.assertEqual(diagnostic.age_group(50.0), "45_64")
        self.assertEqual(diagnostic.age_group(70.0), "65_plus")
        self.assertEqual(diagnostic.parent_status(0, 0), "childless")
        self.assertEqual(diagnostic.parent_status(2, 1), "dependent_children")
        self.assertEqual(diagnostic.parent_status(2, 0), "empty_nest_parent")

    def test_two_factor_fertility_decomposition_adds_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            baseline = tmp_path / "baseline.csv.gz"
            reform = tmp_path / "reform.csv.gz"
            _write_states(baseline, (0.6, 0.4), (0.1, 0.2))
            _write_states(reform, (0.5, 0.5), (0.2, 0.1))
            summary, groups = collector.decompose_fertility_year(
                baseline,
                reform,
                baseline_heads=1.0,
                reform_heads=1.0,
            )
        expected_baseline = (0.6 * 0.1 + 0.4 * 0.2) / 4.0
        expected_reform = (0.5 * 0.2 + 0.5 * 0.1) / 4.0
        self.assertTrue(
            math.isclose(summary["baseline_annual_births_per_head"], expected_baseline)
        )
        self.assertTrue(
            math.isclose(summary["reform_annual_births_per_head"], expected_reform)
        )
        self.assertTrue(
            math.isclose(
                summary["within_state_policy_effect"]
                + summary["distribution_composition_effect"],
                expected_reform - expected_baseline,
                rel_tol=0.0,
                abs_tol=1.0e-15,
            )
        )
        self.assertEqual(
            {row["dimension"] for row in groups},
            {"age_group", "origin_tenure", "parity", "wealth_grid_quartile"},
        )

    def test_state_probability_audit_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            tmp_path = Path(temporary)
            baseline = tmp_path / "baseline.csv.gz"
            reform = tmp_path / "reform.csv.gz"
            _write_states(baseline, (1.0,), (0.2,))
            _write_states(reform, (1.0,), (1.2,))
            with self.assertRaisesRegex(RuntimeError, "Invalid reform"):
                collector.decompose_fertility_year(
                    baseline,
                    reform,
                    baseline_heads=1.0,
                    reform_heads=1.0,
                )


if __name__ == "__main__":
    unittest.main()
