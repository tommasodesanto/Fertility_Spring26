#!/usr/bin/env python3
"""Focused tests for the default-off deterministic-tenure diagnostic."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

import run_e5f_person_policy_tenure_kink_diagnostic as diagnostic


class TenureOptionReconstructionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.grid = np.array([-1.0, 0.0, 1.0, 2.0])
        self.Vd = np.zeros((4, 2, 1, 2, 2))
        self.Vd[:, 0, 0, :, :] = self.grid[:, None, None]
        self.Vd[:, 1, 0, :, :] = 10.0 + 2.0 * self.grid[:, None, None]
        self.heq = np.array([[0.0, 0.8]])
        self.hcost = np.array([[0.0, 1.0]])
        self.dp = np.zeros((1, 2, 2, 2))
        self.dp[:, 1, :, :] = 0.4
        self.bmo = np.zeros((1, 2, 2, 2))
        self.bmo[:, 1, :, :] = -0.6
        self.birth_dp = np.zeros((2, 2, 2, 2), dtype=bool)
        self.grant = np.zeros((1, 2, 2, 2))

    def reconstruct(self) -> np.ndarray:
        return diagnostic.reconstruct_tenure_option_values(
            self.Vd,
            self.grid,
            self.heq,
            self.hcost,
            self.dp,
            self.bmo,
            self.birth_dp,
            self.grant,
        )

    def test_renter_purchase_constraints_and_interpolation(self) -> None:
        values = self.reconstruct()
        self.assertEqual(values[1, 0, 0, 0, 0, 1], diagnostic.NEGATIVE_INFINITY)
        self.assertAlmostEqual(values[2, 0, 0, 0, 0, 1], 10.0)

    def test_owner_sale_and_stay(self) -> None:
        values = self.reconstruct()
        self.assertAlmostEqual(values[1, 1, 0, 0, 0, 0], 0.8)
        self.assertAlmostEqual(values[1, 1, 0, 0, 0, 1], 10.0)

    def test_birth_waiver_uses_collateral_floor(self) -> None:
        self.birth_dp[1, 1, 0, 1] = True
        values = self.reconstruct()
        expected = 10.0 + 2.0 * (-0.6)
        self.assertAlmostEqual(values[0, 0, 0, 1, 1, 1], expected)

    def test_entry_grant_applies_both_feasibility_checks(self) -> None:
        self.grant[0, 1, 1, 1] = 0.75
        values = self.reconstruct()
        self.assertAlmostEqual(values[1, 0, 0, 1, 1, 1], 9.5)
        self.grant[0, 1, 1, 1] = 0.2
        values = self.reconstruct()
        self.assertEqual(values[1, 0, 0, 1, 1, 1], diagnostic.NEGATIVE_INFINITY)

    def test_rebuy_constraint_and_value(self) -> None:
        values = self.reconstruct()
        self.assertAlmostEqual(values[2, 1, 0, 0, 0, 1], 12.0)

    def test_argmax_matches_exact_kernel(self) -> None:
        values = self.reconstruct()
        _, exact = diagnostic.calendar.model.tenure_choice_kernel(
            self.Vd,
            self.grid,
            self.heq,
            self.hcost,
            self.dp,
            self.bmo,
            self.birth_dp,
            self.grant,
        )
        np.testing.assert_array_equal(np.argmax(values, axis=-1), exact)

    def test_bellman_slice_order_descends_in_age(self) -> None:
        self.assertEqual(diagnostic.bellman_slice_indices(0, 17, 15), (16, 0))
        self.assertEqual(diagnostic.bellman_slice_indices(14, 17, 15), (16, 14))
        self.assertEqual(diagnostic.bellman_slice_indices(15, 17, 15), (15, 0))
        self.assertEqual(diagnostic.bellman_slice_indices(254, 17, 15), (0, 14))

    def test_compressed_csv_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "states.csv.gz"
            diagnostic.write_csv(path, [{"state": 3, "value": 1.25}])
            with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(rows, [{"state": "3", "value": "1.25"}])


if __name__ == "__main__":
    unittest.main()
