#!/usr/bin/env python3
"""Focused tests for the tenure-complementarity share-map collector."""

from __future__ import annotations

import unittest

import collect_e5f_person_policy_tenure_complementarity_share_map as collector


class ShareIntervalTests(unittest.TestCase):
    def test_constant_residual_interval(self) -> None:
        self.assertEqual(collector.allowed_interval(1e-5, 1e-5, 2e-5), (0.0, 1.0))
        self.assertIsNone(collector.allowed_interval(3e-5, 3e-5, 2e-5))

    def test_crossing_residual_interval(self) -> None:
        lower, upper = collector.allowed_interval(-2.0, 2.0, 1.0)
        self.assertAlmostEqual(lower, 0.25)
        self.assertAlmostEqual(upper, 0.75)

    def test_signed_market_residual_matches_stored_absolute_value(self) -> None:
        row = {
            "housing_supply": "10.0",
            "housing_demand": "10.002",
            "relative_market_residual": "0.0002",
        }
        self.assertAlmostEqual(collector.signed_market_residual(row), -0.0002)


if __name__ == "__main__":
    unittest.main()
