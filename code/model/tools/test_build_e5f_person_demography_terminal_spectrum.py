"""Tests for the terminal person-law spectral diagnostic."""

from __future__ import annotations

import unittest

import numpy as np

import build_e5f_person_demography_terminal_spectrum as spectrum


class TerminalSpectrumTests(unittest.TestCase):
    def test_four_year_matrix_matches_direct_block_iteration(self) -> None:
        survival = np.array(
            [
                [0.98, 0.97, 0.95, 0.90, 0.70],
                [0.98, 0.98, 0.96, 0.92, 0.75],
            ]
        )
        annual = spectrum.build_existing_cohort_operator(survival)
        newborn = np.zeros_like(survival)
        newborn[:, 0] = np.array([0.51, 0.49]) * survival[:, 0]
        headship = np.array(
            [
                [0.0, 0.1, 0.4, 0.6, 0.5],
                [0.0, 0.1, 0.35, 0.55, 0.45],
            ]
        )
        operator = spectrum.build_four_year_operator(
            annual_existing_operator=annual,
            newborn_survivor_vector=newborn,
            headship_rates=headship,
            annual_births_per_head=0.02,
        )
        gap = spectrum.validate_block_operator(
            four_year_operator=operator,
            annual_existing_operator=annual,
            newborn_survivor_vector=newborn,
            headship_rates=headship,
            annual_births_per_head=0.02,
        )
        self.assertLess(gap, 1e-12)

    def test_minimum_periods_crosses_requested_threshold(self) -> None:
        radius = 0.96
        periods = spectrum.minimum_periods(radius, 0.01)
        self.assertLessEqual(radius**periods, 0.01)
        self.assertGreater(radius ** (periods - 1), 0.01)


if __name__ == "__main__":
    unittest.main()
