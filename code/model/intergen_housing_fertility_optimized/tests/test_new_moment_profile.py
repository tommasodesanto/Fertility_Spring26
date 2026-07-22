from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from intergen_housing_fertility_optimized.calibration_search import validate_contract
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_TARGETS,
    NEW_MOMENT_WEIGHTS,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.solver import (
    add_aggregate_wealth_bequest_flow_moments,
    markov_tenure_residual_variance,
)


class NewMomentProfileTests(unittest.TestCase):
    def test_contract_has_fourteen_identifying_rows_and_relative_weights(self) -> None:
        validate_contract("new-moments")
        system = new_moment_target_system()
        self.assertEqual(system.count, 14)
        self.assertEqual(set(NEW_MOMENT_TARGETS), set(NEW_MOMENT_WEIGHTS))
        for name, target in NEW_MOMENT_TARGETS.items():
            self.assertAlmostEqual(NEW_MOMENT_WEIGHTS[name], 1.0 / target**2)

    def test_aggregate_wealth_and_bequest_flow_units(self) -> None:
        p = SimpleNamespace(
            J=2,
            I=1,
            J_R=2,
            period_years=4.0,
            da=4.0,
            scale_flows_to_period=True,
            use_age_survival=True,
            survival_probs=np.array([0.5]),
            income=np.array([[4.0, 4.0]]),
            property_tax_lump_sum_transfer=0.0,
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
        )
        g = np.zeros((2, 2, 1, 2, 1, 1, 1))
        g[0, 0, 0, 0, 0, 0, 0] = 1.0
        g[1, 1, 0, 1, 0, 0, 0] = 1.0
        stats = SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(
            stats, g, p, np.array([-1.0, 3.0]), np.array([1.0])
        )
        self.assertAlmostEqual(stats.aggregate_wealth, 4.0)
        self.assertAlmostEqual(stats.aggregate_annual_after_tax_earnings, 2.0)
        self.assertAlmostEqual(stats.aggregate_wealth_to_annual_after_tax_earnings, 2.0)
        self.assertAlmostEqual(stats.annual_bequest_flow, 1.25)
        self.assertAlmostEqual(stats.annual_bequest_flow_to_aggregate_wealth, 0.3125)

    def test_tenure_residual_is_expected_bernoulli_variance(self) -> None:
        p = SimpleNamespace(I=1, J=4, age_start=18.0, da=4.0)
        g = np.zeros((1, 1, 1, 4, 1, 1, 1))
        g[0, 0, 0, 2, 0, 0, 0] = 2.0
        choice = np.zeros_like(g, dtype=int)
        probs = np.zeros((*g.shape, 2))
        probs[..., 0] = 1.0
        probs[0, 0, 0, 2, 0, 0, 0, :] = (0.75, 0.25)
        value = markov_tenure_residual_variance(g, choice, probs, p)
        self.assertAlmostEqual(value, 0.25 * 0.75)


if __name__ == "__main__":
    unittest.main()
