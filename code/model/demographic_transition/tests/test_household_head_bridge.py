from __future__ import annotations

import unittest

import numpy as np

from demographic_transition import (
    CohortState,
    aggregate_heads_to_model_age_cells,
    rebalance_household_distribution_by_age,
)


class HouseholdHeadBridgeTest(unittest.TestCase):
    def test_aggregate_heads_to_four_year_cells(self) -> None:
        persons = np.ones((2, 12))
        heads = np.arange(24, dtype=float).reshape(2, 12) / 24.0
        state = CohortState(2025, persons * 2.0, heads)
        result = aggregate_heads_to_model_age_cells(
            state, age_start=2, cell_width=4, number_of_cells=2
        )
        np.testing.assert_allclose(
            result,
            [np.sum(heads[:, 2:6]), np.sum(heads[:, 6:10])],
        )

    def test_rebalance_hits_age_targets_and_preserves_cell_composition(self) -> None:
        distribution = np.zeros((2, 2, 1, 3, 1, 1, 1))
        distribution[:, :, 0, 0, 0, 0, 0] = [[1.0, 2.0], [3.0, 4.0]]
        distribution[:, :, 0, 1, 0, 0, 0] = [[2.0, 1.0], [1.0, 2.0]]
        distribution[:, :, 0, 2, 0, 0, 0] = [[1.0, 1.0], [1.0, 1.0]]
        targets = np.array([20.0, 3.0, 0.0])
        output, ledger = rebalance_household_distribution_by_age(
            distribution, targets
        )
        np.testing.assert_allclose(
            np.sum(output, axis=(0, 1, 2, 4, 5, 6)), targets
        )
        before_share = distribution[:, :, 0, 0, 0, 0, 0] / np.sum(
            distribution[:, :, 0, 0, 0, 0, 0]
        )
        after_share = output[:, :, 0, 0, 0, 0, 0] / np.sum(
            output[:, :, 0, 0, 0, 0, 0]
        )
        np.testing.assert_allclose(after_share, before_share)
        np.testing.assert_allclose(ledger.added_mass, [10.0, 0.0, 0.0])
        np.testing.assert_allclose(ledger.removed_mass, [0.0, 3.0, 4.0])
        self.assertLessEqual(abs(ledger.total_mass_residual), 1e-12)

    def test_empty_age_cell_requires_and_uses_template(self) -> None:
        distribution = np.zeros((2, 1, 1, 2, 1, 1, 1))
        distribution[0, 0, 0, 0, 0, 0, 0] = 1.0
        with self.assertRaisesRegex(RuntimeError, "zero source mass"):
            rebalance_household_distribution_by_age(
                distribution, np.array([1.0, 2.0])
            )
        template = np.zeros((2, 1, 1, 1, 1, 1))
        template[:, 0, 0, 0, 0, 0] = [1.0, 3.0]
        output, ledger = rebalance_household_distribution_by_age(
            distribution,
            np.array([1.0, 2.0]),
            empty_age_templates={1: template},
        )
        np.testing.assert_allclose(output[:, 0, 0, 1, 0, 0, 0], [0.5, 1.5])
        self.assertTrue(np.isinf(ledger.scale_factors[1]))


if __name__ == "__main__":
    unittest.main()
