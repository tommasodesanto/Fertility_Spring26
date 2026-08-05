from __future__ import annotations

import unittest

from intergen_eqscale_seq_optimized.e5_profile import E5_TARGETS
from intergen_eqscale_seq_optimized.e6_plain_weighting import (
    MOMENT_BLOCKS,
    target_relative_block_equal_l1_weights,
    target_relative_block_equal_weights,
)


class E6PlainWeightingTests(unittest.TestCase):
    def test_each_block_has_equal_aggregate_weight_at_unit_relative_gap(self) -> None:
        targets = {name: float(value) for name, value in E5_TARGETS.items()}
        weights = target_relative_block_equal_weights(targets)
        for names in MOMENT_BLOCKS.values():
            contribution = sum(
                weights[name] * (0.10 * targets[name]) ** 2 for name in names
            )
            self.assertAlmostEqual(contribution, 0.01)

    def test_weights_follow_target_order(self) -> None:
        targets = {name: float(value) for name, value in E5_TARGETS.items()}
        weights = target_relative_block_equal_weights(targets)
        self.assertEqual(tuple(weights), tuple(targets))

    def test_l1_each_block_is_mean_absolute_proportional_gap(self) -> None:
        targets = {name: float(value) for name, value in E5_TARGETS.items()}
        weights = target_relative_block_equal_l1_weights(targets)
        for names in MOMENT_BLOCKS.values():
            contribution = sum(
                weights[name] * abs(0.10 * targets[name]) for name in names
            )
            self.assertAlmostEqual(contribution, 0.10)

    def test_l1_weights_follow_target_order(self) -> None:
        targets = {name: float(value) for name, value in E5_TARGETS.items()}
        weights = target_relative_block_equal_l1_weights(targets)
        self.assertEqual(tuple(weights), tuple(targets))

    def test_contract_rejects_missing_row(self) -> None:
        targets = {name: float(value) for name, value in E5_TARGETS.items()}
        del targets["childless_rate"]
        with self.assertRaisesRegex(ValueError, "exact twelve-row"):
            target_relative_block_equal_weights(targets)
        with self.assertRaisesRegex(ValueError, "exact twelve-row"):
            target_relative_block_equal_l1_weights(targets)


if __name__ == "__main__":
    unittest.main()
