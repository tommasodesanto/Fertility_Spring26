"""Synthetic tests for the two-block four-shock acceleration helper."""
from __future__ import annotations

import time
import unittest

import numpy as np

import e5f_four_shock_acceleration as acceleration
import e5f_social_security_root as social_root


class FourShockAccelerationTest(unittest.TestCase):
    def arguments(self, evaluate, **changes):
        values = dict(
            closure="fixed_tax", initial_prices=np.array([1.0]),
            initial_fiscal_values=np.array([0.2]), evaluate=evaluate,
            project_prices=lambda price: price, fiscal_bounds=(0.02, 0.8),
            market_tolerance=2e-4, fiscal_tolerance=1e-8,
            market_slope=1.4, fiscal_slope=0.7, max_log_step=0.2,
            damping=1.0, max_evaluations=6, deadline_monotonic=time.monotonic() + 30,
            max_condition_number=1e9, worsening_factor=3.0,
            final_reproduction_tolerance=2e-10,
        )
        values.update(changes)
        return values

    def test_assemble_preserves_two_block_order_signs_and_zero_filled_lags(self):
        horizon, date = 4, 1
        # First unknown: housing rows are negative; pension rows are positive.
        price_column = np.array([-4., -3., -2., -1., 10., 20., 30., 40.])
        pension_column = np.array([5., 6., 7., 8., -50., -60., -70., -80.])
        jacobian, receipt = acceleration.assemble_jacobian([price_column, pension_column], horizon, date)
        self.assertEqual(jacobian.shape, (8, 8))
        self.assertEqual(receipt["unknown_blocks"], ["log_house_price", "log_period_pension"])
        self.assertEqual(receipt["residual_blocks"], ["housing_imbalance", "pension_imbalance"])
        self.assertFalse(receipt["fake_news_derivatives_constructed"])
        self.assertEqual(receipt["measured_lags"], [-1, 0, 1, 2])
        self.assertEqual(jacobian[1, 1], -3.)       # housing at lag zero
        self.assertEqual(jacobian[2, 1], -2.)       # housing at lag +1
        self.assertEqual(jacobian[1, 5], 6.)        # housing <- pension, lag zero
        self.assertEqual(jacobian[5, 1], 20.)       # pension <- price, lag zero
        self.assertEqual(jacobian[5, 5], -60.)      # pension <- pension, lag zero
        self.assertEqual(jacobian[0, 3], 0.)        # lag -3 was not measured
        self.assertGreater(receipt["zero_filled_entries_per_block"], 0)

    def test_wrapper_scales_physical_jacobian_once_and_uses_uniform_direction(self):
        calls = []
        target_log_step = np.array([0.6, 0.3])
        physical_jacobian = np.array([[-1.4, 0.2], [0.1, -0.7]])

        def evaluate(prices, pension):
            x = np.array([np.log(prices[0]), np.log(pension[0] / 0.2)])
            calls.append(x.copy())
            residual = physical_jacobian @ (x - target_log_step)
            return dict(market_residual=residual[:1], fiscal_residual=residual[1:], mapping_valid=True)

        result = acceleration.solve_joint_with_acceleration(**self.arguments(
            evaluate, initial_jacobian=physical_jacobian, max_evaluations=6))
        self.assertTrue(result["converged"], result["status"])
        np.testing.assert_allclose(result["final_jacobian"], physical_jacobian, rtol=0, atol=1e-12)
        # The first Newton direction [0.6, 0.3] is common-factor scaled to [.2, .1].
        np.testing.assert_allclose(calls[1], [0.2, 0.1], rtol=0, atol=1e-12)
        self.assertTrue(all(result["gates"].values()))

    def test_final_replay_and_separate_physical_gates_are_preserved(self):
        replies = iter([(0.0, 0.0), (0.0, 1e-9)])

        def evaluate(_prices, _pension):
            market, fiscal = next(replies)
            return dict(market_residual=np.array([market]), fiscal_residual=np.array([fiscal]), mapping_valid=True)

        result = acceleration.solve_joint_with_acceleration(**self.arguments(evaluate, max_evaluations=2))
        self.assertFalse(result["converged"])
        self.assertTrue(result["gates"]["housing"])
        self.assertTrue(result["gates"]["social_security"])
        self.assertFalse(result["gates"]["fiscal_replay"])
        self.assertEqual(result["history"][-1]["phase"], "final")

    def test_restores_original_generic_root_after_success_and_exception(self):
        original = social_root.solve_price_path
        success = acceleration.solve_joint_with_acceleration(**self.arguments(
            lambda *_: dict(market_residual=np.array([0.]), fiscal_residual=np.array([0.]), mapping_valid=True),
            max_evaluations=2))
        self.assertTrue(success["converged"])
        self.assertIs(social_root.solve_price_path, original)

        with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
            acceleration.solve_joint_with_acceleration(**self.arguments(
                lambda *_: (_ for _ in ()).throw(RuntimeError("synthetic failure"))))
        self.assertIs(social_root.solve_price_path, original)

    def test_rejects_old_three_block_scaled_jacobian_without_calling_evaluator(self):
        called = []
        with self.assertRaisesRegex(ValueError, "three-block scaled-200"):
            acceleration.solve_joint_with_acceleration(**self.arguments(
                lambda *_: called.append(True), initial_jacobian=np.eye(3)))
        self.assertEqual(called, [])


if __name__ == "__main__":
    unittest.main()
