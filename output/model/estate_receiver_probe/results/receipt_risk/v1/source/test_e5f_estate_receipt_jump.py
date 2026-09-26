"""Focused tests for the standalone estate receipt wealth-jump operator.

These tests are intentionally not run on the Mac.  Run on Torch with:
``python -m unittest code/model/tools/test_e5f_estate_receipt_jump.py``.
"""

import unittest

import numpy as np

from e5f_estate_receipt_jump import (
    OccupiedClippingError,
    adjointness_residual,
    backward_expectation,
    build_estate_receipt_jump_plan,
    forward_transport,
)


class EstateReceiptJumpTests(unittest.TestCase):
    def test_zero_receipt_is_exact_identity(self):
        plan = build_estate_receipt_jump_plan([0.0, 1.0, 2.0], [0.0], [1.0])
        values = np.arange(12, dtype=np.int64).reshape(3, 2, 2)
        mass = np.arange(12, dtype=np.int64).reshape(3, 2, 2)
        backward = backward_expectation(plan, values)
        forward, account = forward_transport(plan, mass)
        self.assertEqual(backward.dtype, values.dtype)
        self.assertEqual(forward.dtype, mass.dtype)
        np.testing.assert_array_equal(backward, values)
        np.testing.assert_array_equal(forward, mass)
        self.assertEqual(account["expected_receipt_flow"], 0.0)

    def test_hand_computable_lottery_and_mean_receipt_difference(self):
        grid = [0.0, 1.0, 2.0, 4.0]
        lottery = build_estate_receipt_jump_plan(grid, [0.0, 1.0], [0.5, 0.5])
        mean = build_estate_receipt_jump_plan(grid, [0.5], [1.0])
        mass = np.array([1.0, 1.0, 0.0, 0.0])
        after, account = forward_transport(lottery, mass)
        np.testing.assert_allclose(after, [0.5, 1.0, 0.5, 0.0])
        self.assertAlmostEqual(after.sum(), mass.sum())
        self.assertAlmostEqual(account["expected_receipt_flow"], 1.0)
        self.assertAlmostEqual(account["wealth_after"] - account["wealth_before"], 1.0)
        # A lottery spanning an interior grid knot differs from its mean.
        # A lottery confined to one interpolation interval cannot differ:
        # the interpolant is affine on that interval.
        spanning_lottery = build_estate_receipt_jump_plan(grid, [0.0, 2.0], [0.5, 0.5])
        same_mean = build_estate_receipt_jump_plan(grid, [1.0], [1.0])
        nonlinear_values = np.asarray(grid) ** 2
        self.assertNotAlmostEqual(
            backward_expectation(spanning_lottery, nonlinear_values)[0],
            backward_expectation(same_mean, nonlinear_values)[0],
        )
        mean_after, mean_account = forward_transport(mean, mass)
        self.assertAlmostEqual(account["expected_receipt_flow"], mean_account["expected_receipt_flow"])
        self.assertAlmostEqual(float(np.dot(grid, mean_after)), float(np.dot(grid, after)))

    def test_adjointness_with_extra_axes_and_nonuniform_grid(self):
        plan = build_estate_receipt_jump_plan([0.0, 0.4, 1.7, 3.0], [0.0, 0.2, 0.7], [0.2, 0.3, 0.5])
        values = np.arange(24.0).reshape(4, 2, 3) / 7.0
        mass = (np.arange(24.0).reshape(4, 2, 3) + 1.0) / 13.0
        self.assertLess(abs(adjointness_residual(plan, values, mass)), 1e-12)

    def test_probability_one_lottery_equals_deterministic_operator(self):
        deterministic = build_estate_receipt_jump_plan([0.0, 1.0, 2.0, 4.0], [0.75], [1.0])
        lottery = build_estate_receipt_jump_plan([0.0, 1.0, 2.0, 4.0], [0.0, 0.75], [0.0, 1.0])
        values = np.arange(8.0).reshape(4, 2)
        mass = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [0.0, 0.0]])
        np.testing.assert_allclose(backward_expectation(deterministic, values), backward_expectation(lottery, values))
        np.testing.assert_allclose(
            forward_transport(deterministic, mass)[0], forward_transport(lottery, mass)[0]
        )

    def test_invalid_inputs_fail(self):
        with self.assertRaises(ValueError):
            build_estate_receipt_jump_plan([0.0, 1.0], [-0.1], [1.0])
        with self.assertRaises(ValueError):
            build_estate_receipt_jump_plan([0.0, 1.0], [np.nan], [1.0])
        with self.assertRaises(ValueError):
            build_estate_receipt_jump_plan([0.0, 1.0], [0.2], [0.8])
        plan = build_estate_receipt_jump_plan([0.0, 1.0, 2.0], [0.5], [1.0])
        with self.assertRaises(ValueError):
            forward_transport(plan, [-1.0, 0.0, 0.0])
        with self.assertRaises(ValueError):
            backward_expectation(plan, np.ones((2, 3)))
        with self.assertRaises(ValueError):
            forward_transport(plan, np.ones((2, 3)))

    def test_occupied_overflow_rejects_and_unoccupied_does_not(self):
        plan = build_estate_receipt_jump_plan([0.0, 1.0, 2.0], [1.0], [1.0])
        with self.assertRaises(OccupiedClippingError) as caught:
            forward_transport(plan, [0.0, 0.0, 2.0])
        account = caught.exception.account
        self.assertEqual(account["clipped_wealth_loss"], 2.0)
        self.assertEqual(account["mass_encountering_clipping"], 2.0)
        after, unoccupied_account = forward_transport(plan, [1.0, 0.0, 0.0])
        np.testing.assert_allclose(after, [0.0, 1.0, 0.0])
        self.assertEqual(unoccupied_account["clipped_wealth_loss"], 0.0)
        after, allowed = forward_transport(plan, [0.0, 0.0, 2.0], max_clipped_wealth=2.0)
        self.assertEqual(allowed["wealth_before"] + allowed["expected_receipt_flow"] - allowed["wealth_after"], allowed["clipped_wealth_loss"])
        np.testing.assert_allclose(after, [0.0, 0.0, 2.0])


if __name__ == "__main__":
    unittest.main()
