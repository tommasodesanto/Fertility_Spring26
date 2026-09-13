"""Pure interface tests; these do not call the native model."""
from types import SimpleNamespace
import unittest

import numpy as np

from e5f_sequence_space_prototype import (
    NativePathEvaluation,
    OriginalQueueSequenceSpaceBridge,
    dense_directional_jacobian,
    pack_unknowns,
    unpack_unknowns,
)


def _affine_evaluator(prices, pensions, rebates):
    # The native boundary receives levels for prices, so make residual affine in
    # log(price) again solely for this interface-only test.
    x = pack_unknowns(np.log(prices), pensions, rebates)
    state = SimpleNamespace(g_pre=np.ones((1,)), scheduled_entries=[2.] * 4,
                            scheduled_raw_entries=[3.] * 4)
    return NativePathEvaluation(residual=2.0 * x - 1.0, state=state,
                                receipt={"native_gates_preserved": True})


class TestSequenceSpacePrototype(unittest.TestCase):
    def test_zero_shock_and_population_contract_are_preserved(self):
        bridge = OriginalQueueSequenceSpaceBridge(_affine_evaluator)
        unknowns = pack_unknowns(np.full(2, 0.5), np.full(2, 0.5), np.full(2, 0.5))
        result = bridge.evaluate(unknowns)
        np.testing.assert_allclose(result.residual, 0.0)
        self.assertEqual(len(result.state.scheduled_entries), 4)
        self.assertEqual(len(result.state.scheduled_raw_entries), 4)

    def test_central_directional_derivative_matches_affine_native_interface(self):
        bridge = OriginalQueueSequenceSpaceBridge(_affine_evaluator)
        unknowns = pack_unknowns(np.full(2, 0.5), np.full(2, 0.5), np.full(2, 0.5))
        direction = np.arange(6.0).reshape(3, 2) / 10.0
        actual = dense_directional_jacobian(bridge, unknowns, direction, 1e-6)
        np.testing.assert_allclose(actual, 2.0 * direction, rtol=0.0, atol=1e-9)

    def test_rejects_lost_queue_state(self):
        def bad(prices, pensions, rebates):
            return NativePathEvaluation(np.zeros((3, 1)), SimpleNamespace(g_pre=np.ones(1)), {})
        with self.assertRaisesRegex(ValueError, "scheduled_entries"):
            OriginalQueueSequenceSpaceBridge(bad).evaluate(pack_unknowns(np.zeros(1), np.zeros(1), np.zeros(1)))

    def test_unpack_rejects_wrong_shape(self):
        with self.assertRaisesRegex(ValueError, "shape"):
            unpack_unknowns(np.zeros((2, 3)))
