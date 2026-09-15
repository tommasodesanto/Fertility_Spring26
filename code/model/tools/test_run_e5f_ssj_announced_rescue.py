"""Pure tests for the announced-rescue Jacobian and warm-start selection; no native model."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "cluster"))
import run_e5f_ssj_announced_rescue as rescue  # noqa: E402


def _receipt(T):
    lags = list(range(-5, 5))
    profiles = {}
    for i in ("housing_relative_imbalance", "paygo_relative_imbalance_scaled_200", "rebate_relative_imbalance_scaled_200"):
        for j in ("log_house_price", "log_pension", "log_rebate"):
            profiles[f"{i}<-{j}"] = [0.0] * 10
    profiles["housing_relative_imbalance<-log_house_price"][5] = -1.9
    profiles["paygo_relative_imbalance_scaled_200<-log_pension"][5] = -200.0
    profiles["rebate_relative_imbalance_scaled_200<-log_rebate"][5] = -200.0
    return dict(horizon=T, measured_lags=lags, lag_profiles=profiles)


class TestRescueSelection(unittest.TestCase):
    def test_toeplitz_mode_builds_full_horizon(self):
        J, info = rescue.select_jacobian("toeplitz", {}, _receipt(10), 104)
        self.assertEqual(J.shape, (312, 312))
        self.assertAlmostEqual(J[0, 0], -1.9)
        self.assertEqual(info["source_horizon"], 10)

    def test_broyden_mode_uses_receipt_matrix(self):
        Jm = -np.eye(6)
        J, info = rescue.select_jacobian("broyden_final", dict(final_jacobian=Jm.tolist()), None, 2)
        np.testing.assert_array_equal(J, Jm)
        with self.assertRaises(ValueError):
            rescue.select_jacobian("broyden_final", dict(final_jacobian=np.eye(5).tolist()), None, 2)

    def test_warm_coordinates_require_exact_replay(self):
        good = dict(best=dict(prices=[1.0] * 6, score=0.1), final_reproduction_max_abs=0.0)
        x, score = rescue.warm_coordinates(good, 2)
        self.assertEqual(x.shape, (3, 2)); self.assertEqual(score, 0.1)
        with self.assertRaisesRegex(ValueError, "exactly reproduced"):
            rescue.warm_coordinates(dict(good, final_reproduction_max_abs=1e-9), 2)


if __name__ == "__main__":
    unittest.main()


class TestLoosenedGatePieces(unittest.TestCase):
    def test_block_tolerance_vector(self):
        v = rescue.block_tolerance_vector(4, 2e-4, 0.5)
        self.assertEqual(v.shape, (12,)); self.assertEqual(v[0], 2e-4); self.assertEqual(v[4], 0.5); self.assertEqual(v[11], 0.5)
        with self.assertRaises(ValueError):
            rescue.block_tolerance_vector(4, 2e-4, 0.0)

    def test_warm_from_checkpoint_bypasses_replay_requirement(self):
        ck = dict(prices=[1.0] * 6, score=0.4)
        x, score = rescue.warm_coordinates(dict(best=None), 2, checkpoint=ck)
        self.assertEqual(x.shape, (3, 2)); self.assertEqual(score, 0.4)
