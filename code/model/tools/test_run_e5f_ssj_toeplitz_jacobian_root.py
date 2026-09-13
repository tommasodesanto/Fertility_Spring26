"""Loop-structure smoke for the derivative stage with a toy mapping; no native model."""
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "cluster"))
import run_e5f_ssj_toeplitz_jacobian_root as driver  # noqa: E402


def _toy(label, u):
    # Linear stationary toy: housing = -2 log q_t - 0.2 log q_{t-1}; paygo = -200 log b_t;
    # rebate = -200 log r_t - 70 log q_t.  Base coordinates are zero so baseline residual is 0.
    T = u.shape[1]
    h = -2.0 * u[0] - 0.2 * np.concatenate([[0.0], u[0][:-1]])
    p = -200.0 * u[1]
    r = -200.0 * u[2] - 70.0 * u[0]
    reply = dict(residual=np.concatenate([h, p, r]), gates=dict(ok=0.0), seconds=0.01)
    if label == "baseline":
        reply["stationary_drift"] = dict(distribution_relative_l1=0.0, population_relative_gap=0.0,
                                         queue_relative_max=0.0, raw_queue_relative_max=0.0)
    return reply


class TestDerivativeStage(unittest.TestCase):
    def test_seven_mappings_and_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            J, receipt = driver.run_derivative_stage(
                _toy, tmp, horizon=10, perturbed_date=5, step=1e-5,
                deadline_monotonic=time.monotonic() + 60, base_coordinates=np.zeros((3, 10)), slope=1.63)
            self.assertEqual(receipt["mappings"], 7)
            self.assertAlmostEqual(J[5, 5], -2.0, places=6)
            self.assertAlmostEqual(J[6, 5], -0.2, places=6)
            self.assertAlmostEqual(J[25, 5], -70.0, places=5)
            self.assertAlmostEqual(J[15, 15], -200.0, places=5)
            self.assertAlmostEqual(J[25, 25], -200.0, places=5)
            self.assertEqual(J[0, 9], 0.0)
            for name in ("latest_completed.json", "derivative_receipt.json", "jacobian.json"):
                self.assertTrue((Path(tmp) / name).exists(), name)
            latest = json.loads((Path(tmp) / "latest_completed.json").read_text())
            self.assertEqual(latest["mappings"], 7)
            self.assertFalse(receipt["fake_news_derivatives_constructed"])

    def test_baseline_gate_failure_stops_before_perturbations(self):
        calls = []

        def bad(label, u):
            calls.append(label)
            reply = _toy(label, u)
            reply["residual"] = reply["residual"] + 1e-3
            return reply
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RuntimeError, "baseline"):
                driver.run_derivative_stage(bad, tmp, horizon=4, perturbed_date=2, step=1e-5,
                                            deadline_monotonic=time.monotonic() + 60,
                                            base_coordinates=np.zeros((3, 4)), slope=1.63)
            self.assertEqual(calls, ["baseline"])
            self.assertTrue((Path(tmp) / "derivative_receipt.json").exists())

    def test_compare_roots_table(self):
        hist = [dict(evaluation=1, phase="initial", residual=[1.0] * 3 + [2.0] * 3 + [3.0] * 3,
                     score=3.0, evaluation_seconds=10.0)]
        out = driver.compare_roots(dict(history=hist, status="x", converged=False, best=dict(score=3.0)),
                                   dict(history=hist, status="y", converged=False, best=dict(score=3.0)), 3)
        self.assertEqual(out["new_history"][0]["max_rebate"], 3.0)
        self.assertEqual(out["reference_history"][0]["max_housing"], 1.0)


if __name__ == "__main__":
    unittest.main()


class TestStepRuleSelection(unittest.TestCase):
    def test_select_step_rule(self):
        from e5f_ssj_scaled_step_root import solve_price_path_scaled
        self.assertIsNone(driver.select_step_rule("clipped"))
        self.assertIs(driver.select_step_rule("scaled"), solve_price_path_scaled)
        with self.assertRaises(ValueError):
            driver.select_step_rule("other")
