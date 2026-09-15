"""Toy tests: scaled-step variant versus componentwise clipping on a coupled linear system."""
import time
import unittest

import numpy as np

from e5f_ssj_scaled_step_root import solve_price_path_scaled
from e5f_ssj_toeplitz_jacobian import toeplitz_block


def _coupled_system(T=10):
    # Near-differencing housing block like the measured one: -1.9 on the diagonal, +1.0 one date ahead.
    J = toeplitz_block(np.array([0, -1]), np.array([-1.9, 1.0]), T)
    root = np.log(np.linspace(0.7, 0.4, T))
    def evaluate(prices):
        x = np.log(prices)
        return dict(residual=J @ (x - root), mapping_valid=True)
    return J, root, evaluate


def _controls(J, budget):
    return dict(slope=1.63, market_tolerance=2e-4, max_log_step=0.2, damping=1.0,
                max_evaluations=budget, deadline_monotonic=time.monotonic() + 30,
                max_condition_number=1e10, worsening_factor=1.5, final_reproduction_tolerance=0.0,
                project=lambda x: x, initial_jacobian=J)


class TestScaledStep(unittest.TestCase):
    def test_scaled_step_preserves_direction_and_converges(self):
        J, root, evaluate = _coupled_system()
        start = np.exp(root + np.linspace(0.0, -0.6, 10))  # far below the root at late dates
        receipt = solve_price_path_scaled(initial_prices=start, evaluate=evaluate, **_controls(J, 8))
        self.assertTrue(receipt["converged"], receipt["status"])
        scores = [h["score"] for h in receipt["history"] if h["phase"] != "final"]
        self.assertTrue(all(b < a for a, b in zip(scores[:-1], scores[1:])), scores)
        first = receipt["history"][1]
        act = np.asarray(first["actual_log_step"])
        self.assertAlmostEqual(float(np.max(np.abs(act))), 0.2, places=12)
        # Exact Newton direction from the start is (root - log start); the scaled step is a
        # common multiple of it, so the max coordinate is 0.2 and the shape is preserved.
        newton = root - np.log(start)
        np.testing.assert_allclose(act, newton * (0.2 / np.max(np.abs(newton))), rtol=0, atol=1e-12)

    def test_clipped_reference_converges_on_linear_toy_but_distorts_direction(self):
        from e5f_ssj_scaled_step_root import np as _np  # noqa: F401  (same numpy)
        import importlib.util, pathlib
        spec = importlib.util.spec_from_file_location(
            "ref_root", pathlib.Path(__file__).resolve().parents[3] / "tmp/e5f_matched_pf/code/model/tools/e5f_matched_pf_path_root.py")
        ref = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref)
        J, root, evaluate = _coupled_system()
        start = np.exp(root + np.linspace(0.0, -0.6, 10))
        receipt = ref.solve_price_path(initial_prices=start, evaluate=evaluate, **_controls(J, 8))
        # On a linear toy both rules converge; the distinction is the direction.
        self.assertTrue(receipt["converged"])
        act = np.asarray(receipt["history"][1]["actual_log_step"])
        newton = root - np.log(start)
        # Componentwise clipping flattens the late-date coordinates to 0.2, distorting the direction.
        self.assertGreater(float(np.max(np.abs(act - newton * (0.2 / np.max(np.abs(newton)))))), 0.05)
        self.assertEqual(int(np.sum(np.isclose(np.abs(act), 0.2))), 7)


if __name__ == "__main__":
    unittest.main()


class TestToleranceVector(unittest.TestCase):
    def test_blockwise_gate_certifies_when_each_block_meets_its_tolerance(self):
        J, root, evaluate = _coupled_system()
        start = np.exp(root + np.linspace(0.0, -0.6, 10))
        tol = np.full(10, 1e-3)
        receipt = solve_price_path_scaled(initial_prices=start, evaluate=evaluate, tolerance_vector=tol,
                                          **_controls(J, 8))
        self.assertTrue(receipt["converged"])
        self.assertEqual(receipt["gate"], 1.0)
        final = receipt["final"]
        self.assertLessEqual(final["score"], 1.0)
        self.assertAlmostEqual(final["score"], final["raw_max_abs"] / 1e-3, places=12)
        self.assertTrue(all("raw_max_abs" in h for h in receipt["history"]))

    def test_rejects_wrong_shape_tolerance(self):
        J, root, evaluate = _coupled_system()
        with self.assertRaisesRegex(ValueError, "tolerance_vector"):
            solve_price_path_scaled(initial_prices=np.exp(root), evaluate=evaluate,
                                    tolerance_vector=np.ones(3), **_controls(J, 4))

    def test_without_tolerance_vector_behaviour_is_unchanged(self):
        J, root, evaluate = _coupled_system()
        start = np.exp(root + np.linspace(0.0, -0.6, 10))
        receipt = solve_price_path_scaled(initial_prices=start, evaluate=evaluate, **_controls(J, 8))
        self.assertIsNone(receipt["tolerance_vector"]); self.assertEqual(receipt["gate"], 2e-4)
        self.assertEqual(receipt["final"]["score"], receipt["final"]["raw_max_abs"])
