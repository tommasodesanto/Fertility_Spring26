from __future__ import annotations

from types import SimpleNamespace
import unittest
import numpy as np

import run_e5f_joint_rebated_initial_probe as j


class JointRootTests(unittest.TestCase):
    def test_three_equation_root(self):
        target = np.array([0.2, -0.3, 0.4])
        matrix = np.array([[1.0, .2, 0.0], [0.1, 1.0, .1], [0.0, .2, 1.0]])
        result = j.solve_three_residual_root(
            lambda x: dict(residual=matrix @ (x-target), payload=x.copy()),
            np.zeros(3), maximum_evaluations=20)
        self.assertEqual(result["status"], "converged")
        self.assertTrue(j.residual_passes(result["result"]["residual"]))
        self.assertLessEqual(result["evaluations"], 20)

    def test_nonfinite_residual_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "nonfinite"):
            j.solve_three_residual_root(lambda _: dict(residual=[0, np.nan, 0]),
                                        np.zeros(3), maximum_evaluations=4)

    def test_rank_deficient_root_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "rank deficient"):
            j.solve_three_residual_root(lambda _: dict(residual=np.ones(3)),
                                        np.zeros(3), maximum_evaluations=4)


class AdapterTests(unittest.TestCase):
    def driver(self, written=None):
        written = [] if written is None else written
        return SimpleNamespace(
            bind_initial_balanced_pension=lambda *a, **k: None,
            certify_initial_pension=lambda *a, **k: None,
            solve_balanced_initial_equilibrium=lambda **k: None,
            calibration=SimpleNamespace(solve_old_steady_state=lambda *a, **k: None),
            primitive=SimpleNamespace(pf=SimpleNamespace(calendar=SimpleNamespace(
                write_json_atomic=lambda path, value: written.append((path, value))))))

    def test_normalizer_calls_joint_closure_once(self):
        calls = []
        class Chain:
            def run_model_cp_dt(self, overrides):
                calls.append(overrides)
                return SimpleNamespace(tfr=2.1), SimpleNamespace(psi_child=.12), np.array([.7])
            def extract_moments(self, solution, _): return {"tfr": solution.tfr}
        driver = self.driver()
        j.install_on_initial_driver(driver)
        answer = driver.calibration.solve_old_steady_state(
            Chain(), {}, initial_psi=.1, completed_fertility_target=2.1,
            completed_fertility_tolerance=5e-4, normalize=True)
        self.assertEqual(len(calls), 1)
        self.assertEqual(answer[-1]["status"], "derived_intercept")
        self.assertEqual(answer[-1]["normalization_method"],
                         "joint_price_fertility_rebate_root")
        self.assertEqual(answer[1].psi_child, .12)

    def test_relaxed_contract_rejected_before_call(self):
        driver = self.driver()
        j.install_on_initial_driver(driver)
        with self.assertRaisesRegex(ValueError, "unchanged"):
            driver.calibration.solve_old_steady_state(
                object(), {}, initial_psi=.1, completed_fertility_target=2.0,
                completed_fertility_tolerance=5e-4, normalize=True)

    def test_receipt_preserves_requested_and_records_solved_psi(self):
        written = []; driver = self.driver(written)
        j.install_on_initial_driver(driver)
        row = dict(psi_child=.1, fiscal=dict(joint_initial_root=dict(psi_child=.12)))
        driver.primitive.pf.calendar.write_json_atomic("stationary_solves.json", [row])
        saved = written[0][1][0]
        self.assertEqual(saved["requested_psi_child"], .1)
        self.assertEqual(saved["solved_psi_child"], .12)
        self.assertEqual(saved["psi_child"], .12)

    def test_tighter_housing_tolerance_is_respected(self):
        residual = np.array([2e-5, 0.0, 0.0])
        self.assertTrue(j.residual_passes(residual))
        self.assertFalse(j.residual_passes(residual, housing_tolerance=1e-5))


if __name__ == "__main__": unittest.main()
