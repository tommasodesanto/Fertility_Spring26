"""Synthetic coupled roots and contract checks; no household/model solves."""
from __future__ import annotations

import time
import unittest
from unittest import mock

import numpy as np

import e5f_social_security_root as root


class SocialSecurityRootTest(unittest.TestCase):
    def arguments(self, evaluate, *, closure="fixed_tax", **changes):
        result = dict(closure=closure, initial_prices=np.array([1.0]),
            initial_fiscal_values=np.array([.2]), evaluate=evaluate,
            project_prices=lambda values: values, fiscal_bounds=(.02, .8),
            market_tolerance=2e-4, fiscal_tolerance=1e-8,
            market_slope=1.4, fiscal_slope=.7, max_log_step=.2, damping=1.,
            max_evaluations=24, deadline_monotonic=time.monotonic()+30,
            max_condition_number=1e8, worsening_factor=3.,
            final_reproduction_tolerance=2e-10)
        result.update(changes)
        return result

    @staticmethod
    def reply(market, fiscal, valid=True, payload=None):
        return dict(market_residual=np.asarray(market), fiscal_residual=np.asarray(fiscal),
                    mapping_valid=valid, payload=payload)

    def coupled_evaluator(self, target_prices, target_fiscal, closure, calls):
        target_prices, target_fiscal = np.asarray(target_prices), np.asarray(target_fiscal)
        def evaluate(prices, fiscal):
            calls.append((prices.copy(), fiscal.copy()))
            x, y = np.log(prices / target_prices), np.log(fiscal / target_fiscal)
            # Both blocks depend on both unknowns. Payroll revenue rises in tax;
            # pension outlays rise in pension, so the fiscal signs differ.
            market = -1.4*x + .2*y
            budget = .1*x + (-.7 if closure == "fixed_tax" else .7)*y
            return self.reply(market, budget, payload={"trial": len(calls)})
        return evaluate

    def test_coupled_two_date_equilibrium_under_each_explicit_closure(self):
        for closure in ("fixed_tax", "fixed_pension"):
            with self.subTest(closure=closure):
                calls, progress = [], []
                target_prices, target_fiscal = np.array([1.2, .9]), np.array([.3, .15])
                args = self.arguments(self.coupled_evaluator(target_prices, target_fiscal, closure, calls),
                    closure=closure, initial_prices=np.ones(2), initial_fiscal_values=np.full(2, .2),
                    callback=progress.append)
                before = (args["initial_prices"].copy(), args["initial_fiscal_values"].copy())
                result = root.solve_social_security_path(**args)
                self.assertTrue(result["converged"], result["status"])
                self.assertTrue(all(result["gates"].values()))
                self.assertEqual(len(calls), result["evaluations"])
                np.testing.assert_array_equal(calls[-1][0], calls[-2][0])
                np.testing.assert_array_equal(calls[-1][1], calls[-2][1])
                np.testing.assert_allclose(result["final"]["prices"], target_prices, atol=1e-6)
                np.testing.assert_allclose(result["final"]["fiscal_values"], target_fiscal, atol=1e-6)
                np.testing.assert_array_equal(args["initial_prices"], before[0])
                np.testing.assert_array_equal(args["initial_fiscal_values"], before[1])
                self.assertEqual(result["contract"]["fiscal_residual_scale"], 20000.)
                self.assertEqual(progress[-1]["event"], "complete")
                self.assertTrue(any(row.get("new_best") for row in progress))
                self.assertTrue(all(len(row["prices"]) == 2 for row in progress if "prices" in row))
                self.assertLessEqual(result["final"]["score"], 1.)

    def test_block_scaled_fallback_retains_physical_step_under_tight_fiscal_gate(self):
        for closure in ("fixed_tax", "fixed_pension"):
            calls = []
            def evaluate(p, fiscal):
                calls.append((p.copy(), fiscal.copy()))
                sign = 1 if closure == "fixed_tax" else -1
                return self.reply([0.], [sign * .7 * np.log(.25/fiscal[0])])
            result = root.solve_social_security_path(**self.arguments(evaluate, closure=closure,
                max_evaluations=3, initial_jacobian=np.zeros((2,2))))
            # A singular supplied Jacobian triggers the stable scaled fallback.
            # The 20000 residual multiplier must not create a huge physical step.
            self.assertAlmostEqual(calls[1][1][0], .2*np.exp(.2))
            self.assertEqual(result["evaluations"], 3)
            self.assertIn("reset_reason", result["history"][1])

    def test_initial_physical_jacobian_round_trip_and_exact_joint_step(self):
        closure, calls = "fixed_pension", []
        jacobian = np.array([[-1.4, .2], [.1, .7]])
        original = jacobian.copy()
        evaluate = self.coupled_evaluator([1.05], [.21], closure, calls)
        result = root.solve_social_security_path(**self.arguments(evaluate, closure=closure,
            initial_jacobian=jacobian, max_evaluations=3))
        self.assertTrue(result["converged"])
        self.assertEqual(result["evaluations"], 3)
        np.testing.assert_allclose(result["final_jacobian"], original, atol=1e-12)
        np.testing.assert_array_equal(jacobian, original)

    def test_separate_fiscal_gate_cannot_be_hidden_by_housing_convergence(self):
        calls = mock.Mock(return_value=self.reply([0.], [1e-5]))
        result = root.solve_social_security_path(**self.arguments(calls, max_evaluations=5))
        self.assertFalse(result["converged"])
        self.assertTrue(result["gates"]["housing"])
        self.assertFalse(result["gates"]["social_security"])
        self.assertLessEqual(calls.call_count, 5)

    def test_boundary_projection_cannot_make_nonzero_tax_residual_pass(self):
        fiscal_values = []
        def evaluate(p, fiscal):
            fiscal_values.append(float(fiscal[0]))
            return self.reply([0.], [fiscal[0] - .9])
        result = root.solve_social_security_path(**self.arguments(evaluate, closure="fixed_pension",
            initial_fiscal_values=np.array([.799]), max_evaluations=8))
        self.assertFalse(result["converged"])
        self.assertFalse(result["gates"]["social_security"])
        self.assertTrue(all(.02 <= value <= .8 for value in fiscal_values))
        self.assertIn(result["status"], ("projection_stalled_without_market_gate", "evaluation_budget"))

    def test_discontinuous_fiscal_residual_never_passes(self):
        def evaluate(p, fiscal):
            return self.reply([0.], [1e-4 if fiscal[0] < .2512345 else -1e-4])
        result = root.solve_social_security_path(**self.arguments(evaluate, max_evaluations=18))
        self.assertFalse(result["converged"])
        self.assertEqual(result["evaluations"], 18)
        self.assertFalse(result["gates"]["social_security"])

    def test_success_at_initial_point_still_requires_uncached_replay(self):
        evaluate = mock.Mock(return_value=self.reply([0.], [0.]))
        result = root.solve_social_security_path(**self.arguments(evaluate, max_evaluations=2))
        self.assertTrue(result["converged"])
        self.assertEqual(evaluate.call_count, 2)
        self.assertEqual(result["history"][-1]["phase"], "final")

    def test_fiscal_only_replay_failure_inside_fiscal_gate_is_rejected(self):
        evaluate = mock.Mock(side_effect=[self.reply([0.], [0.]), self.reply([0.], [1e-9])])
        result = root.solve_social_security_path(**self.arguments(evaluate, max_evaluations=2))
        self.assertFalse(result["converged"])
        self.assertTrue(result["gates"]["housing"])
        self.assertTrue(result["gates"]["social_security"])
        self.assertFalse(result["gates"]["fiscal_replay"])

    def test_replay_has_original_units_even_when_fiscal_scaling_is_below_one(self):
        evaluate = mock.Mock(side_effect=[self.reply([0.], [0.]), self.reply([0.], [1e-8])])
        result = root.solve_social_security_path(**self.arguments(evaluate,
            fiscal_tolerance=.02, max_evaluations=2))
        self.assertFalse(result["converged"])
        self.assertAlmostEqual(result["fiscal_reproduction_max_abs"], 1e-8)

    def test_completion_callback_uses_original_market_replay_gate(self):
        progress = []
        evaluate = mock.Mock(side_effect=[self.reply([0.], [0.]), self.reply([1e-9], [0.])])
        result = root.solve_social_security_path(**self.arguments(evaluate,
            callback=progress.append, max_evaluations=2))
        self.assertFalse(result["converged"])
        self.assertFalse(result["gates"]["market_replay"])
        self.assertEqual(progress[-1]["event"], "complete")
        self.assertFalse(progress[-1]["converged"])

    def test_each_replay_block_can_pass_its_own_original_tolerance(self):
        evaluate = mock.Mock(side_effect=[self.reply([0.], [0.]), self.reply([1e-10], [1e-10])])
        result = root.solve_social_security_path(**self.arguments(evaluate, max_evaluations=2))
        self.assertTrue(result["converged"])
        self.assertTrue(all(result["gates"].values()))

    def test_explicit_domain_rejection_is_counted_and_reported_without_retrying_initial(self):
        evaluate = mock.Mock(side_effect=root.CandidateDomainError("no finite demographic state"))
        result = root.solve_social_security_path(**self.arguments(evaluate))
        self.assertFalse(result["converged"])
        self.assertEqual(evaluate.call_count, 1)
        self.assertEqual(result["evaluations"], 1)
        self.assertIn("no finite", result["history"][0]["candidate_domain_rejection"])
        self.assertEqual(len(result["domain_rejections"]), 1)

    def test_unexpected_economic_exception_is_not_silently_retried(self):
        evaluate = mock.Mock(side_effect=RuntimeError("broken policy source"))
        with self.assertRaisesRegex(RuntimeError, "broken policy"):
            root.solve_social_security_path(**self.arguments(evaluate))
        self.assertEqual(evaluate.call_count, 1)

    def test_invalid_tax_and_closure_fail_before_evaluation(self):
        cases = [dict(closure=""), dict(initial_fiscal_values=np.array([0.])),
            dict(initial_fiscal_values=np.array([1.])), dict(fiscal_bounds=(0., .8)),
            dict(fiscal_bounds=(.02, 1.)), dict(fiscal_bounds=(.02, float("inf"))),
            dict(initial_fiscal_values=np.array([.9]))]
        for changes in cases:
            evaluate = mock.Mock(return_value=self.reply([0.], [0.]))
            args = self.arguments(evaluate, closure="fixed_pension")
            args.update(changes)
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                root.solve_social_security_path(**args)
            evaluate.assert_not_called()

    def test_contract_and_residual_shape_failures(self):
        for reply in (self.reply([0., 0.], [0.]), self.reply([0.], [0.], valid=1), []):
            with self.subTest(reply=reply), self.assertRaises(ValueError):
                root.solve_social_security_path(**self.arguments(lambda *_: reply))
        for changes in (dict(fiscal_tolerance=0.), dict(market_slope=float("nan")),
                        dict(initial_jacobian=np.eye(3)), dict(max_evaluations=1)):
            evaluate = mock.Mock(return_value=self.reply([0.], [0.]))
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                root.solve_social_security_path(**self.arguments(evaluate, **changes))
            evaluate.assert_not_called()

    def test_deadline_never_starts_an_evaluation(self):
        evaluate = mock.Mock(return_value=self.reply([0.], [0.]))
        with self.assertRaises(ValueError):
            root.solve_social_security_path(**self.arguments(evaluate, deadline_monotonic=time.monotonic()-1))
        evaluate.assert_not_called()

    def test_callback_gets_copies_of_physical_latest_records(self):
        def callback(record):
            if "prices" in record:
                record["prices"][:] = -99
                record["fiscal_values"][:] = -99
        result = root.solve_social_security_path(**self.arguments(
            lambda *_: self.reply([0.], [0.]), callback=callback))
        self.assertTrue(result["converged"])
        np.testing.assert_array_equal(result["final"]["prices"], [1.])
        np.testing.assert_array_equal(result["final"]["fiscal_values"], [.2])


if __name__ == "__main__":
    unittest.main()
