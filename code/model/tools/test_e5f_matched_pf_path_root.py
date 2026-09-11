"""Synthetic numerical checks, with no household or equilibrium model calls."""
import time
import unittest
from unittest.mock import patch
import numpy as np
from e5f_matched_pf_path_root import solve_price_path


class PriceRootTests(unittest.TestCase):
    def run_root(self, evaluate, **kwargs):
        arguments = dict(initial_prices=np.ones(2), evaluate=evaluate,
            project=lambda p: p, slope=1., market_tolerance=1e-9,
            max_log_step=.25, damping=1., max_evaluations=30,
            deadline_monotonic=time.monotonic() + 10., max_condition_number=1e8,
            worsening_factor=1.25, final_reproduction_tolerance=1e-12)
        arguments.update(kwargs)
        return solve_price_path(**arguments)

    def test_linear_coupled_market_and_reserved_final_replay(self):
        A = np.array([[-2., .2], [.3, -1.5]])
        target = np.array([.15, -.2])
        evaluated = []
        def evaluate(prices):
            evaluated.append(prices.copy())
            return dict(residual=A @ (np.log(prices) - target), mapping_valid=True)
        result = self.run_root(evaluate)
        self.assertTrue(result['converged'], result['status'])
        np.testing.assert_allclose(np.log(result['final']['prices']), target, atol=1e-8)
        np.testing.assert_array_equal(evaluated[-1], evaluated[-2])
        self.assertEqual(result['history'][-1]['phase'], 'final')

    def test_supplied_coupled_jacobian_gives_exact_first_newton_step(self):
        A = np.array([[-2., .2], [.3, -1.5]])
        untouched = A.copy()
        target = np.array([.15, -.2])
        result = self.run_root(lambda p: dict(residual=A @ (np.log(p) - target), mapping_valid=True),
                               initial_jacobian=A)
        self.assertTrue(result['converged'])
        self.assertEqual(result['evaluations'], 3)
        np.testing.assert_allclose(result['history'][1]['actual_log_step'], target, atol=1e-15)
        np.testing.assert_array_equal(A, untouched)

    def test_supplied_jacobian_shape_and_finiteness_fail_before_evaluation(self):
        calls = []
        def evaluate(p):
            calls.append(p)
            return dict(residual=np.zeros(2), mapping_valid=True)
        for matrix in (np.ones(2), np.eye(3), np.ones((2, 1)),
                       np.array([[np.nan, 0.], [0., -1.]]),
                       np.array([[-1., 0.], [0., np.inf]])):
            with self.subTest(matrix=matrix):
                with self.assertRaisesRegex(ValueError, 'finite NxN'):
                    self.run_root(evaluate, initial_jacobian=matrix)
        self.assertFalse(calls)

    def test_bad_supplied_jacobian_resets_to_slope_not_supplied_matrix(self):
        result = self.run_root(lambda p: dict(residual=np.full(2, .1) - np.log(p), mapping_valid=True),
                               initial_jacobian=np.eye(2), max_evaluations=5)
        first, second = result['history'][1:3]
        np.testing.assert_allclose(first['actual_log_step'], [-.1, -.1])
        self.assertEqual(first['safeguard'], 'restored_best_halved_damping_reset_jacobian')
        np.testing.assert_allclose(second['actual_log_step'], [.05, .05])

    def test_singular_supplied_matrix_uses_positive_residual_fallback(self):
        result = self.run_root(lambda p: dict(residual=np.full(2, .1) - np.log(p), mapping_valid=True),
                               initial_jacobian=np.zeros((2, 2)))
        self.assertTrue(result['converged'])
        self.assertEqual(result['history'][1]['reset_reason'], 'jacobian_reset_positive_residual_price_increase')

    def test_explicit_reset_matrix_respects_different_equation_scales(self):
        matrix = np.diag([-2., -10000.])
        original = matrix.copy()
        target = np.array([.1, -.1])
        result = self.run_root(lambda p: dict(residual=matrix @ (np.log(p)-target),
            mapping_valid=True), initial_jacobian=np.zeros((2, 2)), default_jacobian=matrix)
        self.assertTrue(result['converged'])
        self.assertEqual(result['evaluations'], 3)
        np.testing.assert_allclose(result['history'][1]['actual_log_step'], target, atol=1e-15)
        np.testing.assert_array_equal(matrix, original)

    def test_invalid_reset_matrix_fails_before_evaluation(self):
        calls = []
        def evaluate(p):
            calls.append(p)
            return dict(residual=np.zeros(2), mapping_valid=True)
        for matrix in (np.eye(3), np.zeros((2, 2)), np.full((2, 2), np.nan)):
            with self.subTest(matrix=matrix), self.assertRaisesRegex(ValueError, 'Default Jacobian'):
                self.run_root(evaluate, default_jacobian=matrix)
        self.assertFalse(calls)

    def test_ill_conditioned_reset_raises_price_for_positive_residual(self):
        def evaluate(p):
            return dict(residual=np.array([.1, .1]) - np.log(p), mapping_valid=True)
        with patch('numpy.linalg.cond', return_value=np.inf):
            result = self.run_root(evaluate)
        self.assertTrue(result['converged'])
        trial = result['history'][1]
        self.assertTrue(np.all(trial['actual_log_step'] > 0))
        self.assertIn('positive_residual_price_increase', trial['reset_reason'])

    def test_material_worsening_restores_best_and_halves_damping(self):
        def evaluate(p):
            return dict(residual=np.ones(2) * .1 + np.log(p), mapping_valid=True)
        result = self.run_root(evaluate, max_evaluations=5)
        self.assertFalse(result['converged'])
        self.assertEqual(result['history'][1]['safeguard'], 'restored_best_halved_damping_reset_jacobian')
        np.testing.assert_allclose(result['history'][2]['prices'], np.exp([.05, .05]))
        np.testing.assert_array_equal(result['best']['prices'], np.ones(2))

    def test_budget_exhaustion_never_certifies_unclosed_market(self):
        result = self.run_root(lambda p: dict(residual=np.ones(2), mapping_valid=True), max_evaluations=2)
        self.assertFalse(result['converged'])
        self.assertEqual(result['evaluations'], 2)
        self.assertEqual(result['history'][-1]['phase'], 'final')

    def test_projection_enters_broyden_with_actual_step(self):
        A = np.array([[-2., .4], [.1, -1.]])
        target = np.array([.3, -.2])
        def project(p):
            x = np.log(p)
            return np.exp(np.array([min(x[0], .04), x[1]]))
        result = self.run_root(lambda p: dict(residual=A @ (np.log(p) - target), mapping_valid=True),
            project=project, max_evaluations=3, worsening_factor=10.)
        row = result['history'][1]
        s = row['actual_log_step']
        self.assertNotEqual(s[0], row['requested_log_step'][0])
        y = row['residual'] - result['history'][0]['residual']
        expected = -np.eye(2) + np.outer(y + s, s) / (s @ s)
        np.testing.assert_allclose(result['final_jacobian'], expected)

    def test_discontinuous_jump_cannot_pass_on_small_steps(self):
        result = self.run_root(lambda p: dict(residual=np.where(np.log(p) < .05, .1, -.1), mapping_valid=True),
            max_evaluations=20)
        self.assertFalse(result['converged'])
        self.assertGreater(result['best']['score'], 1e-9)

    def test_fresh_final_drift_rejects_candidate(self):
        calls = []
        def evaluate(p):
            calls.append(1)
            return dict(residual=np.zeros(2) if len(calls) == 1 else np.ones(2) * 1e-7, mapping_valid=True)
        result = self.run_root(evaluate)
        self.assertFalse(result['converged'])
        self.assertEqual(result['evaluations'], 2)
        self.assertEqual(result['status'], 'final_market_gate_failed')

    def test_deadline_cannot_replace_required_final_replay(self):
        now = [0.]
        def evaluate(p):
            now[0] = 3.
            return dict(residual=np.zeros(2), mapping_valid=True)
        with patch('e5f_matched_pf_path_root.time.monotonic', side_effect=lambda: now[0]):
            result = self.run_root(evaluate, deadline_monotonic=2.)
        self.assertFalse(result['converged'])
        self.assertEqual(result['evaluations'], 1)
        self.assertIsNone(result['final'])
        self.assertEqual(result['status'], 'time_or_evaluation_budget')

    def test_zero_residual_with_failed_mapping_does_not_pass(self):
        result = self.run_root(lambda p: dict(residual=np.zeros(2), mapping_valid=False))
        self.assertFalse(result['converged'])
        self.assertIsNone(result['best'])
        self.assertEqual(result['status'], 'invalid_initial_mapping')


if __name__ == '__main__':
    unittest.main()
