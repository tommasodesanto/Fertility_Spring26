"""Small routing/fit tests; no household or equilibrium model solves."""
import unittest
from types import SimpleNamespace as NS

import numpy as np

import run_e5f_long_successive_refit as refit


class RefitContractTests(unittest.TestCase):
    def test_same_terminal_boundary_leaves_matched_policy_tail(self):
        self.assertEqual([refit.horizon(y) for y in (2007, 2011, 2015, 2019, 2023)],
                         [104, 103, 102, 101, 100])
        self.assertEqual(refit.horizon(2019) - 1, refit.horizon(2023))
        with self.assertRaises(ValueError):
            refit.horizon(2008)
        with self.assertRaises(ValueError):
            refit.horizon(2423)

    def test_initial_seed_and_bounds(self):
        self.assertAlmostEqual(refit.next_psi([], .12, (-.05, .17), []), .12)
        self.assertAlmostEqual(refit.next_psi([], .25, (-.05, .17), []), .17)

    def test_bracketed_fit_stays_inside_bracket(self):
        trials = [dict(psi=.08, gap=-.4), dict(psi=.14, gap=.001)]
        x = refit.next_psi(trials, .12, (-.05, .17), [.08, .14])
        self.assertGreaterEqual(x, .08 + .15 * .06 - 1e-12)
        self.assertLessEqual(x, .14 - .15 * .06 + 1e-12)

    def test_failed_candidate_is_not_a_fertility_observation(self):
        trials = [dict(psi=.10, gap=None), dict(psi=.11, gap=float('nan'))]
        x = refit.next_psi(trials, .12, (-.05, .17), [.10, .11, .12])
        self.assertTrue(np.isfinite(x))
        self.assertNotIn(x, [.10, .11, .12])

    def test_scaled_solver_patch_keeps_factory_interface(self):
        original = lambda: 'original'
        actual = lambda **kwargs: kwargs
        c = NS(rebated=NS(_path_root_solver=original))
        with refit.scaled_root_context(c, NS(solve_price_path_scaled=actual)):
            self.assertIs(c.rebated._path_root_solver(), actual)
        self.assertIs(c.rebated._path_root_solver, original)

    def test_unconverged_forecast_cannot_be_accepted(self):
        result = NS(path=NS(), next_state=NS(), root_receipt={
            'finite_horizon_market_fiscal_converged': False})
        ok, gap = refit.accepted(result, [{'period_tfr_topcode_adjusted': 1.974875}], 1.974875)
        self.assertFalse(ok)
        self.assertIsNone(gap)

    def test_valid_wrong_fit_is_available_to_scalar_search(self):
        result = NS(path=NS(), next_state=NS(), root_receipt={
            'finite_horizon_market_fiscal_converged': True})
        ok, gap = refit.accepted(result, [{'period_tfr_topcode_adjusted': 1.9}], 1.8)
        self.assertFalse(ok)
        self.assertAlmostEqual(gap, .1)


if __name__ == '__main__':
    unittest.main()
