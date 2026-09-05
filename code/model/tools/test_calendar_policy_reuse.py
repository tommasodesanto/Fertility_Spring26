"""Regression checks for self-contained calendar policies and bounded searches.

Run with the project interpreter (unittest; no pytest plugin dependency).
"""
from __future__ import annotations

import pickle
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import patch

import numpy as np

MODEL = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(MODEL), str(MODEL / 'tools')]
import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as post
import run_e5f_perfect_foresight_transition as foresight
from intergen_eqscale_seq_optimized import solver


class PolicyReuseTests(unittest.TestCase):
    def setUp(self):
        self.P = NS(J=1, I=1, n_house=0, n_parity=4, n_child_states=4,
                    A_f_start=1, A_f_end=1, sequential_births=True,
                    child_state_mode='independent_count', readiness_gate_enabled=False,
                    fecundity_omega1=0., fecundity_omega2=0., age_start=18, da=4.,
                    user_cost_rate=.1, H0=np.array([1.]), r_bar=.1, xi_supply=np.array([1.]))
        self.g = np.zeros((2, 1, 1, 1, 1, 4, 4))
        self.g[0, 0, 0, 0, 0, 1, 1] = .5
        self.g[1, 0, 0, 0, 0, 0, 0] = .5
        self.fp = np.zeros((2, 1, 1, 1, 1, 4))
        self.fp[..., 1] = 1.
        self.cont = np.zeros((2, 1, 1, 1, 1, 2, 2, 4))
        self.cont[..., 1, :, :] = .2
        self.P._fert2_probs = self.cont
        self.maps = NS(lmm_idx=None, lmm_wt=None, tmx_idx=None, tmx_wt=None)
        self.pol = calendar.PolicyBundle(
            V=np.ones_like(self.g), c_pol=np.ones_like(self.g),
            hR_pol=np.ones_like(self.g), bp_pol=np.zeros_like(self.g),
            tenure_choice=np.zeros_like(self.g, dtype=int), tenure_probs=None,
            loc_probs=np.ones((2, 1, 1, 1, 1, 1, 4, 4)), fert_probs=self.fp,
            fert_value=np.zeros((2, 1, 1, 1, 1)), price=np.array([1.]),
            maps=self.maps, fert2_probs=self.cont)
        self.counter = calendar.SolveCounter()
        for mocked in (
            patch.object(calendar, 'model', solver),
            patch.object(calendar, 'apply_fertility', transition.apply_sequential_fertility),
            patch.object(calendar, 'gate_pre_fertility_distribution',
                         side_effect=lambda g, *args: (g.copy(), 0.)),
            patch.object(solver, 'realize_current_cross_section',
                         side_effect=lambda g, *args, **kw: g.copy()),
            patch.object(calendar, 'build_transition_maps', return_value=self.maps),
        ):
            mocked.start()
            self.addCleanup(mocked.stop)

    def evaluate(self, policy=None):
        return calendar.evaluate_period(np.array([1.]), self.g, self.P,
            np.array([0., 1.]), NS(), self.counter,
            supplied_policy=self.pol if policy is None else policy)

    def test_reuse_after_another_solve_and_pickle_preserves_entire_distribution(self):
        first = self.evaluate()
        self.assertAlmostEqual(first.births, .6)
        # Mutate the original array, then replace the solver side channel too.
        self.cont[...] = 1.
        self.P._fert2_probs = np.zeros_like(self.cont)
        for policy in (self.pol, pickle.loads(pickle.dumps(self.pol))):
            again = self.evaluate(policy)
            self.assertEqual(first.births, again.births)
            np.testing.assert_array_equal(first.g_post_fertility, again.g_post_fertility)
            np.testing.assert_array_equal(first.g_current, again.g_current)
            self.assertEqual(float(again.g_current.sum()), 1.)
            # New first births cannot have a second birth in the same period.
            self.assertEqual(float(again.g_current[1, ..., 2:, :].sum()), 0.)
        self.assertEqual(self.counter.bellman, 0)

    def test_all_factories_capture_the_matching_continuation(self):
        solution = NS(**vars(self.pol))
        self.P._fert2_probs = np.zeros_like(self.cont)
        saved = calendar.policy_from_solution(solution, np.array([1.]), self.P,
                                             np.array([0., 1.]), NS())
        np.testing.assert_array_equal(saved.fert2_probs, self.pol.fert2_probs)
        objects = tuple(getattr(self.pol, name) for name in (
            'V', 'c_pol', 'hR_pol', 'bp_pol', 'tenure_choice', 'tenure_probs',
            'loc_probs', 'fert_probs', 'fert_value')) + (None,)
        def bellman(*args):
            self.P._fert2_probs = self.cont.copy()
            return objects
        with patch.object(solver, 'solve_bellman_full_markov_income', side_effect=bellman):
            solved = calendar.solve_policy(np.array([1.]), self.P,
                                           np.array([0., 1.]), NS(), self.counter)
            forward = foresight.policy_from_objects(objects, 1., self.P,
                                                    np.array([0., 1.]), NS())
        self.P._fert2_probs[...] = 0.
        for policy in (saved, solved, forward):
            np.testing.assert_array_equal(policy.fert2_probs, self.cont)
        self.assertEqual(self.counter.bellman, 1)

    def test_old_incomplete_sequential_checkpoint_fails_without_using_P(self):
        del self.pol.fert2_probs
        with self.assertRaisesRegex(RuntimeError, 'rebuild'):
            self.evaluate(pickle.loads(pickle.dumps(self.pol)))

    def test_wrong_price_is_rejected(self):
        self.pol.price = np.array([2.])
        with self.assertRaisesRegex(ValueError, 'price'):
            self.evaluate()

    def test_failed_market_retry_uses_original_policy(self):
        births = []
        def market(*args, **kwargs):
            evaluation = self.evaluate(kwargs['initial_policy'])
            births.append(evaluation.births)
            if len(births) == 1:
                self.P._fert2_probs[...] = 1.
                raise RuntimeError('Housing market did not clear: forced retry')
            return evaluation
        state = NS(g_pre=self.g, price_guess=1., initial_policy=self.pol)
        with patch.object(calendar, 'clear_scalar_housing_market', side_effect=market):
            _, fallback = post.clear_market(state, self.P, np.array([0., 1.]),
                NS(), self.counter, calendar.HousingSupplyRule('fixed-stock', 1., 1., 0.),
                2e-4, 30)
        self.assertTrue(fallback)
        self.assertEqual(births, [.6, .6])

    def test_legacy_one_shot_does_not_require_continuation(self):
        self.P.sequential_births = False
        self.pol.fert2_probs = None
        # Invoke the original one-shot operator, preserved beneath the patch.
        with patch.object(calendar, 'apply_fertility', ONE_SHOT_FERTILITY):
            self.assertAlmostEqual(self.evaluate().births, .5)


ONE_SHOT_FERTILITY = calendar.apply_fertility


class BracketTests(unittest.TestCase):
    def test_out_of_bounds_guess_fails_before_solving(self):
        with patch.object(calendar, 'evaluate_period') as evaluate:
            with self.assertRaisesRegex(ValueError, 'price bounds'):
                calendar.clear_scalar_housing_market(np.ones(1), .8,
                    NS(I=1, p_min=.9, p_max=1.1), np.array([0., 1.]), NS(),
                    calendar.SolveCounter(), 1e-4, 2,
                    calendar.HousingSupplyRule('fixed-stock', 1., 1., 0.))
        evaluate.assert_not_called()

    def test_each_exhausted_boundary_is_evaluated_only_once(self):
        for excess, bounds, expected in ((1., (.1, 1.1), [1., 1.1]),
                                         (-1., (.9, 10.), [1., .9])):
            calls = []
            def evaluate(price, *args, **kwargs):
                calls.append(float(price[0]))
                return NS(demand_by_loc=np.array([2. + excess]),
                          supply_by_loc=np.array([2.]), relative_market_residual=.5)
            with patch.object(calendar, 'evaluate_period', side_effect=evaluate):
                with self.assertRaisesRegex(RuntimeError, 'Could not bracket'):
                    calendar.clear_scalar_housing_market(np.ones(1), 1.,
                        NS(I=1, p_min=bounds[0], p_max=bounds[1]), np.array([0., 1.]),
                        NS(), calendar.SolveCounter(), 1e-4, 2,
                        calendar.HousingSupplyRule('fixed-stock', 1., 2., 0.))
            self.assertEqual(calls, expected)

    def test_root_at_boundary_still_returns_success(self):
        def evaluate(price, *args, **kwargs):
            residual = 1.1 - float(price[0])
            return NS(demand_by_loc=np.array([1. + residual]),
                      supply_by_loc=np.array([1.]), relative_market_residual=abs(residual))
        with patch.object(calendar, 'evaluate_period', side_effect=evaluate):
            result = calendar.clear_scalar_housing_market(np.ones(1), 1.,
                NS(I=1, p_min=.1, p_max=1.1), np.array([0., 1.]), NS(),
                calendar.SolveCounter(), 1e-4, 2,
                calendar.HousingSupplyRule('fixed-stock', 1., 1., 0.))
        self.assertEqual(result.relative_market_residual, 0.)


if __name__ == '__main__':
    unittest.main()
