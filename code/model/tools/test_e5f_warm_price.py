"""Price-search safeguards and integration routing; run on Torch."""
from contextlib import ExitStack
from types import SimpleNamespace
from unittest.mock import patch
import unittest

import numpy as np
from intergen_eqscale_seq_optimized import solver
from intergen_eqscale_seq_optimized.warm_price import search_warm_price
import e5f_stationary_paygo as paygo


class SearchTests(unittest.TestCase):
    def run_search(self, function, **kwargs):
        calls = []
        def evaluate(price):
            calls.append(price)
            excess = function(price)
            return excess, abs(excess), price
        result = search_warm_price(evaluate, initial_price=1., initial_slope=None,
                                   lower_bound=.5, upper_bound=2., tolerance=1e-10, **kwargs)
        return result, calls

    def test_initial_exact_root_exits_without_probe(self):
        result, calls = self.run_search(lambda p: 1. - p)
        self.assertEqual(calls, [1.])
        self.assertFalse(result[3]['fallback_required'])

    def test_secant_finds_linear_root(self):
        result, calls = self.run_search(lambda p: 1.04 - p)
        self.assertAlmostEqual(result[0], 1.04)
        self.assertLess(result[2], 1e-10)
        self.assertEqual(len(calls), 3)

    def test_flat_demand_cannot_pass_or_repeat_point(self):
        result, calls = self.run_search(lambda p: 1.)
        self.assertTrue(result[3]['fallback_required'])
        self.assertEqual(len(calls), len(set(calls)))
        self.assertLessEqual(len(calls), 9)
        self.assertTrue(all(.5 <= p <= 2 for p in calls))

    def test_stale_slope_and_global_boundary_do_not_escape(self):
        calls = []
        def evaluate(price):
            calls.append(price)
            return 1., 1., None
        result = search_warm_price(evaluate, initial_price=1.99, initial_slope=-.0001,
            lower_bound=.5, upper_bound=2., tolerance=1e-10)
        self.assertEqual(calls[:2], [1.99, 2.])
        self.assertTrue(all(.5 <= p <= 2. for p in calls))
        self.assertTrue(result[3]['fallback_required'])

    def test_discontinuous_market_is_never_certified_by_bracket_width(self):
        result, calls = self.run_search(lambda p: .001 if p < 1.011 else -.001)
        self.assertTrue(result[3]['fallback_required'])
        self.assertEqual(result[2], .001)
        self.assertLessEqual(len(calls), 9)

    def test_invalid_input_and_nonfinite_residual_fail(self):
        for initial, low, high, slope in [(0., .5, 2., None), (1., 2., 3., None),
                                          (1., .5, 2., float('nan'))]:
            with self.assertRaises(ValueError):
                search_warm_price(lambda p: self.fail('must fail before evaluation'),
                    initial_price=initial, initial_slope=slope,
                    lower_bound=low, upper_bound=high, tolerance=1e-10)
        with self.assertRaises(RuntimeError):
            self.run_search(lambda p: float('nan'))


class SolverRoutingTests(unittest.TestCase):
    def solve(self, state, *, root=1.04, fail_refine=False, discontinuous=False):
        P = SimpleNamespace(I=1, tol_eq=2.5e-5, markov_equilibrium_method='direct',
            p_min=.5, p_max=2., scalar_market_refine=True,
            adult_entry_clock='split_birth_vintage')
        prices, refine_calls, upgrades = [], [], []
        def at_prices(p, P, grid, **kwargs):
            value = float(p[0]); prices.append(value)
            excess = (.001 if value < root else -.001) if discontinuous else root - value
            return SimpleNamespace(housing_supply=np.array([1.]), demand=np.array([1.+excess]),
                entry_rate=.04, adult_entry_potential_total=.04, timings={}, _model_payload=object())
        def refine(p, sol, error, P, grid, **kwargs):
            refine_calls.append(float(p[0]))
            if fail_refine:
                return sol, p, error, {'used': True}
            result = at_prices(np.array([root]), P, grid)
            return result, np.array([root]), 0., {'used': True}
        def upgrade(sol, *args):
            upgrades.append(sol)
            return sol
        with ExitStack() as stack:
            stack.enter_context(patch.object(solver, 'income_transition_values', return_value=([], [], np.array([]))))
            stack.enter_context(patch.object(solver, 'precompute_shared', return_value=None))
            stack.enter_context(patch.object(solver, 'solve_markov_income_at_prices', side_effect=at_prices))
            stack.enter_context(patch.object(solver, 'markov_market_housing_demand', side_effect=lambda s,*a:(s.demand,None)))
            stack.enter_context(patch.object(solver, 'refine_one_market_markov_income', side_effect=refine))
            stack.enter_context(patch.object(solver, 'upgrade_fast_markov_solution', side_effect=upgrade))
            attach = stack.enter_context(patch.object(solver, 'attach_markov_market_accounting', side_effect=lambda s,*a:s))
            result = solver.solve_markov_income_equilibrium(np.array([1.]), P, np.array([0.]),
                verbose=False, warm_price_state=state)
            self.assertEqual(attach.call_count, 1)
            self.assertEqual(len(upgrades), 1)
        return result, prices, refine_calls

    def test_absent_state_keeps_cold_route_and_metadata(self):
        (sol,P,p), prices, refine = self.solve(None)
        self.assertEqual(prices, [1., 1.04])
        self.assertEqual(refine, [1.])
        self.assertNotIn('warm_price_search', sol.timings)
        self.assertTrue(sol.converged)
        self.assertEqual(sol.adult_entry_stationary_relative_gap, 0.)

    def test_empty_state_starts_cold_and_carries_only_scalar_summary(self):
        state = {}
        (sol,P,p), prices, refine = self.solve(state)
        self.assertEqual(prices, [1., 1.04])
        self.assertEqual(state['price'], 1.04)
        self.assertEqual(set(state), {'price', 'slope'})
        self.assertFalse(sol.timings['warm_price_search']['used'])

    def test_warm_route_reuses_common_finalizer_and_skips_refinement_on_success(self):
        state = {'price':1., 'slope':-1.}
        (sol,P,p), prices, refine = self.solve(state)
        self.assertEqual(refine, [])
        self.assertEqual(len(prices), 2)
        self.assertTrue(sol.converged)
        self.assertAlmostEqual(state['slope'], -1.)
        self.assertEqual(sol.timings['unique_fast_price_evaluations'], 2)

    def test_initial_warm_root_retains_previous_slope(self):
        state={'price':1.04,'slope':-1.}
        (sol,P,p),prices,refine=self.solve(state)
        self.assertEqual(prices,[1.04])
        self.assertEqual(refine,[])
        self.assertEqual(state['slope'],-1.)

    def test_fallback_and_failed_gate_do_not_overwrite_state(self):
        state = {'price':1., 'slope':None}
        before = dict(state)
        (sol,P,p), prices, refine = self.solve(state, root=1.011, fail_refine=True, discontinuous=True)
        self.assertEqual(len(refine), 1)
        self.assertFalse(sol.converged)
        self.assertEqual(state, before)
        self.assertTrue(sol.timings['warm_price_search']['fallback_required'])

    def test_fallback_success_is_finalized_normally(self):
        state = {'price':1., 'slope':None}
        (sol,P,p), prices, refine = self.solve(state, root=1.011, discontinuous=True)
        self.assertEqual(len(refine), 1)
        self.assertTrue(sol.converged)
        self.assertEqual(state['price'], 1.011)

    def test_pension_failure_does_not_commit_warm_state(self):
        state = {'price':1., 'slope':-1.}
        def solve(*args, warm_price_state, **kwargs):
            warm_price_state.update(price=2., slope=-2.)
            return SimpleNamespace(converged=True, timings={'strict_converged':True},g=None), args[1], np.array([2.])
        P = SimpleNamespace(tol_eq=2.5e-5)
        with patch.object(paygo,'bind_initial_balanced_pension',return_value=(P,{})), \
             patch.object(paygo,'certify_initial_pension',side_effect=RuntimeError('fiscal failure')):
            with self.assertRaisesRegex(RuntimeError,'fiscal failure'):
                paygo.solve_balanced_initial_equilibrium(model=SimpleNamespace(solve_markov_income_equilibrium=solve),
                    parameters=P,b_grid=np.array([0.]),initial_prices=np.array([1.]),payroll_tax=.08,
                    marginal_tolerance=1e-9,fiscal_tolerance=1e-6,warm_price_state=state)
        self.assertEqual(state, {'price':1., 'slope':-1.})


if __name__ == '__main__':
    unittest.main()
