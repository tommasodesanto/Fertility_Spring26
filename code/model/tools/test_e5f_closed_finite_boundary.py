"""Bounded pure adapter checks; no Bellman solve, Numba compile or model run."""
from __future__ import annotations

from dataclasses import make_dataclass
import importlib.util
from pathlib import Path
import sys
import time
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import e5f_closed_finite_boundary as boundary

# Reuse the actual fiscal definitions without importing the numerical runtime.
try:
    import e5f_social_security as social
except ModuleNotFoundError as error:
    if error.name != 'e5f_social_security':
        raise
    ROOT = Path(__file__).resolve().parents[3]
    spec = importlib.util.spec_from_file_location('boundary_test_social',
        ROOT / 'tmp/e5f_matched_pf/code/model/tools/e5f_social_security.py')
    social = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(social)


class FiniteBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.P = NS(I=1, J=2, J_R=1, Nb=2, period_years=4.,
            scale_flows_to_period=True, tau_pay=.2, pension=9.,
            property_tax_lump_sum_transfer=9., z_grid=np.array([1.]),
            w_hat=np.array([2.]), income_age_profile=np.array([1., 0.]),
            retirement_income_z_scale=0., exhaustive_saving_control=True)
        self.g = np.zeros((2, 1, 1, 2, 1, 1, 1))
        self.g[0, 0, 0, 0, 0, 0, 0] = 3.
        self.g[0, 0, 0, 1, 0, 0, 0] = 2.
        self.seen = []
        self.audit_pass = True
        Audit = make_dataclass('Audit', [(n, float) for n in
            ('reconstruction_tolerance', 'feasibility_projection_tolerance',
             'probability_tolerance', 'occupied_mass_tolerance', 'value_drop_tolerance')])

        def bellman(rents, prices, P, grid, shared, *, continuation_V):
            self.assertIsNone(continuation_V)  # Native lifetime/death recursion.
            self.seen.append((P.pension, P.property_tax_lump_sum_transfer,
                              P.income.copy(), continuation_V))
            return np.ones_like(self.g) * (prices[0] + P.pension + P.property_tax_lump_sum_transfer)

        def actual(prices, g, P, grid, shared, counter, *, supply_rule, supplied_policy):
            self.assertIs(supply_rule, self.supply)
            return NS(g_pre=g, g_current=g.copy(), g_post_fertility=g.copy(),
                policy=supplied_policy, demand_by_loc=np.array([2.*g.sum()]),
                supply_by_loc=supply_rule.quantity(prices))

        self.model = NS(precompute_shared=lambda P, grid: NS(),
            solve_bellman_full_markov_income=bellman,
            property_tax_revenue_from_distribution=lambda g, h, prices, P: .1*prices[0]*g.sum())
        self.calendar = NS(evaluate_period=actual, SolveCounter=lambda: None)
        self.balanced = NS(TerminalAuditControls=Audit,
            _household_checks=lambda *args: ({'audit': 'fake'}, {'households': self.audit_pass}))
        self.pf = NS(rents_from_asset_prices=lambda prices, terminal, P: .1*prices,
            policy_from_objects=lambda V, price, P, grid, shared:
                NS(V=V, price=np.array([price]), hR_pol=np.ones_like(V)))
        self.supply = NS(quantity=lambda prices: np.array([10.]))
        self.runtime = (self.model, self.calendar, None, self.balanced, social, self.pf)
        self.args = dict(parameters=self.P, g_pre=self.g, grid=np.array([0., 1.]),
            supply_rule=self.supply, price=2., pension=2.4, transfer=.2)

    def evaluate(self, **changes):
        with patch.object(boundary, '_runtime', return_value=self.runtime):
            return boundary.boundary_evaluation(**dict(self.args, **changes))

    def test_actual_mass_and_two_separate_budgets(self):
        result = self.evaluate()
        self.assertEqual(result.actual_accounts['household_heads'], 5.)
        self.assertAlmostEqual(result.actual_accounts['payroll_tax_revenue'], 4.8)
        self.assertAlmostEqual(result.actual_accounts['pension_outlays'], 4.8)
        self.assertEqual(result.actual_accounts['property_tax_revenue'], 1.)
        self.assertEqual(result.actual_accounts['equal_transfer_outlays'], 1.)
        for name in ('housing_relative', 'pension_relative', 'rebate_relative'):
            self.assertAlmostEqual(result.residuals[name], 0.)
        self.assertTrue(result.mapping_valid)
        self.assertFalse(result.horizon_verified)
        self.assertFalse(result.production_eligible)
        self.assertFalse(hasattr(result, 'fixed_point'))

    def test_population_scaling_is_not_normalized_away(self):
        result = self.evaluate(g_pre=2.*self.g)
        self.assertEqual(result.actual_accounts['household_heads'], 10.)
        self.assertEqual(result.residuals['housing_relative'], 1.)
        self.assertAlmostEqual(result.actual_accounts['payroll_tax_revenue'], 9.6)
        self.assertTrue(result.mapping_valid)  # Mapping-valid is not equilibrium.

    def test_age_composition_changes_actual_paygo_exposure(self):
        altered = self.g.copy()
        altered[0, 0, 0, 0, 0, 0, 0] = 2.
        altered[0, 0, 0, 1, 0, 0, 0] = 3.
        result = self.evaluate(g_pre=altered)
        self.assertAlmostEqual(result.actual_accounts['implied_balanced_pension_period'], 3.2/3.)
        self.assertLess(result.residuals['pension_relative'], 0.)
        self.assertAlmostEqual(result.residuals['rebate_relative'], 0.)

    def test_fiscal_binding_precedes_household_solve_and_preserves_inputs(self):
        before = self.g.copy()
        result = self.evaluate()
        self.assertEqual(self.seen[0][:2], (2.4, .2))
        np.testing.assert_array_equal(self.seen[0][2], [[6.4, 2.4]])
        self.assertEqual(self.P.pension, 9.)
        self.assertEqual(self.P.property_tax_lump_sum_transfer, 9.)
        np.testing.assert_array_equal(self.g, before)
        self.assertFalse(np.shares_memory(result.g_pre, self.g))

    def test_household_failure_blocks_mapping(self):
        self.audit_pass = False
        self.assertFalse(self.evaluate().mapping_valid)

    def test_invalid_population_and_deadline_fail_before_household_solve(self):
        for invalid in (np.zeros_like(self.g), -self.g, self.g*np.nan):
            with self.assertRaises(ValueError):
                self.evaluate(g_pre=invalid)
        with self.assertRaises(TimeoutError):
            self.evaluate(deadline_monotonic=time.monotonic()-1.)
        self.assertFalse(self.seen)

    def test_progress_and_no_fictitious_future_budget_certificate(self):
        records = []
        result = self.evaluate(callback=records.append)
        self.assertEqual(records[-1]['phase'], 'boundary_complete')
        self.assertFalse(result.diagnostics['stationary_population_computed'])
        self.assertFalse(result.diagnostics['future_fiscal_consistency_verified'])


if __name__ == '__main__':
    unittest.main()
