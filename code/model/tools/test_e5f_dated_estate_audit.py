"""Focused synthetic checks for actual next-cohort estate funding; no solves."""
import copy
from types import SimpleNamespace
import unittest
import numpy as np
from e5f_overnight_estate_audit import audit, EstateFundingShortfall


def fixture():
    g = np.zeros((2, 2, 1, 2, 1, 1, 1))
    g[0, 0, 0, 0, 0, 0, 0] = .4
    g[1, 0, 0, 0, 0, 0, 0] = .6
    g[0, 0, 0, 1, 0, 0, 0] = .25
    g[0, 1, 0, 1, 0, 0, 0] = .5
    g[1, 1, 0, 1, 0, 0, 0] = .25
    bp = np.zeros_like(g)
    bp[0, 0, 0, 1, 0, 0, 0] = -2.
    bp[0, 1, 0, 1, 0, 0, 0] = -9.
    bp[1, 1, 0, 1, 0, 0, 0] = 5.
    P = SimpleNamespace(I=1, J=2, H_own=np.array([10.]), period_years=4.,
                        psi=.2, use_age_survival=True, survival_probs=np.array([1.]))
    e = SimpleNamespace(g_current=g, g_pre=g.copy(),
                        policy=SimpleNamespace(bp_pol=bp, price=np.array([1.])))
    return e, P, np.array([-2., 5.])


class DatedEstateTests(unittest.TestCase):
    def test_stationary_equal_cohort_economic_equality(self):
        e, P, grid = fixture()
        stationary = audit(e, P, grid)
        dated = audit(e, P, grid, next_entrant_cohort=e.g_pre[:, :, :, 0].copy())
        self.assertEqual(stationary['audit_id'], 'estate_funded_entry_provisional_net_v1')
        self.assertEqual(stationary['timing'], 'stationary death flow finances the next entrant cohort; no transition implementation')
        self.assertIn('dated', dated['audit_id'])
        self.assertIn('date-t+1', dated['timing'])
        for key in stationary:
            if key not in ('audit_id', 'timing'):
                self.assertEqual(stationary[key], dated[key], key)
        self.assertEqual(stationary['entry']['positive_financial_endowment'], 3.)
        self.assertEqual(stationary['available_estates_period'], 3.25)

    def test_actual_next_cohort_scales_cost_not_estates(self):
        e, P, grid = fixture()
        s = audit(e, P, grid)
        d = audit(e, P, grid, next_entrant_cohort=e.g_pre[:, :, :, 0] * .5)
        self.assertEqual(d['entry']['positive_financial_endowment'], 1.5)
        self.assertEqual(d['residual_sink_period'], 1.75)
        self.assertEqual(s['estate'], d['estate'])
        self.assertEqual(d['entry']['negative_financial_position'], .4)

    def test_zero_next_cohort(self):
        e, P, grid = fixture()
        d = audit(e, P, grid, next_entrant_cohort=np.zeros_like(e.g_pre[:, :, :, 0]))
        self.assertEqual(d['entry']['mass'], 0.)
        self.assertEqual(d['entry']['positive_financial_endowment'], 0.)
        self.assertEqual(d['residual_sink_period'], 3.25)

    def test_shortfall_preserves_dated_failure(self):
        e, P, grid = fixture()
        with self.assertRaises(EstateFundingShortfall) as caught:
            audit(e, P, grid, next_entrant_cohort=e.g_pre[:, :, :, 0] * 2)
        self.assertEqual(caught.exception.audit['funding_shortfall_period'], 2.75)
        self.assertIsNone(caught.exception.audit['residual_sink_period'])
        self.assertIn('dated', caught.exception.audit['audit_id'])

    def test_no_mutation(self):
        e, P, grid = fixture()
        cohort = e.g_pre[:, :, :, 0].copy()
        before = copy.deepcopy((e, P, grid, cohort))
        audit(e, P, grid, next_entrant_cohort=cohort)
        for actual, old in ((e.g_pre, before[0].g_pre), (e.g_current, before[0].g_current),
                            (e.policy.bp_pol, before[0].policy.bp_pol),
                            (grid, before[2]), (cohort, before[3])):
            np.testing.assert_array_equal(actual, old)

    def test_invalid_cohorts_and_original_age_gate(self):
        e, P, grid = fixture()
        base = e.g_pre[:, :, :, 0].copy()
        for invalid in (np.zeros((2,)), base * np.nan, -base):
            with self.assertRaises(ValueError):
                audit(e, P, grid, next_entrant_cohort=invalid)
        owner = base.copy()
        owner[0, 1, 0, 0, 0, 0] = .1
        with self.assertRaises(ValueError):
            audit(e, P, grid, next_entrant_cohort=owner)
        e.g_pre[:, :, :, 1] *= .5
        with self.assertRaisesRegex(ValueError, 'preserve positive age mass'):
            audit(e, P, grid, next_entrant_cohort=base)


if __name__ == '__main__':
    unittest.main()
