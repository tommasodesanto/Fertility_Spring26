"""Conditional PF history tests: real age bridge/queue, stubbed Bellman/KFE."""
from contextlib import ExitStack
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import run_e5f_perfect_foresight_transition as pf


class HistoricalBridgeTests(unittest.TestCase):
    def setUp(self):
        self.P = SimpleNamespace(I=2, J=17, period_years=4., age_start=18., da=4.,
                                 entry_shares=np.array([.25, .75]), R_gross=1.04,
                                 delta=.02, tau_H=.01, entrant_conversion_factor=.5)
        self.shape = (2, 2, 2, 17, 1, 2, 2)
        self.g = np.zeros(self.shape)
        for age in range(17):
            self.g[0, 0, 0, age, 0, 0, 0] = (age + 1) / 153 * .25
            self.g[1, 0, 0, age, 0, 0, 0] = (age + 1) / 153 * .75
        self.state = pf.PFInitialState(self.g, [.1, .2, .3, .4], [.1, .2, .3, .4])
        self.flows = []
        self.events = []
        self.joint_marker = object()
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)

        def backward(**kwargs):
            return [np.zeros(self.shape) for _ in range(len(kwargs['prices']) + 1)], len(kwargs['prices'])

        def solve(**kwargs):
            return SimpleNamespace(V=np.zeros(self.shape), price=np.array([kwargs['price']]),
                                   hR_pol=np.ones(self.shape), joint_choice=self.joint_marker)

        def evaluate(price, g, P, grid, shared, counter, **kwargs):
            self.events.append('evaluate')
            return SimpleNamespace(policy=kwargs['supplied_policy'], g_pre=g.copy(),
                g_post_fertility=g.copy(), g_current=g.copy(), births=.21,
                demand_by_loc=np.array([1., 1.]), supply_by_loc=np.array([1., 1.]),
                relative_market_residual=0., feasibility_projection_mass=0.)

        def advance(e, entries, P, grid, shared):
            self.events.append('advance')
            nxt = np.zeros_like(e.g_post_fertility)
            nxt[:, :, :, 1:] = e.g_post_fertility[:, :, :, :-1]
            deaths = e.g_post_fertility[:, :, :, -1].sum()
            return nxt, np.zeros(P.I), float(deaths), 0.

        def entrants(flows, P, grid):
            self.flows.append(flows.copy())
            out = np.zeros((2, 2, P.I, 1, 2, 2))
            out[0, 0, :, 0, 0, 0] = .25 * flows
            out[1, 0, :, 0, 0, 0] = .75 * flows
            return out

        for target, name, implementation in (
            (pf, 'backward_value_path', backward), (pf, 'solve_date_policy', solve),
            (pf.calendar, 'evaluate_period', evaluate),
            (pf.calendar, 'entrant_cohort', entrants),
            (pf.calendar.model, 'precompute_shared', lambda P, grid: SimpleNamespace()),
            (pf.calendar.model, 'property_tax_revenue_from_distribution', lambda *args: 0.),
            (pf.transition, 'advance_sequential_calendar_distribution', advance),
            (pf.transition, 'calendar_topcode_birth_accounting',
             lambda *args: {'topcode_adjusted_birth_children': .21}),
        ):
            self.stack.enter_context(patch.object(target, name, side_effect=implementation))

    def evaluate(self, periods, conditioning='omit'):
        kwargs = dict(prices=np.ones(periods), psi_path=np.full(periods, .1),
            terminal_price=1., terminal_V=np.zeros(self.shape), base_parameters=self.P,
            b_grid=np.arange(2.), initial_state=self.state, supply_rule=object(),
            birth_to_entry_conversion=1 / 2.1)
        if conditioning != 'omit':
            kwargs['historical_conditioning'] = conditioning
        return pf.evaluate_path_at_prices(**kwargs)

    def conditioning(self, **kwargs):
        result = pf.HistoricalConditioning(2007, 1., {1: 2011, 2: 2015, 3: 2019, 4: 2023},
                                           .08, .5)
        return replace(result, **kwargs)

    def test_real_age_bridge_entry_identity_and_observer(self):
        observed = []

        def observer(period, e, P, grid, shared):
            self.events.append('observe')
            self.assertIs(e.policy.joint_choice, self.joint_marker)
            self.assertEqual(P.psi_child, .1)
            observed.append((period, e.g_pre.copy()))

        result = self.evaluate(5, self.conditioning(observer=observer))
        self.assertEqual([row['calendar_year'] for row in result.rows], [2007, 2011, 2015, 2019, 2023])
        self.assertEqual(self.events, ['evaluate', 'observe', 'advance'] * 5)
        self.assertEqual([p for p, _ in observed], list(range(5)))
        for period, due in enumerate([.1, .2, .3, .4, .1]):
            np.testing.assert_allclose(self.flows[period], [.02 + .5 * due, .06])
            row = result.rows[period]
            self.assertAlmostEqual(row['entrant_flow_next'], .08 + .5 * due)
            self.assertAlmostEqual(row['mass_accounting_residual'], 0., places=13)
            if period < 4:
                audit = row['historical_bridge_audit']
                year = 2011 + 4 * period
                self.assertEqual(audit['year'], year)
                self.assertLess(audit['maximum_absolute_target_gap'], 1e-13)
                self.assertAlmostEqual(audit['model_mass_after_bridge'],
                    row['adult_population'] - row['adult_deaths'] + row['entrant_flow_next']
                    + row['historical_bridge_net_residual'])
                expected_age = np.array([g['target_mass'] for g in audit['groups']])
                actual_age = observed[period + 1][1].sum(axis=(0, 1, 2, 4, 5, 6))
                np.testing.assert_allclose(actual_age, expected_age, atol=1e-14)
                # True bridge preserves all conditional economic states.
                np.testing.assert_allclose(observed[period + 1][1][1],
                                           3 * observed[period + 1][1][0], atol=1e-14)
        self.assertIsNone(result.rows[-1]['historical_next_bridge_year'])
        self.assertIsNone(result.rows[-1]['historical_bridge_audit'])
        self.assertLess(result.maximum_mass_accounting_error, 1e-13)
        np.testing.assert_array_equal(self.state.g_pre, self.g)
        self.assertEqual(self.state.scheduled_entries, [.1, .2, .3, .4])

    def test_default_keeps_2023_clock_closed_entry_and_unbridged_mass(self):
        implicit = self.evaluate(2)
        first_flows = [x.copy() for x in self.flows]
        explicit = self.evaluate(2, None)
        self.assertEqual(implicit.rows, explicit.rows)
        np.testing.assert_array_equal(implicit.terminal_state.g_pre, explicit.terminal_state.g_pre)
        self.assertEqual([r['calendar_year'] for r in implicit.rows], [2023, 2027])
        np.testing.assert_allclose(first_flows, [[.025, .075], [.05, .15]])
        self.assertFalse(any(k.startswith('historical_') for k in implicit.rows[0]))
        for row in implicit.rows:
            self.assertAlmostEqual(row['mass_accounting_residual'], 0., places=13)

    def test_missing_or_fabricated_target_rejected_before_backward_solve(self):
        for conditioning in (
            self.conditioning(next_age_targets={1: 2011}),
            self.conditioning(next_age_targets={1: 2011, 2: 2015, 3: 2019, 4: 2023, 5: 2027}),
            self.conditioning(next_age_targets={1: 2015, 2: 2015, 3: 2019, 4: 2023}),
        ):
            with self.subTest(conditioning=conditioning), self.assertRaises(ValueError):
                self.evaluate(5, conditioning)
        pf.backward_value_path.assert_not_called()

    def test_invalid_conditioning_bounds_and_shapes_fail_closed(self):
        for changes in ({'outside_flow': -1.}, {'retention': 1.1}, {'retention': np.nan},
                        {'initial_mass': 0.}, {'start_year': 2008}, {'observer': 4}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.evaluate(5, self.conditioning(**changes))
        with self.assertRaises(ValueError):
            self.evaluate(6, self.conditioning())
        self.P.entry_shares = np.array([1., np.nan])
        with self.assertRaises(ValueError):
            self.evaluate(5, self.conditioning())
        pf.backward_value_path.assert_not_called()


if __name__ == '__main__':
    unittest.main()
