"""Queue and observer integration checks; no household solve is performed."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import run_e5f_perfect_foresight_transition as pf


class CurrentPFIntegration(unittest.TestCase):
    def setUp(self):
        self.P = SimpleNamespace(adult_entry_clock='split_birth_vintage', I=1,
                                 entry_shares=np.ones(1))
        self.g = np.ones((1, 1, 1, 1, 1, 1, 1))

    def test_actual_entry_prehistory_and_raw_are_separate(self):
        state = pf.stationary_initial_state(self.g, 4., 6.3, self.P, 1/2.1)
        np.testing.assert_array_equal(state.g_pre, self.g)
        self.assertIsNot(state.g_pre, self.g)
        self.assertEqual(state.scheduled_entries.due_in_16, (2.,)*3)
        self.assertEqual(state.scheduled_entries.due_in_20, (2.,)*4)
        self.assertAlmostEqual(state.scheduled_raw_entries.due_in_16[0], 1.5)
        queue = state.scheduled_entries
        for _ in range(10):
            due, queue = pf.transition.advance_adult_entry_clock(queue, 8.4, 1/2.1, 'split-16-20')
            self.assertEqual(due, 4.)

    def test_impulse_adjusted_and_raw(self):
        for births in (2.1, 4.2):
            state = pf.stationary_initial_state(self.g, 0, 0, self.P, 1/2.1)
            queue = state.scheduled_entries
            flows = []
            for t in range(6):
                due, queue = pf.transition.advance_adult_entry_clock(queue, births if t == 0 else 0,
                                                                     1/2.1, pf.entry_clock_timing(self.P))
                flows.append(due)
            np.testing.assert_array_equal(flows, [0, 0, 0, births/4.2, births/4.2, 0])

    def test_copy_and_serialization_preserve_two_lags(self):
        state = pf.stationary_initial_state(self.g, 4, 6.3, self.P, 1/2.1)
        queue = pf.copy_birth_queue(state.scheduled_entries)
        self.assertIsNot(queue, state.scheduled_entries)
        self.assertEqual(pf.jsonable(queue), {'timing':'split-16-20',
                                             'due_in_16':[2.]*3, 'due_in_20':[2.]*4})
        pf.validate_entry_queues(state, self.P)
        state.scheduled_entries = [2.]*7
        with self.assertRaises(TypeError):
            pf.validate_entry_queues(state, self.P)

    def test_legacy_unchanged(self):
        P = SimpleNamespace()
        state = pf.stationary_initial_state(self.g, 3., 8., P, .5)
        self.assertEqual(state.scheduled_entries, [3.]*4)
        self.assertEqual(state.scheduled_raw_entries, [4.]*4)
        due, queue = pf.transition.advance_adult_entry_clock(state.scheduled_entries, 10, .5,
                                                           pf.entry_clock_timing(P))
        self.assertEqual(due, 3.)
        self.assertEqual(queue, [3.,3.,3.,5.])
        with self.assertRaises(ValueError):
            pf.stationary_initial_state(self.g, 3, 8, self.P, .5)

    def test_observer_receives_computed_next_cohort(self):
        state = pf.stationary_initial_state(self.g, 4., 6.3, self.P, 1/2.1)
        policy = SimpleNamespace(V=np.ones(1))
        evaluation = SimpleNamespace(g_pre=self.g, g_post_fertility=self.g, births=6.3)
        received = []
        class Observed(Exception):
            pass
        def observer(t, ev, P, grid, shared, cohort):
            received.append((t, ev, cohort.copy()))
            raise Observed()
        with patch.object(pf.social_security, 'validated_fiscal_paths', return_value=(None,None)), \
             patch.object(pf.social_security, 'apply_fiscal_date'), \
             patch.object(pf, 'rents_from_asset_prices', return_value=np.ones(1)), \
             patch.object(pf, 'backward_value_path', return_value=([np.ones(1)]*2,1)), \
             patch.object(pf.calendar.model, 'precompute_shared', return_value=None), \
             patch.object(pf, 'solve_date_policy', return_value=policy), \
             patch.object(pf.calendar, 'evaluate_period', return_value=evaluation), \
             patch.object(pf.transition, 'calendar_topcode_birth_accounting', return_value={'topcode_adjusted_birth_children':8.4}), \
             patch.object(pf.transition, 'advance_sequential_calendar_distribution', return_value=(self.g.copy(),None,0,None)), \
             patch.object(pf.calendar, 'entrant_cohort', side_effect=lambda flows,P,grid:np.ones((1,)*6)*flows[0]):
            with self.assertRaises(Observed):
                pf.evaluate_path_at_prices(prices=[1],psi_path=[1],terminal_price=1,
                    terminal_V=np.ones(1),base_parameters=self.P,b_grid=np.ones(1),
                    initial_state=state,supply_rule=None,birth_to_entry_conversion=1/2.1,
                    dated_observer=observer)
        self.assertEqual(received[0][0], 0)
        self.assertIs(received[0][1], evaluation)
        self.assertEqual(float(received[0][2].sum()), 4.)
        np.testing.assert_array_equal(state.g_pre, self.g)

if __name__ == '__main__':
    unittest.main()
