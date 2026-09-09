"""Pure joined-history tests: true calendar loops, stubbed economic/demographic kernels."""
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import run_e5f_matched_pf_history as joined
import test_e5f_pf_historical_bridge as bridge_tests

pf = joined.pf
person_pf = joined.person_pf
REAL_BACKWARD = pf.backward_value_path


class JoinedHistoryTests(unittest.TestCase):
    conditioning = bridge_tests.HistoricalBridgeTests.conditioning

    def setUp(self):
        bridge_tests.HistoricalBridgeTests.setUp(self)
        # Reuse the real backward loop. Only individual Bellman calls are
        # replaced by a transparent future-value-sensitive recursion.
        pf.backward_value_path.side_effect = REAL_BACKWARD

        def solve(**kwargs):
            P = kwargs['P']
            return SimpleNamespace(
                V=kwargs['continuation_V'] + kwargs['price'] + P.psi_child
                  + P.property_tax_lump_sum_transfer,
                price=np.array([kwargs['price']]), hR_pol=np.ones(self.shape),
                joint_choice=self.joint_marker,
            )
        pf.solve_date_policy.side_effect = solve
        self.g2023, _ = pf.transition.reweight_distribution_to_observed_age_path(
            self.g, np.arange(18., 86., 4.), year=2023, initial_mass=1.)
        heads = np.zeros((2, 101))
        target = self.g2023.sum(axis=(0, 1, 2, 4, 5, 6))
        heads[0, np.arange(18, 86, 4)] = target
        self.people = person_pf.CohortState(2023, heads * 2., heads)
        self.primitives = SimpleNamespace(initial_person_state=self.people,
            headship_rates=np.full_like(heads, .5), block_inputs=lambda *args: ({}, {}, {}))

        def couple(raw, people, **kwargs):
            next_people = person_pf.CohortState(people.year + 4,
                                                people.persons.copy(), people.heads.copy())
            person_ledger = SimpleNamespace(person_identity_max_abs=0., head_identity_max_abs=0.,
                total_net_migration=0., total_new_heads_from_nonheads=0.,
                total_head_dissolutions=0., total_net_migrant_heads=0.)
            ledger = SimpleNamespace(person=person_ledger, household_person_head_gap=0.,
                household_heads=SimpleNamespace(added_mass=np.zeros(17), removed_mass=np.zeros(17)))
            return self.g2023.copy(), next_people, ledger
        self.couple = self.stack.enter_context(patch.object(person_pf, 'advance_household_person_block', side_effect=couple))

    def run_joined(self, prices=None, observer=None):
        return joined.evaluate_history_and_person_tail(
            years=[2007, 2011, 2015, 2019, 2023, 2027],
            prices=np.ones(6) if prices is None else prices,
            psi_path=np.full(6, .1), transfer_path=np.arange(6) * .01,
            terminal_price=1., terminal_V=np.zeros(self.shape),
            base_parameters=self.P, b_grid=np.arange(2.), initial_state=self.state,
            historical_conditioning=self.conditioning(), initial_2023_persons=self.people,
            demographic_primitives=self.primitives, supply_rule=object(),
            birth_to_entry_conversion=1 / 2.1, observer=observer,
        )

    def test_future_tail_changes_2007_value_and_each_date_is_solved_twice(self):
        observations = []

        def observer(period, evaluation, P, grid, shared):
            observations.append((period, float(evaluation.policy.V.flat[0]),
                                 P.property_tax_lump_sum_transfer))
            self.assertIs(evaluation.policy.joint_choice, self.joint_marker)

        result = self.run_joined(observer=observer)
        self.assertEqual([p for p, _, _ in observations], list(range(6)))
        np.testing.assert_allclose([t for _, _, t in observations], np.arange(6) * .01)
        self.assertEqual([r['calendar_year'] for r in result.rows], [2007, 2011, 2015, 2019, 2023, 2027])
        self.assertEqual(sum(r['calendar_year'] == 2023 for r in result.rows), 1)
        self.assertEqual(result.bellman_solves, 12)
        self.assertEqual(pf.solve_date_policy.call_count, 12)
        self.assertEqual(pf.backward_value_path.call_count, 2)
        self.assertEqual(result.person_tail.bellman_solves, 2)
        self.assertEqual(len(result.values), 7)
        self.assertLess(result.initial_2023_age_head_gap, 1e-13)
        self.assertLess(result.person_tail.maximum_policy_reproduction_error, 1e-13)
        first_value = observations[0][1]
        observations.clear()
        changed = np.ones(6)
        changed[-1] = 1.03
        self.run_joined(prices=changed, observer=observer)
        self.assertAlmostEqual(observations[0][1] - first_value, .03)

    def test_2023_age_mismatch_fails_before_person_forward(self):
        changed_heads = self.people.heads.copy()
        changed_heads[0, 18] += .001
        self.people = person_pf.CohortState(2023, self.people.persons.copy(), changed_heads)
        self.primitives.initial_person_state = self.people
        with self.assertRaisesRegex(RuntimeError, '2023 household/person head-age identity'):
            self.run_joined()
        self.couple.assert_not_called()

    def test_invalid_cached_values_fail_exact_forward_replay(self):
        wrong = [np.zeros(self.shape), np.zeros(self.shape)]
        pf.backward_value_path.reset_mock()
        with self.assertRaisesRegex(RuntimeError, 'fails exact dated replay'):
            person_pf.evaluate_path_at_prices_person_demography(
                prices=[1.], psi_path=[.1], transfer_path=[0.], terminal_price=1.,
                terminal_V=np.zeros(self.shape), base_parameters=self.P,
                b_grid=np.arange(2.), initial_state=person_pf.PersonPFState(self.g2023, self.people),
                demographic_primitives=self.primitives, supply_rule=object(),
                precomputed_value_path=wrong,
            )
        pf.backward_value_path.assert_not_called()
        self.couple.assert_not_called()

    def test_old_person_evaluator_still_performs_backward_and_forward(self):
        result = person_pf.evaluate_path_at_prices_person_demography(
            prices=[1.], psi_path=[.1], transfer_path=[0.], terminal_price=1.,
            terminal_V=np.zeros(self.shape), base_parameters=self.P,
            b_grid=np.arange(2.), initial_state=person_pf.PersonPFState(self.g2023, self.people),
            demographic_primitives=self.primitives, supply_rule=object(),
        )
        self.assertEqual(result.bellman_solves, 2)
        self.assertEqual(pf.backward_value_path.call_count, 1)
        self.assertEqual(pf.solve_date_policy.call_count, 2)
        self.assertEqual(result.rows[0]['calendar_year'], 2023)


if __name__ == '__main__':
    unittest.main()
