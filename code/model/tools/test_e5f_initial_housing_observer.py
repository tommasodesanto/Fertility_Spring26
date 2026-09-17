"""Synthetic-mass observation checks; no equilibrium or regression runs."""
import copy
import json
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np

from e5f_initial_housing_observer import (
    AGE_PROJECTION, MOMENT_NAMES, observe_initial_housing_wealth,
    uniform_age_cell_overlap,
)


def fixture():
    P = SimpleNamespace(
        J=17, da=4.0, period_years=4.0, age_start=18.0, J_R=12,
        n_house=1, I=1, n_parity=4, n_child_states=4,
        H_own=np.array([6.0]), z_grid=np.array([.5, 1.5]),
        child_state_mode="independent_count", child_bin_high_cutoff=3,
        scale_flows_to_period=True, tau_pay=.2,
        income=np.array([[3.2] * 12 + [8.] * 5]),
        retirement_income_z_scale=0., property_tax_lump_sum_transfer=0.,
        use_age_survival=True, survival_probs=np.ones(17),
    )
    g = np.zeros((3, 2, 1, 17, 2, 4, 4))
    policy = SimpleNamespace(hR_pol=np.full_like(g, 4.),
                             bp_pol=np.zeros_like(g), price=np.array([2.]))
    ev = SimpleNamespace(g_current=g, g_post_fertility=g.copy(),
                         g_pre=g.copy(), policy=policy)
    return P, ev, np.array([1., 3., 9.])


def observe(P, ev, bg, **kwargs):
    defaults = dict(diagnostic_enabled=True, age_projection=AGE_PROJECTION,
                    include_wealth=False)
    defaults.update(kwargs)
    return observe_initial_housing_wealth(ev, P, bg, None, **defaults)


def row(result, name):
    return next(r for r in result['rows'] if r['moment'] == name)


class InitialHousingObserverTests(unittest.TestCase):
    def test_default_off_needs_no_model_inputs_or_imports(self):
        result = observe_initial_housing_wealth(None, None, None, None)
        self.assertEqual(result['status'], 'disabled')
        self.assertTrue(all(v is None for v in result['moments'].values()))
        self.assertFalse(result['production_eligible'])
        json.dumps(result, allow_nan=False)

    def test_age_projection_requires_explicit_selection(self):
        P, ev, bg = fixture()
        with self.assertRaisesRegex(ValueError, 'Explicit age_projection'):
            observe(P, ev, bg, age_projection=None)
        with self.assertRaises(TypeError):
            observe(P, ev, bg, diagnostic_allow_family_proxies=1)

    def test_overlap_matches_literal_annual_age_intervals(self):
        P, _, _ = fixture()
        expected = np.zeros(17); expected[14:17] = [.5, 1., .75]
        np.testing.assert_array_equal(uniform_age_cell_overlap(P, 76, 85), expected)
        expected = np.zeros(17); expected[3:9] = 1.; expected[9] = .5
        np.testing.assert_array_equal(uniform_age_cell_overlap(P, 30, 56), expected)
        expected = np.zeros(17); expected[1:5] = [.25, 1., 1., .25]
        np.testing.assert_array_equal(uniform_age_cell_overlap(P, 25, 35), expected)

    def test_caps_full_income_state_before_aggregation(self):
        P, ev, bg = fixture()
        ev.g_current[0, 0, 0, 3, :, 0, 0] = 1.
        ev.policy.hR_pol[0, 0, 0, 3, :, 0, 0] = [4., 12.]
        result = observe(P, ev, bg)
        self.assertEqual(result['moments'][MOMENT_NAMES[0]], 6.5)
        self.assertNotEqual(result['moments'][MOMENT_NAMES[0]], min(8., 9.))
        self.assertTrue(row(result, MOMENT_NAMES[0])['cap_before_income_aggregation'])

    def test_owner_rooms_use_selected_product_and_same_cap(self):
        P, ev, bg = fixture(); P.H_own[0] = 12.
        ev.g_current[0, 1, 0, 3, 0, 0, 0] = 1.
        ev.policy.hR_pol[:] = np.nan  # Irrelevant for this owner-only population.
        self.assertEqual(observe(P, ev, bg)['moments'][MOMENT_NAMES[0]], 9.)

    def test_prime_ownership_half_weights_age54_cell(self):
        P, ev, bg = fixture()
        ev.g_current[0, 0, 0, 3, 0, 0, 0] = 1.
        ev.g_current[0, 1, 0, 9, 0, 0, 0] = 1.
        result = observe(P, ev, bg)
        self.assertAlmostEqual(result['moments']['own_rate_30_55'], 1 / 3)
        self.assertEqual(row(result, 'own_rate_30_55')['denominator'], 1.5)
        self.assertIn('DUE', ' '.join(row(result, 'own_rate_30_55')['approximations']))

    def test_family_proxy_is_disabled_until_explicitly_allowed(self):
        P, ev, bg = fixture()
        ev.g_current[0, 0, 0, 3, 0, 3, 1] = 1.
        ev.g_current[1, 0, 0, 3, 0, 3, 3] = 1.
        ev.policy.hR_pol[0, 0, 0, 3, 0, 3, 1] = 2.
        ev.policy.hR_pol[1, 0, 0, 3, 0, 3, 3] = 12.
        disabled = observe(P, ev, bg)
        self.assertIsNone(disabled['moments'][MOMENT_NAMES[3]])
        result = observe(P, ev, bg, diagnostic_allow_family_proxies=True)
        # Both have lifetime parity 3, but current counts differ: 9 - 2.
        self.assertEqual(result['moments'][MOMENT_NAMES[3]], 7.)
        self.assertTrue(row(result, MOMENT_NAMES[3])['model_dependent_proxy'])
        self.assertIsNone(result['moments'][MOMENT_NAMES[4]])

    def test_shared_clock_does_not_silently_supply_count_proxy(self):
        P, ev, bg = fixture(); P.child_state_mode = 'shared_clock'
        result = observe(P, ev, bg, diagnostic_allow_family_proxies=True)
        self.assertIsNone(result['moments'][MOMENT_NAMES[3]])
        self.assertIn('independent_count', row(result, MOMENT_NAMES[3])['reason'])

    def test_empty_denominators_remain_unavailable_not_zero(self):
        P, ev, bg = fixture()
        result = observe(P, ev, bg, diagnostic_allow_family_proxies=True)
        for name in MOMENT_NAMES[:4]:
            self.assertIsNone(result['moments'][name])
        json.dumps(result, allow_nan=False)

    def test_bad_occupied_policy_or_negative_mass_is_rejected(self):
        P, ev, bg = fixture()
        ev.g_current[0, 0, 0, 3, 0, 0, 0] = 1.
        ev.policy.hR_pol[0, 0, 0, 3, 0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'Occupied renters'):
            observe(P, ev, bg)
        ev.policy.hR_pol[:] = 4.; ev.g_current[0, 0, 0, 3, 0, 0, 0] = -1.
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            observe(P, ev, bg)

    def test_existing_wealth_primitive_uses_correct_two_balance_sheets(self):
        P, ev, bg = fixture()
        ev.g_post_fertility[0, 0, 0, 0, 0, 0, 0] = 2.  # Wealth 2, annual gross earnings 1.
        ev.g_post_fertility[1, 1, 0, 15, 0, 0, 0] = 1.  # Wealth 3+2*6=15.
        # Different post-transaction tenure: old household sold its home.
        ev.g_current[1, 0, 0, 15, 0, 0, 0] = 1.
        ev.g_current[0, 1, 0, 0, 0, 0, 0] = 2.
        ev.policy.bp_pol[1, 0, 0, 15, 0, 0, 0] = 7.
        P.survival_probs[15] = .75
        result = observe(P, ev, bg, include_wealth=True)
        self.assertEqual(result['moments'][MOMENT_NAMES[5]], 17.)
        self.assertAlmostEqual(result['moments'][MOMENT_NAMES[6]], (.25 * 7 / 4) / 17)
        self.assertEqual(result['moments'][MOMENT_NAMES[8]], 7.5)
        self.assertTrue(row(result, MOMENT_NAMES[8])['pension_income_proxy'])
        # Transfer is excluded from gross labor earnings, included in income proxy.
        P.property_tax_lump_sum_transfer = 4.
        transfer = observe(P, ev, bg, include_wealth=True)
        self.assertEqual(transfer['moments'][MOMENT_NAMES[5]], 17.)
        self.assertEqual(transfer['moments'][MOMENT_NAMES[8]], 5.)

    def test_old_quantiles_use_overlap_weights_and_earnings_zero_is_unavailable(self):
        P, ev, bg = fixture(); bg[:] = [2., 4., 18.]
        for index, j in enumerate((14, 15, 16)):
            ev.g_current[index, 0, 0, j, 0, 0, 0] = 1.
        ev.g_post_fertility = ev.g_current.copy()
        result = observe(P, ev, bg, include_wealth=True)
        # Values [1,2,9], weights [.5,1,.75]: median 2, p90 9.
        self.assertEqual(result['moments'][MOMENT_NAMES[8]], 2.)
        self.assertEqual(result['moments'][MOMENT_NAMES[7]], 4.5)
        self.assertIsNone(result['moments'][MOMENT_NAMES[5]])
        self.assertEqual(row(result, MOMENT_NAMES[5])['denominator'], 0.)
        P.income[:, 14:] = 0.
        invalid = observe(P, ev, bg, include_wealth=True)
        self.assertIsNone(invalid['moments'][MOMENT_NAMES[7]])
        self.assertIsNone(invalid['moments'][MOMENT_NAMES[8]])

    def test_birth_adapter_calls_dated_pair_on_same_policy_and_preserves_continuation(self):
        P, ev, bg = fixture()
        fake = SimpleNamespace(begin_dated_first_birth_housing_branch=Mock(return_value={'pending': True}),
                               finish_dated_first_birth_housing_branch=Mock(return_value=dict(
                                   housing_response=.75, treated_mean_housing=10.75,
                                   control_mean_housing=10., treated_continuation_births=.3,
                                   origin_mass=1., destination_mass=1., census_age_bridge_applied=False)))
        with patch.dict('sys.modules', {'run_e5f_transition_calibration': fake}):
            result = observe(P, ev, bg, include_birth_response=True)
        args0 = fake.begin_dated_first_birth_housing_branch.call_args
        args1 = fake.finish_dated_first_birth_housing_branch.call_args
        self.assertIs(args0.args[0], ev); self.assertIs(args1.args[1], ev)
        self.assertEqual(args0.kwargs, {'origin_period': 0})
        self.assertEqual(args1.kwargs, {'destination_period': 1})
        birth = row(result, MOMENT_NAMES[9])
        self.assertEqual(birth['model_value'], .75)
        self.assertEqual(birth['branch']['treated_continuation_births'], .3)
        self.assertTrue(birth['uncapped_rooms'])
        self.assertTrue(birth['destination_continuation_births_allowed'])
        self.assertFalse(result['production_eligible'])

    def test_invalid_working_ages_and_balance_sheet_mass_are_rejected(self):
        P, ev, bg = fixture()
        P.J_R = 11
        with self.assertRaisesRegex(ValueError, '12 working cells'):
            observe(P, ev, bg, include_wealth=True)
        P.J_R = 12
        ev.g_current[0, 0, 0, 0, 0, 0, 0] = 1.
        with self.assertRaisesRegex(ValueError, 'same living mass'):
            observe(P, ev, bg, include_wealth=True)

    def test_observer_does_not_mutate_supplied_arrays(self):
        P, ev, bg = fixture()
        ev.g_current[0, 0, 0, 3, 0, 0, 0] = 1.
        ev.g_post_fertility = ev.g_current.copy()
        original = copy.deepcopy(ev)
        observe(P, ev, bg, include_wealth=True)
        for name in ('g_pre', 'g_post_fertility', 'g_current'):
            np.testing.assert_array_equal(getattr(ev, name), getattr(original, name))
        np.testing.assert_array_equal(ev.policy.hR_pol, original.policy.hR_pol)
        np.testing.assert_array_equal(ev.policy.bp_pol, original.policy.bp_pol)


if __name__ == '__main__':
    unittest.main()
