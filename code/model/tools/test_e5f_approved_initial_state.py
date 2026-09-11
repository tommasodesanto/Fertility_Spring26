"""Small accounting fixtures for the approved bridge; no household/GE solve.

The fixture supplies a stationary age/income marginal and independent birth
ledger. It is deliberately not a solved household equilibrium.
"""
from __future__ import annotations

import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest

import numpy as np

import e5f_approved_initial_state as bridge
from e5f_stationary_paygo import bind_initial_balanced_pension
from run_e5f_matched_pf_baseline import validated_preference_path


def fixture():
    P = SimpleNamespace(
        I=1, J=17, Nb=120, J_R=12, age_start=18., da=4., period_years=4.,
        tau_pay=.179, tau_H=.04, property_tax_lump_sum_transfer=0.,
        child_state_mode='independent_count', n_parity=4, n_child_states=4,
        use_stochastic_aging=True, sequential_births=True,
        joint_nested_choice=False, exhaustive_saving_control=True,
        alpha_cons=.733, preference_spec='eqscale', eqscale_form='power',
        child_room_floor=True, sigma=2., delta_alpha=0., delta_alpha_jump=0.,
        c_bar_0=0., c_bar_n=0., hbar_child_rooms=0., hbar_first_child_jump=.8,
        hR_max=9., beta=.98, psi_child=.3, xi_supply=np.array([.63]),
        H0=np.array([6.]), user_cost_rate=.2, r_bar=np.array([.1]),
        entry_shares=np.array([1.]), income_type_transition='markov',
        z_grid=np.array([.5,1.5]), z_weights=np.array([.5,.5]),
        Pi_z=np.eye(2), use_age_survival=False,
        survival_probs=np.ones(16), scale_flows_to_period=True,
        w_hat=np.array([1.]), income_age_profile=np.ones(17), pension=1.,
        retirement_income_z_scale=0., fertility_units='literal_topcode',
        tfr_top_bin_weight=3.4, entrant_conversion_factor=.5,
        _third_births_by_age=np.r_[.5/17.,np.zeros(16)],
    )
    P, _ = bind_initial_balanced_pension(P, payroll_tax=.179)
    grid = np.linspace(0.,10.,120)
    pre = np.zeros((120,1,1,17,2,4,4))
    # Every age has mass 1/17 and the same earnings marginal, with nontrivial
    # within-age wealth/income dependence that ACS reweighting must preserve.
    pre[0,0,0,:,0,0,0] = .5/17.
    pre[1,0,0,:,1,0,0] = .5/17.
    sol = SimpleNamespace(entry_rate=1/17., total_births_kfe=1.9/17.,
        entrants_mature_total=.95/17., parity_dist=np.array([.1,.4,0.,.5]),
        mean_completed_fertility=1.9)
    policy = SimpleNamespace(price=np.array([.5]))
    supply = bridge.pf.calendar.HousingSupplyRule('static-elastic',.5,6.,.63)
    packet = dict(parameters=P, b_grid=grid, solution=sol,
        evaluation=SimpleNamespace(policy=policy), shared=SimpleNamespace(),
        stationary_g_pre=pre, supply_rule=supply)
    normalization = dict(target=2.1, completed_fertility=2.1, psi_child=.3)
    return packet, normalization


def build(packet, normalization, **overrides):
    args = dict(packet=packet, normalization=normalization,
        outside_origin_entry_share=.2, preference_change_2023=-.4,
        fertility_tolerance=5e-4)
    args.update(overrides)
    return bridge.build_approved_initial_state(**args)


class ApprovedInitialBridgeTests(unittest.TestCase):
    def test_adjusted_queue_raw_queue_and_outside_flow_remain_distinct(self):
        packet, normalization = fixture()
        old = build(packet, normalization)
        np.testing.assert_allclose(old.initial_state.scheduled_entries,[1/17.]*4,
                                   rtol=0,atol=1e-16)
        np.testing.assert_allclose(old.initial_state.scheduled_raw_entries,
                                   [1.9/(2.1*17.)]*4,rtol=0,atol=1e-16)
        self.assertAlmostEqual(old.historical_conditioning.outside_flow,.2/17.)
        self.assertAlmostEqual(old.historical_conditioning.retention,.8)
        self.assertAlmostEqual(old.diagnostics['renewal']['identity_residual'],0.)
        self.assertEqual(old.historical_conditioning.next_age_targets,
                         {1:2011,2:2015,3:2019,4:2023})

    def test_acs_reweight_preserves_conditional_states_and_original_supply(self):
        packet, normalization = fixture()
        before = packet['stationary_g_pre'].copy()
        old = build(packet, normalization)
        after = old.initial_state.g_pre
        self.assertAlmostEqual(after.sum(),before.sum())
        shares = np.asarray(bridge.pf.transition.ACS_NATIONAL_HEAD_AGE_SHARES_18_85[2007])
        np.testing.assert_allclose(after.sum(axis=(0,1,2,4,5,6)),shares,rtol=0,atol=1e-15)
        for age in range(17):
            np.testing.assert_allclose(after[:,:,:,age]/after[:,:,:,age].sum(),
                before[:,:,:,age]/before[:,:,:,age].sum(),rtol=0,atol=1e-15)
        np.testing.assert_array_equal(packet['stationary_g_pre'],before)
        self.assertIs(old.supply_rule,packet['supply_rule'])
        np.testing.assert_array_equal(old.parameters.H0,packet['parameters'].H0)
        self.assertFalse(np.array_equal(after,before))

    def test_announcement_pension_is_reported_without_rebinding_old_income(self):
        packet, normalization = fixture()
        before = packet['parameters'].income.copy()
        old = build(packet, normalization)
        np.testing.assert_array_equal(old.parameters.income,before)
        accounts = old.diagnostics['announcement_state_accounts_at_old_pension']
        self.assertGreater(abs(accounts['scaled_pension_budget_residual']),1e-3)
        self.assertNotAlmostEqual(accounts['implied_balanced_pension_period'],
                                  old.parameters.pension)
        self.assertFalse(old.diagnostics['historical_equilibrium_verified'])

    def test_historical_dates_and_downstream_flat_preference_tail(self):
        packet, normalization = fixture()
        old = build(packet, normalization)
        np.testing.assert_array_equal(old.years,[2007,2011,2015,2019,2023])
        np.testing.assert_allclose(old.psi_path,[.3,.2,.1,0.,-.1],rtol=0,atol=1e-16)
        path = validated_preference_path(old.psi_path,10)
        np.testing.assert_array_equal(path[4:],np.full(6,old.psi_path[-1]))

    def test_changed_supply_or_implicit_outside_share_is_rejected(self):
        packet, normalization = fixture()
        for share in (0.,1.,np.nan):
            with self.subTest(share=share), self.assertRaises(ValueError):
                build(packet,normalization,outside_origin_entry_share=share)
        packet['supply_rule'] = replace(packet['supply_rule'],initial_stock=6.1)
        with self.assertRaisesRegex(ValueError,'supply law differs'):
            build(packet,normalization)

    def test_nonfinite_normalization_receipt_is_rejected(self):
        packet, normalization = fixture()
        for value in (np.nan,np.inf):
            broken = dict(normalization,completed_fertility=value)
            with self.subTest(value=value), self.assertRaises(ValueError):
                build(packet,broken)

    def test_normalization_receipt_must_match_actual_solution_moment(self):
        packet, normalization = fixture()
        packet['solution'] = copy.deepcopy(packet['solution'])
        packet['solution'].parity_dist = np.array([.4,.3,.2,.1])
        with self.assertRaises(ValueError):
            build(packet,normalization)


if __name__=='__main__':
    unittest.main()
