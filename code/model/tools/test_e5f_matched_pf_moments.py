"""Pure target-observer tests: real hazard ledger/target loss, stubbed room statistics."""
import copy
from contextlib import ExitStack
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import e5f_matched_pf_moments as adapter
from intergen_eqscale_seq_optimized import solver as sequential_model

m = adapter.measurement


class HistoricalMomentTests(unittest.TestCase):
    def setUp(self):
        self.targets = m.e5_target_system_for_profile('baseline')
        self.P = SimpleNamespace(J=17, age_start=18., da=4., period_years=4.,
            n_parity=4, n_child_states=4, A_f_start=1, A_f_end=7,
            fecundity_omega1=0., readiness_gate_enabled=False,
            tfr_top_bin_weight=3.602359422009, psi_child=.2,
            tenure_choice_kappa=.005, kappa_fert=.03, kappa_fert_continuation=.04,
            joint_nested_choice=True, fertility_nest_choice=True, two_shock_choice=False,
            exhaustive_saving_control=True, beta=.98, chi=1., H0=12., theta0=.5,
            theta1=.1, hbar_child_rooms=.3, first_birth_fixed_cost=2., hbar_first_child_jump=.2)
        hazard = np.r_[np.full(7, .1), np.zeros(10)]
        self.accounting = dict(at_risk=np.ones(17), flow=hazard.copy(), hazard=hazard.copy())
        self.normalization = dict(target=2.1, completed_fertility=2.1,
                                  psi_child=.2, status='derived_intercept')
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        # The real runner configures this shared calendar module before any
        # measurement. Scope the same model binding to this fixture so the
        # real fecundity/readiness/hazard helpers run without global leakage.
        self.stack.enter_context(patch.object(m.calendar, 'model', sequential_model))
        self.calls = []

        def begin(e, P, grid, shared, *, origin_period):
            self.calls.append(('begin', origin_period, e.marker))
            return dict(origin_period=origin_period, marker=e.marker)

        def finish(branch, e, P, grid, shared, *, destination_period):
            self.calls.append(('finish', destination_period, e.marker))
            self.assertEqual(branch['origin_period'], 3)
            self.assertEqual(branch['marker'], '2019')
            self.assertEqual(e.marker, '2023')
            return dict(housing_response=.55, census_age_bridge_applied=False,
                        origin_period=3, destination_period=4)

        def cross(e, P, grid, shared, names, *, housing_increment_override):
            self.calls.append(('cross', P.psi_child, e.marker))
            self.assertEqual(e.marker, '2023')
            self.assertEqual(tuple(names), self.targets.moment_names)
            result = self.targets.targets_dict()
            result['childless_rate'] = .9**7
            result['housing_increment_0to1'] = housing_increment_override
            # These contemporaneous timing placeholders must be overwritten
            # by the real old-prehistory + dated synthetic-cohort ledger.
            result['mean_age_first_birth'] = 999.
            result['share_first_births_age30plus'] = 999.
            return result
        self.stack.enter_context(patch.object(m, 'begin_dated_first_birth_housing_branch', side_effect=begin))
        self.stack.enter_context(patch.object(m, 'finish_dated_first_birth_housing_branch', side_effect=finish))
        self.stack.enter_context(patch.object(m, 'transition_cross_section_moments', side_effect=cross))

    def observer(self, **overrides):
        args = dict(target_system=self.targets, expected_target_fingerprint=self.targets.fingerprint,
            old_parameters=self.P, old_first_birth_accounting=self.accounting,
            old_normalization=self.normalization, old_normalization_tolerance=5e-4)
        args.update(overrides)
        return adapter.HistoricalMomentObserver(**args)

    def feed(self, observer, periods=range(5)):
        for period in periods:
            P = copy.deepcopy(self.P)
            P.psi_child -= .01 * period
            g = np.zeros((1, 2, 1, 17, 1, 4, 4))
            g[0, 0, 0, :, 0, 0, 0] = 1 + period
            fp = np.zeros((1, 2, 1, 17, 1, 4))
            fp[..., 1] = .1
            e = SimpleNamespace(g_pre=g, policy=SimpleNamespace(fert_probs=fp),
                                marker=str(2007 + 4 * period))
            observer(period, e, P, np.array([0.]), SimpleNamespace())

    def test_complete_rows_real_hazard_timing_existing_loss_and_2019_branch(self):
        observer = self.observer()
        self.feed(observer)
        result = observer.result('test')
        self.assertEqual(len(result['target_fit_rows']), 12)
        self.assertEqual(set(result['moments']), set(self.targets.moment_names))
        self.assertEqual(result['loss'], self.targets.loss(result['moments']))
        self.assertEqual(result['target_fingerprint'], self.targets.fingerprint)
        self.assertEqual(result['moments']['housing_increment_0to1'], .55)
        weights = .1 * .9**np.arange(7)
        expected_mean = np.sum(np.arange(20., 48., 4.) * weights) / weights.sum()
        expected_late = weights[3:].sum() / weights.sum()
        self.assertAlmostEqual(result['moments']['mean_age_first_birth'], expected_mean)
        self.assertAlmostEqual(result['moments']['share_first_births_age30plus'], expected_late)
        self.assertEqual(self.calls[:2], [('begin', 3, '2019'), ('finish', 4, '2023')])
        self.assertEqual(self.calls[2][2], '2023')
        self.assertFalse(result['production_promoted'])
        self.assertTrue(any('parent/control' in issue for issue in result['outstanding']))
        before = copy.deepcopy(result)
        self.feed(observer, range(5, 7))
        self.assertEqual(observer.result('test'), before)
        self.assertEqual(len(self.calls), 3)

    def test_explicit_parameter_domain_and_existing_parameter_helper(self):
        observer = self.observer()
        self.feed(observer)
        theta = {key: getattr(self.P, key) for key in ('beta', 'kappa_fert',
            'kappa_fert_continuation', 'chi', 'H0', 'theta0', 'theta1',
            'hbar_child_rooms', 'first_birth_fixed_cost', 'hbar_first_child_jump')}
        domain, overrides, profile = m.activate_model_profile(m.REPAIRED_MODEL_PROFILE, theta)
        domain, _, _ = m.configure_first_child_room_jump(
            model_profile_name=m.REPAIRED_MODEL_PROFILE, theta=theta, active_domain=domain,
            profile_overrides=overrides, model_profile=profile,
            fixed_jump=None, estimate_jump=True)
        domain = tuple(row for row in domain if row[0] != 'psi_child') + (
            ('psi_child_change_2023', -1.5, .2, 'asinh'),)
        with patch.object(m, 'TRANSITION_SEARCH_DOMAIN', domain):
            result = observer.report('test', theta=theta, parameter_domain=domain,
                                     supply_rule=SimpleNamespace(elasticity=.63))
        self.assertEqual(sum(row['is_free_parameter'] for row in result['parameter_rows']), 11)
        rows = {row['parameter']: row for row in result['parameter_rows']}
        self.assertAlmostEqual(rows['psi_child_change_2023']['value'], -.04)
        self.assertEqual(rows['psi_child_2007']['value'], .2)
        self.assertAlmostEqual(rows['psi_child_2023']['value'], .16)
        self.assertEqual(rows['tenure_choice_kappa']['value'], .005)
        with patch.object(m, 'TRANSITION_SEARCH_DOMAIN', domain[:-1]), self.assertRaisesRegex(ValueError, 'configured parameter_rows'):
            observer.report('test', theta=theta, parameter_domain=domain,
                            supply_rule=SimpleNamespace(elasticity=.63))

    def test_changed_targets_or_unnormalized_seed_rejected(self):
        with self.assertRaisesRegex(ValueError, 'pinned unchanged'):
            self.observer(expected_target_fingerprint='wrong')
        changed = dict(self.normalization, completed_fertility=1.9)
        with self.assertRaisesRegex(ValueError, '2.1 normalization'):
            self.observer(old_normalization=changed)
        bad = copy.deepcopy(self.accounting)
        bad['flow'][0] = .2
        with self.assertRaisesRegex(ValueError, 'inconsistent'):
            self.observer(old_first_birth_accounting=bad)

    def test_missing_reordered_or_mixed_arm_dates_fail(self):
        observer = self.observer()
        with self.assertRaisesRegex(RuntimeError, 'incomplete'):
            observer.result('test')
        with self.assertRaisesRegex(ValueError, 'chronological'):
            self.feed(observer, [1])
        self.feed(observer, [0])
        with self.assertRaisesRegex(ValueError, 'chronological'):
            self.feed(observer, [0])
        self.P.fertility_nest_choice = False
        with self.assertRaisesRegex(ValueError, 'choice arm changed'):
            self.feed(observer, [1])

    def test_terminal_childless_identity_gate_is_not_bypassed(self):
        observer = self.observer()
        self.feed(observer, range(4))
        wrong = self.targets.targets_dict()
        with patch.object(m, 'transition_cross_section_moments', return_value=wrong):
            with self.assertRaisesRegex(RuntimeError, 'childlessness does not match'):
                self.feed(observer, [4])


if __name__ == '__main__':
    unittest.main()
