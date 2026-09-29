"""Small Torch-only checks run AFTER installing the isolated source overlay.

No household solve, hash, source overlay installation, or rendering is performed
here. Synthetic transport is identity so the tests isolate fertility accounting,
policy cache ownership, event-branch selection, and target observation.
"""
from __future__ import annotations

import copy
import sys
import unittest
from contextlib import ExitStack
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


class ReportingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if sys.platform != "linux":
            raise RuntimeError("Reporting numerical checks must run on Torch")
        import run_dynamic_population_transition as calendar
        import run_e5f_open_population_transition as transition
        import run_e5f_transition_calibration as measurement
        import e5f_initial_fertility_observer as fertility_observer
        import e5f_recent_parent_flow_observer as recent_observer
        from intergen_eqscale_seq_optimized import solver as model
        cls.calendar = calendar
        cls.transition = transition
        cls.measurement = measurement
        cls.fertility_observer = fertility_observer
        cls.recent_observer = recent_observer
        cls.model = model

    def setUp(self):
        self.stack = ExitStack()
        self.addCleanup(self.stack.close)
        model = self.model
        self.stack.enter_context(patch.object(self.calendar, "model", model))
        self.stack.enter_context(patch.multiple(
            model,
            independent_child_maturation_active=lambda P: True,
            readiness_gate_active=lambda P: False,
            parent_age_maturation_active=lambda P: False,
            readiness_settled_state=lambda P: 0,
            readiness_childless_states=lambda P: (0,),
            birth_destination_child_state=lambda P, m: m + 1,
            get_fecundity_by_age=lambda P: np.r_[np.full(7, .8), np.zeros(10)],
            realize_current_cross_section=lambda mass, *args, **kwargs: mass.copy(),
            advance_cohort_one_period_markov_income=lambda mass, *args, **kwargs: mass.copy(),
            income_transition_values=lambda P: (np.ones(1), np.ones(1), np.ones((1, 1))),
        ))
        self.P = SimpleNamespace(
            J=17, I=1, n_house=1, n_parity=4, n_child_states=4,
            A_f_start=1, A_f_end=7, age_start=18., da=4., period_years=4.,
            child_state_mode="independent_count", sequential_births=True,
            fertility_units="literal_topcode", tfr_top_bin_weight=3.6,
            two_births_per_period=True, use_stochastic_aging=False,
            use_age_survival=False,
        )
        shape = (2, 2, 1, 17, 1, 4, 4)
        pre = np.zeros(shape)
        for n, m, w in ((0, 0, 3.), (1, 0, 2.), (1, 1, 4.), (2, 0, 5.),
                        (2, 2, 6.), (3, 0, 7.), (3, 3, 8.)):
            pre[..., n, m] = w
        pre /= pre.sum()
        first = np.zeros(shape[:5] + (4,))
        first[..., 0] = .4
        first[..., 1] = .6
        continuation = np.empty(shape[:5] + (2, 2, 4))
        continuation[..., 0, :, :] = .65
        continuation[..., 1, :, :] = .35
        extra = np.empty_like(continuation)
        extra[..., 0, :, :] = .3
        extra[..., 1, :, :] = .7
        maps = self.calendar.TransitionMaps(
            np.zeros((1, 2, 2), dtype=int), np.zeros((1, 2, 2)),
            np.zeros((1, 2, 2, 4, 4, 2), dtype=int), np.zeros((1, 2, 2, 4, 4, 2)))
        self.policy = self.calendar.PolicyBundle(
            V=np.ones(shape), c_pol=np.ones(shape), hR_pol=np.ones(shape),
            bp_pol=np.zeros(shape), tenure_choice=np.zeros(shape, dtype=int),
            tenure_probs=None, loc_probs=np.ones((2, 2, 1, 1, 17, 1, 4, 4)),
            fert_probs=first, fert_value=np.zeros_like(first),
            price=np.ones(1), maps=maps, fert2_probs=continuation,
            fert_extra_probs=extra,
        )
        self.pre = pre

    def evaluate(self, P=None, policy=None):
        P = P or self.P
        policy = policy or self.policy
        post, births, by_location = self.transition.apply_sequential_fertility(
            self.pre, policy.fert_probs, P, policy.fert2_probs,
            fert_extra_probs=self.calendar.policy_extra_birth_probs(policy, P))
        return SimpleNamespace(policy=policy, g_pre=self.pre,
                               g_post_fertility=post, g_current=post.copy(),
                               births=births, births_by_loc=by_location,
                               feasibility_projection_mass=0.)

    def test_owned_cache_cannot_follow_parameter_mutation(self):
        source = self.policy.fert_extra_probs.copy()
        policy = copy.copy(self.policy)
        policy.fert_extra_probs = source
        policy.__post_init__()
        source[...] = 0.
        self.P._fert_extra_probs = source
        np.testing.assert_array_equal(self.calendar.policy_extra_birth_probs(policy, self.P),
                                      self.policy.fert_extra_probs)
        policy.fert_extra_probs = None
        with self.assertRaisesRegex(RuntimeError, "owned"):
            self.calendar.policy_extra_birth_probs(policy, self.P)

    def test_requires_explicit_extra_array(self):
        with self.assertRaisesRegex(RuntimeError, "explicit policy-owned"):
            self.transition.apply_sequential_fertility(
                self.pre, self.policy.fert_probs, self.P, self.policy.fert2_probs)

    def test_birth_counts_and_risk_sets(self):
        ev = self.evaluate()
        d = self.measurement.period_fertility_diagnostics(ev, self.P)
        pre_n = self.pre.sum(axis=(0, 1, 2, 4, 6))
        post_n = ev.g_post_fertility.sum(axis=(0, 1, 2, 4, 6))
        flows = np.column_stack([d[k] for k in
            ("birth_flow_first", "birth_flow_second", "birth_flow_third_bin_entry")])
        np.testing.assert_allclose(np.cumsum(pre_n-post_n, axis=1)[:, :3], flows, atol=1e-15)
        self.assertAlmostEqual(float(flows.sum()), ev.births, places=14)
        self.assertAlmostEqual(float(ev.births_by_loc.sum()), ev.births, places=14)
        self.assertAlmostEqual(float(ev.g_post_fertility.sum()), float(self.pre.sum()), places=14)
        self.assertTrue(np.all(flows <= d["birth_order_risk"] + 1e-15))
        # Exact source risk: third-order extra pool is original n=1 successes,
        # never the new n=2 families that started this period childless.
        expected_second_risk = pre_n[:7, 1] + .8 * .6 * pre_n[:7, 0]
        expected_third_risk = pre_n[:7, 2] + .8 * .35 * pre_n[:7, 1]
        np.testing.assert_allclose(d["birth_order_risk"][:7, 1], expected_second_risk, atol=1e-15)
        np.testing.assert_allclose(d["birth_order_risk"][:7, 2], expected_third_risk, atol=1e-15)

    def test_age25_uses_retained_interpolation_and_reports_approximation(self):
        ev = self.evaluate()
        observed = self.fertility_observer.observe_initial_fertility(
            ev, self.P, age_projection="uniform_birth_time")
        before = self.pre[:, :, :, 1].sum(axis=(0, 1, 2, 3, 5))
        after = ev.g_post_fertility[:, :, :, 1].sum(axis=(0, 1, 2, 3, 5))
        projected = .125*before + .875*after
        expected = float(projected @ np.arange(4) / projected.sum())
        self.assertAlmostEqual(observed["moments"][self.fertility_observer.EARLY_FERTILITY_MOMENT],
                               expected, places=14)
        self.assertIn("common-event", observed["metadata"]["birth_time_assumption"])
        self.assertIn("ordered", observed["metadata"]["early_fertility_age25"]
                      ["two_birth_projection_convention"])

    def test_recent_parent_counts_households_and_tags_original_first_births(self):
        ev = self.evaluate()
        observer = self.recent_observer
        observed = observer.observe_recent_parent_flow(
            ev, self.P, diagnostic_enabled=True, snapshot=observer.SNAPSHOT,
            age_projection=observer.AGE_PROJECTION, diagnostic_allow_residence_proxy=True)
        birth_hh = observed["groups"]["selected_birth"]["unweighted_mass"]
        first_hh = observed["groups"]["first_birth"]["unweighted_mass"]
        later_hh = observed["groups"]["continuation_birth"]["unweighted_mass"]
        self.assertAlmostEqual(first_hh+later_hh, birth_hh, places=14)
        original_childless = float(self.pre[:, :, :, :7, :, 0, 0].sum())
        self.assertAlmostEqual(first_hh, original_childless*.6*.8, places=14)
        self.assertGreater(observed["accounting"]["empty_home_births"], birth_hh)
        self.assertAlmostEqual(observed["accounting"]["empty_home_birth_households"], birth_hh, places=14)

    def test_dated_first_birth_housing_cohort_contains_actual_two_child_destinations(self):
        ev = self.evaluate()
        branch = self.measurement.begin_dated_first_birth_housing_branch(
            ev, self.P, np.arange(2.), SimpleNamespace(), origin_period=0)
        treated = branch["treated_next_pre"]
        original = float(self.pre[:, :, :, :7, :, 0, 0].sum())
        successes = original*.6*.8
        self.assertAlmostEqual(float(treated.sum()), successes, places=14)
        self.assertAlmostEqual(float(treated[..., 2, 2].sum()), successes*.7*.8, places=14)
        self.assertAlmostEqual(float(treated[..., 1, 1].sum()), successes*(1.-.7*.8), places=14)
        self.assertAlmostEqual(float(branch["control_next_pre"][..., 0, 0].sum()), successes, places=14)

    def test_flag_off_ignores_unrelated_extra_cache(self):
        P = copy.copy(self.P)
        P.two_births_per_period = False
        one = self.evaluate(P)
        poisoned = copy.copy(self.policy)
        poisoned.fert_extra_probs = np.full_like(self.policy.fert_extra_probs, np.nan)
        self.assertIsNone(self.calendar.policy_extra_birth_probs(poisoned, P))
        again = self.evaluate(P, poisoned)
        np.testing.assert_array_equal(one.g_post_fertility, again.g_post_fertility)
        self.assertEqual(one.births, again.births)
        # Independent one-birth expected operator, original source pools only.
        expected = self.pre.copy()
        births = 0.
        for j in range(7):
            for n in range(3):
                for m in range(n+1):
                    s = self.pre[:, :, :, j, :, n, m] * .8 * (.6 if n == 0 else .35)
                    expected[:, :, :, j, :, n, m] -= s
                    expected[:, :, :, j, :, n+1, m+1] += s
                    births += float(s.sum())
        np.testing.assert_allclose(one.g_post_fertility, expected, atol=1e-16)
        self.assertAlmostEqual(one.births, births, places=14)


if __name__ == "__main__":
    unittest.main()
