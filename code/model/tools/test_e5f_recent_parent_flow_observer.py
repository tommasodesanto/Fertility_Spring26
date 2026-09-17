"""Small actual-kernel tests; no Bellman, stationary or equilibrium solve."""
import copy
import json
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import e5f_recent_parent_flow_observer as observer
import run_e5f_open_population_transition as transition
from intergen_eqscale_seq_optimized import solver as model


def refresh(ev, P):
    ev.g_post_fertility, ev.births, _ = transition.apply_sequential_fertility(
        ev.g_pre, ev.policy.fert_probs, P, ev.policy.fert2_probs)
    p = ev.policy
    ev.g_current = model.realize_current_cross_section(
        ev.g_post_fertility, p.loc_probs, p.tenure_choice, p.tenure_probs,
        p.maps.lmm_idx, p.maps.lmm_wt, p.maps.tmx_idx, p.maps.tmx_wt)


def fixture(*, readiness=True, locations=1):
    P = SimpleNamespace(J=17, da=4., period_years=4., age_start=18.,
        n_parity=4, n_child_states=4, n_house=1, I=locations,
        A_f_start=1, A_f_end=7, sequential_births=True,
        fertility_units="literal_topcode", child_state_mode="independent_count",
        readiness_gate_enabled=readiness, fecundity_omega1=0.,
        tfr_top_bin_weight=3.602359422009, use_numba_scatter=False)
    shape = (2, 2, locations, 17, 1, 4, 4)
    pre = np.zeros(shape)
    ready = int(readiness)
    pre[0, 0, 0, 3, 0, 0, ready] = .4
    if readiness:
        pre[0, 0, 0, 3, 0, 0, 0] = .2
    pre[0, 0, 0, 3, 0, 1, 0] = .1
    pre[0, 0, 0, 3, 0, 1, 1] = .3
    fp = np.zeros(shape[:-1]); fp[..., 1] = .25
    cp = np.zeros(shape[:5] + (2, 2, 4))
    cp[..., 1, 0, 0] = .5
    cp[..., 1, 0, 1] = 1.  # Existing parent can birth; never joins selected flow.
    cp[..., 1, 1, :] = 1.  # Cannot turn a new second birth into a third now.
    lp = np.zeros((2, 2, locations, locations, 17, 1, 4, 4))
    for i in range(locations):
        lp[:, :, i, i] = 1.
    tp = np.zeros(shape + (2,)); tp[..., 0] = .6; tp[..., 1] = .4
    tp[..., 1, 1, :] = [.2, .8]
    tp[..., 2, 1, :] = [.4, .6]
    lm_shape = (locations, 2, 2)
    tm_shape = (locations, 2, 2, 4, 4, 2)
    maps = SimpleNamespace(lmm_idx=np.zeros(lm_shape, dtype=np.int64),
        lmm_wt=np.broadcast_to([0., 1.], lm_shape).copy(),
        tmx_idx=np.zeros(tm_shape, dtype=np.int64),
        tmx_wt=np.broadcast_to([0., 1.], tm_shape).copy())
    p = SimpleNamespace(fert_probs=fp, fert2_probs=cp, loc_probs=lp,
        tenure_choice=np.zeros(shape, dtype=np.int16), tenure_probs=tp,
        maps=maps, price=np.ones(locations))
    ev = SimpleNamespace(g_pre=pre, policy=p, feasibility_projection_mass=.007)
    refresh(ev, P)
    return ev, P


class RecentParentFlowTests(unittest.TestCase):
    def setUp(self):
        binding = patch.object(transition.calendar, "model", model)
        binding.start()
        self.addCleanup(binding.stop)

    def observe(self, ev, P, **kwargs):
        options = dict(diagnostic_enabled=True, snapshot=observer.SNAPSHOT,
            age_projection=observer.AGE_PROJECTION, diagnostic_allow_residence_proxy=True)
        options.update(kwargs)
        return observer.observe_recent_parent_flow(ev, P, **options)

    def test_default_off_does_not_need_inputs_or_kernels(self):
        with patch.object(transition, "apply_sequential_fertility") as apply:
            result = observer.observe_recent_parent_flow(None, None)
        apply.assert_not_called()
        self.assertEqual(result["status"], "disabled")
        self.assertIsNone(result["model_value"])
        self.assertFalse(result["production_eligible"])
        json.dumps(result, allow_nan=False)

    def test_all_interpretation_opt_ins_are_required(self):
        ev, P = fixture()
        for kwargs, error in ((dict(snapshot=None), "Explicit snapshot"),
                              (dict(snapshot="uniform_birth_time"), "Explicit snapshot"),
                              (dict(age_projection=None), "Explicit age_projection"),
                              (dict(diagnostic_allow_residence_proxy=False), "Requires diagnostic")):
            with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, error):
                self.observe(ev, P, **kwargs)
        with self.assertRaises(TypeError):
            self.observe(ev, P, diagnostic_enabled=1)

    def test_hand_worked_actual_birth_selection_and_control(self):
        ev, P = fixture()
        with patch.object(transition, "apply_sequential_fertility",
                          wraps=transition.apply_sequential_fertility) as apply:
            result = self.observe(ev, P)
        self.assertEqual(apply.call_count, 2)  # Full replay and empty submass.
        groups = result["groups"]
        expected = {"selected_birth": (.11, .15), "current_empty": (.22, .55),
                    "first_birth": (.08, .1), "continuation_birth": (.03, .05),
                    "empty_never_parent": (.2, .5), "empty_former_parent": (.02, .05)}
        for name, (owners, mass) in expected.items():
            self.assertAlmostEqual(groups[name]["owner_numerator"], owners)
            self.assertAlmostEqual(groups[name]["denominator"], mass)
        self.assertAlmostEqual(result["model_value"], 1 / 3)
        self.assertAlmostEqual(result["accounting"]["all_births"], .45)
        self.assertAlmostEqual(result["accounting"]["empty_home_births"], .15)
        self.assertEqual(result["accounting"]["feasibility_projection_mass"], .007)
        for name, value in result["accounting"].items():
            if "error" in name or "excess" in name:
                self.assertLessEqual(value, observer.MASS_ATOL)
        self.assertIsNone(result["actual_weight"])
        self.assertIsNone(result["loss_contribution"])
        self.assertFalse(result["target_contract_activated"])
        json.dumps(result, allow_nan=False)

    def test_readiness_disabled_and_top_parity_former_parent(self):
        ev, P = fixture(readiness=False)
        ev.g_pre[0, 0, 0, 3, 0, 3, 0] = .2
        refresh(ev, P)
        result = self.observe(ev, P)
        self.assertEqual(result["metadata"]["childless_readiness_states"], [0])
        self.assertAlmostEqual(result["groups"]["current_empty"]["denominator"], .55)
        self.assertAlmostEqual(result["groups"]["empty_former_parent"]["denominator"], .25)
        self.assertAlmostEqual(result["groups"]["selected_birth"]["denominator"], .15)
        # Top-bin population representative is irrelevant to household flow.
        P.tfr_top_bin_weight = 20.
        self.assertEqual(self.observe(ev, P)["model_value"], result["model_value"])

    def test_birth_kernel_uses_original_risk_sets_and_fecundity(self):
        ev, P = fixture()
        P.fecundity_omega1 = .5
        refresh(ev, P)
        result = self.observe(ev, P)
        self.assertAlmostEqual(result["groups"]["first_birth"]["denominator"], .05)
        self.assertAlmostEqual(result["groups"]["continuation_birth"]["denominator"], .025)
        self.assertEqual(float(ev.g_post_fertility[..., 3, 1].sum()), 0.)

    def test_moving_owner_uses_destination_renter_wealth_and_tenure_probabilities(self):
        ev, P = fixture(locations=2)
        ev.g_pre[:] = 0.
        # Four origin owners: one realized first birth, three empty controls.
        ev.g_pre[1, 1, 0, 3, 0, 0, 1] = 4.
        ev.policy.loc_probs[:, :, 0, 0] = 0.
        ev.policy.loc_probs[:, :, 0, 1] = 1.
        # Location map moves origin wealth node 1 to destination wealth node 0.
        ev.policy.maps.lmm_wt[0, 1] = 0.
        ev.policy.tenure_probs[..., 0] = 1.
        ev.policy.tenure_probs[..., 1] = 0.
        # Normalize [.1,.3] exactly as the actual tenure kernel: owner share .75.
        ev.policy.tenure_probs[0, 0, 1, 3, 0, 1, 1] = [.1, .3]
        refresh(ev, P)
        result = self.observe(ev, P)
        self.assertAlmostEqual(result["groups"]["selected_birth"]["denominator"], 1.)
        self.assertAlmostEqual(result["groups"]["selected_birth"]["owner_numerator"], .75)
        self.assertEqual(result["groups"]["current_empty"]["owner_numerator"], 0.)
        self.assertAlmostEqual(result["model_value"], .75)

    def test_deterministic_tenure_path_uses_realized_choice(self):
        ev, P = fixture()
        ev.policy.tenure_probs = None
        ev.policy.tenure_choice[..., 1, 1] = 1
        ev.policy.tenure_choice[..., 2, 1] = 1
        refresh(ev, P)
        self.assertAlmostEqual(self.observe(ev, P)["model_value"], 1.)

    def test_age_overlap_and_birth_age_cutoff_preserve_group_composition(self):
        ev, P = fixture()
        ev.g_pre[0, 0, 0, 9, 0, 0, 1] = 2.  # Age54: half weight, no birth.
        ev.policy.tenure_probs[:, :, :, 9, :, 0, 1] = [0., 1.]
        refresh(ev, P)
        result = self.observe(ev, P)
        weights = np.zeros(17); weights[3:9] = 1.; weights[9] = .5
        np.testing.assert_array_equal(result["metadata"]["age_overlap_weights"], weights)
        self.assertAlmostEqual(result["groups"]["selected_birth"]["denominator"], .15)
        self.assertAlmostEqual(result["groups"]["current_empty"]["denominator"], 1.55)
        self.assertAlmostEqual(result["groups"]["current_empty"]["owner_numerator"], 1.22)
        self.assertEqual(result["groups"]["selected_birth"]["mass_by_age"][9], 0.)

    def test_replays_reject_stale_post_fertility_and_current_arrays(self):
        for field, error in (("g_post_fertility", "fertility_replay"),
                             ("g_current", "full current replay")):
            ev, P = fixture()
            getattr(ev, field)[0, 0, 0, 3, 0, 0, 1] += .01
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, error):
                self.observe(ev, P)

    def test_zero_main_denominators_fail_without_floors(self):
        ev, P = fixture()
        ev.policy.fert_probs[:] = 0.; ev.policy.fert2_probs[:] = 0.
        refresh(ev, P)
        with self.assertRaisesRegex(ValueError, "selected_birth denominator"):
            self.observe(ev, P)
        ev, P = fixture(); ev.g_pre[:] = 0.
        ev.g_pre[0, 0, 0, 3, 0, 0, 1] = 1.
        ev.policy.fert_probs[..., 1] = 1.
        refresh(ev, P)
        with self.assertRaisesRegex(ValueError, "current_empty denominator"):
            self.observe(ev, P)

    def test_invalid_family_mass_nan_and_probability_fail(self):
        for field in ("g_pre", "g_post_fertility", "g_current"):
            ev, P = fixture()
            getattr(ev, field)[0, 0, 0, 3, 0, 0, 2] = .01
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "invalid family-state"):
                self.observe(ev, P)
        for field in ("g_pre", "g_post_fertility", "g_current"):
            ev, P = fixture(); getattr(ev, field).flat[0] = np.nan
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.observe(ev, P)
        ev, P = fixture(); ev.policy.fert_probs.flat[0] = 1.1
        with self.assertRaisesRegex(ValueError, "probability above one"):
            self.observe(ev, P)
        ev, P = fixture(); ev.policy.tenure_probs.flat[0] = np.nan
        with self.assertRaises(ValueError):
            self.observe(ev, P)
        ev, P = fixture(); ev.policy.loc_probs[0, 0, 0, 0, 3, 0, 0, 1] = .5
        with self.assertRaisesRegex(ValueError, "location probability sum"):
            self.observe(ev, P)

    def test_missing_owned_continuation_and_unsupported_architecture_fail(self):
        ev, P = fixture(); ev.policy.fert2_probs = None
        with self.assertRaisesRegex(RuntimeError, "lacks continuation"):
            self.observe(ev, P)
        for key in ("joint_nested_choice", "fertility_nest_choice", "two_shock_choice"):
            ev, P = fixture(); setattr(P, key, True)
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Unsupported architecture"):
                self.observe(ev, P)
        ev, P = fixture()
        with patch.object(transition.calendar, "model", None):
            with self.assertRaisesRegex(ValueError, "Caller must configure"):
                self.observe(ev, P)
            self.assertIsNone(transition.calendar.model)

    def test_inputs_unchanged_and_provenance_detached(self):
        ev, P = fixture(); before = copy.deepcopy(ev); p_before = copy.deepcopy(vars(P))
        receipt = {"checkpoint": "synthetic", "source_pin": ["abc"]}
        result = self.observe(ev, P, input_provenance=receipt)
        for name in ("g_pre", "g_post_fertility", "g_current"):
            np.testing.assert_array_equal(getattr(ev, name), getattr(before, name))
        for name in ("fert_probs", "fert2_probs", "loc_probs", "tenure_probs", "tenure_choice", "price"):
            np.testing.assert_array_equal(getattr(ev.policy, name), getattr(before.policy, name))
        for name in ("lmm_idx", "lmm_wt", "tmx_idx", "tmx_wt"):
            np.testing.assert_array_equal(getattr(ev.policy.maps, name), getattr(before.policy.maps, name))
        self.assertEqual(vars(P), p_before)
        self.assertIs(transition.calendar.model, model)
        receipt["source_pin"].append("changed")
        self.assertEqual(result["metadata"]["policy_input_provenance"]["source_pin"], ["abc"])
        self.assertFalse(result["metadata"]["provenance_independently_certified"])
        self.assertTrue(any("annual bridge" in w for w in result["metadata"]["warnings"]))


if __name__ == "__main__":
    unittest.main()
