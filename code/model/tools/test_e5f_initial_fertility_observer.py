"""Pure initial-fertility observation tests; no household or equilibrium solve."""
import copy
import json
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import e5f_initial_fertility_observer as observer
import run_e5f_transition_calibration as measurement
from intergen_eqscale_seq_optimized import solver as sequential_model


def fixture(*, masses=None, parity=None, probabilities=None):
    J = 7
    masses = np.ones(J) if masses is None else np.asarray(masses, dtype=float)
    parity = np.tile([.4, .3, .2, .1], (J, 1)) if parity is None else np.asarray(parity, dtype=float)
    probabilities = np.full((J, 3), .1) if probabilities is None else np.asarray(probabilities, dtype=float)
    P = SimpleNamespace(J=J, age_start=18, da=4., period_years=4., n_parity=4,
        n_child_states=4, sequential_births=True, fertility_units="literal_topcode",
        child_state_mode="independent_count", A_f_start=1, A_f_end=7,
        tfr_top_bin_weight=3.602359422009, fecundity_omega1=0., psi_child=.137)
    shape = (1, 1, 1, J, 1, 4, 4)
    pre = np.zeros(shape)
    pre[0, 0, 0, :, 0, :, 0] = masses[:, None] * parity
    post = pre.copy()
    first = np.zeros((1, 1, 1, J, 1, 4))
    second = np.zeros((1, 1, 1, J, 1, 2, 2, 4))
    first[0, 0, 0, :, 0, 1] = probabilities[:, 0]
    second[0, 0, 0, :, 0, 1, 0, 0] = probabilities[:, 1]
    second[0, 0, 0, :, 0, 1, 1, 0] = probabilities[:, 2]
    flows = masses[:, None] * parity[:, :3] * probabilities
    # Hand-apply each original risk set once; births cannot chain within a cell.
    for j in range(J):
        for n in range(3):
            post[0, 0, 0, j, 0, n, 0] -= flows[j, n]
            post[0, 0, 0, j, 0, n + 1, 1] += flows[j, n]
    ev = SimpleNamespace(g_pre=pre, g_post_fertility=post, births=float(flows.sum()),
        policy=SimpleNamespace(fert_probs=first, fert2_probs=second))
    return ev, P


class InitialFertilityObserverTests(unittest.TestCase):
    def setUp(self):
        # Select the real pure sequential accounting primitives without calling
        # the production configure function, loading a chain, or solving a model.
        self.binding = patch.object(measurement.calendar, "model", sequential_model)
        self.binding.start()
        self.addCleanup(self.binding.stop)

    def observe(self, ev, P, projection="uniform_birth_time"):
        return observer.observe_initial_fertility(ev, P, age_projection=projection)

    def test_projection_is_required_and_has_no_implicit_default(self):
        ev, P = fixture()
        with self.assertRaises(TypeError):
            observer.observe_initial_fertility(ev, P)
        with self.assertRaisesRegex(ValueError, "Explicit age_projection"):
            self.observe(ev, P, "nearest_age42")

    def test_real_flow_helpers_and_nonmutation(self):
        ev, P = fixture(masses=[1., 2., 3., 4., 5., 6., 7.])
        before = copy.deepcopy(ev)
        parameters_before = vars(P).copy()
        with patch.object(measurement, "first_birth_accounting_by_age",
                wraps=measurement.first_birth_accounting_by_age) as first, patch.object(
                measurement, "period_fertility_diagnostics",
                wraps=measurement.period_fertility_diagnostics) as period:
            packet = self.observe(ev, P)
            self.assertGreaterEqual(first.call_count, 1)
            period.assert_called_once_with(ev, P)
        np.testing.assert_array_equal(ev.g_pre, before.g_pre)
        np.testing.assert_array_equal(ev.g_post_fertility, before.g_post_fertility)
        np.testing.assert_array_equal(ev.policy.fert_probs, before.policy.fert_probs)
        np.testing.assert_array_equal(ev.policy.fert2_probs, before.policy.fert2_probs)
        self.assertEqual(vars(P), parameters_before)
        self.assertFalse(packet['metadata']['production_smm_eligible'])
        self.assertFalse(packet['metadata']['weights_or_standard_errors_adopted'])
        self.assertAlmostEqual(sum(packet['parity_shares_40_44'].values()), 1.)
        json.dumps(packet, allow_nan=False)

    def test_constant_parity_and_no_birth_in_window_are_projection_invariant(self):
        probabilities = np.zeros((7, 3)); probabilities[0, 0] = .2
        ev, P = fixture(masses=[1, 2, 3, 4, 5, 6, 7], probabilities=probabilities)
        a = self.observe(ev, P)
        b = self.observe(ev, P, "constant_post_cell")
        self.assertEqual(a['moments'], b['moments'])
        self.assertAlmostEqual(a['moments']['childless_rate_40_44'], .4)
        self.assertAlmostEqual(a['moments']['exactly_one_among_mothers_40_44'], .5)

    def test_unit_first_birth_interpolation_has_hand_computed_stock(self):
        parity = np.tile([1., 0., 0., 0.], (7, 1))
        probabilities = np.zeros((7, 3)); probabilities[5:, 0] = 1.
        ev, P = fixture(masses=[1, 1, 1, 1, 1, 2, 4], parity=parity, probabilities=probabilities)
        a = self.observe(ev, P)
        b = self.observe(ev, P, "constant_post_cell")
        self.assertAlmostEqual(a['moments']['childless_rate_40_44'], .53125)
        self.assertAlmostEqual(a['moments']['exactly_one_among_mothers_40_44'], 1.)
        self.assertEqual(b['moments']['childless_rate_40_44'], 0.)
        self.assertEqual(a['accounting']['overlap_weights'][5:], [.5, .75])
        self.assertEqual(a['accounting']['post_parity_interpolation_share'][5:], [.75, .375])
        self.assertAlmostEqual(a['moments']['period_mean_age_first_birth'], 128./3.)

    def test_continuation_birth_reduces_exactly_one_without_creating_childlessness(self):
        parity = np.tile([0., 1., 0., 0.], (7, 1)); parity[0] = [1., 0., 0., 0.]
        probabilities = np.zeros((7, 3)); probabilities[0, 0] = .2; probabilities[5:, 1] = 1.
        ev, P = fixture(masses=[1, 1, 1, 1, 1, 2, 4], parity=parity, probabilities=probabilities)
        a = self.observe(ev, P)
        b = self.observe(ev, P, "constant_post_cell")
        self.assertEqual(a['moments']['childless_rate_40_44'], 0.)
        self.assertAlmostEqual(a['moments']['exactly_one_among_mothers_40_44'], .53125)
        self.assertEqual(b['moments']['exactly_one_among_mothers_40_44'], 0.)

    def test_unequal_model_age_masses_are_retained_before_forming_ratios(self):
        parity = np.tile([.4, .3, .2, .1], (7, 1))
        parity[5] = [1., 0., 0., 0.]; parity[6] = [0., 1., 0., 0.]
        probabilities = np.zeros((7, 3)); probabilities[0, 0] = .2
        ev, P = fixture(masses=[1, 1, 1, 1, 1, 2, 4], parity=parity, probabilities=probabilities)
        a = self.observe(ev, P)
        self.assertAlmostEqual(a['moments']['childless_rate_40_44'], .25)
        self.assertAlmostEqual(a['moments']['exactly_one_among_mothers_40_44'], 1.)
        self.assertEqual(a['accounting']['window_population_mass'], 4.)

    def test_midpoint_ages_and_age30_classification_use_period_flows(self):
        parity = np.tile([1., 0., 0., 0.], (7, 1))
        probabilities = np.zeros((7, 3)); probabilities[[0, 2, 3, 6], 0] = 1.
        ev, P = fixture(parity=parity, probabilities=probabilities)
        a = self.observe(ev, P)
        self.assertEqual(a['moments']['period_mean_age_first_birth'], 31.)
        self.assertEqual(a['moments']['period_share_first_births_age30plus'], .5)

    def test_zero_first_birth_denominator_fails_instead_of_zero_timing(self):
        ev, P = fixture(probabilities=np.zeros((7, 3)))
        with self.assertRaisesRegex(ValueError, "First-birth timing denominator"):
            self.observe(ev, P)

    def test_zero_projected_mother_or_population_mass_fails(self):
        probabilities = np.zeros((7, 3)); probabilities[0, 0] = .2
        parity = np.tile([1., 0., 0., 0.], (7, 1))
        for masses in ([1] * 7, [1, 1, 1, 1, 1, 0, 0]):
            with self.subTest(masses=masses):
                ev, P = fixture(masses=masses, parity=parity, probabilities=probabilities)
                with self.assertRaisesRegex(ValueError, "total and mother denominators"):
                    self.observe(ev, P)

    def test_geometry_and_model_architecture_fail_closed(self):
        for name, value in [('age_start', 20), ('da', 2.), ('period_years', 1.),
                ('J', 7.5), ('J', float('nan')), ('n_parity', 3), ('A_f_end', 6),
                ('fertility_units', 'parity2x'), ('sequential_births', False),
                ('joint_nested_choice', True), ('fertility_nest_choice', True),
                ('two_shock_choice', True), ('tfr_top_bin_weight', float('nan'))]:
            with self.subTest(name=name, value=value):
                ev, P = fixture(); setattr(P, name, value)
                with self.assertRaises(ValueError):
                    self.observe(ev, P)

    def test_nonfinite_negative_or_wrong_shaped_mass_fails(self):
        for value in [float('nan'), float('inf'), -.1]:
            with self.subTest(value=value):
                ev, P = fixture(); ev.g_pre.flat[0] = value
                with self.assertRaisesRegex(ValueError, 'finite and nonnegative'):
                    self.observe(ev, P)
        ev, P = fixture(); ev.g_post_fertility = ev.g_post_fertility[..., 0]
        with self.assertRaisesRegex(ValueError, 'seven-axis'):
            self.observe(ev, P)

    def test_pre_post_mass_or_parity_mismatch_fails(self):
        ev, P = fixture(); ev.g_post_fertility[0, 0, 0, 5, 0, 0, 0] += .001
        with self.assertRaisesRegex(RuntimeError, 'changes age mass'):
            self.observe(ev, P)
        ev, P = fixture()
        ev.g_post_fertility[0, 0, 0, 5, 0, 0, 0] -= .001
        ev.g_post_fertility[0, 0, 0, 5, 0, 1, 1] += .001
        with self.assertRaisesRegex(RuntimeError, 'stocks disagree'):
            self.observe(ev, P)

    def test_nonfinite_births_and_corrupt_helper_timing_fail(self):
        ev, P = fixture(); ev.births = float('nan')
        with self.assertRaises(RuntimeError):
            self.observe(ev, P)
        ev, P = fixture()
        actual = measurement.period_fertility_diagnostics(ev, P)
        actual['period_first_birth_mean_age'] += 2.
        with patch.object(measurement, 'period_fertility_diagnostics', return_value=actual):
            with self.assertRaisesRegex(RuntimeError, 'timing scalars disagree'):
                self.observe(ev, P)

    def test_invalid_first_birth_probabilities_fail(self):
        ev, P = fixture()
        ev.policy.fert_probs[0, 0, 0, 0, 0, 1] = 1.5
        with self.assertRaisesRegex(RuntimeError, 'hazards'):
            self.observe(ev, P)


if __name__ == '__main__':
    unittest.main()
