"""Pure guards for the matched PF primitive smoke; no model solve."""
import unittest
from types import SimpleNamespace
import numpy as np
import run_e5f_matched_pf_smoke as smoke


class PrimitiveGuards(unittest.TestCase):
    def test_nested_reference_uses_same_pool_and_still_rejects_changed_choices(self):
        from test_e5f_pf_choice_ownership import parameters, choice_object, SHAPE
        P = parameters()
        reference = SimpleNamespace(joint_choice=choice_object(P, .6),
                                    tenure_probs=np.full(SHAPE + (2,), .9))
        pre = np.zeros(SHAPE)
        pre[0, 0, 0, 0, 0, 0, 0] = 1.
        post, effective, births, _, _ = smoke.calendar.factor_joint_distribution(pre, reference, P)
        ev = SimpleNamespace(g_pre=pre, g_post_fertility=post, births=float(births.sum()))
        arrays = smoke.same_population_reference(reference, ev, P)
        np.testing.assert_array_equal(arrays['tenure_probs'], effective)
        self.assertTrue(np.all(reference.tenure_probs == .9))
        reference.joint_choice = choice_object(P, .2)
        with self.assertRaisesRegex(RuntimeError, 'mass/birth'):
            smoke.same_population_reference(reference, ev, P)

    def test_supply_pin_uses_saved_rule_not_legacy_parameter_fields(self):
        P = SimpleNamespace(tenure_choice_kappa=.005, psi_child=-.03,
                            eta_supply=np.array([1.75]), xi_supply=np.array([1.75]))
        selected = {'best_candidate': {'theta': {'tenure_choice_kappa': .005, 'psi_child': 0.},
                                      'new_psi_child': -.03}}
        smoke.verify_selected_parameters(P, selected, SimpleNamespace(elasticity=.63))
        with self.assertRaisesRegex(RuntimeError, 'saved housing-supply'):
            smoke.verify_selected_parameters(P, selected, SimpleNamespace(elasticity=1.75))

    def test_probability_change_is_not_hidden_by_other_arrays(self):
        old = {'joint_probabilities': np.array([.2, .8]), 'V': np.array([1.])}
        new = {'joint_probabilities': np.array([.3, .7]), 'V': np.array([1.])}
        with self.assertRaisesRegex(RuntimeError, 'joint_probabilities'):
            smoke.compare_arrays(old, new)

    def test_missing_joint_arrays_fail(self):
        with self.assertRaisesRegex(RuntimeError, 'field sets'):
            smoke.compare_arrays({'joint_probabilities': np.array([1.])}, {})

    def test_nan_cannot_pass_reproduction(self):
        with self.assertRaises(RuntimeError):
            smoke.compare_arrays({'V': np.array([np.nan])}, {'V': np.array([np.nan])})

    def test_budget_uses_dated_rent(self):
        shape = (1, 1, 1, 1, 1, 1, 1)
        P = SimpleNamespace(J=1, z_grid=[1.], n_parity=1, n_child_states=1,
            R_gross=1., user_cost_rate=.01)
        pol = SimpleNamespace(price=np.array([1.]), hR_pol=np.ones(shape),
            c_pol=np.full(shape, .8), bp_pol=np.zeros(shape))
        ev = SimpleNamespace(policy=pol, g_current=np.ones(shape))
        shared = SimpleNamespace(gb_flat=np.array([0.]))
        from unittest.mock import patch
        with patch.object(smoke.model, 'income_at_state', return_value=1.):
            self.assertEqual(smoke.dated_budget(ev, P, shared, np.array([0.]), .2)['budget_excess_mass'], 0.)
            with self.assertRaisesRegex(RuntimeError, 'Dated budget'):
                smoke.dated_budget(ev, P, shared, np.array([0.]), .3)


if __name__ == '__main__':
    unittest.main()
