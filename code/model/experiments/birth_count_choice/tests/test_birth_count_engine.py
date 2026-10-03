"""No-solve unit checks for the separate direct-count experiment."""
from pathlib import Path
import importlib
import sys
import unittest
from types import SimpleNamespace
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
bc = importlib.import_module('birth_count_choice.model.engine.birth_count')
utils = importlib.import_module('birth_count_choice.model.engine.utils')


class BirthCountKernelTests(unittest.TestCase):
    def test_binomial_normalization_mean_and_endpoints(self):
        for pi in (0., .17, .63, 1.):
            for k in range(4):
                p = bc.binomial_probabilities(k, pi)
                self.assertAlmostEqual(float(p.sum()), 1.)
                self.assertAlmostEqual(float(p @ np.arange(4)), k * pi)
                self.assertTrue(np.all(p[k + 1:] == 0))

    def test_cap1_nests_existing_binary_bellman_at_every_family_state(self):
        rng = np.random.default_rng(771)
        VI = rng.normal(size=(5, 2, 4, 4))
        pi = .61
        for n in range(3):
            for m in range(n + 1):
                kappa = .12 if n == 0 else .36
                cost = .4 if n == 0 else 0.
                old = np.stack((VI[..., n, m], pi * (VI[..., n + 1, m + 1] - cost)
                    + (1 - pi) * VI[..., n, m]), axis=-1)
                old_ls, old_p = utils.logsumexp(old / kappa, axis=-1)
                v, p, r, u = bc.birth_count_menu(VI, n, m, pi, kappa, .4, cap=1)
                np.testing.assert_allclose(v, kappa * old_ls, rtol=1e-14, atol=1e-14)
                np.testing.assert_allclose(p[..., :2], old_p, rtol=1e-14, atol=1e-14)
                np.testing.assert_allclose(r[..., 1], pi * old_p[..., 1])
                np.testing.assert_allclose(r[..., 0], 1 - pi * old_p[..., 1])
                np.testing.assert_array_equal(p[..., 2:], 0.)

    def test_large_negative_values_preserve_strict_probability_mass(self):
        VI = np.full((7, 4, 4), -4e8)
        VI[..., 1, 1] += .3
        VI[..., 2, 2] += .9
        for n in range(3):
            _, action, realized, _ = bc.birth_count_menu(VI, n, n, .37, .021, .4)
            np.testing.assert_allclose(action.sum(axis=-1), 1., rtol=0, atol=1e-15)
            np.testing.assert_allclose(realized.sum(axis=-1), 1., rtol=0, atol=1e-15)

    def test_fixed_cost_once_conditional_on_any_success(self):
        VI = np.zeros((2, 4, 4))
        _, _, _, u = bc.birth_count_menu(VI, 0, 0, .4, .2, 7.)
        for k in range(4):
            np.testing.assert_allclose(u[..., k], -7 * (1 - .6**k))
        _, _, _, parent_u = bc.birth_count_menu(VI, 1, 0, .4, .2, 7.)
        np.testing.assert_array_equal(parent_u[..., :3], 0.)

    def test_menu_remaining_cap_and_noise_scale(self):
        VI = np.zeros((3, 4, 4))
        for n in range(4):
            v, action, realized, _ = bc.birth_count_menu(VI, n, 0, .25, .7)
            np.testing.assert_allclose(v, .7 * np.log(4 - n))
            np.testing.assert_allclose(action[..., :4 - n], 1 / (4 - n))
            np.testing.assert_array_equal(action[..., 4 - n:], 0.)
            np.testing.assert_array_equal(realized[..., 4 - n:], 0.)
            np.testing.assert_allclose(realized.sum(axis=-1), 1.)

    def test_dead_menu_has_no_probability(self):
        VI = np.full((2, 4, 4), -1e10)
        _, action, realized, _ = bc.birth_count_menu(VI, 0, 0, .5, .2)
        np.testing.assert_array_equal(action, 0.)
        np.testing.assert_array_equal(realized, 0.)

    def test_multiple_birth_mixture_mass_orders_and_tagged_destinations(self):
        pre = np.zeros((2, 4, 4))
        pre[0, 0, 0] = 10.
        pre[1, 1, 0] = 6.
        rp = bc.identity_realized_probabilities(pre.shape)
        rp[0, 0, 0] = bc.binomial_probabilities(3, .5)
        rp[1, 1, 0] = bc.binomial_probabilities(2, .5)
        pre_before = pre.copy(); rp_before = rp.copy()
        flow = bc.birth_count_transition(pre, rp)
        self.assertAlmostEqual(flow['post'].sum(), pre.sum())
        np.testing.assert_allclose(flow['births_by_order'], (8.75, 9.5, 2.75))
        self.assertAlmostEqual(flow['expected_births'], 10 * 3 * .5 + 6 * 2 * .5)
        self.assertAlmostEqual(flow['any_birth_mass'], 10 * .875 + 6 * .75)
        np.testing.assert_allclose(flow['first_birth_tagged_post'][0].diagonal(), (0., 3.75, 3.75, 1.25))
        for q in range(3):
            self.assertAlmostEqual(flow['order_tagged_post'][q].sum(), flow['births_by_order'][q])
        np.testing.assert_array_equal(pre, pre_before); np.testing.assert_array_equal(rp, rp_before)
        self.assertFalse(np.shares_memory(flow['post'], pre))
        self.assertFalse(np.shares_memory(flow['first_birth_tagged_post'], flow['post']))

    def test_no_cascade_cap1_nests_binary_forward_from_snapshot(self):
        pre = np.zeros((4, 4))
        pre[0, 0], pre[1, 1], pre[2, 0], pre[3, 0] = 10., 5., 3., 2.
        action = bc.identity_realized_probabilities(pre.shape)
        rp = bc.identity_realized_probabilities(pre.shape)
        pi = .7
        for n in range(3):
            for m in range(n + 1):
                action[n, m, :2] = (.6, .4)
                rp[n, m, :2] = (1 - .4 * pi, .4 * pi)
        old = pre.copy(); born = np.zeros(3)
        for n in range(3):
            for m in range(n + 1):
                flow = pre[n, m] * .4 * pi
                old[n, m] -= flow; old[n + 1, m + 1] += flow; born[n] += flow
        flow = bc.birth_count_transition(pre, rp, action, cap=1)
        np.testing.assert_allclose(flow['post'], old)
        np.testing.assert_allclose(flow['births_by_order'], born)
        np.testing.assert_allclose(flow['at_risk_by_order'], (10., 5., 3.))
        np.testing.assert_allclose(flow['attempts_by_order'], (4., 2., 1.2))

    def test_invalid_mass_and_increment_rejected(self):
        pre = np.zeros((4, 4)); pre[0, 1] = 1.
        rp = bc.identity_realized_probabilities(pre.shape)
        with self.assertRaises(ValueError): bc.birth_count_transition(pre, rp)
        pre[0, 1] = 0.; pre[2, 0] = 1.; rp[2, 0, 2] = 1.
        with self.assertRaises(ValueError): bc.birth_count_transition(pre, rp)
        rp[2, 0] = 0.
        with self.assertRaises(ValueError): bc.birth_count_transition(pre, rp)

    def test_probability_validation_strict_gate_and_action_menu(self):
        pre = np.zeros((4, 4)); pre[0, 0] = 1.
        rp = bc.identity_realized_probabilities(pre.shape)
        action = rp.copy()
        rp[0, 0, 0] += 2e-12
        with self.assertRaises(ValueError): bc.birth_count_transition(pre, rp)
        rp[0, 0, 0] = 1.
        rp[0, 0, :3] = (.5, 0., .5)
        with self.assertRaises(ValueError): bc.birth_count_transition(pre, rp, cap=1)
        rp[0, 0] = (1., 0., 0., 0.)
        for bad in ((-.1, 1.1, 0., 0.), (np.nan, 0., 0., 0.), (.9, 0., 0., 0.)):
            action[0, 0] = bad
            with self.assertRaises(ValueError): bc.birth_count_transition(pre, rp, action)

    def test_fast_solution_probability_snapshot_does_not_alias_parameters(self):
        distribution = importlib.import_module('birth_count_choice.model.engine.distribution')
        action = bc.identity_realized_probabilities((1, 1, 1, 2, 1, 4, 4))
        P = SimpleNamespace(birth_count_choice_enabled=True,
            birth_count_action_probs=action, birth_count_realized_probs=action.copy(),
            birth_count_policy_axes=('test',), user_cost_rate=1., H0=1., r_bar=1., xi_supply=1.)
        st = SimpleNamespace(housing_demand=np.ones(1))
        sol = distribution.pack_fast_solution_markov_income(st, np.ones(1), P)
        before = sol.birth_count_action_probs.copy()
        P.birth_count_action_probs[...] = 0.
        np.testing.assert_array_equal(sol.birth_count_action_probs, before)
        self.assertFalse(np.shares_memory(sol.birth_count_realized_probs, P.birth_count_realized_probs))

    def test_existing_binary_bellman_body_is_ast_identical(self):
        import ast
        original = ast.parse((ROOT.parents[1] / 'production/engine/household.py').read_text())
        experiment = ast.parse((ROOT / 'model/engine/household.py').read_text())
        def binary_body(tree):
            for node in ast.walk(tree):
                if (isinstance(node, ast.If) and node.body
                        and isinstance(node.body[0], ast.Assign)
                        and isinstance(node.body[0].targets[0], ast.Name)
                        and node.body[0].targets[0].id == 'Vfa'
                        and ast.unparse(node.test) == "bool(getattr(P, 'sequential_births', False))"):
                    return ast.dump(ast.Module(body=node.body, type_ignores=[]), include_attributes=False)
            raise AssertionError('No existing binary Bellman body found')
        self.assertEqual(binary_body(original), binary_body(experiment))

    def test_negligible_dead_zero_menu_keeps_origin_mass_without_births(self):
        pre = np.zeros((2, 4, 4)); pre[0, 0, 0] = 1.; pre[1, 1, 0] = 3e-38
        rp = bc.identity_realized_probabilities(pre.shape); rp[1, 1, 0] = 0.
        action = rp.copy()
        values = np.zeros_like(pre); values[1, 1, 0] = -1e10
        result = bc.birth_count_transition(pre, rp, action, state_values=values)
        np.testing.assert_array_equal(result['post'], pre)
        np.testing.assert_array_equal(result['births_by_order'], 0.)
        self.assertEqual(result['zero_menu_identity_mass'], 3e-38)
        self.assertEqual(result['zero_menu_identity_cell_count'], 1)
        # Reporting callbacks lack V but retain the same global mass guard.
        np.testing.assert_array_equal(bc.birth_count_transition(pre, rp)['post'], pre)

    def test_zero_menu_global_mass_guard_and_living_source_rejection(self):
        pre = np.zeros((2, 4, 4)); pre[:, 0, 0] = 6e-13
        rp = bc.identity_realized_probabilities(pre.shape); rp[:, 0, 0] = 0.
        with self.assertRaises(bc.BirthCountProbabilityError) as error:
            bc.birth_count_transition(pre, rp)
        self.assertEqual(error.exception.evidence['kind'], 'zero_menu_mass_exceeds_native_dead_tolerance')
        pre[:, 0, 0] = 1e-30
        with self.assertRaises(bc.BirthCountProbabilityError) as error:
            bc.birth_count_transition(pre, rp, state_values=np.zeros_like(pre))
        self.assertEqual(error.exception.evidence['kind'], 'zero_menu_living_source')
        action = bc.identity_realized_probabilities(pre.shape)
        with self.assertRaises(bc.BirthCountProbabilityError):
            bc.birth_count_transition(pre, rp, action)

    def test_contract_default_off_and_unsupported_variants(self):
        bc.validate_contract(SimpleNamespace())
        P = SimpleNamespace(birth_count_choice_enabled=True, n_parity=4, n_child_states=4,
            child_state_mode='independent_count', sequential_births=True)
        bc.validate_contract(P)
        for name, value in [('joint_nested_choice', True), ('readiness_gate_enabled', True),
                ('child_maturation_mode', 'parent_age'), ('child_state_mode', 'shared_clock')]:
            Q = SimpleNamespace(**vars(P)); setattr(Q, name, value)
            with self.assertRaises(ValueError): bc.validate_contract(Q)


if __name__ == '__main__':
    unittest.main()
