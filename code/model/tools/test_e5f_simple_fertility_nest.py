"""Independent complete-menu checks of the simple fertility-group GEV.

Run with the isolated project Python; no model solve or calibration is performed.
The reference explicitly enumerates K wait plans and K**2 attempted-birth plans.
It does not use the production operator's Cartesian-product factorization.
"""
import itertools
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp, softmax

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code/model"))
from intergen_eqscale_seq_optimized.fertility_nested import choose, bellman_block, factor_age


def brute_force(q0, q1, pi, kappa, sigma, cost=0., available=True):
    """One-state reference, retaining distinct wait/attempt labels at endpoints."""
    q0, q1 = np.asarray(q0), np.asarray(q1)
    k = len(q0)
    groups = [[], []]
    labels = [[], []]
    for h in range(k):
        if np.isfinite(q0[h]):
            groups[0].append(q0[h])
            labels[0].append((h, None))
    if available:
        failed = range(k) if pi < 1 else [None]
        successful = range(k) if pi > 0 else [None]
        for hf, hs in itertools.product(failed, successful):
            if hf is not None and not np.isfinite(q0[hf]):
                continue
            if hs is not None and not np.isfinite(q1[hs]):
                continue
            value = -pi * cost
            if hf is not None:
                value += (1 - pi) * q0[hf]
            if hs is not None:
                value += pi * q1[hs]
            groups[1].append(value)
            labels[1].append((hf, hs))
    inclusive = np.full(2, -np.inf)
    conditional = [np.array([]), np.array([])]
    for a in range(2):
        if groups[a]:
            inclusive[a] = kappa * logsumexp(np.array(groups[a]) / kappa)
            conditional[a] = softmax(np.array(groups[a]) / kappa)
    action = np.zeros(2)
    valid = np.isfinite(inclusive)
    value = -np.inf
    if valid.any():
        # With only wait available, the outer logsum adds no extra entropy.
        value = sigma * logsumexp(inclusive[valid] / sigma)
        action[valid] = softmax(inclusive[valid] / sigma)
    housing = np.zeros((3, k))
    for a in range(2):
        for p, (hf, hs) in zip(conditional[a], labels[a]):
            if a == 0:
                housing[0, hf] += p
            else:
                if hf is not None:
                    housing[1, hf] += p
                if hs is not None:
                    housing[2, hs] += p
    return dict(value=value, action_probability=action,
                wait_housing=housing[0], failure_housing=housing[1],
                success_housing=housing[2]), groups, conditional


class SimpleFertilityNestTests(unittest.TestCase):
    def assert_result(self, actual, expected, atol=4e-12):
        for key in expected:
            np.testing.assert_allclose(actual[key], expected[key], atol=atol,
                                       rtol=2e-13, err_msg=key)

    def test_complete_42_plan_enumeration(self):
        rng = np.random.default_rng(784112)
        for k in (2, 6):
            for kappa, sigma in ((.005, 2.1681), (.13, .8), (.4, .4)):
                for _ in range(30):
                    q0, q1 = rng.normal(0, .1, (2, k))
                    pi = rng.uniform(.00001, .99999)
                    cost = rng.uniform(0, .4)
                    expected, groups, _ = brute_force(q0, q1, pi, kappa, sigma, cost)
                    self.assertEqual(len(groups[0]) + len(groups[1]), k + k*k)
                    self.assert_result(choose(q0, q1, pi, kappa, sigma, cost), expected)

    def test_all_two_product_feasibility_masks(self):
        q = np.array([[.1, .4], [-.2, .3]])
        for mask, pi, available in itertools.product(range(16), (0., .23, 1.), (False, True)):
            qm = q.copy()
            qm.ravel()[[not(mask & (1 << i)) for i in range(4)]] = -np.inf
            expected, _, _ = brute_force(*qm, pi, .13, .8, .2, available)
            self.assert_result(choose(*qm, pi, .13, .8, .2, available), expected)

    def test_unavailable_attempt_has_only_wait_value(self):
        q0 = np.array([.1, -.2, -np.inf, .3, .2, -.3])
        q1 = np.full(6, -np.inf)
        expected, _, _ = brute_force(q0, q1, .7, .005, 2.1681, available=False)
        self.assert_result(choose(q0, q1, .7, .005, 2.1681, available=False), expected)

    def test_endpoint_collapses_only_unused_coordinate(self):
        q0 = np.array([.2, .3, .1, -.2, -np.inf, 0.])
        q1 = np.array([-.1, .4, .2, .3, .5, -np.inf])
        for pi in (0., 1.):
            expected, groups, _ = brute_force(q0, q1, pi, .13, .8, .2)
            self.assertEqual(len(groups[1]), 5)
            self.assert_result(choose(q0, q1, pi, .13, .8, .2), expected)
        # Attempt and wait still have separate labels even with pi=0.
        result = choose(q0, np.full(6, -np.inf), 0., .13, .8)
        np.testing.assert_allclose(result['action_probability'], [.5, .5])
        self.assertAlmostEqual(float(result['value']),
                               .13*logsumexp(q0/.13) + .8*np.log(2))

    def test_flat_scale_matches_one_softmax_over_all_plans(self):
        q0 = np.array([.1, -.2, .4, -np.inf, .2, .5])
        q1 = np.array([.2, .3, -.1, .1, .7, -np.inf])
        for pi in (0., .31, 1.):
            actual = choose(q0, q1, pi, .3, .3, .14)
            _, groups, conditional = brute_force(q0, q1, pi, .3, .3, .14)
            all_values = np.concatenate(groups)
            flat = softmax(all_values/.3)
            self.assertAlmostEqual(float(actual['value']), .3*logsumexp(all_values/.3))
            nwait = len(groups[0])
            np.testing.assert_allclose(actual['action_probability'],
                                       [flat[:nwait].sum(), flat[nwait:].sum()])
            np.testing.assert_allclose(flat[nwait:],
                actual['action_probability'][1]*conditional[1])

    def test_value_envelope_derivatives_include_outcome_probabilities(self):
        rng = np.random.default_rng(66890)
        for kappa, sigma in ((.005, 2.1681), (.13, .8)):
            q0, q1 = rng.normal(0, .05, (2, 6))
            pi, cost, step = .38, .17, 1e-6
            result = choose(q0, q1, pi, kappa, sigma, cost)
            a0, a1 = result['action_probability']
            gradients = (a0*result['wait_housing'] + a1*(1-pi)*result['failure_housing'],
                         a1*pi*result['success_housing'])
            for outcome, h in itertools.product(range(2), range(6)):
                up, down = np.array([q0, q1]), np.array([q0, q1])
                up[outcome, h] += step
                down[outcome, h] -= step
                derivative = (choose(*up, pi, kappa, sigma, cost)['value'] -
                              choose(*down, pi, kappa, sigma, cost)['value'])/(2*step)
                self.assertAlmostEqual(float(derivative), gradients[outcome][h], delta=2e-8)
            dc = (choose(q0, q1, pi, kappa, sigma, cost+step)['value'] -
                  choose(q0, q1, pi, kappa, sigma, cost-step)['value'])/(2*step)
            self.assertAlmostEqual(float(dc), -a1*pi, delta=2e-8)
            dp = (choose(q0, q1, pi+step, kappa, sigma, cost)['value'] -
                  choose(q0, q1, pi-step, kappa, sigma, cost)['value'])/(2*step)
            expected_dp = a1*(np.dot(result['success_housing'], q1) - cost -
                              np.dot(result['failure_housing'], q0))
            self.assertAlmostEqual(float(dp), expected_dp, delta=2e-8)

    def test_realized_joint_mass_and_contingent_conditionals(self):
        q0 = np.array([.01, .015, .03, .01, .02, .025])
        q1 = np.array([.023, .013, .04, .02, .03, .005])
        pi = .24
        r = choose(q0, q1, pi, .005, 2.1681)
        _, _, conditional = brute_force(q0, q1, pi, .005, 2.1681)
        joint_attempt = conditional[1].reshape(6, 6)
        np.testing.assert_allclose(joint_attempt,
            np.outer(r['failure_housing'], r['success_housing']), atol=2e-14)
        failure_mass = r['action_probability'][1]*(1-pi)*r['failure_housing']
        success_mass = r['action_probability'][1]*pi*r['success_housing']
        wait_mass = r['action_probability'][0]*r['wait_housing']
        self.assertAlmostEqual(float((failure_mass + success_mass + wait_mass).sum()), 1.)
        self.assertGreater(np.max(abs(r['failure_housing']-r['wait_housing'])), .01)

    def test_batch_shape_and_state_specific_scales(self):
        rng = np.random.default_rng(7228)
        q0, q1 = rng.normal(0, .03, (2, 2, 3, 6))
        sigma = np.array([[2.1681, .01, .2], [.01, .005, .8]])
        cost = .17
        for pi, available in itertools.product((0., .2, 1.), (False, True)):
            result = choose(q0, q1, pi, .005, sigma, cost, available)
            self.assertEqual(np.shape(result['value']), (2, 3))
            self.assertEqual(np.shape(result['action_probability']), (2, 3, 2))
            for idx in np.ndindex(2, 3):
                expected, _, _ = brute_force(q0[idx], q1[idx], pi, .005,
                                            sigma[idx], cost, available)
                self.assert_result({key: value[idx] for key, value in result.items()
                                    if key in expected}, expected)

    def test_translation_permutation_and_large_value_stability(self):
        q0 = np.array([.1, .11, -np.inf, -.01, -.02, .09])
        q1 = np.array([.12, -.12, .2, -.01, .09, -np.inf])
        r = choose(q0, q1, .34, .005, 2.1681, .2)
        for shift in (-1e6, 1e6):
            shifted = choose(q0+shift, q1+shift, .34, .005, 2.1681, .2)
            shifted['value'] = shifted['value']-shift
            self.assert_result(shifted, r, atol=2e-8)
        perm = np.array([5, 2, 3, 0, 4, 1])
        shuffled = choose(q0[perm], q1[perm], .34, .005, 2.1681, .2)
        for key in ('wait_housing', 'failure_housing', 'success_housing'):
            shuffled[key] = shuffled[key][np.argsort(perm)]
        self.assert_result(shuffled, r)

    def test_near_zero_pi_retains_full_plan_entropy(self):
        q = np.zeros(6)
        positive = choose(q, q, 1e-12, .005, 2.1681)
        endpoint = choose(q, q, 0., .005, 2.1681)
        inclusive_wait = .005*np.log(6)
        inclusive_attempt = .005*np.log(36)
        expected = 2.1681*logsumexp(np.array([inclusive_wait, inclusive_attempt])/2.1681)
        self.assertAlmostEqual(float(positive['value']), expected)
        self.assertGreater(float(positive['value']-endpoint['value']), .004)

    def test_simple_nest_is_not_sequential_recursion(self):
        q0 = np.array([.00, .004])
        q1 = np.array([.001, .012])
        pi, kappa, sigma = .31, .005, 2.1681
        iw = kappa*logsumexp(q0/kappa)
        sequential_attempt = (1-pi)*iw + pi*kappa*logsumexp(q1/kappa)
        sequential_value = sigma*logsumexp(np.array([iw, sequential_attempt])/sigma)
        actual = choose(q0, q1, pi, kappa, sigma)['value']
        self.assertGreater(abs(float(actual)-sequential_value), 1e-4)

    def test_probability_guards_at_extreme_odds(self):
        for offset, pi in itertools.product((-1e5, 0., 1e5), (1e-15, .5, 1-1e-15)):
            q0 = np.array([-20., 0., 20., -np.inf, 19., -19.]) + offset
            q1 = np.array([20., 0., -20., -19., -np.inf, 19.]) + offset
            r = choose(q0, q1, pi, .005, 2.1681, .3)
            self.assertTrue(np.isfinite(r['value']))
            for key in ('action_probability', 'wait_housing', 'failure_housing', 'success_housing'):
                p = np.asarray(r[key])
                self.assertTrue(np.all(np.isfinite(p)))
                self.assertTrue(np.all((p >= 0) & (p <= 1)))
                self.assertAlmostEqual(float(p.sum()), 1., delta=4e-14)

    def test_rejects_invalid_scales_and_probabilities(self):
        q = np.zeros(2)
        for kappa, sigma in ((0., 1.), (-.1, 1.), (.2, .1), (.2, np.nan), (.2, np.inf)):
            with self.assertRaises(ValueError):
                choose(q, q, .4, kappa, sigma)
        for pi in (-.1, 1.1, np.nan):
            with self.assertRaises(ValueError):
                choose(q, q, pi, .1, .2)

    def test_empty_menu_cannot_invent_probability_mass(self):
        q = np.full(6, -np.inf)
        result = choose(q, q, .4, .005, 2.1681)
        self.assertTrue(np.isneginf(result['value']))
        for key in ('action_probability', 'wait_housing', 'failure_housing', 'success_housing'):
            np.testing.assert_array_equal(result[key], 0.)

    def test_bellman_population_and_selected_first_birth_cohorts(self):
        """Six-product Bellman output through every population measurement mode."""
        q = np.zeros((2, 6, 2, 2, 6))
        for wealth, tenure, parity, children in np.ndindex(q.shape[:-1]):
            q[wealth, tenure, parity, children] = (
                np.array([.01, .03, .025, .015, .04, .02]) +
                parity*np.array([.02, -.01, .01, -.02, -.025, .015]) +
                wealth*np.array([.005, -.003, .002, .001, -.004, .003]))
        # Product infeasibility differs by realized family state.
        q[:, :, 0, 0, 2] = -np.inf
        q[:, :, 1, 1, 4] = -np.inf
        q[:, :, 0, 1, :] = -np.inf
        g = np.zeros((2, 6, 1, 1, 2, 2))
        g[0, 0, 0, 0, 0, 0] = .4
        g[1, 2, 0, 0, 0, 0] = .3
        g[1, 5, 0, 0, 1, 1] = .3
        for pi, fertile in itertools.product((0., .24, 1.), (False, True)):
            p = SimpleNamespace(E_loc=np.array([0.]), mu_stay=0.,
                tenure_choice_kappa=.005, kappa_fert=2.1681,
                kappa_fert_continuation=.8, A_f_start=1 if fertile else 2,
                A_f_end=1 if fertile else 2, n_parity=2, n_child_states=2,
                first_birth_fixed_cost=.01, birth_entry_grant=False)
            calls = []
            def kernel(restricted, grants, rental_flag):
                products = np.flatnonzero(np.any(restricted > -1e9, axis=(0, 2, 3)))
                self.assertEqual(len(products), 1)
                self.assertTrue(rental_flag)
                product = int(products[0])
                calls.append(product)
                return np.where(np.isfinite(q[..., product]), q[..., product], -1e10), None
            value, prob, products, wait, failure = bellman_block(
                np.ones((2, 6, 2, 2)), (np.array([False]),), p, 0, np.array([pi]), kernel)
            self.assertEqual(calls, list(range(6)))
            np.testing.assert_array_equal(products, np.broadcast_to(np.arange(6), q.shape))
            joint = SimpleNamespace(
                probabilities=prob[:, :, None, None, None, ...],
                wait_probabilities=wait[:, :, None, None, None, ...],
                failure_probabilities=failure[:, :, None, None, None, ...])
            results = {}
            for mode in ('natural', 'wait', 'first_birth_treated', 'first_birth_control'):
                expected_weighted = np.zeros(g.shape + (6,))
                expected_births, expected_attempts, expected_risk = (np.zeros(2) for _ in range(3))
                for idx in np.argwhere(g > 0):
                    idx = tuple(idx)
                    wealth, tenure, _, _, parity, children = idx
                    mass = g[idx]
                    eligible = fertile and parity == 0
                    r, _, _ = brute_force(q[wealth, tenure, parity, children],
                        q[wealth, tenure, 1, 1], pi, .005, 2.1681,
                        .01, eligible)
                    self.assertAlmostEqual(float(value[wealth, tenure, parity, children]), r['value'])
                    if mode.startswith('first_birth_') and not eligible:
                        continue
                    a0, a1 = r['action_probability']
                    if mode == 'natural':
                        outcomes = [(parity, children, a0*r['wait_housing'])]
                        if eligible:
                            outcomes += [(0, 0, a1*(1-pi)*r['failure_housing']),
                                         (1, 1, a1*pi*r['success_housing'])]
                            expected_births[0] += mass*a1*pi
                            expected_attempts[0] += mass*a1
                            expected_risk[0] += mass
                    elif mode == 'wait':
                        outcomes = [(parity, children, r['wait_housing'])]
                    else:
                        expected_births[0] += mass*a1*pi
                        outcomes = ([(1, 1, a1*pi*r['success_housing'])]
                            if mode == 'first_birth_treated' else
                            [(0, 0, a1*pi*r['wait_housing'])])
                    for dn, dc, weights in outcomes:
                        destination = (wealth, tenure, 0, 0, dn, dc)
                        expected_weighted[destination] += mass*weights
                expected_post = expected_weighted.sum(axis=-1)
                expected_effective = np.divide(expected_weighted, expected_post[..., None],
                    out=np.zeros_like(expected_weighted), where=expected_post[..., None] > 0)
                with patch('intergen_eqscale_seq_optimized.parameters.get_fecundity_by_age',
                           return_value=np.array([pi])):
                    actual = factor_age(g, joint, p, 0, mode)
                for got, expected in zip(actual, (expected_post, expected_effective,
                        expected_births, expected_attempts, expected_risk)):
                    np.testing.assert_allclose(got, expected, atol=3e-14, rtol=2e-13)
                results[mode] = actual
            treated, control = results['first_birth_treated'], results['first_birth_control']
            self.assertAlmostEqual(float(treated[0].sum()), float(control[0].sum()))
            if fertile and 0 < pi < 1:
                # The matched childless control must not use the conditional
                # housing law following a failed attempted birth.
                idx = (0, 0, 0, 0, 0, 0)
                r, _, _ = brute_force(q[0, 0, 0, 0], q[0, 0, 1, 1], pi, .005, 2.1681, .01)
                np.testing.assert_allclose(control[1][idx], r['wait_housing'])
                self.assertGreater(np.max(abs(control[1][idx]-r['failure_housing'])), .01)


if __name__ == '__main__':
    unittest.main(verbosity=2)
