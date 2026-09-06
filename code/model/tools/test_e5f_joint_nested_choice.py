"""Independent mathematical/accounting checks for the experimental operator."""
import unittest
import numpy as np
from e5f_joint_nested_choice import choose, logsum_prob, plan_values, scatter_joint_block


class NestedChoiceTests(unittest.TestCase):
    def test_flat_logit_and_direct_gev_formula(self):
        q = np.array([[.2, -.8], [1.1, -.1]])
        for scale in [.005, .5, 2.5]:
            v, p = choose(q, scale, 1)
            flatv, flatp = logsum_prob(q.ravel(), scale)
            np.testing.assert_allclose(v, flatv, atol=1e-14)
            np.testing.assert_allclose(p.ravel(), flatp, atol=1e-14)
        k, lam = 1.2, .37
        expq = np.exp(q / (k * lam))
        sums = expq.sum(axis=1)
        direct = expq * sums[:, None] ** (lam - 1) / np.sum(sums ** lam)
        np.testing.assert_allclose(choose(q, k, lam)[1], direct, atol=1e-14)

    def test_derivative_of_value_is_joint_probability(self):
        q = np.array([[.4, -.3], [.1, .8]])
        for rule in ['joint', 'sequential']:
            _, p = choose(q, .7, .4, rule)
            for index in np.ndindex(q.shape):
                dq = np.zeros_like(q); dq[index] = 1e-6
                derivative = (choose(q + dq, .7, .4, rule)[0] -
                              choose(q - dq, .7, .4, rule)[0]) / 2e-6
                self.assertAlmostEqual(derivative, p[index], places=8)

    def test_translation_and_conditional_deterministic_limit(self):
        q = np.array([[.4, -.3], [.1, .8]])
        v, p = choose(q, .7, .4)
        vp, pp = choose(q + 127, .7, .4)
        self.assertAlmostEqual(vp - v, 127)
        np.testing.assert_allclose(pp, p, atol=1e-13)
        v0, p0 = choose(q, .7, 1e-6)
        exact, tenure = logsum_prob(q.max(axis=1), .7)
        np.testing.assert_allclose(v0, exact, atol=1e-13)
        np.testing.assert_allclose(p0, [[tenure[0], 0], [0, tenure[1]]], atol=1e-13)

    def test_infeasible_and_singletons(self):
        q = np.array([[.4, -np.inf], [-np.inf, -np.inf]])
        v, p = choose(q, .7, .4)
        self.assertEqual(v, .4)
        np.testing.assert_array_equal(p, [[1, 0], [0, 0]])
        v, p = choose(np.full((2, 2), -np.inf), .7, .4)
        self.assertEqual(v, -np.inf); self.assertEqual(p.sum(), 0)
        for lam in [0, -1, 1.01, np.nan]:
            with self.assertRaises(ValueError): choose(q, .7, lam)

    def test_conception_feasibility_and_first_birth_cost(self):
        q0 = np.array([1., -np.inf]); q1 = np.array([-np.inf, 3.])
        np.testing.assert_array_equal(plan_values(q0, q1, 0)[..., 1], q0)
        np.testing.assert_array_equal(plan_values(q0, q1, 1, .5)[..., 1], q1 - .5)
        self.assertTrue(np.isneginf(plan_values(q0, q1, .5)[..., 1]).all())
        q = plan_values([1, 2], [3, 6], .25, .4)
        np.testing.assert_allclose(q[:, 1], [.75 + .25 * 2.6, 1.5 + .25 * 5.6])
        self.assertTrue(np.isneginf(plan_values(q0, None, .5, available=False)[:, 1]).all())

    def test_joint_scatter_keeps_selection_and_no_double_birth(self):
        nb, nt, npar, ncs = 3, 2, 3, 3
        idx = np.zeros((nt, nt, npar, ncs, nb), dtype=int)
        wt = np.zeros_like(idx, dtype=float)
        idx[...] = [0, 1, 1]; wt[...] = [0, 0, 1]
        mass = np.zeros((nb, nt)); mass[1, 0] = 1
        # Perfect correlation: half wait/rent, half attempt/own.
        prob = np.zeros((nb, nt, 2, 2)); prob[..., 0, 0] = .5; prob[..., 1, 1] = .5
        choices = np.zeros((nb, nt, 2), dtype=int); choices[..., 1] = 1
        g, post, births = scatter_joint_block(mass, prob, choices, choices, .4,
                                            (1, 1), (0, 0), idx, wt)
        self.assertAlmostEqual(g.sum(), 1); self.assertAlmostEqual(post.sum(), 1)
        self.assertAlmostEqual(births, .2)
        self.assertAlmostEqual(g[:, 1, 1, 1].sum(), .2)
        self.assertAlmostEqual(g[:, 1, 0, 0].sum(), .3)
        self.assertAlmostEqual(g[:, 0, 1, 1].sum(), 0)
        self.assertEqual(g[:, :, 2, :].sum(), 0)

    def test_scatter_nonidentity_wealth_maps(self):
        nb, nt, npar, ncs = 3, 2, 2, 2
        idx = np.zeros((nt, nt, npar, ncs, nb), dtype=int)
        wt = np.full_like(idx, .25, dtype=float)
        idx[:, 1] = 1
        m = np.zeros((nb, nt)); m[0, 0] = 2
        p = np.zeros((nb, nt, 2, 2)); p[..., 1, 1] = 1
        c = np.zeros((nb, nt, 2), dtype=int); c[..., 1] = 1
        g, post, births = scatter_joint_block(m, p, c, c, 1, (1, 1), (0, 0), idx, wt)
        np.testing.assert_allclose(g[:, 1, 1, 1], [0, 1.5, .5])
        self.assertAlmostEqual(births, 2); self.assertAlmostEqual(post.sum(), 2)


if __name__ == '__main__':
    unittest.main()
