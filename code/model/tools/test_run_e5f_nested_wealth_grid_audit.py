"""Check that a nested grid preserves the inherited joint probability measure."""
import unittest
import numpy as np

from run_e5f_nested_wealth_grid_audit import nested_grid, embed_measure


class NestedMeasureTests(unittest.TestCase):
    def test_joint_measure_preserves_nonlinear_state_payoffs(self):
        # Unequal interval lengths and mass at both endpoints stress the map.
        bg = np.array([-12., -5., -.17, 0., .09, 1.3, 7., 30.])
        g = np.random.default_rng(20260905).uniform(size=(len(bg), 3, 2, 4))
        g /= g.sum()
        weights = np.arange(24).reshape(3, 2, 4) + 1
        for subdivisions in (1, 2, 4):
            fine, old = nested_grid(bg, subdivisions)
            mapped, receipt = embed_measure(g, bg, fine, old)
            for payoff in (lambda b:b, lambda b:b*b, lambda b:np.exp(.02*b), lambda b:(b >= 0).astype(float)):
                original = np.sum(g * payoff(bg)[:, None, None, None] * weights)
                refined = np.sum(mapped * payoff(fine)[:, None, None, None] * weights)
                self.assertAlmostEqual(original, refined, places=12)
            self.assertEqual(receipt["mass_on_added_nodes"], 0.)
            np.testing.assert_allclose(mapped.sum(axis=0), g.sum(axis=0), atol=1e-15, rtol=0)

    def test_rejects_non_nested_or_unsorted_inputs(self):
        for bg, subdivisions in ((np.array([-1., 0., 1.]), 3),
                                  (np.array([-1., 0., 0., 1.]), 2),
                                  (np.array([0., -1., 1.]), 2)):
            with self.assertRaises(ValueError):
                nested_grid(bg, subdivisions)


if __name__ == "__main__":
    unittest.main()
