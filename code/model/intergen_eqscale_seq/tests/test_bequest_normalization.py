from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from intergen_eqscale_seq.solver import bequest_utility_vec


def _params(*, sigma: float = 2.0, normalized: bool = True) -> SimpleNamespace:
    return SimpleNamespace(
        theta0=0.8,
        theta_n=0.4,
        theta1=0.2,
        bequest_spec="linear_child_scale",
        sigma=float(sigma),
        normalize_bequest_utility=bool(normalized),
        estate_tax_rate=0.0,
        estate_tax_exemption=0.0,
    )


class BequestNormalizationTests(unittest.TestCase):
    def test_zero_estate_is_zero_for_every_family_size(self) -> None:
        p = _params()
        for n_children in (0, 1, 2, 4):
            result = bequest_utility_vec(np.array([-1.0, 0.0]), n_children, p)
            np.testing.assert_allclose(result, np.zeros(2), atol=1e-14, rtol=0.0)

    def test_normalized_bequest_is_strictly_increasing_in_estate_wealth(self) -> None:
        p = _params()
        values = bequest_utility_vec(np.array([0.0, 0.05, 0.5, 2.0]), 1, p)
        self.assertTrue(np.all(np.diff(values) > 0.0))

    def test_family_size_raises_marginal_value_when_theta_n_is_positive(self) -> None:
        p = _params()
        wealth = np.array([0.4, 0.4001])
        marginal_no_children = np.diff(bequest_utility_vec(wealth, 0, p))[0]
        marginal_three_children = np.diff(bequest_utility_vec(wealth, 3, p))[0]
        self.assertGreater(marginal_three_children, marginal_no_children)
        self.assertAlmostEqual(
            marginal_three_children / marginal_no_children,
            1.0 + 3.0 * p.theta_n,
            places=10,
        )

    def test_zero_child_shifter_makes_linear_spec_child_blind(self) -> None:
        p = _params()
        p.theta_n = 0.0
        wealth = np.array([0.0, 0.2, 1.0])
        baseline = bequest_utility_vec(wealth, 0, p)
        for n_children in (1, 2, 4):
            np.testing.assert_allclose(bequest_utility_vec(wealth, n_children, p), baseline)

    def test_parent_gated_luxury_is_child_blind_on_the_intensive_margin(self) -> None:
        p = _params()
        p.bequest_spec = "parent_gated_luxury"
        wealth = np.array([0.0, 0.2, 1.0])
        np.testing.assert_allclose(bequest_utility_vec(wealth, 0, p), np.zeros(wealth.size))
        np.testing.assert_allclose(
            bequest_utility_vec(wealth, 1, p),
            bequest_utility_vec(wealth, 4, p),
        )
        self.assertGreater(float(bequest_utility_vec(np.array([1.0]), 1, p)[0]), 0.0)

    def test_parent_gate_rejects_unnormalized_crra_level(self) -> None:
        p = _params(normalized=False)
        p.bequest_spec = "parent_gated_luxury"
        with self.assertRaisesRegex(ValueError, "zero-estate normalization"):
            bequest_utility_vec(np.array([0.0, 1.0]), 1, p)

    def test_equal_division_matches_normalized_child_level_formula(self) -> None:
        p = _params(sigma=2.5)
        p.bequest_spec = "equal_division_luxury"
        wealth = np.array([0.0, 0.3, 1.2])
        n_children = 2
        expected = p.theta0 * n_children * (
            (p.theta1 + wealth / n_children) ** (1.0 - p.sigma)
            - p.theta1 ** (1.0 - p.sigma)
        ) / (1.0 - p.sigma)
        np.testing.assert_allclose(bequest_utility_vec(wealth, n_children, p), expected)
        np.testing.assert_allclose(bequest_utility_vec(wealth, 0, p), np.zeros(wealth.size))

    def test_equal_division_one_child_equals_parent_gated(self) -> None:
        p = _params()
        wealth = np.array([0.0, 0.2, 1.0])
        p.bequest_spec = "parent_gated_luxury"
        parent_gated = bequest_utility_vec(wealth, 1, p)
        p.bequest_spec = "equal_division_luxury"
        np.testing.assert_allclose(bequest_utility_vec(wealth, 1, p), parent_gated)

    def test_sigma_one_branch_uses_scaled_normalized_log_difference(self) -> None:
        p = _params(sigma=1.0)
        wealth = np.array([0.0, 0.3, 1.0])
        scale = p.theta0 * (1.0 + p.theta_n * 2)
        expected = scale * (np.log(p.theta1 + wealth) - np.log(p.theta1))
        np.testing.assert_allclose(bequest_utility_vec(wealth, 2, p), expected)

    def test_crra_branch_uses_scaled_normalized_power_difference(self) -> None:
        p = _params(sigma=2.5)
        wealth = np.array([0.0, 0.3, 1.0])
        scale = p.theta0 * (1.0 + p.theta_n * 2)
        expected = scale * (
            (p.theta1 + wealth) ** (1.0 - p.sigma)
            - p.theta1 ** (1.0 - p.sigma)
        ) / (1.0 - p.sigma)
        np.testing.assert_allclose(bequest_utility_vec(wealth, 2, p), expected)

    def test_legacy_branch_remains_available_for_fixed_parameter_ab(self) -> None:
        p = _params(sigma=2.5, normalized=False)
        wealth = np.array([0.0, 0.3, 1.0])
        scale = p.theta0 * (1.0 + p.theta_n * 2)
        expected = scale * (p.theta1 + wealth) ** (1.0 - p.sigma) / (1.0 - p.sigma)
        np.testing.assert_allclose(bequest_utility_vec(wealth, 2, p), expected)

    def test_estate_tax_is_applied_before_normalized_utility(self) -> None:
        p = _params(sigma=2.0)
        p.estate_tax_rate = 0.25
        p.estate_tax_exemption = 0.4
        gross_wealth = np.array([0.0, 0.4, 1.2])
        after_tax = np.array([0.0, 0.4, 1.0])
        scale = p.theta0
        expected = scale * (
            (p.theta1 + after_tax) ** (1.0 - p.sigma)
            - p.theta1 ** (1.0 - p.sigma)
        ) / (1.0 - p.sigma)
        np.testing.assert_allclose(bequest_utility_vec(gross_wealth, 0, p), expected)


if __name__ == "__main__":
    unittest.main()
