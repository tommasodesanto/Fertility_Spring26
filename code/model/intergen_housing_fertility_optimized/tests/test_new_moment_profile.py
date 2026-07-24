from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from intergen_housing_fertility_optimized.audit_new_moment_jacobian import (
    MOMENT_PARAMETER_MAP,
)
from intergen_housing_fertility_optimized.calibration_search import DOMAIN, validate_contract
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_TARGETS,
    NEW_MOMENT_WEIGHTS,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.solver import (
    _weighted_ols,
    add_aggregate_wealth_bequest_flow_moments,
    markov_model_feasible_tenure_brier,
    markov_tenure_residual_variance,
)


class NewMomentProfileTests(unittest.TestCase):
    def test_timing_repaired_contract_is_runnable_and_has_fourteen_rows(self) -> None:
        validate_contract("new-moments")
        system = new_moment_target_system()
        self.assertEqual(system.count, 14)
        self.assertEqual(set(NEW_MOMENT_TARGETS), set(NEW_MOMENT_WEIGHTS))
        for name, target in NEW_MOMENT_TARGETS.items():
            self.assertAlmostEqual(NEW_MOMENT_WEIGHTS[name], 1.0 / target**2)

    def test_jacobian_mapping_covers_contract_and_domain_once(self) -> None:
        mapped_moments = [moment for moment, _parameter in MOMENT_PARAMETER_MAP]
        mapped_parameters = [parameter for _moment, parameter in MOMENT_PARAMETER_MAP]
        self.assertEqual(mapped_moments, list(NEW_MOMENT_TARGETS))
        self.assertEqual(set(mapped_parameters), {name for name, *_rest in DOMAIN})
        self.assertEqual(len(mapped_parameters), len(set(mapped_parameters)))

    def test_weighted_auxiliary_regression_recovers_coefficients(self) -> None:
        x = np.arange(1.0, 7.0)
        design = np.column_stack([np.ones(x.size), x])
        weights = np.array([1.0, 2.0, 1.0, 3.0, 2.0, 1.0])
        coefficients = _weighted_ols(2.5 + 0.4 * x, design, weights)
        np.testing.assert_allclose(coefficients, np.array([2.5, 0.4]), atol=1e-12)

    def test_aggregate_wealth_and_bequest_flow_units(self) -> None:
        p = SimpleNamespace(
            J=2,
            I=1,
            J_R=2,
            period_years=4.0,
            da=4.0,
            scale_flows_to_period=True,
            use_age_survival=True,
            survival_probs=np.array([0.5]),
            income=np.array([[4.0, 4.0]]),
            tau_pay=0.2,
            age_start=18.0,
            property_tax_lump_sum_transfer=0.0,
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
        )
        g = np.zeros((2, 2, 1, 2, 1, 1, 1))
        g[0, 0, 0, 0, 0, 0, 0] = 1.0
        g[1, 1, 0, 1, 0, 0, 0] = 1.0
        bp = np.broadcast_to(
            np.array([-1.0, 3.0]).reshape(2, 1, 1, 1, 1, 1, 1),
            g.shape,
        ).copy()
        stats = SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(
            stats,
            g,
            g,
            bp,
            p,
            np.array([-1.0, 3.0]),
            np.array([1.0]),
        )
        self.assertAlmostEqual(stats.aggregate_wealth, 4.0)
        self.assertAlmostEqual(stats.aggregate_liquid_net_worth, 2.0)
        self.assertAlmostEqual(stats.aggregate_positive_liquid_assets, 3.0)
        self.assertAlmostEqual(stats.aggregate_liquid_debt, 1.0)
        self.assertAlmostEqual(stats.aggregate_gross_housing_wealth, 2.0)
        self.assertAlmostEqual(stats.aggregate_annual_after_tax_earnings, 2.0)
        self.assertAlmostEqual(stats.aggregate_wealth_to_annual_after_tax_earnings, 2.0)
        self.assertAlmostEqual(stats.aggregate_annual_gross_labor_earnings, 2.5)
        self.assertAlmostEqual(
            stats.aggregate_wealth_to_annual_gross_labor_earnings,
            1.6,
        )
        self.assertAlmostEqual(stats.annual_bequest_flow, 1.25)
        self.assertAlmostEqual(stats.annual_bequest_flow_to_aggregate_wealth, 0.3125)

    def test_age_profile_prorates_four_year_states_at_bin_boundaries(self) -> None:
        p = SimpleNamespace(
            J=4,
            I=1,
            J_R=4,
            period_years=4.0,
            da=4.0,
            scale_flows_to_period=True,
            use_age_survival=False,
            survival_probs=np.array([]),
            income=np.full((1, 4), 4.0),
            tau_pay=0.2,
            age_start=26.0,
            property_tax_lump_sum_transfer=0.0,
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
        )
        bg = np.array([0.0, 10.0, 20.0, 30.0])
        g = np.zeros((4, 2, 1, 4, 1, 1, 1))
        for j in range(4):
            g[j, 0, 0, j, 0, 0, 0] = 1.0
        bp = np.zeros_like(g)
        stats = SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(
            stats,
            g,
            g,
            bp,
            p,
            bg,
            np.array([1.0]),
        )
        # The 34--37 state contributes one half to the 26--35 bin.
        self.assertAlmostEqual(
            stats.aggregate_wealth_to_annual_gross_labor_earnings_26_35,
            (0.0 + 10.0 + 0.5 * 20.0) / (2.5 * 1.25),
        )

    def test_bequest_flow_uses_post_saving_estate_at_death(self) -> None:
        p = SimpleNamespace(
            J=1,
            I=1,
            J_R=0,
            period_years=4.0,
            da=4.0,
            scale_flows_to_period=True,
            use_age_survival=True,
            survival_probs=np.array([]),
            income=np.array([[0.0]]),
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
        )
        wealth_g = np.zeros((2, 2, 1, 1, 1, 1, 1))
        wealth_g[1, 1, 0, 0, 0, 0, 0] = 1.0
        death_choice_g = np.zeros_like(wealth_g)
        death_choice_g[0, 1, 0, 0, 0, 0, 0] = 1.0
        bp = np.zeros_like(wealth_g)
        bp[0, 1, 0, 0, 0, 0, 0] = 7.0
        stats = SimpleNamespace()
        add_aggregate_wealth_bequest_flow_moments(
            stats,
            wealth_g,
            death_choice_g,
            bp,
            p,
            np.array([-1.0, 3.0]),
            np.array([1.0]),
        )
        self.assertAlmostEqual(stats.aggregate_wealth, 5.0)
        self.assertAlmostEqual(stats.annual_bequest_flow, 9.0 / 4.0)
        self.assertAlmostEqual(
            stats.annual_bequest_flow_to_aggregate_wealth,
            9.0 / 20.0,
        )

    def test_tenure_residual_is_expected_bernoulli_variance(self) -> None:
        p = SimpleNamespace(I=1, J=4, age_start=18.0, da=4.0)
        g = np.zeros((1, 1, 1, 4, 1, 1, 1))
        g[0, 0, 0, 2, 0, 0, 0] = 2.0
        choice = np.zeros_like(g, dtype=int)
        probs = np.zeros((*g.shape, 2))
        probs[..., 0] = 1.0
        probs[0, 0, 0, 2, 0, 0, 0, :] = (0.75, 0.25)
        value = markov_tenure_residual_variance(g, choice, probs, p)
        self.assertAlmostEqual(value, 0.25 * 0.75)

    def test_matched_tenure_brier_includes_unpredictable_choice_dispersion(self) -> None:
        p = SimpleNamespace(
            I=1,
            J=3,
            J_R=3,
            age_start=25.0,
            da=4.0,
            n_child_stages=1,
            n_parity=1,
            n_child_states=1,
            z_grid=np.array([1.0]),
            income=np.ones((1, 3)),
            property_tax_lump_sum_transfer=0.0,
        )
        g = np.ones((2, 2, 1, 3, 1, 1, 1))
        choice = np.zeros_like(g, dtype=int)
        choice[:, 1, ...] = 1
        probabilities = np.full((*g.shape, 2), 0.5)
        value = markov_model_feasible_tenure_brier(
            g,
            choice,
            probabilities,
            p,
            np.array([-0.5, 0.5]),
        )
        self.assertAlmostEqual(value, 0.25, places=10)


if __name__ == "__main__":
    unittest.main()
