"""Independent household fiscal accounting checks; no Bellman or model solve."""

import copy
from types import SimpleNamespace
import unittest

import numpy as np

from e5f_social_security import (
    apply_fiscal_date,
    bind_social_security_income,
    fiscal_accounts,
    validated_fiscal_paths,
)


def parameters():
    """Two locations, two working ages and two retired ages in period units."""
    return SimpleNamespace(
        Nb=2, n_house=1, I=2, J=4, J_R=2, n_parity=2, n_child_states=2,
        period_years=4.0, da=4.0, age_start=18.0,
        scale_flows_to_period=True,
        income_age_breaks=np.array([18.0, 22.0]),
        income_age_values=np.array([0.5, 1.5]),
        income_age_profile=np.array([0.5, 1.5, 1.5, 1.5]),
        normalize_income_profile=False,
        w_hat=np.array([1.0, 2.0]),
        z_grid=np.array([0.5, 1.5]),
        # Deliberately different from actual occupied-state weights.
        z_weights=np.array([0.9, 0.1]),
        Pi_z=np.array([[0.0, 1.0], [1.0, 0.0]]),
        retirement_income_z_scale=0.5,
        tau_pay=0.2, pension=3.0, pension_mode="manual",
        pension_by_loc=np.array([3.0, 3.0]),
        income=np.array([[1.6, 4.8, 3.0, 3.0], [3.2, 9.6, 3.0, 3.0]]),
        property_tax_lump_sum_transfer=7.0,
    )


def distribution():
    """Axes: wealth, tenure, location, age, current z, parity, child state."""
    g = np.zeros((2, 2, 2, 4, 2, 2, 2))
    g[0, 0, 0, 0, 0, 0, 0] = 2.0
    g[1, 1, 1, 1, 1, 1, 1] = 1.0
    g[0, 1, 0, 1, 1, 0, 1] = 0.5
    g[0, 1, 0, 2, 0, 1, 1] = 3.0
    g[1, 0, 1, 3, 1, 0, 0] = 1.0
    return g


class FiscalAccountingTests(unittest.TestCase):
    def test_hand_calculated_current_income_and_retirement_exposure(self):
        result = fiscal_accounts(distribution(), parameters())
        # Annual payroll: .5 + 4.5 + 1.125 = 6.125; four-year payroll = 24.5.
        # Retiree exposure: 3 * .75 + 1 * 1.25 = 3.5 household benefits.
        expected = dict(
            payroll_tax_base_period=24.5,
            retiree_benefit_exposure=3.5,
            payroll_tax_revenue=4.9,
            pension_outlays=10.5,
            pension_budget_residual=-5.6,
            scaled_pension_budget_residual=-5.6 / 10.5,
            implied_balanced_pension_period=1.4,
            implied_balanced_payroll_tax=3.0 / 7.0,
        )
        for key, value in expected.items():
            with self.subTest(key=key):
                self.assertAlmostEqual(result[key], value, places=13)

    def test_flat_retirement_benefits_use_actual_retired_head_mass(self):
        P = parameters()
        P.retirement_income_z_scale = 0.0
        result = fiscal_accounts(distribution(), P)
        self.assertEqual(result["retiree_benefit_exposure"], 4.0)
        self.assertEqual(result["pension_outlays"], 12.0)
        self.assertAlmostEqual(result["implied_balanced_pension_period"], 1.225)

    def test_transfers_and_unoccupied_reference_weights_do_not_enter_payroll(self):
        P = parameters()
        baseline = fiscal_accounts(distribution(), P)
        P.property_tax_lump_sum_transfer = 999.0
        P.z_weights = np.array([0.0, 1.0])
        P.Pi_z = np.eye(2)
        self.assertEqual(fiscal_accounts(distribution(), P), baseline)

    def test_population_scaling_preserves_implied_policy(self):
        g, P = distribution(), parameters()
        baseline = fiscal_accounts(g, P)
        scaled = fiscal_accounts(7.0 * g, P)
        for key in ("payroll_tax_base_period", "retiree_benefit_exposure",
                    "payroll_tax_revenue", "pension_outlays", "pension_budget_residual"):
            with self.subTest(key=key):
                self.assertAlmostEqual(scaled[key], 7.0 * baseline[key])
        for key in ("implied_balanced_pension_period", "implied_balanced_payroll_tax",
                    "scaled_pension_budget_residual"):
            with self.subTest(key=key):
                self.assertAlmostEqual(scaled[key], baseline[key])

    def test_redistribution_of_nonincome_states_preserves_budget(self):
        g, P = distribution(), parameters()
        redistributed = np.zeros_like(g)
        redistributed[1, 1, :, :, :, 1, 1] = g.sum(axis=(0, 1, 5, 6))
        self.assertEqual(fiscal_accounts(redistributed, P), fiscal_accounts(g, P))

    def test_period_scaling_is_applied_once(self):
        g, P = distribution(), parameters()
        period_result = fiscal_accounts(g, P)
        P.scale_flows_to_period = False
        P.pension /= 4.0
        unit_result = fiscal_accounts(g, P)
        for key in ("payroll_tax_base_period", "payroll_tax_revenue", "pension_outlays",
                    "pension_budget_residual", "implied_balanced_pension_period"):
            with self.subTest(key=key):
                self.assertAlmostEqual(period_result[key], 4.0 * unit_result[key])
        self.assertEqual(period_result["retiree_benefit_exposure"],
                         unit_result["retiree_benefit_exposure"])
        self.assertEqual(period_result["implied_balanced_payroll_tax"],
                         unit_result["implied_balanced_payroll_tax"])

    def test_either_implied_instrument_closes_the_hand_calculated_budget(self):
        g, original = distribution(), parameters()
        accounts = fiscal_accounts(g, original)
        for kwargs in (
            dict(pension_period=accounts["implied_balanced_pension_period"]),
            dict(payroll_tax=accounts["implied_balanced_payroll_tax"]),
        ):
            with self.subTest(kwargs=kwargs):
                P = copy.deepcopy(original)
                bind_social_security_income(P, **kwargs)
                self.assertAlmostEqual(fiscal_accounts(g, P)["pension_budget_residual"], 0.0)

    def test_zero_exposure_does_not_manufacture_a_balancing_pension(self):
        g = distribution()
        g[:, :, :, 2:, :, :, :] = 0.0
        result = fiscal_accounts(g, parameters())
        self.assertIsNone(result["implied_balanced_pension_period"])
        self.assertEqual(result["implied_balanced_payroll_tax"], 0.0)
        self.assertGreater(result["pension_budget_residual"], 0.0)

    def test_zero_payroll_does_not_manufacture_a_balancing_tax(self):
        g = distribution()
        g[:, :, :, :2, :, :, :] = 0.0
        result = fiscal_accounts(g, parameters())
        self.assertEqual(result["implied_balanced_pension_period"], 0.0)
        self.assertIsNone(result["implied_balanced_payroll_tax"])
        self.assertLess(result["pension_budget_residual"], 0.0)

    def test_empty_population_has_zero_budget_and_undefined_ratios(self):
        result = fiscal_accounts(np.zeros_like(distribution()), parameters())
        self.assertEqual(result["pension_budget_residual"], 0.0)
        self.assertEqual(result["scaled_pension_budget_residual"], 0.0)
        self.assertIsNone(result["implied_balanced_pension_period"])
        self.assertIsNone(result["implied_balanced_payroll_tax"])

    def test_invalid_distribution_is_rejected(self):
        for value in (float("nan"), float("inf"), -0.1):
            with self.subTest(value=value):
                g = distribution()
                g[0, 0, 0, 0, 0, 0, 0] = value
                with self.assertRaises(ValueError):
                    fiscal_accounts(g, parameters())
        with self.assertRaises(ValueError):
            fiscal_accounts(distribution().sum(axis=4), parameters())

    def test_retirement_boundaries_allow_all_workers_or_all_retirees(self):
        P = parameters()
        P.J_R = P.J
        result = fiscal_accounts(distribution(), P)
        self.assertEqual(result["payroll_tax_base_period"], 51.5)
        self.assertEqual(result["retiree_benefit_exposure"], 0.0)
        P.J_R = 0
        result = fiscal_accounts(distribution(), P)
        self.assertEqual(result["payroll_tax_base_period"], 0.0)
        self.assertEqual(result["retiree_benefit_exposure"], 6.875)

    def test_mismatched_income_dimensions_and_invalid_multipliers_fail(self):
        for changes in (
            dict(I=3), dict(J=5), dict(J_R=5),
            dict(z_grid=np.array([0.5])),
            dict(z_grid=np.array([0.0, 1.5])),
            dict(w_hat=np.array([1.0, float("nan")])),
            dict(income_age_profile=np.array([0.5, 1.5])),
            dict(retirement_income_z_scale=-0.1),
            dict(retirement_income_z_scale=2.0),
            dict(retirement_income_z_scale=float("nan")),
            dict(period_years=0.0),
        ):
            with self.subTest(changes=changes):
                P = parameters()
                for key, value in changes.items():
                    setattr(P, key, value)
                with self.assertRaises(ValueError):
                    fiscal_accounts(distribution(), P)

    def test_finite_inputs_cannot_silently_overflow_fiscal_aggregates(self):
        g = np.full_like(distribution(), 1e308)
        with np.errstate(over="ignore", invalid="ignore"), self.assertRaises(ValueError):
            fiscal_accounts(g, parameters())


class DatedIncomeTests(unittest.TestCase):
    def test_binder_is_idempotent_with_existing_period_manual_pension(self):
        for mode in ("manual", "custom", "balanced_stationary"):
            with self.subTest(mode=mode):
                P = parameters()
                P.pension_mode = mode
                returned = bind_social_security_income(P, pension_period=1.4, payroll_tax=0.25)
                self.assertIs(returned, P)
                expected_income = np.array([[1.5, 4.5, 1.4, 1.4],
                                            [3.0, 9.0, 1.4, 1.4]])
                np.testing.assert_allclose(P.income, expected_income, rtol=0, atol=1e-14)
                np.testing.assert_array_equal(P.pension_by_loc, [1.4, 1.4])
                self.assertEqual(P.pension, 1.4)
                self.assertEqual(P.tau_pay, 0.25)
                self.assertEqual(P.property_tax_lump_sum_transfer, 7.0)
                bind_social_security_income(P, pension_period=1.4, payroll_tax=0.25)
                np.testing.assert_allclose(P.income, expected_income, rtol=0, atol=1e-14)
                self.assertEqual(P.pension, 1.4)

    def test_partial_updates_preserve_the_other_instrument(self):
        P = parameters()
        bind_social_security_income(P, pension_period=1.4)
        self.assertEqual(P.tau_pay, 0.2)
        np.testing.assert_allclose(P.income[:, :2], [[1.6, 4.8], [3.2, 9.6]])
        bind_social_security_income(P, payroll_tax=0.25)
        self.assertEqual(P.pension, 1.4)
        np.testing.assert_array_equal(P.income[:, 2:], np.full((2, 2), 1.4))

    def test_binder_respects_unscaled_flows_and_zero_benefits(self):
        P = parameters()
        P.scale_flows_to_period = False
        bind_social_security_income(P, pension_period=0.0, payroll_tax=0.0)
        np.testing.assert_array_equal(P.income, [[0.5, 1.5, 0.0, 0.0],
                                                [1.0, 3.0, 0.0, 0.0]])

    def test_invalid_instrument_updates_are_rejected(self):
        for value in (-0.1, float("nan"), float("inf")):
            with self.subTest(pension=value), self.assertRaises(ValueError):
                bind_social_security_income(parameters(), pension_period=value)
        for value in (-0.1, 1.0, 1.1, float("nan"), float("inf")):
            with self.subTest(tax=value), self.assertRaises(ValueError):
                bind_social_security_income(parameters(), payroll_tax=value)


class FiscalPathTests(unittest.TestCase):
    def test_optional_paths_and_no_path_do_not_change_income(self):
        pensions, taxes = validated_fiscal_paths(3)
        self.assertIsNone(pensions)
        self.assertIsNone(taxes)
        P = parameters()
        before = copy.deepcopy(vars(P))
        apply_fiscal_date(P, 1, pensions, taxes)
        self.assertEqual(set(vars(P)), set(before))
        for key, value in before.items():
            with self.subTest(key=key):
                if isinstance(value, np.ndarray):
                    np.testing.assert_array_equal(getattr(P, key), value)
                else:
                    self.assertEqual(getattr(P, key), value)

    def test_dated_binding_selects_the_same_explicit_units(self):
        pensions, taxes = validated_fiscal_paths(3, [1.0, 1.4, 0.0], [0.1, 0.25, 0.0])
        P = parameters()
        apply_fiscal_date(P, 1, pensions, taxes)
        self.assertEqual(P.pension, 1.4)
        self.assertEqual(P.tau_pay, 0.25)
        np.testing.assert_allclose(P.income, [[1.5, 4.5, 1.4, 1.4],
                                             [3.0, 9.0, 1.4, 1.4]])

    def test_single_instrument_paths_preserve_the_other_instrument(self):
        for kwargs, expected_pension, expected_tax in (
            (dict(pension_path=[1.4, 2.0]), 2.0, 0.2),
            (dict(payroll_tax_path=[0.1, 0.25]), 3.0, 0.25),
        ):
            with self.subTest(kwargs=kwargs):
                paths = validated_fiscal_paths(2, **kwargs)
                P = parameters()
                apply_fiscal_date(P, 1, *paths)
                self.assertEqual(P.pension, expected_pension)
                self.assertEqual(P.tau_pay, expected_tax)

    def test_invalid_path_lengths_dimensions_and_values_are_rejected(self):
        invalid = ([1.0], [[1.0, 1.0]], 1.0, [float("nan"), 1.0],
                   [1.0, float("inf")], [-0.1, 0.1])
        for path in invalid:
            with self.subTest(pension_path=path), self.assertRaises(ValueError):
                validated_fiscal_paths(2, pension_path=path)
        for path in invalid + ([0.0, 1.0], [0.0, 1.1]):
            with self.subTest(payroll_tax_path=path), self.assertRaises(ValueError):
                validated_fiscal_paths(2, payroll_tax_path=path)

    def test_validated_paths_are_owned_copies(self):
        pension_input = np.array([1.0, 1.4])
        tax_input = np.array([0.1, 0.2])
        pensions, taxes = validated_fiscal_paths(2, pension_input, tax_input)
        pension_input[:] = 999.0
        tax_input[:] = 0.99
        np.testing.assert_array_equal(pensions, [1.0, 1.4])
        np.testing.assert_array_equal(taxes, [0.1, 0.2])


if __name__ == "__main__":
    unittest.main()
