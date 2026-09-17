"""No model solves: contract, residual, and failure-gate checks."""
from __future__ import annotations

import copy
from dataclasses import replace
from types import SimpleNamespace as NS
import unittest
from unittest import mock

import numpy as np

import e5f_matched_pf_endpoint as endpoint


class ExplicitEndpointTest(unittest.TestCase):
    def setUp(self):
        self.controls = endpoint.EndpointControls(
            maximum_inner_iterations=17, inner_damping=0.4,
            distribution_tolerance=1e-9, birth_rate_tolerance=1e-10,
            one_step_tolerance=1e-8, market_tolerance=2e-4,
            fiscal_absolute_tolerance=2.5e-5,
            accounting_absolute_tolerance=2e-10)
        self.g = np.ones((2, 1, 1, 1, 1, 1, 1))
        self.params = NS(I=1, J=1, Nb=2, psi_child=-0.1,
                         property_tax_lump_sum_transfer=0.2,
                         user_cost_rate=0.1, tau_H=0.04,
                         exhaustive_saving_control=True)
        self.policy = NS(price=np.array([2.0]))
        self.supply = NS(mode="static-elastic", initial_price=2.0,
                         initial_stock=3.0, elasticity=0.63)
        self.primitives = NS(start_year=2023, last_empirical_year=2100)
        self.args = dict(parameters=self.params, b_grid=np.array([0., 1.]),
                         policy=self.policy, asset_price=2.0, transfer=0.2,
                         psi_child=-0.1, demographic_primitives=self.primitives,
                         initial_g_pre=self.g, supply_rule=self.supply,
                         fiscal_regime="equal_rebate", controls=self.controls)
        self.inner = NS(
            converged=True, g_pre=self.g.copy(),
            persons=NS(persons=np.array([[3., 2.]]), heads=np.array([[1., 1.]])),
            housing_demand=3., housing_supply=3., government_budget_residual=0.,
            equal_transfer_gap=0., annual_births_per_head=0.02, renewal_ratio=0.8,
            distribution_mapping_relative_l1=0., annual_birth_rate_relative_gap=0.,
            person_one_step_relative_l1=0., head_one_step_relative_l1=0.,
            age_head_one_step_max_abs=0., household_person_head_gap=0.)

    def run_endpoint(self, **changes):
        args = dict(self.args, **changes)
        with mock.patch.object(endpoint, "_solve_inner", return_value=self.inner) as call:
            result = endpoint.evaluate_endpoint(**args)
        return result, call.call_args.kwargs

    def test_exact_inputs_and_budgets_pass_without_reconstruction_or_reanchor(self):
        before = copy.deepcopy(self.params.__dict__)
        result, kwargs = self.run_endpoint()
        self.assertTrue(result.accepted)
        self.assertTrue(result.mapping_valid)
        self.assertIs(kwargs["policy"], self.policy)
        self.assertIs(kwargs["parameters"], self.params)
        self.assertIs(kwargs["demographic_primitives"], self.primitives)
        self.assertIs(kwargs["supply_rule"], self.supply)
        self.assertEqual(kwargs["maximum_iterations"], 17)
        self.assertEqual(kwargs["damping"], .4)
        self.assertEqual(kwargs["birth_rate_tolerance"], 1e-10)
        self.assertEqual(self.params.__dict__, before)
        self.assertEqual(result.contract["choice_flags"]["exhaustive_saving_control"], True)
        self.assertEqual(result.contract["supply"]["initial_stock"], 3.)

    def test_signed_residuals_have_no_hidden_scaling(self):
        self.inner.housing_demand = 3.3
        self.inner.government_budget_residual = -.03
        self.inner.equal_transfer_gap = -.015
        result, _ = self.run_endpoint()
        np.testing.assert_allclose(result.root_residuals, [.1, -.03])
        self.assertEqual(result.root_coordinates, ("log_asset_price", "transfer"))
        self.assertAlmostEqual(result.residuals["tax_revenue"], .37)
        self.assertTrue(result.mapping_valid)
        self.assertFalse(result.accepted)

    def test_fixed_transfer_reports_surplus_without_imposing_rebate(self):
        self.inner.government_budget_residual = .4
        result, _ = self.run_endpoint(fiscal_regime="fixed_transfer")
        self.assertTrue(result.accepted)
        np.testing.assert_array_equal(result.root_residuals, [0.])
        self.assertEqual(result.residuals["fiscal_absolute"], .4)

    def test_bad_declared_inputs_fail_before_inner_work(self):
        cases = [dict(asset_price=float("nan")), dict(transfer=-1),
                 dict(psi_child=float("inf")), dict(fiscal_regime=""),
                 dict(asset_price=2.1), dict(transfer=.21), dict(psi_child=-.2),
                 dict(b_grid=np.array([1., 0.])),
                 dict(initial_g_pre=-self.g)]
        for change in cases:
            with self.subTest(change=change), mock.patch.object(endpoint, "_solve_inner") as inner:
                with self.assertRaises(ValueError):
                    endpoint.evaluate_endpoint(**dict(self.args, **change))
                inner.assert_not_called()

    def test_budget_and_gate_validation_precedes_inner_work(self):
        changes = [dict(maximum_inner_iterations=0), dict(maximum_inner_iterations=1.5),
                   dict(inner_damping=1.1), dict(market_tolerance=float("nan")),
                   dict(fiscal_absolute_tolerance=0.)]
        for change in changes:
            with self.subTest(change=change), mock.patch.object(endpoint, "_solve_inner") as inner:
                with self.assertRaises(ValueError):
                    endpoint.evaluate_endpoint(**dict(self.args, controls=replace(self.controls, **change)))
                inner.assert_not_called()

    def test_no_default_supply_or_demographic_anchor(self):
        for change in [dict(supply_rule=NS(mode="static-elastic", initial_price=2., initial_stock=3., elasticity=0.)),
                       dict(demographic_primitives=NS(start_year=2007))]:
            with self.subTest(change=change):
                with self.assertRaises(ValueError):
                    self.run_endpoint(**change)

    def test_nonfinite_output_and_failed_inner_are_not_valid_root_evaluations(self):
        for field, value in [("government_budget_residual", float("nan")),
                             ("annual_birth_rate_relative_gap", float("nan")),
                             ("converged", False),
                             ("person_one_step_relative_l1", 1e-5),
                             ("annual_birth_rate_relative_gap", 1e-9)]:
            original = getattr(self.inner, field)
            with self.subTest(field=field):
                setattr(self.inner, field, value)
                result, _ = self.run_endpoint(fiscal_regime="fixed_transfer")
                self.assertFalse(result.mapping_valid)
                self.assertFalse(result.accepted)
                setattr(self.inner, field, original)

    def test_actual_returned_population_and_supply_identities_are_gated(self):
        self.inner.persons.heads[0, 0] += 1e-5
        result, _ = self.run_endpoint()
        self.assertFalse(result.gates["returned_head_identity"])
        self.inner.persons.heads[0, 0] -= 1e-5
        self.inner.housing_supply += 1e-5
        result, _ = self.run_endpoint()
        self.assertFalse(result.gates["inherited_supply"])

    def test_inner_failure_is_propagated_without_retry(self):
        with mock.patch.object(endpoint, "_solve_inner", side_effect=RuntimeError("finite demographic root unavailable")) as call:
            with self.assertRaisesRegex(RuntimeError, "finite demographic root"):
                endpoint.evaluate_endpoint(**self.args)
        self.assertEqual(call.call_count, 1)

    def test_negative_birth_flow_is_rejected(self):
        self.inner.annual_births_per_head = -0.01
        result, _ = self.run_endpoint()
        self.assertFalse(result.mapping_valid)
        self.assertFalse(result.gates["nonnegative_flows"])


if __name__ == "__main__":
    unittest.main()
