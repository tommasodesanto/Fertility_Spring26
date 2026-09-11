"""Compiled two-date Social Security smoke; conditional paths, not equilibria.

Run explicitly on the cluster with NUMBA_DISABLE_JIT unset or zero. One existing
tiny-grid fixture is constructed, then six two-date paths are evaluated: twelve
forward dates and twenty-four backward/forward Bellman calls. No price root,
historical data, calibration, person-law endpoint or figure is involved.
"""
from __future__ import annotations

import copy
import unittest
from unittest.mock import patch

import numpy as np
from numba import config as numba_config

# Import the module, not its TestCase into this module's namespace: otherwise
# unittest would discover and execute all of the fixture's unrelated tests.
import test_run_e5f_perfect_foresight_transition as fixture
import run_e5f_matched_pf_smoke as primitive
from e5f_social_security import fiscal_accounts


class CompiledSocialSecurityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if numba_config.DISABLE_JIT or not fixture.model.NUMBA_AVAILABLE:
            raise RuntimeError("This numerical smoke requires enabled Numba compilation")
        fixture.TinyPerfectForesightTests.setUpClass()
        source = fixture.TinyPerfectForesightTests
        cls.parameters = source.parameters
        cls.initial_g = source.stationary_g_pre.copy()
        cls.original_income = source.parameters.income.copy()
        cls.original_pension = float(source.parameters.pension)
        cls.original_tax = float(source.parameters.tau_pay)
        cls.initial_accounts = fiscal_accounts(cls.initial_g, cls.parameters)
        balanced_pension = cls.initial_accounts["implied_balanced_pension_period"]
        balanced_tax = cls.initial_accounts["implied_balanced_payroll_tax"]
        if (balanced_pension is None or not np.isfinite(balanced_pension)
                or balanced_pension <= 0 or balanced_tax is None
                or not np.isfinite(balanced_tax) or not 0 < balanced_tax < 0.9):
            raise RuntimeError("Tiny fixture cannot support both fiscal smoke instruments")

        # The last two cases change only future income relative to their own
        # balanced first-date control. Their terminal continuation stays fixed.
        cases = {
            "baseline": {},
            "explicit_baseline": dict(
                pension_path=[cls.original_pension] * 2,
                payroll_tax_path=[cls.original_tax] * 2),
            "fixed_tax": dict(
                pension_path=[balanced_pension] * 2,
                payroll_tax_path=[cls.original_tax] * 2),
            "fixed_pension": dict(
                pension_path=[cls.original_pension] * 2,
                payroll_tax_path=[balanced_tax] * 2),
            "future_pension": dict(
                pension_path=[balanced_pension, balanced_pension * 1.10],
                payroll_tax_path=[cls.original_tax] * 2),
            "future_tax": dict(
                pension_path=[cls.original_pension] * 2,
                payroll_tax_path=[balanced_tax, balanced_tax * 1.10]),
        }
        cls.results = {}
        actual_evaluate_period = fixture.calendar.evaluate_period
        rent = float(source.parameters.user_cost_rate) * source.price

        for name, fiscal_paths in cases.items():
            dated = []

            def inspect_actual_date(*args, **kwargs):
                # This wrapper calls the real Bellman-supplied/KFE evaluator;
                # no model result, policy or accounting function is mocked.
                evaluation = actual_evaluate_period(*args, **kwargs)
                P, grid, shared = args[2:5]
                budget = primitive.dated_budget(evaluation, P, shared, grid, rent)
                dated.append(dict(
                    evaluation=evaluation,
                    income=P.income.copy(),
                    pension=float(P.pension),
                    tax=float(P.tau_pay),
                    accounts=fiscal_accounts(evaluation.g_current, P),
                    pre_accounts=fiscal_accounts(evaluation.g_pre, P),
                    budget=budget,
                ))
                return evaluation

            with patch.object(fixture.calendar, "evaluate_period", side_effect=inspect_actual_date):
                path = fixture.driver.evaluate_path_at_prices(
                    prices=np.full(2, source.price),
                    psi_path=np.full(2, float(source.parameters.psi_child)),
                    terminal_price=source.price,
                    terminal_V=source.policy.V,
                    base_parameters=source.parameters,
                    b_grid=source.b_grid,
                    initial_state=copy.deepcopy(source.initial_state),
                    supply_rule=source.supply_rule,
                    birth_to_entry_conversion=source.conversion,
                    **fiscal_paths,
                )
            if len(dated) != 2:
                raise RuntimeError(f"{name} did not execute exactly two real forward dates")
            cls.results[name] = dict(path=path, dated=dated, fiscal_paths=fiscal_paths)

    def test_explicit_constant_baseline_reproduces_the_unmodified_path(self):
        baseline = self.results["baseline"]
        explicit = self.results["explicit_baseline"]
        for index in range(3):
            np.testing.assert_allclose(explicit["path"].values[index],
                                       baseline["path"].values[index], rtol=0, atol=2e-10)
        np.testing.assert_allclose(explicit["path"].terminal_state.g_pre,
                                   baseline["path"].terminal_state.g_pre, rtol=0, atol=2e-10)
        for left, right in zip(explicit["dated"], baseline["dated"], strict=True):
            np.testing.assert_allclose(left["income"], right["income"], rtol=0, atol=1e-14)
            primitive.compare_arrays(primitive.policy_arrays(left["evaluation"].policy),
                                     primitive.policy_arrays(right["evaluation"].policy))
            np.testing.assert_allclose(left["evaluation"].g_current,
                                       right["evaluation"].g_current, rtol=0, atol=2e-10)
        np.testing.assert_array_equal(self.parameters.income, self.original_income)
        self.assertEqual(self.parameters.pension, self.original_pension)
        self.assertEqual(self.parameters.tau_pay, self.original_tax)

    def test_both_instruments_balance_the_actual_first_date(self):
        for name in ("fixed_tax", "fixed_pension"):
            with self.subTest(case=name):
                result = self.results[name]
                first = result["dated"][0]
                accounts = first["accounts"]
                scale = max(accounts["payroll_tax_revenue"], accounts["pension_outlays"], 1.0)
                self.assertLessEqual(abs(accounts["pension_budget_residual"]), 2e-10 * scale)
                self.assertLessEqual(abs(accounts["scaled_pension_budget_residual"]), 2e-10)
                for key in ("payroll_tax_base_period", "retiree_benefit_exposure"):
                    self.assertAlmostEqual(accounts[key], self.initial_accounts[key], delta=2e-10)
                for key, value in accounts.items():
                    self.assertEqual(result["path"].rows[0][key], value)
                for date in result["dated"]:
                    if name == "fixed_tax":
                        self.assertEqual(date["tax"], self.original_tax)
                    else:
                        self.assertEqual(date["pension"], self.original_pension)

    def test_future_fiscal_income_changes_current_values_and_choices(self):
        for changed_name, control_name in (("future_pension", "fixed_tax"),
                                           ("future_tax", "fixed_pension")):
            with self.subTest(case=changed_name):
                changed, control = self.results[changed_name], self.results[control_name]
                np.testing.assert_array_equal(changed["dated"][0]["income"],
                                              control["dated"][0]["income"])
                self.assertGreater(float(np.max(np.abs(changed["dated"][1]["income"]
                                                      - control["dated"][1]["income"]))), 1e-8)
                occupied_pre = self.initial_g > 1e-12
                value_change = (changed["path"].values[0] - control["path"].values[0])[occupied_pre]
                self.assertTrue(np.isfinite(value_change).all())
                self.assertGreater(float(np.max(np.abs(value_change))), 1e-8)
                left, right = changed["dated"][0]["evaluation"], control["dated"][0]["evaluation"]
                occupied_current = (left.g_current + right.g_current) > 1e-12
                choice_change = max(float(np.max(np.abs(
                    getattr(left.policy, field)[occupied_current]
                    - getattr(right.policy, field)[occupied_current])))
                    for field in ("c_pol", "bp_pol"))
                self.assertGreater(choice_change, 1e-9)

    def test_every_date_passes_actual_budget_mass_and_replay_checks(self):
        for name, result in self.results.items():
            with self.subTest(case=name):
                path = result["path"]
                self.assertEqual(path.bellman_solves, 4)
                self.assertLess(path.maximum_policy_reproduction_error, 1e-12)
                self.assertLess(path.maximum_mass_accounting_error, 1e-10)
                self.assertLessEqual(path.maximum_feasibility_projection_mass, 2e-10)
                for date in result["dated"]:
                    self.assertLessEqual(date["budget"]["budget_excess_mass"], 2e-10)
                    for key in ("payroll_tax_base_period", "retiree_benefit_exposure"):
                        self.assertAlmostEqual(date["accounts"][key], date["pre_accounts"][key],
                                               delta=2e-10)
                # Changed fiscal paths are conditional fixed-price experiments;
                # no housing-market-clearing requirement is imposed on them.

    def test_entry_queue_preserves_due_vintages_and_appends_each_new_birth_flow(self):
        source = fixture.TinyPerfectForesightTests
        for name, result in self.results.items():
            with self.subTest(case=name):
                path = result["path"]
                self.assertEqual(len(path.rows), 2)
                for row in path.rows:
                    self.assertAlmostEqual(row["effective_mature_entrant_flow_B"], source.entry_flow)
                    self.assertAlmostEqual(row["raw_state_scheduled_mature_entrant_flow_B"], source.entry_flow)
                    self.assertAlmostEqual(row["entrant_flow_next"], source.entry_flow)
                for queue, field in (
                    (path.terminal_state.scheduled_entries, "birth_children_topcode_adjusted"),
                    (path.terminal_state.scheduled_raw_entries, "birth_children"),
                ):
                    expected = [source.entry_flow, source.entry_flow] + [
                        source.conversion * row[field] for row in path.rows]
                    np.testing.assert_allclose(queue, expected, rtol=0, atol=1e-14)
                self.assertAlmostEqual(float(path.terminal_state.g_pre[:, :, :, 0].sum()),
                                       source.entry_flow, delta=1e-12)

    def test_real_compiled_kernels_were_used(self):
        self.assertFalse(numba_config.DISABLE_JIT)
        self.assertTrue(fixture.model.NUMBA_AVAILABLE)
        self.assertTrue(fixture.model.full_renter_block_kernel.signatures)
        self.assertTrue(fixture.model.full_owner_block_kernel.signatures)


if __name__ == "__main__":
    unittest.main()
