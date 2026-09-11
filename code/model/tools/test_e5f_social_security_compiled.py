"""Compiled two-date Social Security smoke; conditional paths, not equilibria.

Run explicitly on the cluster with NUMBA_DISABLE_JIT unset or zero. One existing
tiny-grid fixture is constructed, then six two-date paths are evaluated: twelve
forward dates and twenty-four backward/forward Bellman calls. No price root,
historical data, calibration, person-law endpoint or figure is involved.

Historical evidence retained, not relabelled as passed: ``17352615`` failed
its universal occupied-control-response assertion. Diagnostic ``17353361``
then found off-support control responses. But both inherited the tiny fixture's
LEGACY nonexhaustive saving method. Stronger check ``17357939`` found wrong-signed
values and a pension action inferior to the old action under its own new
objective (gain -0.0700371545). Thus those legacy responses do not establish
correct anticipation optimization; unchanged occupied controls alone were
neither proof of failure nor a diagnosis of borrowing corners.

Diagnostic ``17358206`` set exhaustive_saving_control=True BEFORE fixture
initialization, matching the active model without changing economic primitives.
Wrong-signed values and negative crossed-policy gains disappeared. The old
pension state instead keeps saving exactly 0.9375 at a continuation knot while
its value increases. A separate knot test preserves that evidence. The new
pension probe (19,0,0,4,0,3,0) was identified in that saved active diagnostic;
it has interior saving and strict crossed-policy gains on both paths. The tax
probe is unchanged and releases a zero-saving corner into the interior. These
states are fixed BEFORE the next run, not selected from its output. Actual kernel
inputs are checked against next-date/next-age Markov expectations; conditional
Bellman values and crossed-policy incentive gains are reconstructed in NumPy.
This remains a conditional fixed-price test, not an equilibrium, occupied
response, or general optimizer-optimality certificate. Original source is
retained at a654219c and 1a5e88b7; original failure logs,
anticipation_diagnostic.json and active_optimizer_diagnostic.json remain under
output/model/e5f_matched_pf_20260909a/social_security_repair/ in the main checkout.
"""
from __future__ import annotations

import copy
import inspect
import unittest
from unittest.mock import patch

import numpy as np
from numba import config as numba_config

# Import the module, not its TestCase into this module's namespace: otherwise
# unittest would discover and execute all of the fixture's unrelated tests.
import test_run_e5f_perfect_foresight_transition as fixture
import run_e5f_matched_pf_smoke as primitive
from e5f_social_security import fiscal_accounts


# Axes: wealth index, conditional tenure, location, age index, income state,
# lifetime parity, dependent-child count. Both have zero dependent children.
ANTICIPATION_PROBES = {
    "future_pension": (19, 0, 0, 4, 0, 3, 0),
    "future_tax": (10, 1, 0, 0, 1, 0, 0),
}
PENSION_KNOT_PROBE = (1, 0, 0, 4, 0, 1, 0)
KERNEL_PROBES = dict(ANTICIPATION_PROBES,
                    pension_knot=PENSION_KNOT_PROBE,
                    pension_recipient=(1, 0, 0, 5, 0, 1, 0))


def conditional_probe_objective(inputs, state, saving):
    """Independent scalar u(c,h) + beta E[V] at a recorded feasible control.

    This evaluates two saved controls only. It never calls a model kernel or
    solves a household problem. Kernel output floors must be inactive at the
    probes; assertions below check that the reported bundles match this formula.
    """
    wealth, tenure, _, _, _, parity, children = state
    if children != 0:
        raise ValueError("The declared probes have no dependent children")
    column = parity  # flattened family order: parity + n_parity * children
    resource = float(inputs["Rv1d"][wealth]) + float(np.clip(
        inputs["gb_v"][column] - inputs["Rvt1d"][wealth],
        0.0, inputs["gb_v"][column]))
    cb, hb = float(inputs["cb_v"][column]), float(inputs["hb_v"][column])
    alpha, scale = float(inputs["alpha_v"][column]), float(inputs["esc_v"][column])
    if tenure == 0:
        rent = float(inputs["ri"])
        surplus = resource - cb - rent * hb - saving
        h_net = min((1.0 - alpha) * surplus / rent, float(inputs["hR_max"]) - hb)
        c_net = resource - cb - rent * (hb + h_net) - saving
        flow_housing = h_net
        continuation = inputs["Vc_flat"][:, column]
        spending_floor = cb + rent * hb
        collateral_floor = 0.0
    else:
        c_net = resource - float(inputs["oc"]) - cb - saving
        h_net = float(inputs["hsv"]) - float(inputs["owner_h_bar_scale"]) * hb
        flow_housing = float(inputs["owner_service_premium"]) * h_net
        continuation = inputs["Vco_flat"][:, column]
        spending_floor = cb + float(inputs["oc"])
        collateral_floor = float(inputs["bf_v"][column])
    if c_net <= 0 or flow_housing <= 0:
        raise ValueError("Declared probe control is not economically feasible")
    oms = float(inputs["oms"])
    utility = scale * (c_net**alpha * flow_housing**(1.0 - alpha))**oms / oms
    utility += float(inputs["psi_v"][column])
    expected = float(np.interp(saving, inputs["b_grid"], continuation))
    current_unsecured = float(inputs["b_grid"][wealth]) - collateral_floor
    unsecured_floor = min(float(inputs["s_next"]) * min(current_unsecured, 0.0),
                          -float(inputs["D_next"]))
    lower = max(float(inputs["b_grid"][0]), collateral_floor + unsecured_floor)
    upper = resource - spending_floor - 1e-6
    return dict(value=utility + float(inputs["beta"]) * expected, utility=utility,
                expected=expected, consumption=cb + c_net, net_consumption=c_net,
                renter_housing=hb + h_net, net_housing=h_net, lower=lower, upper=upper)


class CompiledSocialSecurityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if numba_config.DISABLE_JIT or not fixture.model.NUMBA_AVAILABLE:
            raise RuntimeError("This numerical smoke requires enabled Numba compilation")
        # Match active saving from the initial solve onward, including the
        # terminal policy. Only this numerical flag differs from the legacy
        # fixture; its failed runs remain preserved in the main output folder.
        original_overrides = fixture._tiny_overrides

        def active_overrides():
            return dict(original_overrides(), exhaustive_saving_control=True)

        with patch.object(fixture, "_tiny_overrides", side_effect=active_overrides):
            fixture.TinyPerfectForesightTests.setUpClass()
        source = fixture.TinyPerfectForesightTests
        if not bool(getattr(source.parameters, "exhaustive_saving_control", False)):
            raise RuntimeError("Compiled fiscal fixture requires the active exhaustive-saving method")
        cls.parameters = source.parameters
        cls.b_grid = source.b_grid
        cls.initial_g = source.stationary_g_pre.copy()
        cls.original_income = source.parameters.income.copy()
        cls.original_pension = float(source.parameters.pension)
        cls.original_tax = float(source.parameters.tau_pay)
        cls.initial_accounts = fiscal_accounts(cls.initial_g, cls.parameters)
        # These are properties of this unchanged diagnostic fixture, not model
        # restrictions. A changed fixture needs a newly declared probe design.
        if (cls.parameters.J != 6 or cls.parameters.J_R != 5
                or cls.parameters.I != 1 or cls.parameters.n_house != 2
                or bool(getattr(cls.parameters, "use_age_survival", False))
                or bool(getattr(cls.parameters, "readiness_gate_enabled", False))):
            raise RuntimeError("Anticipation probe contract no longer matches the tiny fixture")
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
        actual_solve_date_policy = fixture.driver.solve_date_policy
        actual_kernels = {0: fixture.model.full_renter_block_kernel,
                          1: fixture.model.full_owner_block_kernel}
        kernel_signatures = {key: inspect.signature(kernel.py_func)
                             for key, kernel in actual_kernels.items()}
        rent = float(source.parameters.user_cost_rate) * source.price

        for name, fiscal_paths in cases.items():
            dated = []
            bellman_calls = []
            active_call = {}

            def observe_kernel(kind, *args, **kwargs):
                # The Python Bellman loops age descending, then income,
                # location, and (for owners) product. Count every real call;
                # continuation/resource checks below also validate the mapping.
                ordinal = active_call["counts"][kind]
                active_call["counts"][kind] += 1
                output = actual_kernels[kind](*args, **kwargs)
                P = active_call["parameters"]
                for probe_name, state in KERNEL_PROBES.items():
                    _, tenure, location, age, income_state, _, _ = state
                    if (tenure == 0) != (kind == 0):
                        continue
                    expected_ordinal = ((P.J - 1 - age) * len(P.z_grid) + income_state) * P.I + location
                    if kind == 1:
                        expected_ordinal = expected_ordinal * P.n_house + tenure - 1
                    if ordinal == expected_ordinal:
                        bound = kernel_signatures[kind].bind(*args, **kwargs)
                        bound.apply_defaults()
                        active_call["probes"][probe_name] = dict(
                            inputs={key: value.copy() if isinstance(value, np.ndarray) else value
                                    for key, value in bound.arguments.items()},
                            output=tuple(np.asarray(value).copy() for value in output),
                        )
                return output

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

            def inspect_bellman_call(*args, **kwargs):
                # evaluate_path_at_prices calls this first in backward order
                # (date 1, date 0) and then in forward/replay order (date 0,
                # date 1).  Retain copies: the production path creates a fresh
                # dated P for each call, and no solver object is mocked.
                parameters = kwargs["P"]
                active_call.clear()
                active_call.update(parameters=parameters, counts={0: 0, 1: 0}, probes={})
                policy = actual_solve_date_policy(*args, **kwargs)
                count = parameters.J * len(parameters.z_grid) * parameters.I
                if (active_call["counts"] != {0: count, 1: count * parameters.n_house}
                        or set(active_call["probes"]) != set(KERNEL_PROBES)):
                    raise RuntimeError("Compiled kernel loop did not match declared probe mapping")
                bellman_calls.append(dict(
                    income=np.asarray(parameters.income, dtype=float).copy(),
                    pension=float(parameters.pension),
                    tax=float(parameters.tau_pay),
                    continuation=np.asarray(kwargs["continuation_V"], dtype=float).copy(),
                    value=np.asarray(policy.V, dtype=float).copy(),
                    probes=active_call["probes"],
                ))
                return policy

            with patch.object(fixture.calendar, "evaluate_period", side_effect=inspect_actual_date), \
                    patch.object(fixture.driver, "solve_date_policy", side_effect=inspect_bellman_call), \
                    patch.object(fixture.model, "full_renter_block_kernel",
                                 side_effect=lambda *a, **k: observe_kernel(0, *a, **k)), \
                    patch.object(fixture.model, "full_owner_block_kernel",
                                 side_effect=lambda *a, **k: observe_kernel(1, *a, **k)):
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
            if len(bellman_calls) != 4:
                raise RuntimeError(f"{name} did not execute exactly two backward and two forward Bellman calls")
            cls.results[name] = dict(
                path=path, dated=dated, fiscal_paths=fiscal_paths,
                bellman_calls=bellman_calls,
            )

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

    def test_dated_fiscal_income_and_age_support_reach_both_bellman_passes(self):
        """Check the dated chain and the exact one-step age support of shocks."""
        for changed_name, control_name in (("future_pension", "fixed_tax"),
                                           ("future_tax", "fixed_pension")):
            with self.subTest(case=changed_name):
                changed, control = self.results[changed_name], self.results[control_name]
                changed_calls = changed["bellman_calls"]
                control_calls = control["bellman_calls"]
                # Calls are backward date 1/date 0, then forward date 0/date
                # 1.  This proves the bound income reaches each evaluator,
                # rather than only the KFE/accounting call observed below.
                for result, calls in ((changed, changed_calls), (control, control_calls)):
                    for call_index, date_index in enumerate((1, 0, 0, 1)):
                        np.testing.assert_array_equal(
                            calls[call_index]["income"], result["dated"][date_index]["income"])
                        self.assertEqual(calls[call_index]["pension"], result["dated"][date_index]["pension"])
                        self.assertEqual(calls[call_index]["tax"], result["dated"][date_index]["tax"])
                    # The backward date-zero and forward date-zero solves use
                    # exactly the date-one backward value as continuation.
                    np.testing.assert_array_equal(calls[1]["continuation"], calls[0]["value"])
                    np.testing.assert_array_equal(calls[2]["continuation"], calls[0]["value"])
                    np.testing.assert_array_equal(calls[1]["value"], calls[2]["value"])
                    np.testing.assert_array_equal(calls[3]["value"], calls[0]["value"])
                    np.testing.assert_array_equal(calls[3]["continuation"], calls[0]["continuation"])
                    # Independent pension/payroll reconstruction: these are
                    # period-unit incomes before household income-state scaling.
                    for call in calls:
                        P = self.parameters
                        period_scale = (float(P.period_years) if P.scale_flows_to_period else 1.0)
                        expected_income = np.full((P.I, P.J), call["pension"])
                        expected_income[:, :P.J_R] = (period_scale * (1.0 - call["tax"])
                            * np.asarray(P.w_hat)[:, None]
                            * np.asarray(P.income_age_profile)[None, :P.J_R])
                        np.testing.assert_allclose(call["income"], expected_income, rtol=0, atol=1e-14)
                        # Inspect a retiree as well as the working anticipation
                        # states in ALL four real calls. Thus both pensions and
                        # taxes are checked inside the compiled resource inputs.
                        for probe_name, state in KERNEL_PROBES.items():
                            _, _, location, age, income_state, _, _ = state
                            z = float(P.z_grid[income_state])
                            multiplier = (z if age < P.J_R else
                                1.0 + float(P.retirement_income_z_scale) * (z - 1.0))
                            inputs = call["probes"][probe_name]["inputs"]
                            np.testing.assert_allclose(inputs["Rv1d"], P.R_gross * self.b_grid
                                + expected_income[location, age] * multiplier, rtol=0, atol=1e-14)
                            self.assertEqual(inputs["beta"], P.beta)
                np.testing.assert_array_equal(changed["dated"][0]["income"],
                                              control["dated"][0]["income"])
                np.testing.assert_array_equal(changed_calls[0]["continuation"],
                                              control_calls[0]["continuation"])
                self.assertGreater(float(np.max(np.abs(changed["dated"][1]["income"]
                                                      - control["dated"][1]["income"]))), 1e-8)
                date_one_value_change = np.max(np.abs(
                    changed_calls[0]["value"] - control_calls[0]["value"]
                ))
                self.assertGreater(float(date_one_value_change), 1e-8)
                # The sole date-zero difference is the anticipated date-one
                # value: terminal values and current dated income are shared.
                continuation_change = np.max(np.abs(
                    changed_calls[1]["continuation"] - control_calls[1]["continuation"]
                ))
                self.assertGreater(float(continuation_change), 1e-8)
                occupied_pre = self.initial_g > 1e-12
                value_change = (changed["path"].values[0] - control["path"].values[0])[occupied_pre]
                self.assertTrue(np.isfinite(value_change).all())
                self.assertGreater(float(np.max(np.abs(value_change))), 1e-8)
                np.testing.assert_array_equal(changed_calls[1]["value"], changed["path"].values[0])
                # At t=1 only current recipients/payers can respond because
                # the terminal continuation is fixed. At t=0 only households
                # whose NEXT age receives/pays at t=1 can anticipate the shock.
                future_affected = {5} if changed_name == "future_pension" else set(range(5))
                current_affected = {4} if changed_name == "future_pension" else set(range(4))
                sign = 1.0 if changed_name == "future_pension" else -1.0
                for date, affected in ((1, future_affected), (0, current_affected)):
                    delta = changed["path"].values[date] - control["path"].values[date]
                    for age in range(self.parameters.J):
                        if age not in affected:
                            np.testing.assert_allclose(delta[:, :, :, age], 0.0, rtol=0, atol=2e-12)
                    self.assertGreater(float(np.max(sign * delta)), 1e-8)
                    self.assertGreaterEqual(float(np.min(sign * delta)), -2e-10)

    def test_predeclared_controls_satisfy_exact_continuation_and_crossed_policy_incentives(self):
        """Verify the optimization channel at two declared, feasible states.

        Active diagnostic 17358206 pension saving: 7.6287956948 -> 7.5109834889,
        interior on both paths, with crossed gains 4.62546e-5 and 4.53927e-5.
        Tax saving: 0 -> 0.8376618080, a release from its lower bound, with
        crossed gains 0.00113803 and 0.00169014. The pension renter's housing
        cap binds; its saving control remains interior. No universal interior
        or occupied response is assumed.
        """
        P = self.parameters
        for changed_name, control_name in (("future_pension", "fixed_tax"),
                                           ("future_tax", "fixed_pension")):
            state = ANTICIPATION_PROBES[changed_name]
            wealth, tenure, location, age, income_state, parity, children = state
            with self.subTest(case=changed_name, state=state):
                records = []
                for name in (control_name, changed_name):
                    result = self.results[name]
                    # Check both actual current-date calls, not only that
                    # the enclosing function received the right full array.
                    for call_index in (1, 2):
                        call = result["bellman_calls"][call_index]
                        probe = call["probes"][changed_name]
                        inputs, output = probe["inputs"], probe["output"]
                        self.assertEqual(inputs["has_prev"], 0)
                        # E[V] = sum_z' Pi_z[z,z'] sum_m' Pi_child[m,m';n]
                        #             V_{t+1}(b',tenure,location,age+1,z',n,m').
                        # No survival/bequest mixture or readiness gate is
                        # active in this explicit fixture contract.
                        z_weights = np.asarray(P.Pi_z[income_state], dtype=float)
                        z_weights = z_weights / z_weights.sum()
                        child_weights = np.asarray(P.Pi_child[children, :, parity], dtype=float)
                        expected = np.zeros(len(self.b_grid))
                        for znext, probability_z in enumerate(z_weights):
                            for child_next, probability_child in enumerate(child_weights):
                                expected += probability_z * probability_child * call["continuation"][
                                    :, tenure, location, age + 1, znext, parity, child_next]
                        continuation_key = "Vc_flat" if tenure == 0 else "Vco_flat"
                        np.testing.assert_allclose(inputs[continuation_key][:, parity], expected,
                                                   rtol=0, atol=2e-12)
                        # Prove fiscal income reaches the compiled optimizer's
                        # resources, including its income-state multiplier.
                        income = call["income"][location, age] * float(P.z_grid[income_state])
                        np.testing.assert_allclose(inputs["Rv1d"], P.R_gross * self.b_grid + income,
                                                   rtol=0, atol=1e-14)
                        saving = float(output[1][wealth, parity])
                        objective = conditional_probe_objective(inputs, state, saving)
                        np.testing.assert_allclose(output[0][wealth, parity], objective["value"],
                                                   rtol=0, atol=2e-12)
                        self.assertGreater(objective["net_consumption"], inputs["c_min"])
                        np.testing.assert_allclose(output[2][wealth, parity], objective["consumption"],
                                                   rtol=0, atol=2e-12)
                        if tenure == 0:
                            self.assertGreater(objective["net_housing"], 0.01)
                            np.testing.assert_allclose(output[3][wealth, parity], objective["renter_housing"],
                                                       rtol=0, atol=2e-12)
                        self.assertAlmostEqual(saving, result["dated"][0]["evaluation"].policy.bp_pol[state],
                                               delta=2e-12)
                        if call_index == 1:
                            records.append((inputs, saving, objective))

                before_inputs, before, before_value = records[0]
                after_inputs, after, after_value = records[1]
                continuation_key = "Vc_flat" if tenure == 0 else "Vco_flat"
                # Every current primitive and feasible-set argument is equal;
                # only the date-one continuation presented to the kernel differs.
                for key in before_inputs:
                    if key != continuation_key:
                        np.testing.assert_array_equal(before_inputs[key], after_inputs[key])
                self.assertGreater(float(np.max(np.abs(
                    after_inputs[continuation_key][:, parity] - before_inputs[continuation_key][:, parity]))), 1e-8)
                self.assertGreater(abs(after - before), 1e-8)
                for saving, objective in ((after, after_value), (before, before_value)):
                    self.assertGreaterEqual(saving, objective["lower"] - 1e-12)
                    self.assertLess(saving, min(objective["upper"], self.b_grid[-1]) - 1e-6)
                self.assertGreater(after, after_value["lower"] + 1e-6)
                if changed_name == "future_pension":
                    self.assertGreater(before, before_value["lower"] + 1e-6)
                    self.assertLess(after, before)
                else:
                    self.assertAlmostEqual(before, before_value["lower"], delta=1e-10)
                    self.assertGreater(after, before)

                old_control_new_income = conditional_probe_objective(after_inputs, state, before)
                new_control_old_income = conditional_probe_objective(before_inputs, state, after)
                # Strict revealed-preference reversals make these economically
                # informative control tests, not arbitrary nonzero-array tests.
                self.assertGreater(after_value["value"] - old_control_new_income["value"], 1e-10)
                self.assertGreater(before_value["value"] - new_control_old_income["value"], 1e-10)
                # At fixed saving, the ENTIRE value change is beta times the
                # anticipated continuation change. Current utility is identical.
                self.assertAlmostEqual(old_control_new_income["utility"], before_value["utility"], delta=2e-12)
                self.assertAlmostEqual(old_control_new_income["value"] - before_value["value"],
                    P.beta * (old_control_new_income["expected"] - before_value["expected"]), delta=2e-12)
                self.assertAlmostEqual(after_value["value"] - before_value["value"],
                    after_value["utility"] - before_value["utility"]
                    + P.beta * (after_value["expected"] - before_value["expected"]), delta=2e-12)

    def test_preserved_pension_knot_has_exact_continuation_value_channel(self):
        """A control can stay at a continuation knot while its value responds.

        In active diagnostic 17358206 both choices equal 0.9375. The unchanged
        choice is admissible here; there is no strict saving-response assertion.
        The strict interior response remains required at the separate declared
        pension state above. This probe keeps exact continuation/value checks
        and weak optimality against the crossed control, even if it later moves.
        """
        P = self.parameters
        state = PENSION_KNOT_PROBE
        wealth, tenure, location, age, income_state, parity, children = state
        for call_index in (1, 2):
            with self.subTest(bellman_call=call_index):
                records = []
                for name in ("fixed_tax", "future_pension"):
                    result = self.results[name]
                    call = result["bellman_calls"][call_index]
                    probe = call["probes"]["pension_knot"]
                    inputs, output = probe["inputs"], probe["output"]
                    expected = np.zeros(len(self.b_grid))
                    probabilities = np.asarray(P.Pi_z[income_state], dtype=float)
                    probabilities = probabilities / probabilities.sum()
                    for znext, probability_z in enumerate(probabilities):
                        for child_next, probability_child in enumerate(P.Pi_child[children, :, parity]):
                            expected += probability_z * probability_child * call["continuation"][
                                :, tenure, location, age + 1, znext, parity, child_next]
                    np.testing.assert_allclose(inputs["Vc_flat"][:, parity], expected,
                                               rtol=0, atol=2e-12)
                    saving = float(output[1][wealth, parity])
                    objective = conditional_probe_objective(inputs, state, saving)
                    self.assertAlmostEqual(float(output[0][wealth, parity]), objective["value"], delta=2e-12)
                    self.assertAlmostEqual(saving, result["dated"][0]["evaluation"].policy.bp_pol[state],
                                           delta=2e-12)
                    self.assertGreaterEqual(saving, objective["lower"] - 1e-12)
                    self.assertLess(saving, min(objective["upper"], self.b_grid[-1]) - 1e-6)
                    records.append((inputs, saving, objective))
                before_inputs, before, before_value = records[0]
                after_inputs, after, after_value = records[1]
                for key in before_inputs:
                    if key != "Vc_flat":
                        np.testing.assert_array_equal(before_inputs[key], after_inputs[key])
                self.assertLessEqual(float(np.min(np.abs(self.b_grid - before))), 1e-12)
                fixed_control = conditional_probe_objective(after_inputs, state, before)
                self.assertAlmostEqual(fixed_control["utility"], before_value["utility"], delta=2e-12)
                self.assertGreater(fixed_control["expected"] - before_value["expected"], 1e-8)
                self.assertAlmostEqual(fixed_control["value"] - before_value["value"],
                    P.beta * (fixed_control["expected"] - before_value["expected"]), delta=2e-12)
                crossed_control = conditional_probe_objective(before_inputs, state, after)
                self.assertGreaterEqual(after_value["value"] - fixed_control["value"], -2e-12)
                self.assertGreaterEqual(before_value["value"] - crossed_control["value"], -2e-12)
                self.assertAlmostEqual(after_value["value"] - before_value["value"],
                    after_value["utility"] - before_value["utility"]
                    + P.beta * (after_value["expected"] - before_value["expected"]), delta=2e-12)
                if before == after:
                    self.assertAlmostEqual(after_value["value"] - before_value["value"],
                        P.beta * (after_value["expected"] - before_value["expected"]), delta=2e-12)

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
        self.assertTrue(self.parameters.exhaustive_saving_control)
        self.assertTrue(fixture.model.full_renter_block_kernel.signatures)
        self.assertTrue(fixture.model.full_owner_block_kernel.signatures)
        for result in self.results.values():
            for call in result["bellman_calls"]:
                for probe in call["probes"].values():
                    self.assertEqual(probe["inputs"]["exhaustive_saving"], 1)


if __name__ == "__main__":
    unittest.main()
