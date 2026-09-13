"""Pure and routing tests for the isolated rebated-surprise adapter."""
from dataclasses import dataclass
from types import SimpleNamespace as NS
import unittest
from unittest.mock import patch

import numpy as np

import e5f_rebated_surprises as r


@dataclass
class Conditioning:
    start_year: int = 2007
    observer: object = None
    next_age_targets: object = None
    initial_mass: float = 1.


class AccountingTest(unittest.TestCase):
    def test_all_property_tax_revenue_is_rebated_over_all_heads(self):
        ledger = r.rebated_tax_accounts(
            property_tax_revenue=12., transfer_per_head=2., head_mass=6.)
        self.assertEqual(ledger["equal_transfer_outlays"], 12.)
        self.assertEqual(ledger["government_budget_residual"], 0.)
        self.assertEqual(ledger["scaled_government_budget_residual"], 0.)
        self.assertEqual(ledger["implied_equal_transfer"], 2.)

    def test_three_signed_residuals_have_policy_scaling(self):
        residual = r.dated_residual(
            demand=11., supply=10.,
            payroll_accounts={"payroll_tax_revenue": 8., "pension_outlays": 10.},
            tax_accounts={"property_tax_revenue": 12., "equal_transfer_outlays": 8.})
        np.testing.assert_allclose(residual, [.1, -40., 200. * 4. / 12.])

    def test_invalid_head_mass_cannot_fake_a_balanced_rebate(self):
        with self.assertRaisesRegex(ValueError, "positive head mass"):
            r.rebated_tax_accounts(
                property_tax_revenue=0., transfer_per_head=0., head_mass=0.)

    def test_root_residual_order_matches_coordinate_blocks(self):
        stacked = r.stack_dated_residuals([[1., 10., 100.], [2., 20., 200.]])
        np.testing.assert_array_equal(stacked, [1., 2., 10., 20., 100., 200.])

    def test_closed_boundary_needs_no_stationary_fixed_point(self):
        terminal = NS(parameters=NS(), policy=NS(price=np.array([2.]), V="V"),
                      g_pre=np.array([1.]))
        parameters, policy, price, state = r._terminal_parts(terminal)
        self.assertIs(parameters, terminal.parameters)
        self.assertIs(policy.V, terminal.policy.V)
        self.assertEqual(price, 2.)
        self.assertIsNone(state)


class ReplayTest(unittest.TestCase):
    def test_first_period_replays_accepted_transfer_next_price_and_value(self):
        inherited = r.InheritedState(2007, NS(g_pre=np.array([1.])))
        next_house = NS(g_pre=np.array([2.]))
        row = {key: 1. for key in (
            "asset_price", "renter_price", "housing_demand", "owner_rate",
            "birth_children_topcode_adjusted", "pension_period_units",
            "payroll_tax_revenue", "pension_outlays", "property_tax_revenue",
            "equal_transfer_outlays")}
        path = NS(values=[np.array([3.]), np.array([4.])], rows=[row])
        seen = []

        def replay(**kwargs):
            seen.append(kwargs)
            return NS(values=[path.values[0]], rows=[row], terminal_state=next_house,
                maximum_mass_accounting_error=0.,
                maximum_policy_reproduction_error=0.,
                maximum_feasibility_projection_mass=0.)

        joined = NS(pf=NS(evaluate_path_at_prices=replay), person_pf=NS())
        old = NS(parameters=NS(), b_grid=np.array([0.]), supply_rule=object(),
                 historical_conditioning=Conditioning(next_age_targets={}))
        with patch.object(r, "_runtime", return_value=(None, joined, None, None, None)):
            state = r.first_period_state(
                inherited=inherited, old_state=old, demographics=None, path=path,
                prices=[2., 5.], pensions=[3., 4.], transfers=[.7, .8], psi=.2)
        self.assertEqual(seen[0]["transfer_path"], [.7])
        self.assertEqual(seen[0]["terminal_price"], 5.)
        self.assertIs(seen[0]["terminal_V"], path.values[1])
        self.assertEqual(seen[0]["payroll_tax_path"], [.179])
        self.assertIs(state.households, next_house)


class ForecastRoutingTest(unittest.TestCase):
    def test_explicit_transfers_reach_history_backward_tail_and_hook(self):
        n, h = 6, 4
        g = np.ones((1, 1, 1, 1, 1, 1, 1))
        people = NS(year=2023, persons=np.array([1.]), heads=np.array([1.]))
        people.validated = lambda: people
        calls = {}

        def backward(**kwargs):
            calls["backward"] = kwargs
            return [np.array([float(i)]) for i in range(n - h + 1)], n - h

        def history(**kwargs):
            calls["history"] = kwargs
            return NS(rows=[dict(calendar_year=2007 + 4 * i, period=i) for i in range(h)],
                      values=[np.array([float(i)]) for i in range(h + 1)],
                      terminal_state=NS(g_pre=g), bellman_solves=h)

        def hook(**kwargs):
            calls["tail"] = kwargs
            return NS(rows=[dict(calendar_year=2023 + 4 * i, period=i)
                            for i in range(n - h)],
                      values=kwargs["precomputed_value_path"], bellman_solves=n - h)

        joined = NS(
            pf=NS(backward_value_path=backward,
                  rents_from_asset_prices=lambda p, *a: np.asarray(p) * .1,
                  evaluate_path_at_prices=history),
            person_pf=NS(
                aggregate_heads_to_model_age_cells=lambda *a, **k: np.array([1.]),
                PersonPFState=lambda **kwargs: NS(**kwargs)),
            ConditionalHistoryEvaluation=lambda **kwargs: NS(**kwargs),
            check_smoke_gates=lambda *a, **k: None)
        old = NS(parameters=NS(age_start=18, da=4, J=1), b_grid=np.array([0.]),
                 supply_rule=None,
                 historical_conditioning=Conditioning(next_age_targets={}))
        terminal = NS(parameters=NS(psi_child=.2),
                      policy=NS(price=[2.], V=np.array([9.])),
                      fixed_point=NS(g_pre=g, persons=people))
        transfers = np.arange(n, dtype=float) + .25
        with patch.object(r, "_runtime", return_value=(None, joined, None, None, None)):
            result = r.evaluate_forecast(
                inherited=r.InheritedState(2007, NS(g_pre=g)), old_state=old,
                demographics=NS(initial_person_state=people),
                prices=np.arange(n, dtype=float) + 2., pensions=np.ones(n),
                transfers=transfers, psi=.2, terminal=terminal,
                demographic_evaluator=hook)
        np.testing.assert_array_equal(calls["history"]["transfer_path"], transfers[:h])
        np.testing.assert_array_equal(calls["backward"]["transfer_path"], transfers[h:])
        np.testing.assert_array_equal(calls["tail"]["transfer_path"], transfers[h:])
        np.testing.assert_array_equal(calls["tail"]["psi_path"], np.full(n - h, .2))
        self.assertEqual(len(result.rows), n)


if __name__ == "__main__":
    unittest.main()
