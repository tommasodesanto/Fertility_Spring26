from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

import run_e5f_perfect_foresight_person_demography_policy as policy


def state(*, persons: np.ndarray, heads: np.ndarray, g_pre: np.ndarray, year: int):
    return SimpleNamespace(
        persons=SimpleNamespace(persons=persons, heads=heads, year=year),
        g_pre=g_pre,
    )


class TerminalConvergenceDiagnosticsTest(unittest.TestCase):
    def test_terminal_supply_elasticity_is_explicit_and_positive(self) -> None:
        self.assertEqual(
            policy.terminal_supply_elasticity(
                {"supply_contract": {"supply_elasticity": 0.63}}
            ),
            0.63,
        )
        invalid_contracts = (
            {},
            {"supply_contract": {}},
            {"supply_contract": {"supply_elasticity": 0.0}},
        )
        for invalid in invalid_contracts:
            with self.assertRaises(RuntimeError):
                policy.terminal_supply_elasticity(invalid)

    def test_normalized_l1_is_invariant_to_total_mass(self) -> None:
        reference = np.array([1.0, 2.0, 3.0])
        self.assertAlmostEqual(policy.normalized_l1(7.0 * reference, reference), 0.0)

    def test_all_declared_terminal_checks_are_jointly_enforced(self) -> None:
        reference = state(
            persons=np.array([[2.0, 3.0], [4.0, 1.0]]),
            heads=np.array([[0.5, 1.5], [2.0, 0.5]]),
            g_pre=np.array([[[1.0, 2.0], [3.0, 4.0]]]),
            year=2100,
        )
        terminal = SimpleNamespace(
            state=reference,
            asset_price=2.0,
            renter_price=0.5,
            equal_transfer=0.25,
            psi_child=-0.07,
        )
        passing = SimpleNamespace(
            terminal_state=state(
                persons=reference.persons.persons * 1.005,
                heads=reference.persons.heads * 1.005,
                g_pre=reference.g_pre * 1.005,
                year=2535,
            ),
            prices=np.array([2.0 * 1.005]),
            rents=np.array([0.5 * 1.005]),
            rows=[{"equal_transfer_period_units": 0.25 * 1.005}],
        )
        result = policy.terminal_convergence_diagnostics(
            passing, terminal=terminal, psi_path=[-0.0695]
        )
        self.assertTrue(result["all_checks_pass"])
        self.assertTrue(all(result["checks"].values()))

        failing = SimpleNamespace(
            terminal_state=passing.terminal_state,
            prices=np.array([2.0 * 1.011]),
            rents=passing.rents,
            rows=passing.rows,
        )
        result = policy.terminal_convergence_diagnostics(
            failing, terminal=terminal, psi_path=[-0.0695]
        )
        self.assertFalse(result["all_checks_pass"])
        self.assertFalse(result["checks"]["asset_price_relative_gap"])
        self.assertTrue(
            all(
                passed
                for name, passed in result["checks"].items()
                if name != "asset_price_relative_gap"
            )
        )


class AdaptiveDampingTest(unittest.TestCase):
    def test_soft_time_stop_returns_a_restartable_best_path(self) -> None:
        terminal_state = state(
            persons=np.array([[1.0]]),
            heads=np.array([[1.0]]),
            g_pre=np.array([1.0]),
            year=2031,
        )
        evaluation = SimpleNamespace(
            prices=np.array([1.0]),
            rents=np.array([0.1]),
            rows=[
                {
                    "housing_demand": 2.0,
                    "housing_supply": 1.0,
                    "government_budget_residual": 0.1,
                    "property_tax_revenue": 0.2,
                    "household_heads": 1.0,
                    "equal_transfer_period_units": 0.1,
                    "resident_persons": 1.0,
                }
            ],
            terminal_state=terminal_state,
            maximum_market_residual=0.5,
            maximum_policy_reproduction_error=0.0,
            maximum_person_identity_error=0.0,
            maximum_head_identity_error=0.0,
            maximum_household_person_head_gap=0.0,
            maximum_age_head_gap=0.0,
            maximum_feasibility_projection_mass=0.0,
            bellman_solves=1,
        )
        terminal = SimpleNamespace(
            case=SimpleNamespace(
                name="rebated-tax1-baseline",
                label="1%",
                annual_tax_rate=0.01,
            ),
            asset_price=1.0,
            renter_price=0.1,
            equal_transfer=0.1,
            psi_child=0.0,
            state=terminal_state,
            parameters=SimpleNamespace(),
            policy=SimpleNamespace(V=None),
        )
        supply_rule = SimpleNamespace(
            mode="static-elastic",
            elasticity=1.0,
            initial_price=1.0,
            initial_stock=1.0,
        )
        with mock.patch.object(
            policy.person_pf,
            "evaluate_path_at_prices_person_demography",
            return_value=evaluation,
        ):
            result = policy.solve_person_funded_path(
                initial_prices=[1.0],
                initial_transfers=[0.1],
                psi_path=[0.0],
                terminal=terminal,
                b_grid=np.array([0.0]),
                initial_state=terminal_state,
                primitives=SimpleNamespace(scale_model_units_per_person=1.0),
                supply_rule=supply_rule,
                maximum_iterations=5,
                price_damping=0.25,
                transfer_damping=0.50,
                maximum_log_price_step=0.12,
                maximum_transfer_step=0.08,
                stop_after_iteration=lambda _record: True,
            )
        self.assertEqual(result.status, "soft_time_budget_reached")
        self.assertFalse(result.path_converged)
        self.assertEqual(result.iterations, 1)

    def test_deliberately_low_declared_damping_is_not_raised_after_worsening(self) -> None:
        price, transfer = policy.adapt_path_damping(
            current_price_damping=0.02,
            current_transfer_damping=0.04,
            declared_price_damping=0.02,
            declared_transfer_damping=0.04,
            score=1.03,
            previous_score=1.0,
        )
        self.assertEqual(price, 0.02)
        self.assertEqual(transfer, 0.04)

    def test_standard_declared_damping_retains_existing_floors(self) -> None:
        price, transfer = policy.adapt_path_damping(
            current_price_damping=0.0625,
            current_transfer_damping=0.125,
            declared_price_damping=0.25,
            declared_transfer_damping=0.50,
            score=1.03,
            previous_score=1.0,
        )
        self.assertEqual(price, 0.04)
        self.assertEqual(transfer, 0.08)


if __name__ == "__main__":
    unittest.main()
