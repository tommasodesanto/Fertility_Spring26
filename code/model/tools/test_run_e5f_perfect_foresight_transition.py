"""Tiny-grid tests for the isolated perfect-foresight transition driver."""
from __future__ import annotations

import sys
import tempfile
import types
import unittest
import weakref
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

# The local lightweight Python runtime used by this regression test does not
# bundle plotting libraries.  The tested solver path never creates figures, so
# provide only the import surface required by the production driver modules.
if "matplotlib" not in sys.modules:
    matplotlib_stub = types.ModuleType("matplotlib")
    matplotlib_stub.use = lambda *_args, **_kwargs: None
    pyplot_stub = types.ModuleType("matplotlib.pyplot")
    matplotlib_stub.pyplot = pyplot_stub
    sys.modules["matplotlib"] = matplotlib_stub
    sys.modules["matplotlib.pyplot"] = pyplot_stub

from intergen_eqscale_seq_optimized import solver as model

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_perfect_foresight_transition as driver
import run_e5f_perfect_foresight_rebated_property_tax as policy_driver


def _tiny_overrides() -> dict[str, object]:
    return {
        "J": 6,
        "J_R": 5,
        "A_f_start": 1,
        "A_f_end": 4,
        "Nb": 20,
        "b_min": 0.0,
        "b_core_lo": 0.0,
        "b_core_hi": 3.0,
        "b_mid_hi": 6.0,
        "b_max": 10.0,
        "n_house": 2,
        "H_own": np.array([2.0, 4.0]),
        "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]),
        "solve_mode": "pe",
        "p_fixed": np.array([1.0]),
        "w_fixed": np.array([1.0]),
        "entry_shares_fixed": np.array([1.0]),
        "use_income_types": True,
        "income_type_transition": "persistent",
        "z_grid": np.array([0.8, 1.2]),
        "z_weights": np.array([0.5, 0.5]),
        "Pi_z": np.array([[0.9, 0.1], [0.1, 0.9]]),
        "entry_wealth_mode": "scalar",
        "b_entry_fixed": 0.5,
        "entry_wealth_spread_nodes": 1,
        "c_bar_0": 0.05,
        "c_bar_n": 0.02,
        "h_bar_0": 0.25,
        "h_bar_jump": 0.10,
        "h_bar_n": 0.05,
        "lambda_d": 0.0,
        "use_full_kernel": True,
        "use_tenure_kernel": True,
        "use_loc_kernel": True,
        "tenure_choice_kappa": 0.0,
        "sequential_births": True,
        "n_parity": 4,
        "fertility_units": "literal_topcode",
        "tfr_top_bin_weight": 3.4,
        "child_state_mode": "independent_count",
        "A_m": 18.0,
        "use_stochastic_aging": True,
        "entrant_conversion_factor": 0.5,
    }


class TinyPerfectForesightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        calendar.model = model
        calendar.apply_fertility = transition.apply_sequential_fertility
        calendar.advance_calendar_distribution = (
            transition.advance_sequential_calendar_distribution
        )
        solution, parameters, price = model.run_model_cp_dt(
            _tiny_overrides(), verbose=False
        )
        parameters._fert2_probs = np.asarray(solution.fert2_probs, dtype=float).copy()
        b_grid = np.asarray(solution.b_grid, dtype=float)
        shared = model.precompute_shared(parameters, b_grid)
        policy = calendar.policy_from_solution(
            solution, price, parameters, b_grid, shared
        )
        stationary_g_pre, diagnostics = calendar.reconstruct_stationary_pre_fertility(
            solution, policy, parameters, b_grid, shared
        )
        supply_rule, _ = calendar.normalize_date0_housing_supply(
            stationary_g_pre,
            policy,
            parameters,
            b_grid,
            shared,
            "static-elastic",
        )
        evaluation = calendar.evaluate_period(
            price,
            stationary_g_pre,
            parameters,
            b_grid,
            shared,
            calendar.SolveCounter(),
            supply_rule=supply_rule,
            supplied_policy=policy,
        )
        accounting = transition.calendar_topcode_birth_accounting(
            evaluation.g_pre,
            evaluation.g_post_fertility,
            float(evaluation.births),
            parameters,
        )
        entry_flow = float(np.sum(stationary_g_pre[:, :, :, 0, :, :, :]))
        conversion = entry_flow / float(
            accounting["topcode_adjusted_birth_children"]
        )
        cls.solution = solution
        cls.parameters = parameters
        cls.price = float(np.asarray(price).reshape(-1)[0])
        cls.b_grid = b_grid
        cls.policy = policy
        cls.stationary_g_pre = stationary_g_pre
        cls.supply_rule = supply_rule
        cls.entry_flow = entry_flow
        cls.conversion = conversion
        cls.initial_state = driver.PFInitialState(
            g_pre=stationary_g_pre.copy(),
            scheduled_entries=[entry_flow] * 4,
            scheduled_raw_entries=[entry_flow] * 4,
        )
        assert diagnostics["stationary_post_fertility_nesting_l1"] < 1e-9

    def test_stationary_prices_nest_zero_shock_path(self) -> None:
        periods = 3
        path = driver.evaluate_path_at_prices(
            prices=np.full(periods, self.price),
            psi_path=np.full(periods, float(self.parameters.psi_child)),
            terminal_price=self.price,
            terminal_V=self.policy.V,
            base_parameters=self.parameters,
            b_grid=self.b_grid,
            initial_state=self.initial_state,
            supply_rule=self.supply_rule,
            birth_to_entry_conversion=self.conversion,
        )

        expected_rent = float(self.parameters.user_cost_rate) * self.price
        np.testing.assert_allclose(path.rents, expected_rent, atol=1e-14, rtol=0.0)
        self.assertLess(path.maximum_market_residual, 1e-9)
        self.assertLess(path.maximum_policy_reproduction_error, 1e-12)
        self.assertLess(path.maximum_mass_accounting_error, 1e-10)
        self.assertLess(
            float(
                np.sum(
                    np.abs(
                        path.terminal_state.g_pre - self.stationary_g_pre
                    )
                )
            ),
            1e-8,
        )

    def test_equal_transfer_enters_backward_problem_and_fiscal_ledger(self) -> None:
        baseline_path = driver.evaluate_path_at_prices(
            prices=np.array([self.price]),
            psi_path=np.array([float(self.parameters.psi_child)]),
            terminal_price=self.price,
            terminal_V=self.policy.V,
            base_parameters=self.parameters,
            b_grid=self.b_grid,
            initial_state=self.initial_state,
            supply_rule=self.supply_rule,
            birth_to_entry_conversion=self.conversion,
        )
        transfer = 0.10
        funded_path = driver.evaluate_path_at_prices(
            prices=np.array([self.price]),
            psi_path=np.array([float(self.parameters.psi_child)]),
            terminal_price=self.price,
            terminal_V=self.policy.V,
            base_parameters=self.parameters,
            b_grid=self.b_grid,
            initial_state=self.initial_state,
            supply_rule=self.supply_rule,
            birth_to_entry_conversion=self.conversion,
            transfer_path=np.array([transfer]),
        )
        row = funded_path.rows[0]
        self.assertAlmostEqual(row["equal_transfer_period_units"], transfer)
        self.assertAlmostEqual(
            row["equal_transfer_outlays"],
            transfer * row["adult_population"],
        )
        self.assertAlmostEqual(
            row["government_budget_residual"],
            row["property_tax_revenue"] - row["equal_transfer_outlays"],
        )
        self.assertGreater(
            float(np.max(np.abs(funded_path.values[0] - baseline_path.values[0]))),
            1e-8,
        )

    def test_equal_transfer_path_rejects_negative_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            driver.evaluate_path_at_prices(
                prices=np.array([self.price]),
                psi_path=np.array([float(self.parameters.psi_child)]),
                terminal_price=self.price,
                terminal_V=self.policy.V,
                base_parameters=self.parameters,
                b_grid=self.b_grid,
                initial_state=self.initial_state,
                supply_rule=self.supply_rule,
                birth_to_entry_conversion=self.conversion,
                transfer_path=np.array([-0.01]),
            )

    def test_endpoint_solver_rejects_looser_root_tolerance(self) -> None:
        with self.assertRaises(ValueError):
            driver.continuation.solve_closed_stationary_endpoint(
                chain=None,
                base_overrides={},
                old_price=1.0,
                new_psi_child=0.0,
                supply_rule=self.supply_rule,
                price_min_ratio=0.1,
                price_max_ratio=2.0,
                grid_points=5,
                outdir=Path("unused_endpoint_test_output"),
                root_residual_tolerance=1.0,
            )

    def test_collector_contract_loader_is_complete_and_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.json"
            source.write_text("source\n", encoding="utf-8")
            source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            transition_path = root / "best_transition_path.csv"
            transition_path.write_text(
                "period,relative_market_residual,mass_accounting_residual,"
                "population_target_gap\n"
                + "".join(f"{period},0.00003,0,0\n" for period in range(5)),
                encoding="utf-8",
            )
            shared = {
                "code_fingerprints": {"bundle_sha256": "live-code"},
                "external_closure_contract": {
                    "housing_supply_elasticity": 0.63,
                    "housing_supply_elasticity_status": (
                        "externally_fixed_profile_not_estimated"
                    ),
                    "tenure_choice_kappa": 0.005,
                    "tenure_choice_kappa_status": (
                        "externally_fixed_profile_not_estimated"
                    ),
                },
                "model_profile": {"name": "e5f-income-entry"},
                "renewal_accounting_contract": {"test": "mocked"},
                "renewal_accounting_old_state": {"test": "mocked"},
                "outside_origin_entry_share": 0.169,
                "population_validation_status": "matched",
            }
            best = {
                **shared,
                "valid": True,
                "policy_case": "none",
                "post_2023_periods": 0,
                "target_count": 12,
                "target_fingerprint": "live-target",
                "source_sha256": source_hash,
                "theta": {"tenure_choice_kappa": 0.005},
                "old_psi_child": 0.2,
                "new_psi_child": -0.1,
                "max_market_residual": 0.00003,
            }
            report = root / "summary.json"
            report.write_text(
                json.dumps(
                    {
                        **shared,
                        "status": "complete",
                        "expected_tasks": 23,
                        "completed_tasks": 23,
                        "valid_tasks": 23,
                        "failed_or_missing_tasks": [],
                        "best_candidate": best,
                    }
                ),
                encoding="utf-8",
            )
            with (
                mock.patch.object(
                    driver.continuation,
                    "e5_target_system",
                    return_value=SimpleNamespace(fingerprint="live-target"),
                ),
                mock.patch.object(
                    driver.transition,
                    "configure_sequential_model",
                    return_value=(object(), object()),
                ),
                mock.patch.object(
                    driver.continuation.calibration,
                    "code_fingerprint_contract",
                    return_value={"bundle_sha256": "live-code"},
                ),
                mock.patch.object(
                    driver.continuation, "validate_renewal_contract"
                ),
            ):
                contracts, provenance = driver.load_diagnostic_contracts(
                    report, transition_path, source
                )
                self.assertEqual(
                    contracts["report_best"]["theta"]["tenure_choice_kappa"],
                    0.005,
                )
                self.assertEqual(
                    provenance["selected_report_sha256"],
                    hashlib.sha256(report.read_bytes()).hexdigest(),
                )
                broken = json.loads(report.read_text(encoding="utf-8"))
                broken["valid_tasks"] = 22
                report.write_text(json.dumps(broken), encoding="utf-8")
                with self.assertRaisesRegex(RuntimeError, "incomplete or invalid"):
                    driver.load_diagnostic_contracts(
                        report, transition_path, source
                    )

    def test_recorded_price_seed_extends_at_terminal_price(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "seed.csv"
            path.write_text(
                "period,asset_price\n0,1.20\n1,1.10\n",
                encoding="utf-8",
            )
            prices, contract = driver.load_initial_price_path(
                path,
                horizon=4,
                terminal_price=0.90,
            )
        np.testing.assert_allclose(
            prices,
            np.array([1.20, 1.10, 0.90, 0.90]),
            rtol=0.0,
            atol=0.0,
        )
        self.assertEqual(contract["source_rows"], 2)
        self.assertEqual(contract["terminal_extension_periods"], 2)
        self.assertEqual(len(contract["source_sha256"]), 64)

    def test_recorded_price_seed_rejects_noncontiguous_periods(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad_seed.csv"
            path.write_text(
                "period,asset_price\n0,1.20\n2,1.10\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "contiguous"):
                driver.load_initial_price_path(
                    path,
                    horizon=4,
                    terminal_price=0.90,
                )

    def test_policy_price_seed_projection_preserves_positive_rents(self) -> None:
        terminal = SimpleNamespace(
            asset_price=1.0,
            parameters=SimpleNamespace(R_gross=1.10, delta=0.05, tau_H=0.05),
        )
        projected, contract = policy_driver.project_price_seed_to_positive_rents(
            np.array([0.50, 0.60]), terminal=terminal
        )
        rents = driver.rents_from_asset_prices(
            projected, terminal.asset_price, terminal.parameters
        )
        self.assertTrue(np.all(rents > 0.0))
        self.assertEqual(contract["adjusted_period_count"], 2)
        self.assertFalse(contract["equilibrium_contract_changed"])

    def test_policy_price_path_projection_is_date_specific(self) -> None:
        terminal = SimpleNamespace(
            asset_price=1.0,
            parameters=SimpleNamespace(R_gross=1.10, delta=0.05, tau_H=0.05),
        )
        original = np.array([0.80, 0.60])
        projected, contract = policy_driver.project_price_path_to_positive_rents(
            original, terminal=terminal, minimum_rent_share=1.0e-6
        )
        self.assertAlmostEqual(projected[0], original[0])
        self.assertGreater(projected[1], original[1])
        self.assertEqual(contract["adjusted_periods"], [1])
        rents = driver.rents_from_asset_prices(
            projected, terminal.asset_price, terminal.parameters
        )
        self.assertTrue(np.all(rents > 0.0))

    def test_anticipated_preference_path_runs_backward_and_forward(self) -> None:
        periods = 3
        old_psi = float(self.parameters.psi_child)
        psi_path = driver.preference_path(old_psi - 0.10, old_psi, 0.5, periods)
        solution = driver.solve_perfect_foresight(
            initial_prices=np.full(periods, self.price),
            psi_path=psi_path,
            terminal_price=self.price,
            terminal_V=self.policy.V,
            base_parameters=self.parameters,
            b_grid=self.b_grid,
            initial_state=self.initial_state,
            supply_rule=self.supply_rule,
            birth_to_entry_conversion=self.conversion,
            market_tolerance=2e-2,
            max_iterations=8,
            damping=0.40,
            maximum_log_price_step=0.15,
            terminal_reference_state=self.initial_state,
            terminal_reference_entry_flow=self.entry_flow,
            terminal_reference_psi=old_psi,
        )
        self.assertTrue(np.all(np.isfinite(solution.prices)))
        self.assertTrue(np.all(np.isfinite(solution.rents)))
        self.assertGreater(float(np.min(solution.rents)), 0.0)
        self.assertLess(solution.maximum_policy_reproduction_error, 1e-12)
        self.assertLess(solution.maximum_mass_accounting_error, 1e-10)
        self.assertEqual(
            solution.terminal_convergence["status"], "not_converged"
        )
        self.assertFalse(solution.terminal_convergence["all_checks_pass"])
        self.assertLessEqual(
            solution.maximum_market_residual,
            float(solution.iteration_history[0]["maximum_market_residual"]),
        )

    def test_zero_shock_terminal_gate_passes(self) -> None:
        periods = 3
        solution = driver.solve_perfect_foresight(
            initial_prices=np.full(periods, self.price),
            psi_path=np.full(periods, float(self.parameters.psi_child)),
            terminal_price=self.price,
            terminal_V=self.policy.V,
            base_parameters=self.parameters,
            b_grid=self.b_grid,
            initial_state=self.initial_state,
            supply_rule=self.supply_rule,
            birth_to_entry_conversion=self.conversion,
            market_tolerance=1e-8,
            max_iterations=1,
            damping=0.40,
            maximum_log_price_step=0.15,
            terminal_reference_state=self.initial_state,
            terminal_reference_entry_flow=self.entry_flow,
            terminal_reference_psi=float(self.parameters.psi_child),
        )
        self.assertTrue(solution.converged)
        self.assertEqual(solution.terminal_convergence["status"], "passed")
        self.assertTrue(solution.terminal_convergence["all_checks_pass"])
        self.assertTrue(all(solution.terminal_convergence["checks"].values()))

    def test_nonunit_population_terminal_state_nests(self) -> None:
        scale = 0.73
        scaled_state = driver.PFInitialState(
            g_pre=scale * self.stationary_g_pre,
            scheduled_entries=[scale * self.entry_flow] * 4,
            scheduled_raw_entries=[scale * self.entry_flow] * 4,
        )
        scaled_supply = calendar.HousingSupplyRule(
            mode=self.supply_rule.mode,
            initial_price=self.supply_rule.initial_price,
            initial_stock=scale * self.supply_rule.initial_stock,
            elasticity=self.supply_rule.elasticity,
        )
        solution = driver.solve_perfect_foresight(
            initial_prices=np.full(3, self.price),
            psi_path=np.full(3, float(self.parameters.psi_child)),
            terminal_price=self.price,
            terminal_V=self.policy.V,
            base_parameters=self.parameters,
            b_grid=self.b_grid,
            initial_state=scaled_state,
            supply_rule=scaled_supply,
            birth_to_entry_conversion=self.conversion,
            market_tolerance=1e-8,
            max_iterations=1,
            damping=0.40,
            maximum_log_price_step=0.15,
            terminal_reference_state=scaled_state,
            terminal_reference_entry_flow=scale * self.entry_flow,
            terminal_reference_psi=float(self.parameters.psi_child),
        )
        self.assertTrue(solution.converged)
        self.assertEqual(solution.terminal_convergence["status"], "passed")
        self.assertAlmostEqual(
            solution.terminal_convergence["stationary_population"],
            scale,
            places=12,
        )

    def test_joint_policy_iteration_clears_equal_transfer_fixed_point(self) -> None:
        state = driver.PFInitialState(
            g_pre=np.ones((1, 1, 1, 1, 1, 1, 1)),
            scheduled_entries=[1.0] * 4,
            scheduled_raw_entries=[1.0] * 4,
        )
        parameters = SimpleNamespace(
            R_gross=1.1,
            delta=0.05,
            tau_H=0.05,
            user_cost_rate=0.20,
        )
        case = policy_driver.impact.CASES["rebated-tax1-baseline"]
        terminal = policy_driver.PolicyTerminalSteadyState(
            case=case,
            psi_child=0.0,
            asset_price=1.1,
            renter_price=0.22,
            equal_transfer=0.2,
            population_scale=1.0,
            entry_flow=1.0,
            raw_queue_flow=1.0,
            policy=SimpleNamespace(V=np.zeros(1)),
            parameters=parameters,
            state=state,
            endpoint={},
            reconstruction={},
            fixed_point_gates={"status": "passed"},
        )
        supply_rule = calendar.HousingSupplyRule(
            mode="static-elastic",
            initial_price=1.0,
            initial_stock=1.0,
            elasticity=1.0,
        )

        evaluation_refs: list[weakref.ReferenceType[driver.PathEvaluation]] = []

        def fake_evaluation(**kwargs):
            if evaluation_refs:
                self.assertIsNone(
                    evaluation_refs[-1](),
                    "The prior path evaluation must be released before the next solve.",
                )
            prices = np.asarray(kwargs["prices"], dtype=float)
            transfers = np.asarray(kwargs["transfer_path"], dtype=float)
            rows = []
            for period, (price, transfer) in enumerate(zip(prices, transfers)):
                supply = float(price)
                rows.append(
                    {
                        "period": period,
                        "calendar_year": 2023 + 4 * period,
                        "housing_demand": 1.1,
                        "housing_supply": supply,
                        "relative_market_residual": abs(1.1 - supply) / supply,
                        "government_budget_residual": 0.2 - float(transfer),
                        "implied_equal_transfer": 0.2,
                        "equal_transfer_period_units": float(transfer),
                        "birth_children_topcode_adjusted": 1.0,
                        "adult_population": 1.0,
                        "entry_flow_E": 1.0,
                        "queue_B_over_current_E": 1.0,
                        "owner_rate": 0.5,
                        "asset_price": float(price),
                        "renter_price": 0.22,
                    }
                )
            result = driver.PathEvaluation(
                prices=prices,
                rents=np.full_like(prices, 0.22),
                values=[np.zeros(1)] * (len(prices) + 1),
                rows=rows,
                terminal_state=state,
                maximum_market_residual=max(
                    float(row["relative_market_residual"]) for row in rows
                ),
                maximum_policy_reproduction_error=0.0,
                maximum_mass_accounting_error=0.0,
                maximum_feasibility_projection_mass=0.0,
                bellman_solves=2 * len(prices),
                elapsed_seconds=0.0,
            )
            evaluation_refs.append(weakref.ref(result))
            return result

        with mock.patch.object(
            policy_driver.pf,
            "evaluate_path_at_prices",
            side_effect=fake_evaluation,
        ):
            solution = policy_driver.solve_funded_perfect_foresight(
                case=case,
                initial_prices=np.full(2, 1.1),
                initial_transfers=np.zeros(2),
                psi_path=np.zeros(2),
                terminal=terminal,
                b_grid=np.array([0.0]),
                initial_state=state,
                supply_rule=supply_rule,
                market_tolerance=1e-8,
                fiscal_tolerance=1e-8,
                max_iterations=3,
                price_damping=1.0,
                transfer_damping=1.0,
                maximum_log_price_step=0.2,
                maximum_transfer_step=0.5,
            )
        self.assertTrue(solution.converged)
        np.testing.assert_allclose(solution.transfers, 0.2)
        self.assertLessEqual(solution.maximum_fiscal_residual, 1e-8)

    def test_stationary_equal_transfer_root_does_not_need_solution_cache(self) -> None:
        case = policy_driver.impact.CASES["rebated-tax1-baseline"]

        def fake_stationary_evaluation(**kwargs):
            transfer = float(kwargs["transfer"])
            return SimpleNamespace(
                ledger={"implied_equal_transfer": 0.2},
                parameters=SimpleNamespace(
                    property_tax_lump_sum_transfer=transfer
                ),
            )

        with mock.patch.object(
            policy_driver,
            "_stationary_evaluation",
            side_effect=fake_stationary_evaluation,
        ) as mocked:
            result, evaluations = policy_driver.solve_stationary_equal_transfer(
                prepared=SimpleNamespace(),
                case=case,
                asset_price=1.0,
                psi_child=0.0,
                tolerance=1e-12,
            )
        self.assertAlmostEqual(
            result.parameters.property_tax_lump_sum_transfer, 0.2
        )
        self.assertEqual(evaluations, 3)
        self.assertEqual(mocked.call_count, evaluations)


if __name__ == "__main__":
    unittest.main()
