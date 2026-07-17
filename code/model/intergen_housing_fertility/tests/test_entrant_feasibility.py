from __future__ import annotations

import hashlib
import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intergen_housing_fertility.kernels import (
    full_owner_block_kernel,
    full_renter_block_kernel,
)
from intergen_housing_fertility.local_panel import (
    income_process_overrides,
    run_local_panel_case,
)
from intergen_housing_fertility.parameters import (
    apply_overrides,
    setup_parameters,
    unsecured_debt_floor,
)
from intergen_housing_fertility.solver import (
    DEAD_VALUE_CUTOFF,
    InfeasibleThetaError,
    build_forward_tenure_transition_maps,
    eval_owner,
    get_phi_state_matrix,
    owner_borrowing_floor,
    precompute_shared,
    renter_borrowing_floor,
    run_model_cp_dt,
    solve_bellman_full,
    solve_bellman_full_markov_income,
)
from intergen_housing_fertility.utils import make_grid


def _array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def _python_kernel(function):
    """Use the kernel body directly; compiled and fallback paths share it."""

    return getattr(function, "py_func", function)


def _mapped_values(grid: np.ndarray, idx: np.ndarray, wt: np.ndarray) -> np.ndarray:
    return (1.0 - wt) * grid[idx] + wt * grid[idx + 1]


def _inert_overrides() -> dict[str, object]:
    # This deliberately keeps the entire wealth grid non-negative. With
    # lambda_d=0 the debt repair must then be exactly dormant, byte for byte.
    return {
        "J": 17,
        "J_R": 12,
        "Nb": 24,
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
        "use_income_types": False,
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
    }


def _indebted_solution_overrides() -> dict[str, object]:
    overrides = _inert_overrides()
    overrides.update(
        {
            "Nb": 28,
            "b_min": -3.0,
            "b_core_lo": -2.0,
            "b_entry_fixed": -1.0,
            "entry_wealth_spread_nodes": 1,
        }
    )
    return overrides


def _all_dead_overrides() -> dict[str, object]:
    return {
        "J": 3,
        "J_R": 3,
        "Nb": 16,
        "b_min": -3.0,
        "b_core_lo": -2.0,
        "b_core_hi": 2.0,
        "b_mid_hi": 4.0,
        "b_max": 6.0,
        "n_house": 1,
        "H_own": np.array([2.0]),
        "H0": np.array([4.0]),
        "eta_supply": np.array([1.0]),
        "solve_mode": "pe",
        "p_fixed": np.array([1.0]),
        "w_fixed": np.array([1.0]),
        "entry_shares_fixed": np.array([1.0]),
        "use_income_types": False,
        "entry_wealth_mode": "scalar",
        "b_entry_fixed": 0.0,
        "entry_wealth_spread_nodes": 1,
        "c_bar_0": 100.0,
        "c_bar_n": 0.0,
        "h_bar_0": 0.25,
        "h_bar_jump": 0.0,
        "h_bar_n": 0.0,
        "lambda_d": 0.0,
        "tenure_choice_kappa": 0.01,
        "use_full_kernel": True,
        "use_tenure_kernel": True,
        "use_loc_kernel": True,
    }


class DebtCapFormulaTests(unittest.TestCase):
    def test_debt_cap_construction_uses_income_and_exact_age_taper(self) -> None:
        params = apply_overrides(
            setup_parameters(),
            {"lambda_d": 0.4, "use_income_types": False},
        )
        ages = params.age_start + np.arange(params.J) * params.da
        expected_taper = np.where(
            ages <= 42.0,
            1.0,
            np.where(ages >= 62.0, 0.0, (62.0 - ages) / 20.0),
        )
        np.testing.assert_array_equal(
            params.debt_taper_weights,
            np.append(expected_taper, 0.0),
        )
        np.testing.assert_allclose(
            params.debt_caps[:-1],
            0.4 * params.mean_labor_income_by_age * expected_taper,
            atol=0.0,
            rtol=0.0,
        )
        self.assertEqual(float(params.debt_caps[-1]), 0.0)
        self.assertTrue(np.all(params.debt_caps >= 0.0))
        self.assertTrue(np.all(params.debt_caps[:-1][ages >= 62.0] == 0.0))

    def test_debt_cap_income_mean_accepts_all_markov_transition_aliases(self) -> None:
        for alias in ("markov", "stochastic", "persistent"):
            with self.subTest(alias=alias):
                params = apply_overrides(
                    setup_parameters(),
                    {
                        "lambda_d": 0.4,
                        "use_income_types": True,
                        "income_type_transition": alias,
                        "z_grid": np.array([0.5, 1.5]),
                        "z_weights": np.array([0.25, 0.75]),
                        "Pi_z": np.array([[0.8, 0.2], [0.1, 0.9]]),
                    },
                )
                expected_z_mean = 1.25
                expected_income = expected_z_mean * (
                    params.entry_shares @ params.income[:, 0]
                )
                self.assertAlmostEqual(
                    float(params.mean_labor_income_by_age[0]),
                    float(expected_income),
                    places=12,
                )

    def test_unsecured_floor_rolls_legacy_debt_and_opens_standing_line(self) -> None:
        unsecured = np.array([-2.0, -0.5, 0.0, 1.0])
        expected = np.array([-1.0, -0.4, -0.4, -0.4])
        np.testing.assert_allclose(
            unsecured_debt_floor(unsecured, s_next=0.5, D_next=0.4),
            expected,
            atol=0.0,
            rtol=0.0,
        )

    def test_owner_floor_does_not_count_the_mortgage_twice(self) -> None:
        params = SimpleNamespace(
            J=2,
            debt_taper_weights=np.array([1.0, 1.0, 0.0]),
            debt_caps=np.zeros(3),
        )
        collateral_floor = -8.0
        self.assertEqual(
            float(owner_borrowing_floor(params, -8.0, collateral_floor, 0)),
            -8.0,
        )
        self.assertEqual(
            float(owner_borrowing_floor(params, -9.0, collateral_floor, 0)),
            -9.0,
        )
        self.assertEqual(
            float(renter_borrowing_floor(params, -1.5, 0)),
            -1.5,
        )

    def test_owner_ltv_taper_reduces_next_period_collateral_capacity(self) -> None:
        params = SimpleNamespace(
            J=2,
            debt_taper_weights=np.array([1.0, 0.0, 0.0]),
            debt_caps=np.zeros(3),
            owner_ltv_multipliers=np.array([1.0, 0.5, 0.0]),
        )
        self.assertEqual(
            float(owner_borrowing_floor(params, -8.0, -8.0, 0)),
            -4.0,
        )
        self.assertEqual(
            float(owner_borrowing_floor(params, -4.0, -8.0, 1)),
            0.0,
        )


class LowLevelKernelTests(unittest.TestCase):
    def test_indebted_renter_is_feasible_and_budget_identity_holds(self) -> None:
        grid = np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0])
        gross_return = 1.08
        income = 3.0
        resources = gross_return * grid + income
        continuation = np.zeros((grid.size, 1))
        previous = np.minimum(grid, 0.0)[:, None]
        kernel = _python_kernel(full_renter_block_kernel)
        values, bp, consumption, housing = kernel(
            resources,
            continuation,
            previous,
            1,
            grid,
            np.array([0.2]),
            np.array([0.5]),
            np.array([0.0]),
            0.3,
            10.0,
            0.04,
            0.2,
            0.5,
            0.7,
            -1.0,
            0.9,
            1.0,
            0.0,
            (3.0 - math.sqrt(5.0)) / 2.0,
            (math.sqrt(5.0) - 1.0) / 2.0,
            1e-6,
        )
        indebted = int(np.flatnonzero(grid == -1.0)[0])
        self.assertGreater(float(values[indebted, 0]), DEAD_VALUE_CUTOFF)
        self.assertLess(float(bp[indebted, 0]), 0.0)  # warm start must not restore b'>=0
        self.assertGreaterEqual(float(bp[indebted, 0]), -1.0)
        slack = resources[indebted] - (
            consumption[indebted, 0] + 0.3 * housing[indebted, 0] + bp[indebted, 0]
        )
        self.assertAlmostEqual(float(slack), 0.0, places=12)

    def test_owner_kernel_floor_sweep_and_budget_identity(self) -> None:
        grid = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0])
        collateral_floor = -2.0
        expected_floor = collateral_floor + np.minimum(
            np.minimum(grid - collateral_floor, 0.0),
            0.0,
        )
        resources = 1.08 * grid + 5.0
        kernel = _python_kernel(full_owner_block_kernel)
        values, bp, consumption = kernel(
            resources,
            np.zeros((grid.size, 1)),
            expected_floor[:, None],
            1,
            grid,
            np.array([0.2]),
            np.array([0.5]),
            np.array([0.0]),
            np.array([collateral_floor]),
            0.4,
            2.0,
            1.0,
            1.0,
            0.04,
            0.7,
            -1.0,
            0.9,
            1.0,
            0.0,
            (3.0 - math.sqrt(5.0)) / 2.0,
            (math.sqrt(5.0) - 1.0) / 2.0,
            1e-6,
        )
        self.assertTrue(np.all(values > DEAD_VALUE_CUTOFF))
        self.assertTrue(np.all(bp[:, 0] >= expected_floor - 1e-12))
        fully_leveraged = int(np.flatnonzero(grid == collateral_floor)[0])
        self.assertEqual(float(bp[fully_leveraged, 0]), collateral_floor)
        np.testing.assert_allclose(
            consumption[:, 0] + 0.4 + bp[:, 0],
            resources,
            atol=1e-12,
            rtol=0.0,
        )

    def test_seller_map_preserves_negative_post_liquidation_wealth(self) -> None:
        grid = np.arange(-4.0, 3.0)
        params = SimpleNamespace(n_house=1, I=1, propagate_birth_entry_grant=True)
        house_cost = np.array([[0.0, 4.0]])
        sale_proceeds = np.array([[0.0, 1.0]])
        phi_choice = np.ones((1, 2, 1, 1))
        phi_choice[:, 1, :, :] = 0.8
        birth_dp = np.zeros((1, 1, 2, 2), dtype=bool)
        grant = np.zeros((1, 2, 1, 1))
        idx, wt = build_forward_tenure_transition_maps(
            params,
            grid,
            house_cost,
            sale_proceeds,
            phi_choice,
            birth_dp,
            grant,
        )
        mapped = _mapped_values(grid, idx[0, 1, 0, 0, 0], wt[0, 1, 0, 0, 0])
        np.testing.assert_allclose(mapped, np.clip(grid + 1.0, grid[0], grid[-1]))
        source = int(np.flatnonzero(grid == -3.0)[0])
        self.assertEqual(float(mapped[source]), -2.0)


class InertBaselineRegressionTests(unittest.TestCase):
    def test_nonnegative_grid_default_credit_rule_is_bit_identical(self) -> None:
        solution, _, _ = run_model_cp_dt(_inert_overrides(), verbose=False)
        expected_hashes = {
            "V": "b038ec25f8a3ab374da1a7cdf7046085354afe7aa06fe331b480f1d8b3c27e01",
            "c_pol": "47f7949c6b019fa4d2436d11b15f6fdf6275908f8c38608577b823b082a70711",
            "hR_pol": "44598b095bc6f2977228f95d80368da47ddc5d7c78c0be8d2a1e93c8fd6b0e99",
            "bp_pol": "bf7c917a915636a821c46160af404626c719b81e20b90d99dbc3b2cc5d4f737a",
            "tenure_choice": "cbf72ebcc017afcfa2e1f372b922130873ca2422056b2d05c83a977afd8c2220",
            "loc_probs": "e6477a97667f53f8852dc349e20c6eeb9bf01a528ca935dd072174b3ad33dd09",
            "fert_probs": "05d48b4673afa13f53cd74df1a941f52fb6ebb021df9a9e0d5ba0e88248629fa",
            "fert_value": "74f1ec2749cb156721cd7f150a036088cb54192b06926335db0cd806859d78ff",
            "g": "26a0bc3d969bd2e4b16ac46205106d7d9a714e0dcca6cff2e086a2a061dd4fe2",
        }
        actual_hashes = {
            name: _array_sha256(getattr(solution, name)) for name in expected_hashes
        }
        self.assertEqual(actual_hashes, expected_hashes)
        expected_scalars = {
            "mean_parity": 1.4688890253116398,
            "own_rate": 0.5294121199224653,
            "childless_rate": 0.02764776006503948,
            "total_mass": 1.0,
            "aggregate_housing_demand": 4.644972335106958,
            "best_max_abs_rel_excess": 0.3537213993779335,
        }
        self.assertEqual(
            {name: float(getattr(solution, name)) for name in expected_scalars},
            expected_scalars,
        )


class IndebtedLifecycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.solution, cls.params, cls.prices = run_model_cp_dt(
            _indebted_solution_overrides(),
            verbose=False,
        )
        cls.grid = make_grid(cls.params)
        cls.phi_state = get_phi_state_matrix(cls.params)

    def test_negative_entrant_is_feasible_and_all_policies_respect_floors(self) -> None:
        solution = self.solution
        params = self.params
        grid = self.grid
        self.assertGreater(float(np.sum(solution.g[grid < 0.0, 0, 0, 0, :, :])), 0.0)

        for j in range(params.J):
            renter_floor = np.maximum(renter_borrowing_floor(params, grid, j), grid[0])
            renter_gap = solution.bp_pol[:, 0, :, j, :, :] - renter_floor[:, None, None, None]
            self.assertGreaterEqual(float(np.min(renter_gap)), -1e-12)
            for tenure in range(1, 1 + params.n_house):
                collateral = -self.phi_state * self.prices[0] * params.H_own[tenure - 1]
                owner_floor = np.maximum(
                    owner_borrowing_floor(
                        params,
                        grid[:, None, None],
                        collateral[None, :, :],
                        j,
                    ),
                    grid[0],
                )
                gap = solution.bp_pol[:, tenure, 0, j, :, :] - owner_floor
                self.assertGreaterEqual(float(np.min(gap)), -1e-12)

    def test_taper_leaves_retired_renters_debt_free_and_owners_secured_only(self) -> None:
        solution = self.solution
        params = self.params
        grid = self.grid
        for j in range(params.J):
            age = params.age_start + j * params.da
            if age < params.debt_taper_end_age:
                continue
            renter_debt_mass = float(np.sum(solution.g[grid < -1e-12, 0, 0, j, :, :]))
            self.assertLessEqual(renter_debt_mass, 1e-14)
            for tenure in range(1, 1 + params.n_house):
                collateral = -self.phi_state * self.prices[0] * params.H_own[tenure - 1]
                unsecured_debt = grid[:, None, None] < collateral[None, :, :] - 1e-12
                bad_mass = float(
                    np.sum(
                        np.where(
                            unsecured_debt,
                            solution.g[:, tenure, 0, j, :, :],
                            0.0,
                        )
                    )
                )
                self.assertLessEqual(bad_mass, 1e-14)


class FeasibilityGateTests(unittest.TestCase):
    def test_python_owner_fallback_marks_raw_nonpositive_consumption_dead(self) -> None:
        values = eval_owner(
            np.array([-1.0, 0.0, 1.0]),
            np.array([-2.0, -1.0, 0.0]),
            np.zeros(3),
            np.array([-1.0, 0.0, 1.0]),
            0.5,
            0.5,
            0.0,
            1.0,
            0.3,
            -1.0,
            0.95,
        )
        np.testing.assert_array_equal(values, np.full(3, -1e10))

        overrides = _all_dead_overrides()
        overrides.update(
            {
                "interp_method": "pchip",
                "use_full_kernel": False,
                "use_tenure_kernel": False,
                "use_loc_kernel": False,
            }
        )
        params = apply_overrides(setup_parameters(), overrides)
        grid = make_grid(params)
        result = solve_bellman_full(
            np.array([0.16]),
            np.array([1.0]),
            params,
            grid,
            precompute_shared(params, grid),
        )
        self.assertTrue(np.all(result[0][:, 1, ...] <= DEAD_VALUE_CUTOFF))
        self.assertTrue(np.all(result[6] == 0.0))
        self.assertTrue(np.all(result[7] == 0.0))

    def test_markov_all_dead_masks_and_full_gate(self) -> None:
        overrides = _all_dead_overrides()
        overrides.update(
            {
                "use_income_types": True,
                "income_type_transition": "persistent",
                "z_grid": np.array([0.8, 1.2]),
                "z_weights": np.array([0.5, 0.5]),
                "Pi_z": np.array([[0.9, 0.1], [0.1, 0.9]]),
            }
        )
        params = apply_overrides(setup_parameters(), overrides)
        grid = make_grid(params)
        result = solve_bellman_full_markov_income(
            np.array([0.16]),
            np.array([1.0]),
            params,
            grid,
            precompute_shared(params, grid),
        )
        self.assertTrue(np.all(result[0] <= DEAD_VALUE_CUTOFF))
        self.assertIsNotNone(result[5])
        self.assertTrue(np.all(result[5] == 0.0))
        self.assertTrue(np.all(result[6] == 0.0))
        self.assertTrue(np.all(result[7] == 0.0))
        with self.assertRaises(InfeasibleThetaError):
            run_model_cp_dt(overrides, verbose=False)

    def test_all_dead_discrete_choice_rows_have_exactly_zero_probabilities(self) -> None:
        params = apply_overrides(setup_parameters(), _all_dead_overrides())
        grid = make_grid(params)
        state_details = precompute_shared(params, grid)
        result = solve_bellman_full(
            np.array([0.16]),
            np.array([1.0]),
            params,
            grid,
            state_details,
        )
        values, tenure_probs, location_probs, fertility_probs = (
            result[0],
            result[5],
            result[6],
            result[7],
        )
        self.assertTrue(np.all(values <= DEAD_VALUE_CUTOFF))
        self.assertIsNotNone(tenure_probs)
        self.assertTrue(np.all(tenure_probs == 0.0))
        self.assertTrue(np.all(location_probs == 0.0))
        self.assertTrue(np.all(fertility_probs == 0.0))

    def test_full_solve_rejects_positive_entry_mass_on_dead_state(self) -> None:
        with self.assertRaises(InfeasibleThetaError) as caught:
            run_model_cp_dt(_all_dead_overrides(), verbose=False)
        self.assertEqual(caught.exception.stage, "entry")
        self.assertGreater(caught.exception.dead_mass, 1e-12)
        self.assertTrue(caught.exception.census)

    def test_local_panel_records_infeasible_theta_separately(self) -> None:
        error = InfeasibleThetaError(
            "entry",
            0.2,
            [{"age": 18.0, "b": -1.0, "mass": 0.2, "slack": -0.1}],
        )
        income = income_process_overrides(1)
        with patch(
            "intergen_housing_fertility.local_panel.run_model_cp_dt",
            side_effect=error,
        ):
            record = run_local_panel_case(
                7,
                {"label": "dead", "theta": {}},
                J=3,
                Nb=16,
                n_house=1,
                max_iter_eq=1,
                income=income,
                rank_targets={},
                rank_weights={},
            )
        self.assertEqual(record["status"], "infeasible_theta")
        self.assertTrue(math.isinf(record["rank_loss"]))
        self.assertTrue(math.isinf(record["full_old_nonlocation_loss"]))
        self.assertFalse(record["strict_converged"])
        self.assertEqual(record["feasibility_census"], error.census)
        self.assertIn("dead-node mass", record["feasibility_error"])


if __name__ == "__main__":
    unittest.main()
