"""Regression tests for the optional calendar-time Bellman continuation."""
from __future__ import annotations

import unittest

import numpy as np

from intergen_eqscale_seq_optimized.solver import (
    precompute_shared,
    run_model_cp_dt,
    solve_bellman_full_markov_income,
)


def _tiny_markov() -> dict[str, object]:
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
    }


def _stationary_objects():
    solution, parameters, price = run_model_cp_dt(_tiny_markov(), verbose=False)
    b_grid = np.asarray(solution.b_grid, dtype=float)
    shared = precompute_shared(parameters, b_grid)
    rent = parameters.user_cost_rate * np.asarray(price, dtype=float)
    return solution, parameters, np.asarray(price, dtype=float), b_grid, shared, rent


class CalendarTimeContinuationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        (
            cls.solution,
            cls.parameters,
            cls.price,
            cls.b_grid,
            cls.shared,
            cls.rent,
        ) = _stationary_objects()

    def solve(self, continuation_V=None):
        return solve_bellman_full_markov_income(
            self.rent,
            self.price,
            self.parameters,
            self.b_grid,
            self.shared,
            continuation_V=continuation_V,
        )

    def test_absent_calendar_continuation_is_bitwise_inert(self) -> None:
        implicit = solve_bellman_full_markov_income(
            self.rent, self.price, self.parameters, self.b_grid, self.shared
        )
        explicit_none = self.solve(continuation_V=None)
        for left, right in zip(implicit[:-1], explicit_none[:-1], strict=True):
            if left is None:
                self.assertIsNone(right)
            else:
                self.assertTrue(np.array_equal(left, right))

    def test_stationary_value_is_calendar_time_fixed_point(self) -> None:
        stationary = self.solve()
        calendar = self.solve(
            continuation_V=np.asarray(self.solution.V, dtype=float)
        )
        for left, right in zip(stationary[:-1], calendar[:-1], strict=True):
            if left is None:
                self.assertIsNone(right)
            else:
                np.testing.assert_array_equal(left, right)

    def test_calendar_continuation_changes_only_nonterminal_recursion(self) -> None:
        baseline = self.solve(
            continuation_V=np.asarray(self.solution.V, dtype=float)
        )
        shifted = self.solve(
            continuation_V=np.asarray(self.solution.V, dtype=float) + 0.1
        )
        np.testing.assert_array_equal(
            baseline[0][:, :, :, -1], shifted[0][:, :, :, -1]
        )
        self.assertFalse(
            np.array_equal(baseline[0][:, :, :, :-1], shifted[0][:, :, :, :-1])
        )

    def test_calendar_continuation_rejects_wrong_shape(self) -> None:
        with self.assertRaisesRegex(ValueError, "continuation value has shape"):
            self.solve(
                continuation_V=np.asarray(self.solution.V, dtype=float)[..., :-1]
            )


if __name__ == "__main__":
    unittest.main()
