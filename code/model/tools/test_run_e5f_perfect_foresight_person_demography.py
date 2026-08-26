"""Tiny-grid tests for the coherent-person perfect-foresight bridge."""

from __future__ import annotations

import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from intergen_eqscale_seq_optimized import solver as model

import run_dynamic_population_transition as calendar
import run_e5f_open_population_transition as transition
import run_e5f_perfect_foresight_person_demography as person_pf


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


class TinyPersonPerfectForesightTests(unittest.TestCase):
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
        unit_g_pre, _ = calendar.reconstruct_stationary_pre_fertility(
            solution, policy, parameters, b_grid, shared
        )
        supply_rule, _ = calendar.normalize_date0_housing_supply(
            unit_g_pre,
            policy,
            parameters,
            b_grid,
            shared,
            "static-elastic",
        )
        maximum_age = 60
        shape = (2, maximum_age + 1)
        survival = np.full(shape, 0.90)
        survival[:, -1] = 0.50
        migration = np.zeros(shape)
        migration[:, 18:30] = 2.5e-4
        headship = np.zeros(shape)
        headship[:, 18:42] = 0.03
        primitives = person_pf.AnnualDemographicPrimitives(
            start_year=2023,
            last_empirical_year=2100,
            scale_model_units_per_person=1.0,
            initial_person_state=person_pf.CohortState(
                2023, np.zeros(shape), np.zeros(shape)
            ),
            birth_sex_shares={2100: np.array([0.51, 0.49])},
            survival={2100: survival},
            net_migration={2100: migration},
            headship_rates=headship,
            raw_acs_headship_rates=headship,
            model_age_headship_alignment_factors=np.ones(int(parameters.J)),
            initial_household_age_mass=np.zeros(int(parameters.J)),
            initial_person_head_age_mass=np.zeros(int(parameters.J)),
            source_paths={"synthetic": Path("synthetic")},
        )
        cls.parameters = parameters
        cls.b_grid = b_grid
        cls.policy = policy
        cls.unit_g_pre = unit_g_pre
        cls.supply_rule = supply_rule
        cls.primitives = primitives

    def test_terminal_household_person_mapping_closes(self) -> None:
        result = person_pf.solve_terminal_household_person_fixed_point(
            policy=self.policy,
            parameters=self.parameters,
            b_grid=self.b_grid,
            initial_g_pre=self.unit_g_pre,
            demographic_primitives=self.primitives,
            supply_rule=self.supply_rule,
            distribution_tolerance=2e-7,
            birth_rate_tolerance=2e-7,
            one_step_tolerance=2e-7,
            maximum_iterations=180,
            damping=0.65,
        )
        diagnostics = {
            "iterations": result.iterations,
            "distribution_mapping_relative_l1": (
                result.distribution_mapping_relative_l1
            ),
            "annual_birth_rate_relative_gap": (
                result.annual_birth_rate_relative_gap
            ),
            "person_one_step_relative_l1": result.person_one_step_relative_l1,
            "head_one_step_relative_l1": result.head_one_step_relative_l1,
            "age_head_one_step_max_abs": result.age_head_one_step_max_abs,
            "household_person_head_gap": result.household_person_head_gap,
            "renewal_ratio": result.renewal_ratio,
        }
        self.assertTrue(result.converged, diagnostics)
        self.assertLess(result.renewal_ratio, 1.0)
        self.assertLess(result.distribution_mapping_relative_l1, 2e-7)
        self.assertLess(result.person_one_step_relative_l1, 2e-7)
        self.assertLess(result.head_one_step_relative_l1, 2e-7)
        self.assertLess(abs(result.household_person_head_gap), 2e-7)

    def test_terminal_mapping_rejects_heads_outside_model_age_support(self) -> None:
        invalid_headship = self.primitives.headship_rates.copy()
        invalid_headship[:, 42] = 0.01
        invalid = replace(self.primitives, headship_rates=invalid_headship)
        with self.assertRaisesRegex(ValueError, "outside the household model"):
            person_pf.solve_terminal_household_person_fixed_point(
                policy=self.policy,
                parameters=self.parameters,
                b_grid=self.b_grid,
                initial_g_pre=self.unit_g_pre,
                demographic_primitives=invalid,
                supply_rule=self.supply_rule,
                maximum_iterations=1,
            )


if __name__ == "__main__":
    unittest.main()
