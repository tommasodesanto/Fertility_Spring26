from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intergen_housing_fertility_optimized.parameters import apply_overrides, setup_parameters
from intergen_housing_fertility_optimized.m5_profile import M5_TARGET_SET, m5_target_system
from intergen_housing_fertility_optimized.promotion_contract import (
    CALIBRATION_CASES,
    DEMAND_PRICES,
    FREE_PARAMETER_BOUNDS,
    ROOT_CASES,
    bound_cases,
)
from intergen_housing_fertility_optimized.promotion_worker import compare_failures
from intergen_housing_fertility_optimized.state_layout import StateLayout
from intergen_housing_fertility_optimized.solver import refine_one_market_markov_income
from intergen_housing_fertility_optimized.utils import decode_flat_family_state, flat_nc


class IncomeProcessOverrideTests(unittest.TestCase):
    def test_persistence_override_rebuilds_same_shape_transition(self) -> None:
        baseline = setup_parameters()
        old_transition = baseline.Pi_z.copy()
        updated = apply_overrides(baseline, {"income_shock_persistence": 0.10})
        self.assertFalse(np.array_equal(old_transition, updated.Pi_z))
        np.testing.assert_allclose(updated.z_weights @ updated.Pi_z, updated.z_weights, atol=1e-12)

    def test_weight_override_rebuilds_stationary_transition(self) -> None:
        updated = apply_overrides(setup_parameters(), {"z_weights": [0.8, 0.1, 0.1]})
        np.testing.assert_allclose(updated.z_weights @ updated.Pi_z, updated.z_weights, atol=1e-12)

    def test_inconsistent_explicit_transition_is_rejected(self) -> None:
        baseline = setup_parameters()
        with self.assertRaisesRegex(ValueError, "inconsistent with z_weights"):
            apply_overrides(
                baseline,
                {"z_weights": [0.8, 0.1, 0.1], "Pi_z": baseline.Pi_z.copy()},
            )

    def test_unknown_override_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown parameter override"):
            apply_overrides(setup_parameters(), {"income_shock_persistnce": 0.5})


class FamilyStateLayoutTests(unittest.TestCase):
    def test_decoder_inverts_fortran_flattening(self) -> None:
        n_parity = 3
        n_child_states = 4
        values = np.arange(n_parity * n_child_states).reshape(
            1, n_parity, n_child_states, order="F"
        )
        flattened = flat_nc(values, 1, n_parity * n_child_states)
        for column in range(flattened.shape[1]):
            parity, child_state = decode_flat_family_state(column, n_parity)
            self.assertEqual(flattened[0, column], values[0, parity, child_state])

    def test_named_layout_round_trips_every_family_state(self) -> None:
        layout = StateLayout(120, 6, 1, 17, 5, 3, 4)
        for parity in range(layout.parity):
            for child_state in range(layout.child_state):
                column = layout.encode_family_state(parity, child_state)
                self.assertEqual(layout.decode_family_state(column), (parity, child_state))


class ProductionTargetContractTests(unittest.TestCase):
    def test_m5_target_system_is_atomic_identified_and_stable(self) -> None:
        system = m5_target_system()
        self.assertEqual(system.name, M5_TARGET_SET)
        self.assertEqual(system.count, 15)
        self.assertIn("aggregate_mean_occupied_rooms_18_85", system.moment_names)
        self.assertEqual(len(system.fingerprint), 64)
        system.require_identified(14)
        with self.assertRaisesRegex(ValueError, "underidentified"):
            system.require_identified(16)

    def test_target_system_loss_uses_declared_rows_and_weights(self) -> None:
        system = m5_target_system()
        exact = system.targets_dict()
        self.assertEqual(system.loss(exact), 0.0)
        changed = dict(exact)
        changed[system.moment_names[0]] += 0.25
        self.assertAlmostEqual(system.loss(changed), system.weights[0] * 0.25**2)


class DirectRootFallbackTests(unittest.TestCase):
    def test_wrong_direction_uses_bilateral_fallback_and_finds_root(self) -> None:
        parameters = SimpleNamespace(
            p_min=0.1,
            p_max=3.0,
            scalar_market_refine_expand=1.5,
            scalar_market_refine_max_expand=6,
            markov_equilibrium_method="direct_brent",
            scalar_market_refine_iter=16,
            scalar_market_direct_iter=24,
            scalar_market_direct_expand=1.1,
            scalar_market_refine_method="brent",
            tol_eq=1e-8,
        )

        def increasing_excess(price, *_args, **_kwargs):
            p = float(np.asarray(price).reshape(-1)[0])
            return SimpleNamespace(
                housing_demand=np.array([1.0 + p - 1.2]),
                housing_supply=np.array([1.0]),
            )

        initial = increasing_excess(np.array([1.0]))
        with patch(
            "intergen_housing_fertility_optimized.solver.solve_markov_income_at_prices",
            side_effect=increasing_excess,
        ):
            _, price, residual, info = refine_one_market_markov_income(
                np.array([1.0]),
                initial,
                0.2,
                parameters,
                np.array([0.0, 1.0]),
                verbose=False,
            )
        self.assertTrue(info["directional_fallback_used"])
        self.assertTrue(info["bracket_found"])
        self.assertLessEqual(residual, parameters.tol_eq)
        self.assertAlmostEqual(float(price[0]), 1.2, places=7)


class PromotionContractTests(unittest.TestCase):
    def test_battery_dimensions_are_fixed_and_cover_every_live_bound(self) -> None:
        self.assertEqual(len(FREE_PARAMETER_BOUNDS), 14)
        self.assertEqual(len(bound_cases()), 1 + 2 * len(FREE_PARAMETER_BOUNDS))
        self.assertEqual(len(DEMAND_PRICES), 17)
        self.assertEqual((min(DEMAND_PRICES), max(DEMAND_PRICES)), (0.01, 30.0))
        self.assertEqual(len(ROOT_CASES), 6)
        self.assertEqual(len(CALIBRATION_CASES), 10)

    def test_matching_infeasibility_passes_but_matching_program_error_fails(self) -> None:
        infeasible = {
            "status": "infeasible",
            "exception_type": "InfeasibleThetaError",
            "stage": "forward_age_22",
            "dead_mass": 1e-5,
        }
        self.assertTrue(compare_failures(infeasible, dict(infeasible))["pass"])
        error = {
            "status": "error",
            "exception_type": "ValueError",
            "stage": None,
            "dead_mass": None,
        }
        self.assertFalse(compare_failures(error, dict(error))["pass"])


if __name__ == "__main__":
    unittest.main()
