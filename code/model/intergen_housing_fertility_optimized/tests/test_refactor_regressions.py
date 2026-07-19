from __future__ import annotations

import unittest

import numpy as np

from intergen_housing_fertility_optimized.parameters import apply_overrides, setup_parameters
from intergen_housing_fertility_optimized.m5_profile import M5_TARGET_SET, m5_target_system
from intergen_housing_fertility_optimized.state_layout import StateLayout
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


if __name__ == "__main__":
    unittest.main()
