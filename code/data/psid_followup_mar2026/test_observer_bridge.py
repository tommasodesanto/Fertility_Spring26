"""Focused tests for the four-year first-birth observation bridge."""

import unittest

from observer_bridge import (
    ObserverBridgeContract,
    model_age_cell_lower_bound,
    model_period_for_event_year,
    observe_model_rooms,
)


class ObserverBridgeTests(unittest.TestCase):
    def test_contract_pins_auxiliary_measurement_without_adoption(self) -> None:
        contract = ObserverBridgeContract()
        contract.validate()
        self.assertEqual(contract.baseline_event_years, (-3, -2))
        self.assertEqual(contract.outcome_event_years, (3, 4))
        self.assertEqual(contract.main_birth_offset_years, 0.0)
        self.assertIn("not adopted", contract.adoption_status)

    def test_decision_date_mapping_uses_floor_for_prebirth_interviews(self) -> None:
        # U=0: baseline -3/-2 is before the birth decision; +3/+4 straddles
        # the destination boundary and must not be collapsed into one period.
        self.assertEqual(model_period_for_event_year(7, -3, 0.0), 6)
        self.assertEqual(model_period_for_event_year(7, -2, 0.0), 6)
        self.assertEqual(model_period_for_event_year(7, 3, 0.0), 7)
        self.assertEqual(model_period_for_event_year(7, 4, 0.0), 8)

    def test_midpoint_sensitivity_changes_observation_periods(self) -> None:
        self.assertEqual(model_period_for_event_year(7, -3, 2.0), 6)
        self.assertEqual(model_period_for_event_year(7, -2, 2.0), 7)
        self.assertEqual(model_period_for_event_year(7, 3, 2.0), 8)
        self.assertEqual(model_period_for_event_year(7, 4, 2.0), 8)

    def test_exact_period_boundary_and_near_boundary(self) -> None:
        self.assertEqual(model_period_for_event_year(4, 4, 0.0), 5)
        self.assertEqual(model_period_for_event_year(4, 1, 3.0), 5)
        self.assertEqual(model_period_for_event_year(4, -3, 3.999), 4)
        self.assertEqual(model_period_for_event_year(4, -4, 3.999), 3)

    def test_age_cells_match_four_year_model_grid(self) -> None:
        self.assertEqual(model_age_cell_lower_bound(18), 18)
        self.assertEqual(model_age_cell_lower_bound(21), 18)
        self.assertEqual(model_age_cell_lower_bound(22), 22)
        self.assertEqual(model_age_cell_lower_bound(82), 82)
        with self.assertRaises(ValueError):
            model_age_cell_lower_bound(83)

    def test_room_lookup_uses_mapped_date_and_fails_on_missing_period(self) -> None:
        history = {6: 4.0, 7: 5.0, 8: 5.5}
        self.assertEqual(observe_model_rooms(history, 7, 3, 0.0), 5.0)
        self.assertEqual(observe_model_rooms(history, 7, 4, 0.0), 5.5)
        with self.assertRaisesRegex(KeyError, "period 9"):
            observe_model_rooms(history, 7, 8, 0.0)

    def test_invalid_period_offset_and_age_fail_loudly(self) -> None:
        with self.assertRaises(ValueError):
            model_period_for_event_year(0, 0, 4.0)
        with self.assertRaises(ValueError):
            model_period_for_event_year(0, 0, -0.1)
        with self.assertRaises(TypeError):
            model_period_for_event_year(True, 0, 0.0)
        with self.assertRaises(ValueError):
            model_age_cell_lower_bound(17)


if __name__ == "__main__":
    unittest.main()
