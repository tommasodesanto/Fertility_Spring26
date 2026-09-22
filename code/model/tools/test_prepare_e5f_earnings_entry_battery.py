#!/usr/bin/env python3
"""Lightweight no-solve checks for the four-cell battery plan generator."""
import importlib.util
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "prepare_e5f_earnings_entry_battery", HERE / "prepare_e5f_earnings_entry_battery.py")
BATTERY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BATTERY)
PLAN = HERE.parents[2] / "tmp/earnings_wealth_direct_period_20260922_v5/plan.json"


class BatteryPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = json.loads(PLAN.read_text())

    def test_four_cell_factorial_and_fixed_entry_receipt(self):
        BATTERY.validate_cells()
        smoke = BATTERY.smoke_parameters(self.base)
        self.assertEqual(smoke["beta_annual"], 0.985)
        for cell_id, spec in BATTERY.CELLS.items():
            plan = BATTERY._cell_plan(self.base, cell_id, smoke, f"{cell_id}_smoke")
            self.assertEqual(plan["income_specification"]["mapping"], "direct_period")
            self.assertEqual(plan["income_specification"]["max_relative_discrete_level_covariance_error"], 0.15)
            self.assertEqual(plan["income_entry_battery"]["persistent_states"], 7)
            self.assertEqual(plan["income_entry_battery"]["iid_states"], spec["income"]["n_iid"])
            self.assertEqual(plan["entry_specification"]["rule"], spec["entry_rule"])
            self.assertNotIn("entry_wealth", plan["entry_specification"])
            if spec["entry_rule"] == "fixed_reference_marginal":
                receipt = plan["entry_specification"]["reference_marginal_receipt"]
                self.assertEqual(receipt["baseline_mean"], 0.18651967924681825)
                self.assertIn("must be measured", receipt["frontier_censored_mass"])
            self.assertEqual(plan["initial_psi"], self.base["initial_psi"])

    def test_sixty_matched_native_proposals_are_deterministic(self):
        first = BATTERY.candidate_pool(self.base)
        second = BATTERY.candidate_pool(self.base)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 60)
        self.assertEqual(len({BATTERY.fingerprint(point) for point in first}), 60)
        self.assertEqual(first[0], {k: float(v) for k, v in self.base["starting_structural_parameters"].items()})


if __name__ == "__main__":
    unittest.main()
