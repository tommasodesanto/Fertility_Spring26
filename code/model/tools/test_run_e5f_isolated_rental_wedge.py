import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import run_e5f_isolated_rental_wedge as driver


class TestIsolatedRentalWedge(unittest.TestCase):
    def test_case_contract_is_six_arms_and_fourteen_total_bound(self):
        self.assertEqual([case.name for case in driver.CASES], [
            "cap6zero", "cap10zero", "cap10s005", "cap10s02", "cap10s1", "cap10s02_phi1",
        ])
        self.assertEqual(driver.MAX_LIFECYCLE_SOLVES, 14)
        self.assertEqual(driver.CASE_BY_NAME["cap6zero"].slope, 0.0)
        self.assertTrue(all(case.cap == 10.0 for case in driver.CASES[2:]))
        self.assertEqual(driver.CASE_BY_NAME["cap10s02_phi1"].financed_share, 1.0)

    def test_cost_has_zero_intercept_and_six_room_knee(self):
        self.assertAlmostEqual(driver.rental_cost(.2, 4, .2), .8)
        self.assertAlmostEqual(driver.rental_cost(.2, 8, .2), 8 * (.2 + .2 * 2))

    def test_budget_gate_charges_wedge_and_fails_if_omitted(self):
        row = {
            "tenure": "renter", "rent": .2, "rooms": 8,
            "consumption": 1.0, "saving": 0.0,
            "resources": 2.55, "mass": 1.0,
        }
        with self.assertRaisesRegex(RuntimeError, "budget gate"):
            driver.dated_budget_gate([row], slope=.2)
        row["resources"] = 5.8
        dated = driver.dated_budget_gate([row], slope=.2)
        independent = driver.independent_budget_gate([row], slope=.2)
        self.assertEqual(dated, independent)
        self.assertAlmostEqual(dated["weighted_positive_budget_gap"], 0.0)

    def test_missing_renter_cost_inputs_fail_closed(self):
        row = {"tenure": "renter", "rooms": 8, "consumption": 1, "saving": 0, "resources": 3}
        with self.assertRaisesRegex(driver.ContractError, "rent and realized rooms"):
            driver.audit_budget_rows([row], slope=.2)

    def test_wedge_cases_refuse_unreviewed_port(self):
        source = {"wedge_port_reviewed": False}
        with self.assertRaisesRegex(driver.ReviewRequiredError, "not marked reviewed"):
            driver.arm_parameters({"chi": 1.0, "hR_max": 6.0}, "cap10s005", source=source)

    def test_control_preserves_checkpoint_parameters(self):
        source = {"wedge_port_reviewed": False}
        base = {"chi": 1.049, "hR_max": 6.0, "psi_child": .1, "entry": .2}
        result = driver.arm_parameters(base, "cap6zero", source=source)
        self.assertEqual(result, base)

    def test_financing_contrast_updates_only_phi_in_mapping_contract(self):
        source = {"wedge_port_reviewed": True}
        base = {"chi": 1.049, "hR_max": 6.0, "phi": np.array([0.8, 0.8]), "psi_child": .1, "entry": .2}
        result = driver.arm_parameters(base, "cap10s02_phi1", source=source, control_receipt={"status": "passed_control"})
        self.assertTrue(np.array_equal(result["phi"], [1.0, 1.0]))
        self.assertEqual(result["psi_child"], base["psi_child"])
        self.assertEqual(result["entry"], base["entry"])

    def test_build_arm_financing_contrast_rejects_debt_derivative_changes(self):
        import types
        base = types.SimpleNamespace(hR_max=6.0, phi=np.array([0.8, 0.8]), lambda_d=0.0,
                                     income_candidate_fingerprint=None, chi=1.0,
                                     debt_caps=np.array([0.0]), debt_taper_weights=np.array([1.0]),
                                     mean_labor_income_by_age=np.array([2.0]), owner_ltv_multipliers=np.array([1.0]))
        class Params:
            @staticmethod
            def build_debt_caps(p):
                p.debt_caps = np.asarray(p.phi) * 2.0
            @staticmethod
            def rental_wedge_active(p):
                return getattr(p, "rental_wedge_slope", 0.0) > 0
        case = driver.CASE_BY_NAME["cap10s02_phi1"]
        with self.assertRaisesRegex(driver.ContractError, "unexpected parameter changes"):
            driver.build_arm(base, case, Params)

    def test_build_arm_financing_contrast_accepts_unchanged_debt_derivatives(self):
        import types
        base = types.SimpleNamespace(hR_max=6.0, phi=np.array([0.8, 0.8]), lambda_d=0.0,
                                     income_candidate_fingerprint=None, chi=1.0,
                                     debt_caps=np.array([0.0]), debt_taper_weights=np.array([1.0]),
                                     mean_labor_income_by_age=np.array([2.0]), owner_ltv_multipliers=np.array([1.0]))
        class Params:
            @staticmethod
            def build_debt_caps(p):
                return p
            @staticmethod
            def rental_wedge_active(p):
                return getattr(p, "rental_wedge_slope", 0.0) > 0
        result, record = driver.build_arm(base, driver.CASE_BY_NAME["cap10s02_phi1"], Params)
        self.assertTrue(np.array_equal(result.phi, [1.0, 1.0]))
        self.assertEqual(record["financed_share"], 1.0)
        self.assertIn("phi", record["changed_fields"])

    def test_manifest_rejects_active_checkout(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.json"
            manifest.write_text(json.dumps({"source_files": {"code/model/solver.py": "x"}}))
            with self.assertRaisesRegex(driver.ContractError, "active checkout"):
                driver.load_source_manifest(driver.ROOT / "code/model", manifest)


if __name__ == "__main__":
    unittest.main()
