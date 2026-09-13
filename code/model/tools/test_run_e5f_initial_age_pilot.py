from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import e5f_initial_age_profile as age
import run_e5f_initial_age_pilot as pilot


class InitialAgePilotTests(unittest.TestCase):
    def test_all_fixed_source_and_target_pins_match_local_artifacts(self) -> None:
        tools = Path(__file__).resolve().parent
        self.assertEqual(
            pilot.digest(tools / "run_e5f_rebated_initial_overnight.py"),
            pilot.REBATED_HELPER_SHA256,
        )
        self.assertEqual(
            pilot.digest(tools / "run_e5f_joint_rebated_initial_probe.py"),
            pilot.JOINT_ADAPTER_SHA256,
        )
        self.assertEqual(
            pilot.digest(tools / "run_e5f_joint_rebated_initial_scored.py"),
            pilot.SCORED_EVALUATOR_SHA256,
        )
        self.assertEqual(pilot.digest(Path(age.__file__)), pilot.AGE_HELPER_SHA256)
        self.assertEqual(pilot.digest(Path(age.EMPIRICAL_PACKET)), pilot.CPS_AGE_DATA_SHA256)

    def test_evaluator_command_uses_rebated_joint_cli_and_fixed_pins(self) -> None:
        command = pilot.evaluator_command(
            Path("scored.py"), Path("helper.py"), Path("joint.py"),
            Path("template"), Path("proposal.json"), Path("output"),
        )
        self.assertEqual(command[1:3], ["-B", str(Path("scored.py").resolve())])
        self.assertEqual(command[command.index("--helper-sha256") + 1], pilot.REBATED_HELPER_SHA256)
        self.assertEqual(command[command.index("--joint-sha256") + 1], pilot.JOINT_ADAPTER_SHA256)
        self.assertEqual(command[command.index("--proposal") + 1], str(Path("proposal.json").resolve()))

    def test_augmented_table_excludes_separate_normalization(self) -> None:
        original = {"target_fit": [
            *[{"moment": f"m{i}", "scored": True} for i in range(12)],
            {"moment": "fertility_normalization", "scored": False},
        ]}
        extra = {"rows": [{"moment": f"x{i}"} for i in range(6)]}
        rows = pilot.augmented_target_rows(original, extra)
        self.assertEqual(len(rows), 18)
        self.assertEqual(sum(row["system"] == "original_12" for row in rows), 12)
        self.assertNotIn("fertility_normalization", {row["moment"] for row in rows})

    def test_exact_six_proposals_use_half_percent_of_existing_span(self) -> None:
        center = {
            "first_birth_fixed_cost": 4.0,
            "kappa_fert": 1.0,
            "kappa_fert_continuation": 49.9,
            "other": 7.0,
        }
        initial = {"structural_candidate": center}
        objective = {
            "parameter_restrictions": [
                {"parameter": "first_birth_fixed_cost", "lower": 0.0, "upper": 8.0},
                {"parameter": "kappa_fert", "lower": 0.02, "upper": 50.0},
                {"parameter": "kappa_fert_continuation", "lower": 0.02, "upper": 50.0},
            ]
        }
        cases = pilot.proposals(initial, objective)
        self.assertEqual(len(cases), 6)
        self.assertEqual(cases[0]["parameters"], center)
        self.assertAlmostEqual(cases[1]["parameters"]["first_birth_fixed_cost"], 3.96)
        self.assertAlmostEqual(cases[2]["parameters"]["first_birth_fixed_cost"], 4.04)
        self.assertAlmostEqual(cases[3]["parameters"]["kappa_fert"], 1.0 - 0.005 * 49.98)
        self.assertAlmostEqual(cases[4]["parameters"]["kappa_fert"], 1.0 + 0.005 * 49.98)
        self.assertEqual(cases[5]["parameters"]["kappa_fert_continuation"], 50.0)
        self.assertTrue(cases[5]["clipped_to_existing_bounds"])
        for case in cases:
            self.assertEqual(set(case["parameters"]), set(center))

    def test_explicit_data_override_changes_metadata_and_uses_requested_file(self) -> None:
        packet_path = Path(
            "output/model/e5f_matched_pf_20260909a/initial_calibration_contract/"
            "exact_loop_smoke/completed_17370427/raw/repetition_01/early_measurement.json"
        )
        packet = json.loads(packet_path.read_text())["fertility"]["uniform_birth_time"]
        source = Path(age.EMPIRICAL_PACKET)
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory) / "age.json"
            copied.write_bytes(source.read_bytes())
            score = pilot.extra_score(packet, copied)
            self.assertEqual(score["metadata"]["empirical_source"], str(copied.resolve()))
            self.assertEqual(score["metadata"]["empirical_source_sha256"], pilot.digest(copied))
            self.assertEqual(len(score["rows"]), 6)
            self.assertEqual(len(score["age_profiles"]), 5)


if __name__ == "__main__":
    unittest.main()
