from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

from e5f_initial_age_profile import ROOT, TOP_BIN_REPRESENTATIVE, score_extra


def fixture_packet() -> dict:
    ages = np.arange(18.0, 46.0, 4.0)
    pre = np.zeros((ages.size, 4))
    post = np.zeros_like(pre)
    masses = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 4.0])
    pre[:, 0] = masses
    post[:, 1] = masses
    # Reproduce the existing observer's original 40-44 calculation explicitly.
    original_40_44 = {"0": 0.53125, "1": 0.46875, "2": 0.0, "3plus": 0.0}
    return {
        "accounting": {
            "age_cell_start": ages.tolist(),
            "pre_parity_mass_by_age": pre.tolist(),
            "post_parity_mass_by_age": post.tolist(),
        },
        "parity_shares_40_44": original_40_44,
        "metadata": {"age_projection": "uniform_birth_time"},
    }


class InitialAgeProfileTests(unittest.TestCase):
    def test_numeric_age_interpolation_and_six_score_rows(self) -> None:
        scored = score_extra(fixture_packet())
        self.assertEqual(len(scored["rows"]), 6)
        self.assertEqual(len(scored["age_profiles"]), 5)
        profile_25 = next(row for row in scored["age_profiles"] if row["age_lower"] == 25)
        # [25,26) in the age-22 cell has post share .875; [26,30) has
        # post share .5. Their mass-weighted average is .575.
        self.assertAlmostEqual(profile_25["model_share_1"], 0.575)
        self.assertAlmostEqual(profile_25["model_share_0"], 0.425)
        self.assertEqual(profile_25["overlap_weights"][1:3], [0.25, 1.0])
        self.assertEqual(
            profile_25["post_fertility_interpolation_shares"][1:3],
            [0.875, 0.5],
        )
        mean_row = next(row for row in scored["rows"] if row["moment"] == "mean_model_coded_CEB_25_29")
        self.assertAlmostEqual(mean_row["model"], 0.575)
        self.assertAlmostEqual(mean_row["weight"], 1.0 / (0.05 * mean_row["target"]) ** 2)
        self.assertFalse(scored["metadata"]["existing_twelve_rows_changed"])
        self.assertFalse(scored["metadata"]["weight_is_empirical_standard_error"])
        json.dumps(scored, allow_nan=False)

    def test_saved_observer_packet_reproduces_original_40_44_exactly(self) -> None:
        path = ROOT / (
            "output/model/e5f_matched_pf_20260909a/initial_calibration_contract/"
            "exact_loop_smoke/completed_17370427/raw/repetition_01/early_measurement.json"
        )
        with path.open(encoding="utf-8") as handle:
            packet = json.load(handle)["fertility"]["uniform_birth_time"]
        scored = score_extra(packet)
        profile = next(row for row in scored["age_profiles"] if row["age_lower"] == 40)
        for key, original_key in (
            ("model_share_0", "0"),
            ("model_share_1", "1"),
            ("model_share_2", "2"),
            ("model_share_3plus", "3plus"),
        ):
            self.assertAlmostEqual(profile[key], packet["parity_shares_40_44"][original_key], places=14)
        reconstructed = (
            profile["target_share_1"]
            + 2.0 * profile["target_share_2"]
            + TOP_BIN_REPRESENTATIVE * profile["target_share_3plus"]
        )
        self.assertAlmostEqual(reconstructed, profile["target_mean_model_coded_CEB"], places=14)

    def test_wrong_projection_or_original_gate_fails_closed(self) -> None:
        packet = fixture_packet()
        packet["metadata"]["age_projection"] = "constant_post_cell"
        with self.assertRaisesRegex(ValueError, "uniform_birth_time"):
            score_extra(packet)
        packet = fixture_packet()
        packet["parity_shares_40_44"]["0"] += 0.01
        with self.assertRaisesRegex(RuntimeError, "does not reproduce"):
            score_extra(packet)


if __name__ == "__main__":
    unittest.main()
