from __future__ import annotations

import unittest
import json
import tempfile
from pathlib import Path

import numpy as np

from intergen_housing_fertility_optimized.calibration_search import (
    DOMAIN,
    load_start_theta,
    parameter_rows,
    theta_from_unit,
    unit_from_theta,
    validate_contract,
)
from intergen_housing_fertility_optimized.m5_profile import M5_THETA
from intergen_housing_fertility_optimized.new_moment_profile import (
    NEW_MOMENT_SEED,
    new_moment_target_system,
)
from intergen_housing_fertility_optimized.promotion_contract import FREE_PARAMETER_BOUNDS


class CalibrationSearchContractTests(unittest.TestCase):
    def test_domain_matches_promoted_live_bounds_and_is_identified(self) -> None:
        validate_contract()
        self.assertEqual(len(DOMAIN), 14)
        self.assertEqual(
            tuple((name, lower, upper) for name, lower, upper, _ in DOMAIN),
            FREE_PARAMETER_BOUNDS,
        )

    def test_m5_theta_round_trips_through_transformed_search_coordinates(self) -> None:
        recovered = theta_from_unit(unit_from_theta(M5_THETA))
        for name, expected in M5_THETA.items():
            np.testing.assert_allclose(recovered[name], expected, rtol=1e-12, atol=1e-12)

    def test_parameter_table_reports_every_free_coordinate_and_restriction(self) -> None:
        rows = parameter_rows(M5_THETA)
        self.assertEqual(len(rows), 15)
        self.assertEqual(sum(row["role"] == "estimated" for row in rows), 14)
        restriction = next(row for row in rows if row["parameter"] == "theta_n")
        self.assertEqual(restriction["role"], "fixed")
        self.assertEqual(restriction["estimate"], 0.0)

    def test_collector_result_can_seed_a_same_target_continuation(self) -> None:
        payload = {
            "target_fingerprint": new_moment_target_system().fingerprint,
            "selected": {"theta": NEW_MOMENT_SEED},
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "results.json"
            path.write_text(json.dumps(payload))
            loaded = load_start_theta(
                path,
                expected_target_fingerprint=new_moment_target_system().fingerprint,
            )
        self.assertEqual(loaded, NEW_MOMENT_SEED)

    def test_continuation_seed_rejects_wrong_target_or_external_restriction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "results.json"
            path.write_text(json.dumps({"theta": NEW_MOMENT_SEED}))
            with self.assertRaisesRegex(ValueError, "declare a target fingerprint"):
                load_start_theta(
                    path,
                    expected_target_fingerprint=new_moment_target_system().fingerprint,
                )
            path.write_text(json.dumps({"target_fingerprint": "wrong", "theta": NEW_MOMENT_SEED}))
            with self.assertRaisesRegex(ValueError, "fingerprint"):
                load_start_theta(
                    path,
                    expected_target_fingerprint=new_moment_target_system().fingerprint,
                )
            bad_theta = {**NEW_MOMENT_SEED, "theta_n": 0.1}
            path.write_text(json.dumps({"theta": bad_theta}))
            with self.assertRaisesRegex(ValueError, "theta_n=0"):
                load_start_theta(path)


if __name__ == "__main__":
    unittest.main()
