from __future__ import annotations

import unittest
import json
import tempfile
from pathlib import Path

import numpy as np

from intergen_housing_fertility_optimized.calibration_search import (
    DOMAIN,
    build_search_target_system,
    load_start_theta,
    parameter_rows,
    theta_from_unit,
    theta_from_unit_with_fixed_beta,
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

    def test_fixed_beta_profile_overrides_domain_and_marks_parameter_fixed(self) -> None:
        unit = unit_from_theta(M5_THETA)
        theta = theta_from_unit_with_fixed_beta(unit, 0.9999)
        self.assertAlmostEqual(theta["beta"], 0.9999**4)
        beta_row = next(
            row
            for row in parameter_rows(theta, fixed_beta_annual=0.9999)
            if row["parameter"] == "beta_annual"
        )
        self.assertEqual(beta_row["role"], "profile_fixed")
        self.assertAlmostEqual(beta_row["estimate"], 0.9999)

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

    def test_search_weight_tilt_changes_navigation_not_canonical_contract(self) -> None:
        canonical = new_moment_target_system()
        moment = "aggregate_wealth_to_annual_after_tax_earnings"
        search, multipliers = build_search_target_system(canonical, [f"{moment}=4"])
        self.assertEqual(multipliers, {moment: 4.0})
        self.assertEqual(search.moment_names, canonical.moment_names)
        self.assertEqual(search.target_values, canonical.target_values)
        self.assertEqual(
            search.weights[canonical.moment_names.index(moment)],
            4.0 * canonical.weights[canonical.moment_names.index(moment)],
        )
        self.assertEqual(canonical, new_moment_target_system())

    def test_search_weight_tilt_rejects_unknown_duplicate_or_extreme_specs(self) -> None:
        canonical = new_moment_target_system()
        moment = canonical.moment_names[0]
        for specifications in (
            ["not_a_moment=2"],
            [f"{moment}=2", f"{moment}=3"],
            [f"{moment}=100"],
            ["malformed"],
        ):
            with self.subTest(specifications=specifications):
                with self.assertRaises(ValueError):
                    build_search_target_system(canonical, specifications)


if __name__ == "__main__":
    unittest.main()
