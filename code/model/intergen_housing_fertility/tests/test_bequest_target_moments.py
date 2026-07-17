from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from intergen_housing_fertility.calibration import get_target_set
from intergen_housing_fertility.solver import add_annual_gross_estate_wealth_moments


class BequestTargetMomentTests(unittest.TestCase):
    def test_estate_targets_use_gross_home_value_and_live_parity_bins(self) -> None:
        p = SimpleNamespace(
            J=5,
            I=1,
            age_start=66.0,
            da=4.0,
            period_years=1.0,
            J_R=0,
            income=np.ones((1, 5)),
            n_parity=3,
            n_child_states=1,
            n_house=1,
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
            psi=0.5,
        )
        bg = np.array([0.0, 10.0])
        ph = np.array([10.0])
        g = np.zeros((2, 2, 1, 5, 1, 3, 1))

        # Ages 66 and 70 identify the 2+ minus 1-child median gap.
        g[0, 0, 0, 0, 0, 1, 0] = 1.0  # one child, estate 0
        g[0, 1, 0, 0, 0, 2, 0] = 1.0  # two plus, estate 20
        g[1, 0, 0, 1, 0, 1, 0] = 1.0  # one child, estate 10
        g[1, 1, 0, 1, 0, 2, 0] = 1.0  # two plus, estate 30

        # Ages 78 and 82 yield old-tail values 0, 10, 20, and 30.
        g[0, 0, 0, 3, 0, 1, 0] = 1.0
        g[0, 1, 0, 3, 0, 2, 0] = 1.0
        g[1, 1, 0, 4, 0, 1, 0] = 1.0
        g[1, 0, 0, 4, 0, 2, 0] = 1.0

        stats = SimpleNamespace()
        add_annual_gross_estate_wealth_moments(stats, g, p, bg, ph)

        self.assertEqual(stats.old_total_estate_wealth_to_annual_income_median_7684, 10.0)
        self.assertEqual(stats.old_total_estate_wealth_to_annual_income_p90_7684, 30.0)
        self.assertEqual(stats.old_total_estate_wealth_to_annual_income_p90_p50_7684, 3.0)
        self.assertEqual(
            stats.old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575,
            20.0,
        )

    def test_internal_target_set_replaces_three_legacy_old_age_moments(self) -> None:
        targets, weights = get_target_set("candidate_replacement_bequest_internal_v1")
        self.assertEqual(len(targets), 14)
        self.assertEqual(set(targets), set(weights))
        self.assertNotIn("old_nonhousing_wealth_to_income_median_6575", targets)
        self.assertNotIn("old_parent_childless_nonhousing_wealth_to_income_gap_6575", targets)
        self.assertNotIn("old_age_own_rate", targets)
        self.assertIn("old_total_estate_wealth_to_annual_income_median_7684", targets)
        self.assertIn("old_total_estate_wealth_to_annual_income_p90_p50_7684", targets)
        self.assertIn(
            "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575",
            targets,
        )

    def test_median_composition_target_set_keeps_estate_median_only(self) -> None:
        targets, weights = get_target_set("candidate_replacement_bequest_median_composition_v1")
        self.assertEqual(len(targets), 13)
        self.assertEqual(set(targets), set(weights))
        self.assertIn("old_nonhousing_wealth_to_income_median_6575", targets)
        self.assertIn("old_total_estate_wealth_to_annual_income_median_7684", targets)
        self.assertNotIn("old_total_estate_wealth_to_annual_income_p90_p50_7684", targets)
        self.assertNotIn(
            "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575",
            targets,
        )
        # Estate median target/weight are unchanged; the composition median
        # uses the fresh reference-person bootstrap target and inverse variance.
        self.assertEqual(targets["old_total_estate_wealth_to_annual_income_median_7684"], 6.50131577436537)
        self.assertEqual(weights["old_total_estate_wealth_to_annual_income_median_7684"], 18.585767349158665)
        self.assertEqual(targets["old_nonhousing_wealth_to_income_median_6575"], 1.90821154211154)
        self.assertEqual(weights["old_nonhousing_wealth_to_income_median_6575"], 83.74916751466371)


class StandardBequestArmContractTests(unittest.TestCase):
    def _args(self, arm: str) -> SimpleNamespace:
        return SimpleNamespace(
            arm=arm,
            ltv_terminal=0.4,
            theta1=0.25,
            seed_theta0=0.30,
            fixed_theta0=None,
        )

    def test_m4_estimates_both_remaining_de_nardi_parameters(self) -> None:
        from tools.run_intergen_bequest_exit_chain import arm_contract

        active, fixed, mechanism = arm_contract(self._args("M4"))
        self.assertEqual(len(active), 13)
        self.assertEqual(
            active[-2:],
            [
                ("theta0", 0.0, 8.0, "softzero"),
                ("theta1", 0.02, 16.0, "log"),
            ],
        )
        self.assertEqual(fixed, {"theta_n": 0.0, "tenure_choice_kappa": 0.0})
        self.assertEqual(mechanism["bequest_spec"], "linear_child_scale")
        self.assertTrue(mechanism["normalize_bequest_utility"])
        self.assertFalse(mechanism["owner_ltv_taper"])
        self.assertTrue(mechanism["use_age_survival"])
        active_m2, fixed_m2, mechanism_m2 = arm_contract(self._args("M2"))
        self.assertNotEqual(active, active_m2)
        self.assertNotEqual(fixed, fixed_m2)
        self.assertEqual(mechanism, mechanism_m2)

    def test_m4_target_system_has_14_moments_including_rooms(self) -> None:
        from tools.run_intergen_bequest_exit_chain import target_system

        targets, weights = target_system("candidate_replacement_bequest_median_composition_v1")
        self.assertEqual(len(targets), 14)
        self.assertEqual(set(targets), set(weights))
        self.assertIn("aggregate_mean_occupied_rooms_18_85", targets)
        # The other arms keep the audited 15-moment system.
        targets_m3, _ = target_system("candidate_replacement_bequest_internal_v1")
        self.assertEqual(len(targets_m3), 15)

    def test_m4_uses_ssa_survival_schedule(self) -> None:
        from tools.run_intergen_bequest_exit_chain import survival_schedule

        schedule = survival_schedule(SimpleNamespace(arm="M4", J=17))
        self.assertEqual(int(np.count_nonzero(schedule < 1.0)), 4)

    def test_m4_rejects_silent_start_clipping(self) -> None:
        from tools.run_intergen_bequest_exit_chain import arm_contract

        args = self._args("M4")
        args.theta1 = 16.01
        with self.assertRaisesRegex(ValueError, "theta1 start"):
            arm_contract(args)

    def test_m4_profile_fixes_theta1_and_reestimates_other_coordinates(self) -> None:
        from tools.run_intergen_bequest_exit_chain import arm_contract, survival_schedule

        args = self._args("M4_PROFILE")
        args.theta1 = 0.75
        active, fixed, mechanism = arm_contract(args)
        self.assertEqual(len(active), 12)
        self.assertEqual(active[-1], ("theta0", 0.0, 8.0, "softzero"))
        self.assertEqual(
            fixed,
            {"theta_n": 0.0, "tenure_choice_kappa": 0.0, "theta1": 0.75},
        )
        self.assertEqual(mechanism["bequest_spec"], "linear_child_scale")
        self.assertTrue(mechanism["normalize_bequest_utility"])
        self.assertFalse(mechanism["owner_ltv_taper"])
        self.assertTrue(mechanism["use_age_survival"])
        schedule = survival_schedule(SimpleNamespace(arm="M4_PROFILE", J=17))
        self.assertEqual(int(np.count_nonzero(schedule < 1.0)), 4)


class StandardBequestCollectorContractTests(unittest.TestCase):
    def test_tight_repeat_must_be_bit_identical(self) -> None:
        from tools.collect_intergen_internal_bequest_recalibration import eligible_tight

        summary = {
            "best_tight": {"strict_converged": True},
            "tight_repeat_check": {
                "both_strict": True,
                "loss_abs_difference": 0.0,
                "max_abs_moment_difference": 0.0,
            },
        }
        self.assertIsNotNone(eligible_tight(summary))
        summary["tight_repeat_check"]["loss_abs_difference"] = 1e-15
        self.assertIsNone(eligible_tight(summary))

    def test_nested_reference_proves_exact_strict_m1_seed(self) -> None:
        from tools.collect_intergen_internal_bequest_recalibration import (
            ESTATE_MEDIAN,
            ESTATE_MEDIAN_TARGET,
            ESTATE_MEDIAN_WEIGHT,
            EXPECTED_M4_ACTIVE_DOMAIN,
            EXPECTED_M4_FIXED,
            EXPECTED_M4_MECHANISM,
            EXPECTED_TIGHT_EVALUATOR,
            M4_SHARED_THETA,
            NONHOUSING_MEDIAN,
            NONHOUSING_MEDIAN_TARGET,
            NONHOUSING_MEDIAN_WEIGHT,
            validate_m4_chain_metadata,
            validate_m4_nested_reference,
        )

        theta = {name: float(ii + 1) for ii, name in enumerate(M4_SHARED_THETA)}
        theta.update(theta0=0.0, theta1=0.25, theta_n=0.0, tenure_choice_kappa=0.0)
        with tempfile.TemporaryDirectory() as tmp:
            seed_path = Path(tmp) / "m1_results.json"
            seed_path.write_text(
                json.dumps(
                    {
                        "winners": {
                            "M1": {
                                "strict_converged": True,
                                "theta": theta,
                            }
                        }
                    }
                )
            )
            summary = {
                "metadata": {
                    "status": "proper_joint_smm_chain",
                    "arm": "M4",
                    "target_set": "candidate_replacement_bequest_median_composition_v1",
                    "free_parameter_count": 13,
                    "target_count": 14,
                    "max_evals": 1,
                    "start_mix": 0.0,
                    "seed_arm": "M1",
                    "seed_record": str(seed_path),
                    "J": 17,
                    "Nb": 120,
                    "max_iter_eq": 10,
                    "tol_eq": 1e-4,
                    "tight_winner_evaluator": EXPECTED_TIGHT_EVALUATOR,
                    "active_domain": EXPECTED_M4_ACTIVE_DOMAIN,
                    "fixed_parameters": EXPECTED_M4_FIXED,
                    "mechanism": EXPECTED_M4_MECHANISM,
                    "targets": {
                        ESTATE_MEDIAN: ESTATE_MEDIAN_TARGET,
                        NONHOUSING_MEDIAN: NONHOUSING_MEDIAN_TARGET,
                    },
                    "weights": {
                        ESTATE_MEDIAN: ESTATE_MEDIAN_WEIGHT,
                        NONHOUSING_MEDIAN: NONHOUSING_MEDIAN_WEIGHT,
                    },
                },
                "best_tight": {
                    "strict_converged": True,
                    "theta": theta,
                    "rank_loss": 1.0,
                },
                "tight_repeat_check": {
                    "both_strict": True,
                    "loss_abs_difference": 0.0,
                    "max_abs_moment_difference": 0.0,
                },
            }
            validate_m4_chain_metadata(summary["metadata"])
            bad_metadata = dict(summary["metadata"])
            bad_metadata["active_domain"] = []
            with self.assertRaisesRegex(RuntimeError, "exact production contract"):
                validate_m4_chain_metadata(bad_metadata)
            self.assertIs(
                validate_m4_nested_reference(
                    summary,
                    "candidate_replacement_bequest_median_composition_v1",
                ),
                summary["best_tight"],
            )
            # A log-domain encode/decode maps 0.25 to the immediately adjacent
            # representable float. The nested validator should admit only this
            # machine-roundoff discrepancy, not an economically distinct theta.
            summary["best_tight"]["theta"] = dict(theta)
            summary["best_tight"]["theta"]["theta1"] = math.nextafter(0.25, math.inf)
            self.assertIs(
                validate_m4_nested_reference(
                    summary,
                    "candidate_replacement_bequest_median_composition_v1",
                ),
                summary["best_tight"],
            )
            summary["best_tight"]["theta"]["theta1"] = 0.25 + 1e-12
            with self.assertRaisesRegex(RuntimeError, "exact strict M1 theta"):
                validate_m4_nested_reference(
                    summary,
                    "candidate_replacement_bequest_median_composition_v1",
                )
            summary["best_tight"]["theta"] = theta
            summary["metadata"]["free_parameter_count"] = 12
            with self.assertRaisesRegex(RuntimeError, "metadata"):
                validate_m4_nested_reference(
                    summary,
                    "candidate_replacement_bequest_median_composition_v1",
                )


if __name__ == "__main__":
    unittest.main()
