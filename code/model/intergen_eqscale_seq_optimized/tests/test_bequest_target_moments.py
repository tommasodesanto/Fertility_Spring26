from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from intergen_eqscale_seq_optimized.calibration import (
    PSID_ENTRY_WEALTH_RATIO_NODES_1824,
    PSID_ENTRY_WEALTH_RATIO_NODES_2535,
    PSID_ENTRY_WEALTH_RATIO_WEIGHTS_1824,
    get_target_set,
)
from intergen_eqscale_seq_optimized.solver import (
    add_aggregate_wealth_gross_labor_diagnostics,
    add_annual_gross_estate_wealth_moments,
    add_annual_gross_old_wealth_moments,
    add_old_nonhousing_income_share_moments,
)


class BequestTargetMomentTests(unittest.TestCase):
    def test_aggregate_gross_ratio_uses_all_wealth_and_worker_earnings(self) -> None:
        p = SimpleNamespace(
            J=2,
            I=1,
            J_R=1,
            age_start=62.0,
            da=4.0,
            period_years=4.0,
            tau_pay=0.2,
            income=np.array([[4.0, 2.0]]),
            n_parity=1,
            n_child_states=1,
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
        )
        bg = np.array([-1.0, 3.0])
        ph = np.array([1.0])
        g = np.zeros((2, 2, 1, 2, 1, 1, 1))
        g[0, 0, 0, 0, 0, 0, 0] = 1.0
        g[1, 1, 0, 1, 0, 0, 0] = 1.0
        stats = SimpleNamespace()
        add_aggregate_wealth_gross_labor_diagnostics(stats, g, p, bg, ph)
        self.assertAlmostEqual(stats.aggregate_wealth, 4.0)
        self.assertAlmostEqual(stats.aggregate_annual_gross_labor_earnings, 1.25)
        self.assertAlmostEqual(
            stats.aggregate_wealth_to_annual_gross_labor_earnings,
            3.2,
        )

    def test_age_profile_prorates_four_year_states_at_bin_boundaries(self) -> None:
        p = SimpleNamespace(
            J=4,
            I=1,
            J_R=4,
            age_start=26.0,
            da=4.0,
            period_years=4.0,
            tau_pay=0.2,
            income=np.full((1, 4), 4.0),
            n_parity=1,
            n_child_states=1,
            H_own=np.array([2.0]),
            z_grid=np.array([1.0]),
        )
        bg = np.array([0.0, 10.0, 20.0, 30.0])
        g = np.zeros((4, 2, 1, 4, 1, 1, 1))
        for j in range(4):
            g[j, 0, 0, j, 0, 0, 0] = 1.0
        stats = SimpleNamespace()
        add_aggregate_wealth_gross_labor_diagnostics(
            stats,
            g,
            p,
            bg,
            np.array([1.0]),
        )
        # The 34--37 state contributes one half to the 26--35 bin.
        self.assertAlmostEqual(
            stats.aggregate_wealth_to_annual_gross_labor_earnings_26_35,
            (0.0 + 10.0 + 0.5 * 20.0) / (2.5 * 1.25),
        )

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
        self.assertEqual(stats.old_total_wealth_to_annual_income_median_7684, 10.0)
        self.assertEqual(stats.old_total_wealth_to_annual_income_p90_7684, 30.0)
        self.assertEqual(stats.old_total_wealth_to_annual_income_p90_p50_7684, 3.0)
        self.assertEqual(
            stats.old_total_wealth_to_annual_income_median_7684,
            stats.old_total_estate_wealth_to_annual_income_median_7684,
        )
        self.assertEqual(
            stats.old_total_wealth_to_annual_income_p90_7684,
            stats.old_total_estate_wealth_to_annual_income_p90_7684,
        )
        self.assertEqual(
            stats.old_total_wealth_to_annual_income_p90_p50_7684,
            stats.old_total_estate_wealth_to_annual_income_p90_p50_7684,
        )
        self.assertEqual(
            stats.old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575,
            20.0,
        )

    def test_old_wealth_moment_name_is_the_live_implementation(self) -> None:
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
        g[0, 0, 0, 3, 0, 1, 0] = 1.0
        g[0, 1, 0, 3, 0, 2, 0] = 1.0
        g[1, 1, 0, 4, 0, 1, 0] = 1.0
        g[1, 0, 0, 4, 0, 2, 0] = 1.0

        stats = SimpleNamespace()
        add_annual_gross_old_wealth_moments(stats, g, p, bg, ph)

        self.assertEqual(stats.old_total_wealth_to_annual_income_median_7684, 10.0)
        self.assertEqual(stats.old_total_wealth_to_annual_income_p90_7684, 30.0)
        self.assertEqual(stats.old_total_wealth_to_annual_income_p90_p50_7684, 3.0)

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

    def test_income_disciplined_target_set_swaps_composition_median_for_share(self) -> None:
        targets, weights = get_target_set("candidate_replacement_income_disciplined_v1")
        self.assertEqual(len(targets), 15)
        self.assertEqual(targets["aggregate_mean_occupied_rooms_18_85"], 5.779970481941968)
        self.assertEqual(weights["aggregate_mean_occupied_rooms_18_85"], 6.0)
        self.assertEqual(set(targets), set(weights))
        self.assertIn("old_nonhousing_ge_1x_income_share_6575", targets)
        self.assertIn("old_age_own_rate", targets)
        self.assertIn("old_total_estate_wealth_to_annual_income_median_7684", targets)
        self.assertNotIn("old_nonhousing_wealth_to_income_median_6575", targets)
        self.assertNotIn("old_nonhousing_wealth_to_income_6575", targets)
        self.assertNotIn("old_total_estate_wealth_to_annual_income_p90_p50_7684", targets)
        self.assertNotIn(
            "old_2plus_minus_1_total_estate_wealth_to_annual_income_median_gap_6575",
            targets,
        )
        self.assertNotIn("old_parent_childless_nonhousing_wealth_to_income_gap_6575", targets)
        # Share target/weight come from the block3 person bootstrap; the
        # old-age ownership rate keeps its legacy ACS value and weight.
        self.assertEqual(targets["old_nonhousing_ge_1x_income_share_6575"], 0.608333139649131)
        self.assertEqual(weights["old_nonhousing_ge_1x_income_share_6575"], 9435.18732291246)
        self.assertEqual(targets["old_age_own_rate"], 0.76426097)
        self.assertEqual(weights["old_age_own_rate"], 160.0)
        self.assertEqual(targets["old_total_estate_wealth_to_annual_income_median_7684"], 6.50131577436537)
        self.assertEqual(weights["old_total_estate_wealth_to_annual_income_median_7684"], 18.585767349158665)

    def test_old_nonhousing_ge_1x_income_share_uses_raw_b_and_estate_window(self) -> None:
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
        bg = np.array([-5.0, 0.5, 1.0, 10.0])
        ph = np.array([10.0])
        g = np.zeros((4, 2, 1, 5, 1, 3, 1))

        # Ages 66-74 lie inside the 65-75 window; income is 1 everywhere.
        g[0, 0, 0, 0, 0, 0, 0] = 1.0  # b=-5: raw negative b stays in the denominator
        g[2, 0, 0, 0, 0, 1, 0] = 1.0  # b=1: ratio exactly 1 counts (>= is inclusive)
        g[1, 1, 0, 1, 0, 2, 0] = 2.0  # owner b=0.5: housing is excluded, does not count
        g[3, 0, 0, 2, 0, 1, 0] = 1.0  # b=10 counts
        g[3, 0, 0, 3, 0, 1, 0] = 5.0  # age 78 lies outside the window

        stats = SimpleNamespace()
        add_old_nonhousing_income_share_moments(stats, g, p, bg)
        self.assertAlmostEqual(stats.old_nonhousing_ge_1x_income_share_6575, 2.0 / 5.0)


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

    def test_m5_frees_tenure_choice_kappa_with_wide_bequest_domains(self) -> None:
        from tools.run_intergen_bequest_exit_chain import arm_contract

        args = self._args("M5")
        args.seed_kappa = 0.01
        active, fixed, mechanism = arm_contract(args)
        self.assertEqual(len(active), 14)
        self.assertEqual(
            active[-3:],
            [
                ("theta0", 0.0, 8.0, "softzero"),
                ("theta1", 0.02, 16.0, "log"),
                ("tenure_choice_kappa", 0.0, 0.12, "softzero"),
            ],
        )
        self.assertEqual(fixed, {"theta_n": 0.0})
        args_m4 = self._args("M4")
        _, _, mechanism_m4 = arm_contract(args_m4)
        self.assertTrue(mechanism["entry_wealth_censor_to_frontier"])
        self.assertFalse(mechanism_m4["entry_wealth_censor_to_frontier"])
        self.assertEqual(
            {k: v for k, v in mechanism.items() if k != "entry_wealth_censor_to_frontier"},
            {k: v for k, v in mechanism_m4.items() if k != "entry_wealth_censor_to_frontier"},
        )
        self.assertTrue(mechanism["use_age_survival"])

    def test_m5_uses_ssa_survival_schedule(self) -> None:
        from tools.run_intergen_bequest_exit_chain import survival_schedule

        schedule = survival_schedule(SimpleNamespace(arm="M5", J=17))
        self.assertEqual(int(np.count_nonzero(schedule < 1.0)), 4)

    def test_m5_rejects_silent_start_clipping(self) -> None:
        from tools.run_intergen_bequest_exit_chain import arm_contract

        args = self._args("M5")
        args.seed_kappa = 0.13
        with self.assertRaisesRegex(ValueError, "tenure_choice_kappa start"):
            arm_contract(args)
        args.seed_kappa = 0.0
        args.theta1 = 16.01
        with self.assertRaisesRegex(ValueError, "theta1 start"):
            arm_contract(args)

    def test_m5_target_system_has_15_moments_including_rooms(self) -> None:
        from tools.run_intergen_bequest_exit_chain import target_system

        targets, weights = target_system("candidate_replacement_income_disciplined_v1")
        self.assertEqual(len(targets), 15)
        self.assertEqual(set(targets), set(weights))
        self.assertIn("aggregate_mean_occupied_rooms_18_85", targets)
        self.assertIn("old_nonhousing_ge_1x_income_share_6575", targets)

    def test_m5_income_and_entry_overrides(self) -> None:
        from intergen_eqscale_seq_optimized.local_panel import income_process_overrides
        from tools.run_intergen_bequest_exit_chain import (
            SS_ANNUAL_INNOVATION_SD,
            SS_ANNUAL_RHO,
            arm_contract,
            common_overrides,
        )

        self.assertEqual(SS_ANNUAL_RHO, 0.9601845894041878)
        self.assertEqual(SS_ANNUAL_INNOVATION_SD, 0.06453733259357768)
        args = self._args("M5")
        args.seed_kappa = 0.0
        args.J = 17
        args.Nb = 60
        args.max_iter_eq = 2
        args.tol_eq = 0.25
        _, _, mechanism = arm_contract(args)
        overrides = common_overrides(args, mechanism)
        expected_income = income_process_overrides(
            5, "rouwenhorst", SS_ANNUAL_INNOVATION_SD, SS_ANNUAL_RHO
        )
        np.testing.assert_allclose(overrides["z_grid"], expected_income["z_grid"])
        np.testing.assert_allclose(overrides["Pi_z"], expected_income["Pi_z"])
        self.assertEqual(
            overrides["income_shock_persistence"], expected_income["income_shock_persistence"]
        )
        np.testing.assert_array_equal(
            overrides["entry_wealth_ratio_nodes"], PSID_ENTRY_WEALTH_RATIO_NODES_1824
        )
        np.testing.assert_array_equal(
            overrides["entry_wealth_ratio_weights"], PSID_ENTRY_WEALTH_RATIO_WEIGHTS_1824
        )
        args_m4 = self._args("M4")
        args_m4.J = 17
        args_m4.Nb = 60
        args_m4.max_iter_eq = 2
        args_m4.tol_eq = 0.25
        _, _, mechanism_m4 = arm_contract(args_m4)
        overrides_m4 = common_overrides(args_m4, mechanism_m4)
        np.testing.assert_array_equal(
            overrides_m4["entry_wealth_ratio_nodes"], PSID_ENTRY_WEALTH_RATIO_NODES_2535
        )
        # July-16 feasibility probes forced M5 back onto the matched income
        # process (see the m5 contract); M5 and M4 share the z-grid until the
        # M6 forbearance/default margin unlocks higher income risk.
        self.assertTrue(
            np.allclose(overrides["z_grid"], overrides_m4["z_grid"])
        )

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

    def test_m5_nested_reference_and_acceptance_rows(self) -> None:
        from tools.collect_intergen_internal_bequest_recalibration import (
            ESTATE_MEDIAN,
            ESTATE_MEDIAN_TARGET,
            ESTATE_MEDIAN_WEIGHT,
            EXPECTED_M5_ACTIVE_DOMAIN,
            EXPECTED_M5_FIXED,
            EXPECTED_M5_MECHANISM,
            EXPECTED_TIGHT_EVALUATOR,
            M4_SHARED_THETA,
            NONHOUSING_SHARE,
            NONHOUSING_SHARE_TARGET,
            NONHOUSING_SHARE_WEIGHT,
            OLD_AGE_OWN_RATE,
            OLD_AGE_OWN_RATE_TARGET,
            OLD_AGE_OWN_RATE_WEIGHT,
            YOUNG_LIQUID,
            m5_acceptance_rows,
            validate_m5_chain_metadata,
            validate_m5_nested_reference,
        )

        self.assertEqual(len(EXPECTED_M5_ACTIVE_DOMAIN), 14)
        self.assertEqual(
            EXPECTED_M5_ACTIVE_DOMAIN[-1],
            {"name": "tenure_choice_kappa", "lower": 0.0, "upper": 0.12, "transform": "softzero"},
        )
        self.assertEqual(EXPECTED_M5_FIXED, {"theta_n": 0.0})

        theta = {name: float(ii + 1) for ii, name in enumerate(M4_SHARED_THETA)}
        theta.update(theta0=0.0, theta1=0.25, theta_n=0.0, tenure_choice_kappa=0.0)
        metadata = {
            "status": "proper_joint_smm_chain",
            "arm": "M5",
            "target_set": "candidate_replacement_income_disciplined_v1",
            "free_parameter_count": 14,
            "target_count": 15,
            "max_evals": 1,
            "start_mix": 0.0,
            "seed_arm": "M1",
            "J": 17,
            "Nb": 120,
            "max_iter_eq": 10,
            "tol_eq": 1e-4,
            "tight_winner_evaluator": EXPECTED_TIGHT_EVALUATOR,
            "active_domain": EXPECTED_M5_ACTIVE_DOMAIN,
            "fixed_parameters": EXPECTED_M5_FIXED,
            "mechanism": EXPECTED_M5_MECHANISM,
            "income_process": {
                "states": 5,
                "process": "rouwenhorst",
                "annual_rho": 0.9601845894041878,
                "annual_innovation_sd": 0.06453733259357768,
            },
            "entry_wealth_ages": "18_24",
            "targets": {
                ESTATE_MEDIAN: ESTATE_MEDIAN_TARGET,
                NONHOUSING_SHARE: NONHOUSING_SHARE_TARGET,
                OLD_AGE_OWN_RATE: OLD_AGE_OWN_RATE_TARGET,
            },
            "weights": {
                ESTATE_MEDIAN: ESTATE_MEDIAN_WEIGHT,
                NONHOUSING_SHARE: NONHOUSING_SHARE_WEIGHT,
                OLD_AGE_OWN_RATE: OLD_AGE_OWN_RATE_WEIGHT,
            },
        }
        validate_m5_chain_metadata(metadata)
        bad_metadata = dict(metadata)
        bad_metadata["entry_wealth_ages"] = "25_35"
        with self.assertRaisesRegex(RuntimeError, "exact production contract"):
            validate_m5_chain_metadata(bad_metadata)

        with tempfile.TemporaryDirectory() as tmp:
            seed_path = Path(tmp) / "m1_results.json"
            seed_path.write_text(
                json.dumps({"winners": {"M1": {"strict_converged": True, "theta": theta}}})
            )
            summary = {
                "metadata": {**metadata, "seed_record": str(seed_path)},
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
            self.assertIs(
                validate_m5_nested_reference(
                    summary, "candidate_replacement_income_disciplined_v1"
                ),
                summary["best_tight"],
            )
            summary["best_tight"]["theta"] = dict(theta)
            summary["best_tight"]["theta"]["tenure_choice_kappa"] = 0.01
            with self.assertRaisesRegex(RuntimeError, "tenure_choice_kappa=0"):
                validate_m5_nested_reference(
                    summary, "candidate_replacement_income_disciplined_v1"
                )

        winner = {
            "rank_loss": 3.0,
            "target_fit": [
                {"moment": ESTATE_MEDIAN, "gap": 0.1, "loss_contribution": 0.2},
                {"moment": NONHOUSING_SHARE, "gap": 0.01, "loss_contribution": 0.9},
                {"moment": OLD_AGE_OWN_RATE, "gap": -0.02, "loss_contribution": 0.06},
                {"moment": YOUNG_LIQUID, "gap": -0.05, "loss_contribution": 0.03},
                {"moment": "tfr", "gap": 0.1, "loss_contribution": 0.2},
            ],
        }
        rows, established_loss = m5_acceptance_rows(winner, 4.0, None)
        self.assertEqual(
            [row["criterion"] for row in rows],
            [
                "strict_exact_tight_repeat",
                "estate_median_absolute_gap",
                "nonhousing_ge_1x_share_absolute_gap",
                "established_12_moment_loss",
                "young_liquid_absolute_gap",
                "old_age_own_rate_absolute_gap",
                "free_winner_loss_minus_nested_zero",
                "identification_reported",
            ],
        )
        self.assertAlmostEqual(established_loss, 0.23)
        self.assertEqual(rows[-1]["pass"], "pending")
        self.assertTrue(all(row["pass"] is True for row in rows[:-1]))


if __name__ == "__main__":
    unittest.main()
