#!/usr/bin/env python3
"""Focused contract tests for the fixed-parameter slope-frontier orchestrator."""

from __future__ import annotations

import copy
import tempfile
import unittest
from pathlib import Path

import run_e5f_final_closed_root_slope_frontier as frontier


class FrontierContractTests(unittest.TestCase):
    def selected_candidate(self) -> dict:
        return {
            "theta": {
                "kappa_fert": 2.0,
                "kappa_fert_continuation": 1.5,
                "chi": 1.0,
                "H0": 10.0,
            },
            "old_psi_child": 0.2,
            "new_psi_child": -0.15,
        }

    def test_factor_grid_is_exact_and_ordered(self) -> None:
        self.assertEqual(frontier.FACTORS, (1.0, 0.5, 0.25, 0.1, 0.05))
        self.assertEqual(
            [frontier.factor_for_index(i) for i in range(1, 6)],
            list(frontier.FACTORS),
        )
        with self.assertRaises(frontier.ContractError):
            frontier.factor_for_index(0)
        with self.assertRaises(frontier.ContractError):
            frontier.factor_for_index(6)

    def test_center_changes_only_kappas_and_preserves_preference_change(self) -> None:
        selected = self.selected_candidate()
        original = copy.deepcopy(selected)
        center, scaled, change = frontier.build_factor_center(selected, 0.25)
        self.assertEqual(selected, original)
        self.assertAlmostEqual(scaled["kappa_fert"], 0.5)
        self.assertAlmostEqual(scaled["kappa_fert_continuation"], 0.375)
        self.assertAlmostEqual(scaled["chi"], selected["theta"]["chi"])
        self.assertAlmostEqual(scaled["H0"], selected["theta"]["H0"])
        self.assertAlmostEqual(change, -0.35)
        encoded = center["best_candidate"]
        self.assertAlmostEqual(
            encoded["new_psi_child"] - encoded["old_psi_child"], change
        )
        self.assertEqual(
            center["classification"], "fixed_parameter_diagnostic_not_calibration"
        )

    def test_renewal_contract_is_exactly_2p1_birth_accounting(self) -> None:
        contract = frontier.required_renewal_contract()
        frontier.validate_exact_renewal_contract(contract, "test")
        bad = copy.deepcopy(contract)
        bad["replacement_fertility"] = 2.12
        with self.assertRaises(frontier.ContractError):
            frontier.validate_exact_renewal_contract(bad, "test")
        bad = copy.deepcopy(contract)
        bad["birth_measure"] = "raw_birth_children"
        with self.assertRaises(frontier.ContractError):
            frontier.validate_exact_renewal_contract(bad, "test")

    def test_smallest_factor_stays_inside_live_kappa_bounds(self) -> None:
        selected = self.selected_candidate()
        selected["theta"]["kappa_fert"] = 0.3
        selected["theta"]["kappa_fert_continuation"] = 0.3
        with self.assertRaises(frontier.ContractError):
            frontier.build_factor_center(selected, 0.05)

    def test_diagnostic_bundle_hashes_every_executable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = [root / name for name in ("frontier.py", "continuation.py", "launch.sh")]
            for index, path in enumerate(paths):
                path.write_text(f"artifact-{index}\n", encoding="utf-8")
            contract = frontier.diagnostic_code_contract(*paths)
            self.assertEqual(
                set(contract),
                {
                    "frontier_driver_sha256",
                    "continuation_driver_sha256",
                    "launcher_sha256",
                },
            )
            original = frontier.canonical_json_sha256(contract)
            paths[1].write_text("changed\n", encoding="utf-8")
            changed = frontier.diagnostic_code_contract(*paths)
            self.assertNotEqual(original, frontier.canonical_json_sha256(changed))

    def test_smoke_receipt_requires_matching_full_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            hashes = {
                "ridge_report_sha256": "a" * 64,
                "frontier_driver_sha256": "b" * 64,
                "continuation_driver_sha256": "c" * 64,
                "launcher_sha256": "d" * 64,
                "diagnostic_bundle_sha256": "e" * 64,
            }
            summary = {
                "schema": frontier.SCHEMA,
                "status": "complete_fixed_parameter_slope_diagnostic",
                "run_stage": "smoke",
                "promotion_eligible": False,
                "factor_index": 2,
                "kappa_factor": 0.5,
                "target_count": frontier.TARGET_COUNT,
                "old_completed_fertility": 2.1,
                "old_queue_B_over_E": 1.0,
                "transition_loss": 50.0,
                "maximum_B_over_E": 0.9,
                "price_grid_contract": {
                    "minimum_ratio": frontier.PRICE_MIN_RATIO,
                    "maximum_ratio": frontier.PRICE_MAX_RATIO,
                    "declared_grid_points": frontier.PRICE_GRID_POINTS,
                },
                "input_provenance": hashes,
            }
            receipt = root / "summary.json"
            frontier.atomic_write_json(receipt, summary)
            receipt_hash = frontier.file_sha256(receipt)
            manifest = {
                "schema": "e5f_final_closed_root_slope_frontier_factor_manifest_v1",
                "status": "complete_fixed_parameter_slope_diagnostic_manifest",
                "run_stage": "smoke",
                "factor_index": 2,
                "kappa_factor": 0.5,
                "diagnostic_code_contract": {
                    name: hashes[name]
                    for name in (
                        "frontier_driver_sha256",
                        "continuation_driver_sha256",
                        "launcher_sha256",
                    )
                },
                "diagnostic_bundle_sha256": hashes["diagnostic_bundle_sha256"],
                "artifacts": {"summary.json": receipt_hash},
            }
            frontier.atomic_write_json(root / "manifest.json", manifest)
            validated = frontier.validate_smoke_receipt(receipt, receipt_hash, hashes)
            self.assertEqual(validated["smoke_receipt_sha256"], receipt_hash)
            summary["input_provenance"]["launcher_sha256"] = "f" * 64
            frontier.atomic_write_json(receipt, summary)
            bad_hash = frontier.file_sha256(receipt)
            with self.assertRaises(frontier.ContractError):
                frontier.validate_smoke_receipt(receipt, bad_hash, hashes)


if __name__ == "__main__":
    unittest.main()
