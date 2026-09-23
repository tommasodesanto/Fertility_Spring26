"""Contract tests for the utility collector's critical selection and audit rules."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).with_name("collect_e5f_utility_overnight.py")
SPEC = importlib.util.spec_from_file_location("collect_e5f_utility_overnight", MODULE)
collector = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(collector)


class UtilityCollectorTests(unittest.TestCase):
    def test_timeout_is_incomplete_and_scientific_rejection_is_distinct(self):
        with tempfile.TemporaryDirectory() as temp:
            case = Path(temp)
            evaluation = case / "result/evaluation"
            (evaluation / "raw").mkdir(parents=True)
            (case / "status.json").write_text(json.dumps({"status": "started"}))
            running_status, reason = collector._attempt_status(case, evaluation, False)
            self.assertEqual(running_status, "running")
            self.assertIn("liveness is unknown", reason)
            failure = evaluation / "raw/failure.json"
            failure.write_text(json.dumps({"error_type": "TimeoutExpired", "error": "native timeout"}))
            self.assertEqual(collector._attempt_status(case, evaluation, False)[0], "incomplete")
            failure.write_text(json.dumps({"error_type": "InfeasibleThetaError", "error": "dead mass"}))
            self.assertEqual(collector._attempt_status(case, evaluation, False)[0], "rejected")
            failure.write_text(json.dumps({"error_type": "KeyError", "error": "missing contract"}))
            self.assertEqual(collector._attempt_status(case, evaluation, False)[0], "failed")

    def test_full_parameter_dimensions_are_cell_specific(self):
        targets = [{}] * 13
        collector.validate_row_counts(targets, [{}] * 17, "B_floor")
        collector.validate_row_counts(targets, [{}] * 19, "D_shares")
        with self.assertRaises(collector.CollectionError):
            collector.validate_row_counts(targets, [{}] * 17, "B_shares")
        with self.assertRaises(collector.CollectionError):
            collector.validate_row_counts(targets[:-1], [{}] * 19, "D_shares")

    def test_mixed_target_or_weight_fingerprint_fails_closed(self):
        collector.assert_one_target_system(["same", "same", "same", "same"], "same")
        with self.assertRaises(collector.CollectionError):
            collector.assert_one_target_system(["same", "changed"], "same")

    def test_science_rows_reconcile_twelve_weighted_gaps_and_separate_normalization(self):
        targets, objective = [], []
        for index in range(12):
            restriction = f"moment_{index}"
            target, model, weight = 0.2 + index, 0.3 + index, 2.0 + index
            targets.append({"restriction_id": restriction, "target": target, "model": model,
                            "gap": model-target, "actual_weight": weight,
                            "loss_contribution": weight*(model-target)**2,
                            "role": "proposed_scored_restriction", "scored": True})
            objective.append({"restriction_id": restriction, "target": target,
                              "actual_weight": weight, "role": "proposed_scored_restriction"})
        targets.append({"restriction_id": "initial_normalization", "target": 2.1,
                        "model": 2.1, "gap": 0.0, "actual_weight": None,
                        "loss_contribution": None, "role": "normalization_separate_from_scored_objective",
                        "scored": False})
        objective.append({"restriction_id": "initial_normalization", "target": 2.1,
                          "actual_weight": None, "role": "normalization_separate_from_scored_objective"})
        loss = sum(row["loss_contribution"] for row in targets if row["scored"])
        collector.validate_science_rows(targets, objective, loss)
        targets[0]["loss_contribution"] += 0.01
        with self.assertRaises(collector.CollectionError):
            collector.validate_science_rows(targets, objective, loss)

    def test_smoke_can_win_but_verification_never_enters_selection_pool(self):
        verified = {
            ("smoke", "B_floor", "smoke_01"): {"loss": 5.0},
            ("production", "B_floor", "worker01_proposal01"): {"loss": 6.0},
            ("verification", "B_floor", "selected_repeat_01"): {"loss": 1.0},
        }
        candidates = [row for row in collector.selection_candidates([], verified)
                      if row[1] == "B_floor"]
        self.assertEqual({row[0] for row in candidates}, {"smoke", "production"})
        self.assertEqual(min(row[3]["loss"] for row in candidates), 5.0)

    def test_case_counts_keep_unrun_and_all_failure_classes_visible(self):
        inventory = [
            {"cell": "D_shares", "status": "unrun"},
            {"cell": "D_shares", "status": "verified_scored"},
            {"cell": "D_shares", "status": "rejected"},
            {"cell": "D_shares", "status": "failed"},
            {"cell": "D_shares", "status": "incomplete"},
            {"cell": "D_shares", "status": "collection_rejected"},
        ]
        counts = collector.count_statuses(inventory, "D_shares")
        self.assertEqual(counts["unrun"], 1)
        self.assertEqual(counts["verified_scored"], 1)
        self.assertEqual(counts["rejected"], 1)
        self.assertEqual(counts["failed"], 1)
        self.assertEqual(counts["incomplete"], 1)
        self.assertEqual(counts["collection_rejected"], 1)

    def test_corrupt_selected_candidate_is_removed_before_score_table_export(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            bad_checkpoint, good_checkpoint = root / "bad.pkl", root / "good.pkl"
            bad_checkpoint.write_bytes(b"corrupted")
            good_checkpoint.write_bytes(b"verified")
            bad = {"loss": 1.0, "checkpoint_path": bad_checkpoint,
                   "checkpoint_sha256": "0" * 64,
                   "score": {"target_fit": [{"restriction_id": "bad"}], "parameters": []},
                   "plan": {"parameter_bounds": {}}}
            import hashlib
            good_hash = hashlib.sha256(good_checkpoint.read_bytes()).hexdigest()
            good = {"loss": 2.0, "checkpoint_path": good_checkpoint,
                    "checkpoint_sha256": good_hash,
                    "score": {"target_fit": [{"restriction_id": "good"}], "parameters": []},
                    "plan": {"parameter_bounds": {}}}
            verified = {("production", "B_floor", "bad_case"): bad,
                        ("smoke", "B_floor", "good_case"): good}
            inventory = [
                {"stage": "production", "cell": "B_floor", "case_id": "bad_case", "status": "verified_scored"},
                {"stage": "smoke", "cell": "B_floor", "case_id": "good_case", "status": "verified_scored"},
            ]
            candidates = [("production", "B_floor", "bad_case", bad),
                          ("smoke", "B_floor", "good_case", good)]
            selected = collector.selected_checkpoint_candidates(candidates, inventory, verified, "B_floor")
            targets, _ = collector.score_table_rows(verified)
            self.assertEqual([item[1] for item in selected], ["good_case"])
            self.assertNotIn(("production", "B_floor", "bad_case"), verified)
            self.assertEqual({row["restriction_id"] for row in targets}, {"good"})
            self.assertEqual(inventory[0]["status"], "collection_rejected")

    def test_scientific_receipt_mismatch_requires_review(self):
        receipt = {"status": "verified_exact_twice", "selected": {
            "case_id": "old_case", "checkpoint_sha256": "a" * 64}}
        selected = {"case_id": "selected_case", "checkpoint_sha256": "b" * 64}
        self.assertEqual(
            collector.scientific_receipt_status(receipt, selected, "a" * 64),
            "not_applicable_selected_binding_mismatch_requires_review")

    def test_raw_completed_repeat_count_is_not_verified_repeat_count(self):
        repeats = [{"status": "completed"}, {"status": "completed"}]
        verified = {("verification", "B_floor", "selected_repeat_01"): {"score": {}}}
        self.assertEqual(collector.verification_counts(repeats, verified, "B_floor"), (2, 1))


if __name__ == "__main__":
    unittest.main()
