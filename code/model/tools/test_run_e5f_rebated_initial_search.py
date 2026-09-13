from __future__ import annotations

import threading
import time
import unittest
from pathlib import Path

import run_e5f_rebated_initial_search as s


class SearchTests(unittest.TestCase):
    def test_failure_classification_keeps_numerical_rejections_distinct(self):
        self.assertEqual(s.failure_status("housing equilibrium did not converge"),
                         "rejected_equilibrium")
        self.assertEqual(s.failure_status("stationary mass feasibility gate"),
                         "rejected_mass_gate")
        self.assertEqual(s.failure_status("source fingerprint changed"), "failed")
        self.assertEqual(s.failure_status("joint root evaluation budget exhausted"),
                         "rejected_numerical")

    def test_six_worker_batch_can_process_eighteen_cases(self):
        lock = threading.Lock(); live = 0; maximum = 0; completed = []
        items = [dict(case_id=str(i)) for i in range(18)]
        def worker(item):
            nonlocal live, maximum
            with lock: live += 1; maximum = max(maximum, live)
            time.sleep(0.002)
            with lock: live -= 1
            return dict(item, status="verified", loss=float(item["case_id"]))
        results, hard = s.bounded_batch(items, worker, 6, completed.append)
        self.assertEqual(len(results), 18); self.assertFalse(hard)
        self.assertLessEqual(maximum, 6)

    def test_hard_error_stops_new_submissions_for_this_branch(self):
        calls = []
        items = [dict(case_id=str(i)) for i in range(18)]
        def worker(item):
            calls.append(item["case_id"])
            return dict(item, status="failed" if item["case_id"] == "0" else "verified")
        _, hard = s.bounded_batch(items, worker, 2, lambda _: None)
        self.assertTrue(hard)
        self.assertLess(len(calls), len(items))

    def test_one_candidate_numerical_failure_does_not_stop_search(self):
        calls = []; items = [dict(case_id=str(i)) for i in range(8)]
        def worker(item):
            calls.append(item["case_id"])
            return dict(item, status="rejected_numerical", failure_detail="root failed") \
                if item["case_id"] == "0" else dict(item, status="verified")
        results, hard = s.bounded_batch(items, worker, 2, lambda _: None)
        self.assertFalse(hard); self.assertEqual(len(results), len(items))

    def test_three_matching_candidate_failures_stop_as_systemic(self):
        calls = []; items = [dict(case_id=str(i)) for i in range(12)]
        def worker(item):
            calls.append(item["case_id"])
            return dict(item, status="rejected_numerical",
                        failure_detail="same root failure at evaluation 12")
        _, hard = s.bounded_batch(items, worker, 1, lambda _: None)
        self.assertTrue(hard); self.assertEqual(len(calls), 3)

    def test_budget_constants_match_contract(self):
        self.assertEqual(s.DEFAULT_WORKERS, 6)
        self.assertEqual(s.MAXIMUM_WORKERS, 18)
        self.assertEqual(s.MAXIMUM_COORDINATES, 18)
        self.assertEqual(s.MAXIMUM_JOINT, 18)
        self.assertEqual(s.MAXIMUM_WALL_SECONDS, 10800)
        self.assertEqual(s.FINAL_RESERVE_SECONDS, 2100)

    def test_candidate_command_uses_joint_scored_invoker_and_all_pins(self):
        command = s.scored_command(Path("scored.py"), Path("helper.py"), "h",
            Path("joint.py"), "j", Path("template"), Path("proposal.json"), Path("out"))
        self.assertEqual(command[2], "scored.py")
        self.assertEqual(command[command.index("--helper-sha256") + 1], "h")
        self.assertEqual(command[command.index("--joint-sha256") + 1], "j")
        self.assertEqual(command[command.index("--proposal") + 1], "proposal.json")

    def test_unapproved_source_pin_fails_before_candidate(self):
        with self.assertRaisesRegex(ValueError, "helper"):
            s.require_source_pins(Path("missing-helper"), "wrong",
                Path("missing-joint"), s.JOINT_SHA256,
                Path("missing-scored"), s.SCORED_SHA256)


if __name__ == "__main__": unittest.main()
