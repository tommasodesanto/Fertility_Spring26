#!/usr/bin/env python3
"""Focused synthetic subprocess tests. Run on Torch, never a native model solve."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import Mock

import run_e5f_utility_comparison_search as controller


class SchedulerTests(unittest.TestCase):
    def test_failed_smoke_cannot_trigger_automatic_selected_retries(self):
        instance=controller.Controller.__new__(controller.Controller)
        instance.best={"point":{"fixture":1.}}
        instance.clock={"start_epoch":1.}
        instance.smoke_verified=False
        instance.summary=Mock()
        instance.batch=Mock()
        instance.repeat_and_collect()
        instance.batch.assert_not_called()
        instance.summary.assert_called_once()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.started = []
        self.completed = []
        self.heartbeats = []

    def tearDown(self):
        self.temporary.cleanup()

    def candidate(self, name, code="pass", cap=3.):
        return dict(id=name, code=code, cap=cap)

    def launch(self, candidate, deadline):
        self.started.append(candidate["id"])
        return controller.ManagedProcess([sys.executable, "-c", candidate["code"]],
            self.root / (candidate["id"] + ".log"),
            min(deadline, time.time() + candidate["cap"]), os.environ.copy())

    def finish(self, candidate, process, code):
        row = dict(candidate_id=candidate["id"], loss=1., returncode=code,
                   status="timed_out" if process.timed_out else "success" if code == 0 else "failed")
        self.completed.append(row)
        return row

    def batch(self, candidates, workers=2, duration=5.):
        return controller.run_batch(candidates, workers=workers, deadline=time.time() + duration,
            launch=self.launch, finish=self.finish, heartbeat=lambda **row: self.heartbeats.append(row),
            poll_seconds=.01)

    def test_success_uses_real_subprocesses_and_completes_finite_list(self):
        result = self.batch([self.candidate(str(i)) for i in range(4)])
        self.assertTrue(result["complete"])
        self.assertIsNone(result["stop_reason"])
        self.assertEqual(result["unrun_ids"], [])
        self.assertEqual(len(self.completed), 4)
        self.assertLessEqual(max(row["active"] for row in self.heartbeats), 2)

    def test_error_stops_dispatch_but_preserves_running_sibling(self):
        candidates = [self.candidate("error", "raise SystemExit(9)"),
                      self.candidate("running", "import time; time.sleep(.2)"),
                      self.candidate("unrun")]
        result = self.batch(candidates)
        self.assertEqual(self.started, ["error", "running"])
        self.assertEqual(result["unrun_ids"], ["unrun"])
        self.assertEqual(result["stop_reason"], "case_failed_no_retry")
        self.assertEqual({row["candidate_id"]: row["status"] for row in result["results"]},
                         {"error": "failed", "running": "success"})

    def test_timeout_kills_owned_group_without_retry(self):
        before = time.monotonic()
        result = self.batch([self.candidate("timeout", "import time; time.sleep(30)", cap=.15),
                             self.candidate("unrun")], workers=1)
        self.assertLess(time.monotonic() - before, 3.)
        self.assertEqual(self.started, ["timeout"])
        self.assertEqual(result["results"][0]["status"], "timed_out")
        self.assertEqual(result["unrun_ids"], ["unrun"])

    def test_exhausted_budget_never_dispatches(self):
        result = self.batch([self.candidate("unrun")], duration=-1.)
        self.assertEqual(self.started, [])
        self.assertEqual(result["stop_reason"], "stage_deadline_exhausted")
        self.assertEqual(result["unrun_ids"], ["unrun"])

    def test_stage_cutoff_bounds_running_child(self):
        before = time.monotonic()
        result = self.batch([self.candidate("cutoff", "import time; time.sleep(30)")], duration=.15)
        self.assertLess(time.monotonic() - before, 3.)
        self.assertEqual(result["results"][0]["status"], "timed_out")

    def test_interruption_stops_new_dispatch_and_closes_own_children(self):
        before = time.monotonic()
        result = controller.run_batch([self.candidate("running", "import time; time.sleep(30)"),
            self.candidate("unrun")], workers=1, deadline=time.time() + 20,
            launch=self.launch, finish=self.finish, heartbeat=lambda **_: None,
            interrupted=lambda: time.monotonic() - before > .15, poll_seconds=.01)
        self.assertLess(time.monotonic() - before, 3.)
        self.assertEqual(result["stop_reason"], "controller_interrupted")
        self.assertEqual(result["unrun_ids"], ["unrun"])

    def test_process_group_timeout_does_not_kill_unrelated_process(self):
        unrelated = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"],
                                     start_new_session=True)
        try:
            self.batch([self.candidate("owned", "import time; time.sleep(30)", cap=.1)])
            self.assertIsNone(unrelated.poll())
        finally:
            unrelated.kill()
            unrelated.wait(timeout=3)

    def test_child_process_is_killed_with_its_owned_parent_group(self):
        marker = self.root / "descendant_completed"
        child = f"import pathlib,time; time.sleep(.7); pathlib.Path({str(marker)!r}).write_text('bad')"
        parent = f"import subprocess,sys,time; subprocess.Popen([sys.executable,'-c',{child!r}]); time.sleep(30)"
        self.batch([self.candidate("tree", parent, cap=.15)])
        time.sleep(.8)
        self.assertFalse(marker.exists())

    def test_missing_and_failed_de_results_cannot_become_parents(self):
        points = [dict(id="a"), dict(id="b")]
        a = dict(candidate_id="a", status="success", loss=2.)
        with self.assertRaisesRegex(RuntimeError, "incomplete"):
            controller.complete_scores(points, [a])
        with self.assertRaisesRegex(RuntimeError, "incomplete"):
            controller.complete_scores(points, [a, dict(candidate_id="b", status="failed", loss=None)])
        scores = controller.complete_scores(points, [a, dict(candidate_id="b", status="success", loss=3.)])
        self.assertEqual(scores, {"a": 2., "b": 3.})

    def test_wrong_generation_ids_are_not_treated_as_rejections(self):
        with self.assertRaisesRegex(RuntimeError, "IDs"):
            controller.complete_scores([dict(id="a")], [dict(candidate_id="x", status="success", loss=1.)])

    def test_only_exact_classification_in_search_consumes_inadmissible_slot(self):
        failure = dict(status="inadmissible_parameter_proposal", error_type="RuntimeError")
        self.assertTrue(controller.classified_inadmissible(failure, "initial"))
        self.assertTrue(controller.classified_inadmissible(failure, "de"))
        for stage in ("smoke", "repeat"):
            self.assertFalse(controller.classified_inadmissible(failure, stage))
        self.assertFalse(controller.classified_inadmissible({**failure, "status": "other"}, "de"))
        self.assertFalse(controller.classified_inadmissible({**failure, "error_type": "TimeoutError"}, "de"))

    def test_inadmissible_slot_is_consumed_once_then_next_distinct_point_runs(self):
        def classified_finish(candidate, process, code):
            row = self.finish(candidate, process, code)
            if candidate["id"] == "inadmissible":
                row.update(status="inadmissible", loss=None)
            return row
        result = controller.run_batch([self.candidate("inadmissible", "raise SystemExit(1)"),
            self.candidate("distinct_next")], workers=1, deadline=time.time() + 5,
            launch=self.launch, finish=classified_finish, heartbeat=lambda **_: None, poll_seconds=.01)
        self.assertEqual(self.started, ["inadmissible", "distinct_next"])
        self.assertTrue(result["complete"])
        self.assertIsNone(result["stop_reason"])

    def test_completed_barrier_stops_above_half_inadmissible(self):
        points = [dict(id=str(i)) for i in range(4)]
        rows = [dict(candidate_id=str(i), status="success", loss=float(i)) for i in range(4)]
        for i in (0, 1):
            rows[i].update(status="inadmissible", loss=None)
        self.assertEqual(controller.complete_scores(points, rows), {"0": None, "1": None, "2": 2., "3": 3.})
        rows[2].update(status="inadmissible", loss=None)
        with self.assertRaisesRegex(controller.BarrierStop, "more_than_half"):
            controller.complete_scores(points, rows)

    def test_exclusive_artifact_creation_does_not_overwrite_existing_work(self):
        path = self.root / "selected.json"
        controller.create(path, {"original": True})
        with self.assertRaises(FileExistsError):
            controller.create(path, {"original": False})
        self.assertEqual(json.loads(path.read_text()), {"original": True})
        self.assertFalse(list(self.root.glob("*.new.*")))


if __name__ == "__main__":
    unittest.main()
