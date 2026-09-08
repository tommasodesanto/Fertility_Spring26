#!/usr/bin/env python3
"""Orchestration regressions: temporary receipts and fake work, never model solves."""
from __future__ import annotations

import copy
import concurrent.futures
import json
import os
from pathlib import Path
import random
import tempfile
import threading
import unittest
from unittest.mock import patch

import run_e5f_simple_fertility_overnight as overnight


class OvernightTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, path, payload):
        overnight._write_json(path, payload)
        return path

    def controller(self):
        controller = object.__new__(overnight.Controller)
        controller.c = {"max_evaluations": 300, "case_timeout_seconds": 5400,
                        "seed": 20260908, "workers": 23}
        controller.root = self.root
        controller.start = overnight._now()
        controller.hard_deadline = controller.start + 43200
        controller.search_deadline = controller.start + 32400
        controller.active, controller.lock = set(), threading.Lock()
        controller.rows, controller.valid, controller.rejections = [], [], []
        controller.launched, controller.cancelled = 0, False
        controller.search_rejection_waves, controller.search_stopped = [], False
        return controller

    def failure(self, message, kind="RuntimeError", *, log=""):
        out = self.root / "failure"
        self.write(out / "adapter_failure.json", {"error": message, "type": kind})
        log_path = self.root / "case.log"
        log_path.write_text(log)
        return overnight.classify_failure(1, False, log_path, out)

    def test_declared_numerical_failures_are_rejections(self):
        examples = [
            ("Housing market did not clear: residual=0.01", "RuntimeError", "market_nonconvergence"),
            ("Old-steady-state fertility normalization is not bracketed: bounds", "RuntimeError", "old_fertility_no_bracket"),
            ("The dated first-birth branch has zero treated mass", "RuntimeError", "undefined_firstbirth_support"),
            ("bad theta", "InfeasibleThetaError", "infeasible_theta"),
        ]
        for message, kind, expected in examples:
            with self.subTest(message=message):
                self.assertEqual(self.failure(message, kind), expected)

    def test_unrelated_no_bracket_is_fatal(self):
        self.assertIsNone(self.failure("Unexpected policy interpolation: no bracket"))

    def test_logs_do_not_downgrade_unknown_exceptions(self):
        self.assertIsNone(self.failure("broken receipt", "KeyError",
                                       log="market gate failed; InfeasibleThetaError; no bracket"))

    def test_embedded_known_words_do_not_downgrade_unknown_exceptions(self):
        self.assertIsNone(self.failure("Missing file named Housing market did not clear.json", "FileNotFoundError"))
        self.assertIsNone(self.failure("Unexpected serializer field InfeasibleThetaError", "KeyError"))

    def test_missing_failure_receipt_is_fatal(self):
        log = self.root / "case.log"
        log.write_text("Housing market did not clear")
        self.assertIsNone(overnight.classify_failure(1, False, log, self.root))

    def test_declared_timeout_is_rejected_without_failure_receipt(self):
        self.assertEqual(overnight.classify_failure(124, True, self.root / "absent", self.root), "timeout")

    def test_reflection_is_finite_and_preserves_interior(self):
        for value, expected in [(-.2, .2), (0., 0.), (.3, .3), (1., 1.), (1.2, .8), (2.2, .2)]:
            self.assertAlmostEqual(overnight.repair_unit(value), expected)
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(ValueError):
                overnight.repair_unit(value)

    def test_inward_probe_moves_bounds_into_domain(self):
        self.assertEqual(overnight.inward_probe([0., 1.]), [.002, .998])

    def test_population_is_reproducible_and_covers_all_coordinates(self):
        anchor = [.5] * 11
        design = overnight.initial_population(anchor, 20260908)
        self.assertEqual(design, overnight.initial_population(anchor, 20260908))
        self.assertEqual(len(design), 23)
        self.assertEqual(design[0][1], anchor)
        self.assertEqual(len({tuple(u) for _, u in design}), 23)
        self.assertTrue(all(0 <= v <= 1 for _, u in design for v in u))
        self.assertTrue(all(any(u[j] != anchor[j] for _, u in design) for j in range(11)))

    def test_prepared_proposal_is_exact_driver_center(self):
        template = {"source_sha256": "source-pin", "target_fingerprint": "target-pin",
                    "code_bundle_sha256": "bundle-pin", "numerical_gates": {"market": .0002}}
        original = copy.deepcopy(template)
        center = {"proposal": {"unit_vector": [.2345678901234567] * 11}}
        path, sha = overnight._new_plan(template, self.root / "prepared", "test",
                                       [("one", center, "coordinate_polish")], 20260908)
        plan = overnight.adapter.read_json(path)
        case = plan["cases"][0]
        self.assertEqual((case["panel_task_id"], case["panel_size"], case["panel_design"]), (1, 1, "mixed"))
        self.assertEqual(overnight.adapter.read_json(path.parent / case["center"]), center)
        self.assertEqual(template, original)
        for key in original:
            self.assertEqual(plan[key], original[key])
        overnight.adapter.verify(path, sha)
        overnight.adapter.verify(path.parent / case["center"], case["center_sha256"])

    def test_final_probe_cannot_replace_frozen_best(self):
        controller = self.controller()
        selected = {"status": "valid", "loss": 26., "label": "selected"}
        controller._record(selected)
        frozen = (self.root / "best_so_far.json").read_bytes()
        controller._record({"status": "valid", "loss": 1., "label": "final_probe"}, update_best=False)
        self.assertEqual((self.root / "best_so_far.json").read_bytes(), frozen)
        self.assertEqual(len(controller.valid), 2)

    def stage_fixture(self, mode, *, update_best=True, allow_rejections=True):
        controller = self.controller()
        path = self.root / "stage" / "plan.json"
        plan = {"cases": [{"id": 1, "label": "fake", "output": "task_001"}]}
        self.write(path, plan)
        sha = overnight._digest(path)
        out = path.parent / "task_001"
        out.mkdir()
        self.write(out / "summary.json", {"fake": True})
        receipt = {"status": "complete", "plan_sha256": sha, "loss": 4.,
                   "artifact_sha256": {"summary.json": overnight._digest(out / "summary.json")}}
        if mode != "missing":
            self.write(out / "case_receipt.json", receipt)
        if mode == "tampered":
            self.write(out / "summary.json", {"tampered": True})
        if mode == "mismatched":
            receipt["plan_sha256"] = "different-plan"
            self.write(out / "case_receipt.json", receipt)
        if mode == "reject":
            self.write(out / "adapter_failure.json", {"error": "bad theta", "type": "InfeasibleThetaError"})
        if mode == "fatal":
            self.write(out / "adapter_failure.json", {"error": "unexpected", "type": "KeyError"})
        result = {"returncode": 1 if mode in ("reject", "fatal") else 0, "timeout": False}
        with patch.object(overnight.adapter, "load_plan", return_value=plan), \
             patch.object(overnight.adapter, "validate_result") as validate, \
             patch.object(overnight, "_run_process", return_value=result), \
             patch.object(controller, "_reap") as reap:
            rows = controller._stage(path, sha, "fake", deadline=controller.search_deadline,
                                     update_best=update_best, allow_rejections=allow_rejections)
        return controller, rows, validate, reap

    def test_stage_pass_requires_receipt_and_validates_result(self):
        controller, rows, validate, _ = self.stage_fixture("pass")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["loss"], 4.)
        self.assertEqual(controller.launched, 1)
        validate.assert_called_once()

    def test_rejected_proposal_uses_slot_but_is_not_ranked(self):
        controller, rows, validate, _ = self.stage_fixture("reject")
        self.assertEqual(rows, [])
        self.assertEqual(controller.launched, 1)
        self.assertEqual(len(controller.rejections), 1)
        self.assertEqual(controller.valid, [])
        self.assertFalse((self.root / "best_so_far.json").exists())
        validate.assert_not_called()

    def test_unknown_failure_is_fatal(self):
        with self.assertRaisesRegex(RuntimeError, "Fatal case failure"):
            self.stage_fixture("fatal")

    def test_smoke_rejection_is_fatal(self):
        with self.assertRaisesRegex(RuntimeError, "Fatal case failure"):
            self.stage_fixture("reject", allow_rejections=False)

    def test_success_without_receipt_is_fatal(self):
        with self.assertRaises((RuntimeError, FileNotFoundError)):
            self.stage_fixture("missing")

    def test_success_with_wrong_plan_receipt_is_fatal(self):
        with self.assertRaises(RuntimeError):
            self.stage_fixture("mismatched")

    def test_success_with_tampered_artifact_is_fatal(self):
        with self.assertRaisesRegex(RuntimeError, "Hash mismatch"):
            self.stage_fixture("tampered")

    def test_no_launch_when_full_case_cannot_finish(self):
        controller = self.controller()
        with patch.object(overnight.adapter, "load_plan", return_value={"cases": [{}]}), \
             patch.object(overnight, "_run_process") as run:
            self.assertEqual(controller._stage(self.root / "plan.json", "sha", "search",
                                              deadline=overnight._now() + 100), [])
            run.assert_not_called()
        self.assertEqual(controller.launched, 0)

    def test_search_preserves_twenty_four_final_launch_slots(self):
        controller = self.controller()
        controller.launched = 276
        with patch.object(overnight.adapter, "load_plan", return_value={"cases": [{}]}), \
             patch.object(overnight, "_run_process") as run:
            try:
                result = controller._stage(self.root / "plan.json", "sha", "search",
                                           deadline=controller.search_deadline, search=True)
            except RuntimeError as error:
                self.assertIn("budget", str(error))
            else:
                self.assertEqual(result, [])
            run.assert_not_called()
        self.assertEqual(controller.launched, 276)

    def test_de_rand_one_has_three_separate_donors(self):
        class DeterministicRng:
            def random(self):
                return 0.
            def randrange(self, _):
                return 0
        target, best, a, b, c = [.2] * 11, [.4] * 11, [.5] * 11, [.3] * 11, [.6] * 11
        trial = overnight.de_proposal(target, best, a, b, DeterministicRng(), rule="rand_1", c=c)
        for value in trial:
            self.assertAlmostEqual(value, .6 + .55 * (.5 - .3))
        trial = overnight.de_proposal(target, best, a, b, DeterministicRng(), rule="current_to_best_1")
        for value in trial:
            self.assertAlmostEqual(value, .2 + .55 * (.4 - .2) + .55 * (.5 - .3))

    def test_final_repeats_replay_original_generator_and_freeze_selection(self):
        controller = self.controller()
        selected_dir = self.root / "original" / "task_001"
        selected = {"panel_design": {"unit_vector": [.5] * 11},
                    "best_candidate": {"candidate": "chosen"}}
        self.write(selected_dir / "summary.json", selected)
        (selected_dir / "target_fit_long.csv").write_text("moment,target,model\n")
        history = selected_dir / "cases" / "chosen" / "transition_path.csv"
        history.parent.mkdir(parents=True)
        history.write_text("year\n2023\n")
        center = {"exact_original_float": .12345678901234568}
        self.write(selected_dir.parent / "center.json", center)
        generator = {"panel_task_id": 7, "panel_size": 23, "panel_design": "local",
                     "radius": .005, "panel_seed": 918}
        original_case = {"id": 1, "center": "center.json", **generator}
        original_plan = {"cases": [original_case]}
        path = self.write(selected_dir.parent / "plan.json", original_plan)
        best = {"status": "valid", "case_id": 1, "loss": 26., "summary": str(selected_dir / "summary.json"),
                "plan": str(path), "plan_sha256": overnight._digest(path)}
        self.write(self.root / "best_so_far.json", best)
        frozen = (self.root / "best_so_far.json").read_bytes()
        plans = []

        def fake_stage(plan_path, plan_sha, stage, **kwargs):
            self.assertFalse(kwargs["update_best"])
            plan = overnight.adapter.read_json(plan_path)
            plans.append(plan)
            rows = []
            for case in plan["cases"]:
                repeat = case["label"].startswith("final_exact_repeat")
                if repeat:
                    for key, value in generator.items():
                        self.assertEqual(case[key], value, key)
                    self.assertEqual(overnight.adapter.read_json(plan_path.parent / case["center"]), center)
                receipt_path = self.write(plan_path.parent / case["output"] / "case_receipt.json",
                                          {"reference": {"exact_twelve_row_fit": True}})
                row = {"label": case["label"], "loss": 26. if repeat else 1., "status": "valid",
                       "receipt": str(receipt_path)}
                controller._record(row, update_best=kwargs["update_best"])
                rows.append(row)
            return rows

        with patch.object(overnight.adapter, "load_plan", return_value=original_plan), \
             patch.object(overnight, "_theta_center", side_effect=lambda summary, units, label: {"unit": units}), \
             patch.object(controller, "_stage", side_effect=fake_stage), \
             patch.object(controller, "_write_jacobian"):
            returned_best, verified, probes = controller._final_checks({})
        self.assertEqual([len(plan["cases"]) for plan in plans], [23, 1])
        self.assertEqual(returned_best, best)
        self.assertTrue(verified)
        self.assertEqual(len(probes), 22)
        self.assertEqual((self.root / "best_so_far.json").read_bytes(), frozen)

    def test_linear_jacobian_and_missing_column(self):
        import numpy as np
        controller = self.controller()
        matrix = np.vstack((np.eye(11), np.ones((1, 11))))
        center = [.5] * 11
        selected = {"panel_design": {"unit_vector": center,
                    "domain": [{"name": f"p{j}"} for j in range(11)]}}
        def save(label, units):
            path = self.write(self.root / label / "summary.json",
                              {"panel_design": {"unit_vector": units}})
            residuals = matrix @ (np.array(units) - .5) + np.arange(12)
            overnight.planner.write_csv(path.parent / "target_fit_long.csv", [
                dict(moment=f"m{k}", target=0, weight=1, standardized_gap=float(value))
                for k, value in enumerate(residuals)])
            return str(path)
        base = save("base", center)
        self.write(self.root / "best_so_far.json", {"summary": base})
        probes = {}
        for j in range(11):
            for sign, delta in (("minus", -.0025), ("plus", .0025)):
                units = center.copy(); units[j] += delta
                label = f"final_jacobian_{j}_{sign}"
                probes[label] = {"status": "valid", "summary": save(label, units)}
        controller._write_jacobian(selected, probes)
        actual = overnight.adapter.read_csv(self.root / "weighted_moment_jacobian.csv")
        actual = np.array([float(row["derivative"]) for row in actual]).reshape(11, 12).T
        np.testing.assert_allclose(actual, matrix, atol=1e-10)
        self.assertEqual(overnight.adapter.read_json(self.root / "jacobian_diagnostic.json")["rank"], 11)
        probes["final_jacobian_0_minus"] = {"status": "not_completed"}
        controller._write_jacobian(selected, probes)
        self.assertEqual(overnight.adapter.read_json(self.root / "jacobian_diagnostic.json")["columns"][0]["method"], "one_sided")
        probes["final_jacobian_0_plus"] = {"status": "not_completed"}
        controller._write_jacobian(selected, probes)
        diagnostic = overnight.adapter.read_json(self.root / "jacobian_diagnostic.json")
        self.assertFalse(diagnostic["complete_matrix"])
        self.assertIsNone(diagnostic["rank"])

    def test_bad_receipt_reaps_before_executor_shutdown(self):
        controller = self.controller()
        events = []
        class FakeExecutor:
            def __init__(self, **kwargs):
                pass
            def __enter__(self):
                return self
            def submit(self, *args):
                future = concurrent.futures.Future()
                future.set_result({"returncode": 0, "timeout": False})
                return future
            def __exit__(self, *args):
                events.append("executor_shutdown")
        plan = {"cases": [{"id": 1, "label": "missing_receipt", "output": "absent"}]}
        with patch.object(overnight.adapter, "load_plan", return_value=plan), \
             patch.object(overnight.futures, "ThreadPoolExecutor", FakeExecutor), \
             patch.object(controller, "_reap", side_effect=lambda: events.append("reap")):
            with self.assertRaises((RuntimeError, FileNotFoundError)):
                controller._stage(self.root / "plan.json", "sha", "fake", deadline=controller.search_deadline)
        self.assertIn("reap", events, "unexpected receipt failures must cancel other running cases")
        self.assertLess(events.index("reap"), events.index("executor_shutdown"))

    def test_fresh_heartbeat_allows_more_than_thirty_minutes_without_completion(self):
        controller = self.controller()
        virtual_time = [1000.]
        plan = {"cases": [{"id": 1, "label": "healthy", "output": "healthy"}]}
        out = self.root / "healthy"
        heartbeat = self.write(out / "heartbeat.json", {"alive": True})
        os.utime(heartbeat, (3001., 3001.))
        self.write(out / "case_receipt.json", {"status": "complete", "plan_sha256": "sha",
                                               "loss": 4., "artifact_sha256": {}})
        waits = []
        def fake_wait(pending, **kwargs):
            waits.append(True)
            if len(waits) == 1:
                virtual_time[0] = 3001.
                return set(), pending
            return pending, set()
        with patch.object(overnight.adapter, "load_plan", return_value=plan), \
             patch.object(overnight.adapter, "validate_result"), \
             patch.object(overnight, "_run_process", return_value={"returncode": 0, "timeout": False}), \
             patch.object(overnight.futures, "wait", side_effect=fake_wait), \
             patch.object(overnight.time, "time", side_effect=lambda: virtual_time[0]), \
             patch.object(controller, "_reap") as reap:
            rows = controller._stage(self.root / "plan.json", "sha", "fake", deadline=controller.search_deadline)
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(waits), 2)
        reap.assert_not_called()

    def test_contract_rejects_changed_pinned_code_before_run(self):
        controller = self.controller()
        source = self.root / "pinned_source.py"
        source.write_text("original source\n")
        original_digest = overnight._digest(source)
        controller.c.update(schema=overnight.SCHEMA, output=str(self.root / "output"),
                            imported_repeat_plan=str(self.root / "repeat.json"),
                            imported_repeat_plan_sha256="repeat-sha", comparison_reference=str(self.root),
                            source=str(source), source_sha256=overnight.adapter.SOURCE,
                            code_sha256={str(source): original_digest},
                            reference_sha256={}, smoke_plan=str(self.root / "smoke.json"),
                            smoke_plan_sha256="smoke-sha", total_seconds=43200, search_seconds=32400,
                            final_reserve_seconds=10800)
        source.write_text("changed source\n")
        with self.assertRaisesRegex(RuntimeError, "Hash mismatch"):
            controller._validate_contract()

    def test_sigterm_cancels_processes_before_unwinding(self):
        controller = self.controller()
        handlers, events = {}, []
        def fake_signal(signum, handler):
            handlers[signum] = handler
        def fake_run():
            try:
                handlers[overnight.signal.SIGTERM](overnight.signal.SIGTERM, None)
            finally:
                events.append("run_unwound")
        with patch.object(overnight, "Controller", return_value=controller), \
             patch.object(overnight.sys, "argv", ["overnight", "--contract", str(self.root / "contract.json"),
                                                   "--contract-sha256", "sha"]), \
             patch.object(overnight.signal, "signal", side_effect=fake_signal), \
             patch.object(controller, "run", side_effect=fake_run), \
             patch.object(controller, "_reap", side_effect=lambda: events.append("reap")):
            with self.assertRaises(KeyboardInterrupt):
                overnight.main()
        self.assertLess(events.index("reap"), events.index("run_unwound"))


if __name__ == "__main__":
    unittest.main()
