"""Focused guards for the bounded calibration search's measurement and selection."""
import copy
import json
import tempfile
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import build_e5f_bounded_refinement_plan as planner
import run_e5f_bounded_calibration_refinement as adapter


ROOT = Path(__file__).resolve().parents[3]
REFERENCE = ROOT / "output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_youngown_overid_coord_20260904a/task_003/summary.json"


class BoundedRefinementTests(unittest.TestCase):
    def test_generation_rejects_planner_drift_before_reading_cases(self):
        with patch.object(sys, "argv", [planner.__file__, "joint", "--plan", "must_not_read.json",
                 "--plan-sha256", "unused", "--outdir", "must_not_create", "--expected-planner-sha256", "changed"]):
            with self.assertRaisesRegex(RuntimeError, "Hash mismatch"):
                planner.main()

    def test_plan_rejects_changed_target_fingerprint(self):
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "plan.json"
            adapter.write_json(path, dict(schema="e5f_bounded_refinement_v1", source_sha256=adapter.SOURCE,
                target_fingerprint="different", code_bundle_sha256=adapter.BUNDLE, cases=[]))
            with self.assertRaisesRegex(RuntimeError, "target system"):
                adapter.load_plan(path, adapter.digest(path))

    def test_plan_accepts_only_the_two_verified_bundles(self):
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "plan.json"
            for bundle in adapter.SUPPORTED_BUNDLES:
                adapter.write_json(path, dict(schema="e5f_bounded_refinement_v1",
                    source_sha256=adapter.SOURCE, target_fingerprint=adapter.TARGET,
                    code_bundle_sha256=bundle, cases=[]))
                self.assertEqual(adapter.load_plan(path, adapter.digest(path))["code_bundle_sha256"], bundle)
            changed = adapter.read_json(path)
            changed["code_bundle_sha256"] = "unreviewed"
            adapter.write_json(path, changed)
            with self.assertRaisesRegex(RuntimeError, "unverified scientific bundle"):
                adapter.load_plan(path, adapter.digest(path))

    def test_first_child_pattern_preserves_three_child_floor(self):
        reference = adapter.read_json(REFERENCE)
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            rows = []
            for i in range(1, 24):
                path = root / f"case_{i}.json"
                adapter.write_json(path, reference)
                # Distinct improvements in four parameters exercise both joint
                # coordinate families; the floor parameters alone worsen fit.
                loss = 30.0 - i/20 if i in (3, 5, 9, 11) else 30.0 + i/100
                if i == 1:
                    loss = 30.0
                rows.append(dict(id=i, loss=loss, summary=str(path), summary_sha256=adapter.digest(path)))
            rows.sort(key=lambda r: r["loss"])
            out = root / "joint"
            with patch("builtins.print"):
                planner.joint(dict(stage="round1_small_coordinate"), {}, rows, out)
            plan = adapter.read_json(out / "plan.json")
            self.assertLessEqual(len(plan["cases"]), 12)
            theta = reference["best_candidate"]["theta"]
            initial_three = theta["hbar_first_child_jump"] + 3*theta["hbar_child_rooms"]
            initial_one = theta["hbar_first_child_jump"] + theta["hbar_child_rooms"]
            rebalanced = 0
            for case in plan["cases"]:
                payload = adapter.read_json(out / case["center"])
                delta = payload["proposal"]["delta_first_child_rooms"]
                if delta is None:
                    continue
                rebalanced += 1
                new = payload["best_candidate"]["theta"]
                self.assertAlmostEqual(new["hbar_first_child_jump"] + 3*new["hbar_child_rooms"], initial_three, places=14)
                self.assertAlmostEqual(new["hbar_first_child_jump"] + new["hbar_child_rooms"] - initial_one, 2*delta/3, places=14)
            self.assertEqual(rebalanced, 6)

    def test_repeat_uses_cross_round_best_and_original_generator(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw); stage = root / "round1"; stage.mkdir()
            center = adapter.read_json(REFERENCE)
            adapter.write_json(stage / "center.json", center)
            chosen = dict(id=7, center="center.json", panel_task_id=7, panel_size=23,
                          panel_design="coordinate", panel_seed=2026090506, radius=.005)
            plan = dict(schema="e5f_bounded_refinement_v1", source_sha256=adapter.SOURCE,
                target_fingerprint=adapter.TARGET, code_bundle_sha256=adapter.BUNDLE, cases=[chosen])
            adapter.write_json(stage / "plan.json", plan)
            summary = stage / "summary.json"; adapter.write_json(summary, center)
            (stage / "target_fit_long.csv").write_text("test fixture\n")
            case_dir = stage / "cases" / center["best_candidate"]["candidate"]
            case_dir.mkdir(parents=True); (case_dir / "transition_path.csv").write_text("test fixture\n")
            best = dict(id=7, summary=str(summary), plan=str(stage / "plan.json"), plan_sha256=adapter.digest(stage / "plan.json"))
            adapter.write_json(root / "best_so_far.json", dict(best=best))
            last = root / "round2"; last.mkdir()
            with patch("builtins.print"):
                planner.repeats(plan, {}, [dict(plan=str(last / "plan.json"), summary="must not select this worse point")], root / "repeat")
            result = adapter.read_json(root / "repeat/plan.json")
            self.assertEqual(result["selected_summary_sha256"], adapter.digest(summary))
            self.assertTrue(all(c["panel_task_id"] == 7 and c["panel_size"] == 23 for c in result["cases"]))


if __name__ == "__main__":
    unittest.main()
