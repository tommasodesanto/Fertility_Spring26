from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import run_e5f_earnings_entry_battery as battery


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class EarningsEntryBatteryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.cells = {
            "A": {"earnings_specification": "single_ar1", "income_mapping": "direct_period", "persistent_states": 7, "iid_states": 1, "entry_rule": "reference_marginal", "source_inventory_sha256": "pending", "source_file_count": 1},
            "B": {"earnings_specification": "single_ar1", "income_mapping": "direct_period", "persistent_states": 7, "iid_states": 1, "entry_rule": "zero_assets", "source_inventory_sha256": "pending", "source_file_count": 1},
            "C": {"earnings_specification": "persistent_iid", "income_mapping": "direct_period", "persistent_states": 7, "iid_states": 3, "entry_rule": "reference_marginal", "source_inventory_sha256": "pending", "source_file_count": 1},
            "D": {"earnings_specification": "persistent_iid", "income_mapping": "direct_period", "persistent_states": 7, "iid_states": 3, "entry_rule": "zero_assets", "source_inventory_sha256": "pending", "source_file_count": 1},
        }
        self.manifest = {"schema": battery.SCHEMA, "status": "ready", "cells": self.cells,
            "workers_per_cell": 10, "case_seconds": 3200, "total_seconds": 3600,
            "minimum_next_proposal_seconds": 2500, "objective_canonical_sha256": "objective",
            "target_system_sha256": "targets", "scored_moment_count": 13, "parameter_row_count": 17,
            "smoke": {"cases": []}, "production": {"cases": []}}
        for cell_id in self.cells:
            self.manifest["smoke"]["cases"].append({"cell_id": cell_id, "task_id": 1,
                "proposals": [{"proposal_id": "smoke", "case_id": f"smoke-{cell_id}", "plan_path": "plan", "plan_sha256": "hash"}]})
            for task_id in range(1, 11):
                self.manifest["production"]["cases"].append({"cell_id": cell_id, "task_id": task_id,
                    "proposals": [{"proposal_id": f"p{task_id}", "case_id": f"case-{cell_id}-{task_id}",
                                   "plan_path": "plan", "plan_sha256": "hash"}]})

    def tearDown(self):
        self.tmp.cleanup()

    def test_four_cell_factorial_and_40_worker_shape(self):
        battery.validate_manifest(self.manifest)
        self.assertEqual(len(self.manifest["production"]["cases"]), 40)

    def test_rejects_wrong_worker_counts_and_more_than_six_proposals(self):
        bad = json.loads(json.dumps(self.manifest))
        bad["production"]["cases"].pop()
        with self.assertRaisesRegex(battery.ContractError, "40 tasks"):
            battery.validate_manifest(bad)
        bad = json.loads(json.dumps(self.manifest))
        row = bad["production"]["cases"][0]
        row["proposals"] = [{"proposal_id": str(i), "case_id": f"c{i}", "plan_path": "p", "plan_sha256": "h"} for i in range(7)]
        with self.assertRaisesRegex(battery.ContractError, "one to six"):
            battery.validate_manifest(bad)

    def make_plan(self, entry_rule="reference_marginal"):
        source_root = self.root / "source"
        source_root.mkdir(exist_ok=True)
        source = source_root / "model.py"
        source.write_text("source pin\n")
        inventory = {"model.py": sha(source)}
        cell = {"earnings_specification": "persistent_iid", "income_mapping": "direct_period",
                "persistent_states": 7, "iid_states": 3, "entry_rule": entry_rule,
                "source_inventory_sha256": battery.canonical_sha256(inventory), "source_file_count": 1}
        initial = self.root / "initial.json"
        initial.write_text(json.dumps({"source_sha256": inventory}))
        adapter = self.root / "adapter.py"
        adapter.write_text("# pinned adapter\n")
        params = {name: .2 for name in battery.PARAMETERS}
        plan = {"objective_canonical_sha256": "objective", "target_system_sha256": "targets",
            "structural_parameters": params, "source_root": str(source_root),
            "income_specification": {"mapping": "direct_period", "constructor_arguments": {"n_persistent": 7, "n_iid": 3}},
            "entry_specification": {"rule": entry_rule},
            "files": {"initial_contract": {"path": str(initial), "sha256": sha(initial)}},
            "adapter_path": str(adapter), "adapter_sha256": sha(adapter),
            "cases": [{"id": "one", "arm": "literature_income_purchase", "repetitions": 1}]}
        path = self.root / "plan.json"
        path.write_text(json.dumps(plan))
        proposal = {"proposal_id": "p1", "case_id": "one", "plan_path": str(path), "plan_sha256": sha(path)}
        return cell, plan, proposal

    def test_full_source_and_plan_contract_accepts_nonzero_entry_rule(self):
        cell, plan, proposal = self.make_plan("reference_marginal")
        self.assertEqual(battery.verify_source_inventory(plan, cell), {"model.py": sha(self.root / "source/model.py")})
        manifest = {"cells": {"C": cell}, "objective_canonical_sha256": "objective", "target_system_sha256": "targets"}
        path, validated = battery.validate_plan(manifest, "C", proposal)
        self.assertEqual(path.name, "plan.json")
        self.assertEqual(validated["entry_specification"]["rule"], "reference_marginal")

    def test_source_content_drift_fails_before_scoring(self):
        cell, plan, _ = self.make_plan()
        (self.root / "source/model.py").write_text("changed\n")
        with self.assertRaisesRegex(battery.ContractError, "source pin mismatch"):
            battery.verify_source_inventory(plan, cell)

    def test_wrong_objective_or_target_hash_fails(self):
        cell, plan, proposal = self.make_plan()
        manifest = {"cells": {"C": cell}, "objective_canonical_sha256": "wrong", "target_system_sha256": "targets"}
        with self.assertRaisesRegex(battery.ContractError, "objective fingerprint"):
            battery.validate_plan(manifest, "C", proposal)

    def test_execute_task_runs_two_then_budget_stops_and_keeps_latest_best(self):
        cell, plan, _ = self.make_plan()
        proposals = []
        for i in range(3):
            plan_i = dict(plan)
            plan_i["cases"] = [{"id": f"case-{i}", "arm": "literature_income_purchase", "repetitions": 1}]
            path = self.root / f"plan-{i}.json"
            path.write_text(json.dumps(plan_i))
            proposals.append({"proposal_id": f"p{i}", "case_id": f"case-{i}",
                              "plan_path": str(path), "plan_sha256": sha(path)})
        cell.update({"source_file_count": 1})
        manifest = {"cells": {"C": cell}, "objective_canonical_sha256": "objective",
                    "target_system_sha256": "targets", "minimum_next_proposal_seconds": 2500}
        row = {"task_id": 1, "proposals": proposals}
        output = self.root / "output"
        smoke_dir = output / "smoke/C/worker_01"
        smoke_dir.mkdir(parents=True)
        (smoke_dir / "status.json").write_text(json.dumps({"status": "completed", "elapsed_seconds": 100}))
        clock = [0.0]
        class Proc:
            calls = 0
            def __init__(self, cmd, **kwargs):
                type(self).calls += 1
                self.pid = 1000 + type(self).calls
                self.returncode = 0
                self._output = Path(cmd[cmd.index("--output") + 1])
                self._output.mkdir(parents=True, exist_ok=True)
                score_dir = self._output / "evaluation/scored_repetition_01"
                score_dir.mkdir(parents=True)
                structural = [{"parameter": n, "structural_coordinate": True} for n in battery.PARAMETERS]
                structural += [{"parameter": f"fixed_{j}", "structural_coordinate": False} for j in range(8)]
                score = {"schema": battery.OBJECTIVE_SCHEMA, "loss": 5.0 - type(self).calls,
                         "target_fit": [{} for _ in range(13)], "parameters": structural}
                (score_dir / "score.json").write_text(json.dumps(score))
                (self._output / "evaluation/summary.json").write_text(json.dumps({
                    "status": "verified_scored_candidate", "objective_canonical_sha256": "objective", "repetitions": 1}))
                clock[0] += 600
            def poll(self): return 0
        result = battery.execute_task(manifest, "C", row, "production", output,
            popen=Proc, monotonic=lambda: clock[0], sleep=lambda _: None, killer=lambda _: None)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["proposal_count"], 2)
        self.assertEqual(result["stop_reason"], "insufficient_budget_for_next_proposal")
        worker = output / "production/C/worker_01"
        self.assertTrue((worker / "latest.json").is_file())
        self.assertEqual(json.loads((worker / "latest.json").read_text())["proposal_id"], "p1")
        self.assertEqual(json.loads((worker / "best_so_far.json").read_text())["proposal_id"], "p1")
        self.assertEqual(json.loads((worker / "status.json").read_text())["proposal_count"], 2)
        self.assertFalse((worker / "proposal_03_p2/result").exists())


if __name__ == "__main__":
    unittest.main()
