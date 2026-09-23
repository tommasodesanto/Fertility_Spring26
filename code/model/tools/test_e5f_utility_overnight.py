"""Zero-solve contract checks for the inherited-entry utility comparison."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import prepare_e5f_utility_overnight as prepare
import run_e5f_utility_overnight as worker

ROOT = Path(__file__).resolve().parents[3]
V5 = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_wealth_v1/smoke_v5/frozen_v5_plan.json"
SELECTED_B = ROOT / "output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/final_readout/B/selected/plan.json"


class OvernightContracts(unittest.TestCase):
    def test_four_templates_and_zero_solve_preflight(self):
        with tempfile.TemporaryDirectory(prefix="e5f_utility_test_") as scratch:
            bundle = Path(scratch) / "bundle"
            receipt = prepare.build(V5, SELECTED_B, bundle, bundle, sys.executable,
                ROOT / "code/model/tools/build_period_earnings_process.py",
                ROOT / "code/model/tools/run_e5f_preference_share_candidate.py",
                ROOT / "code/model/tools/run_e5f_utility_overnight.py",
                ROOT / "code/cluster/submit_e5f_utility_overnight.sh")
            manifest = worker.read(receipt["manifest"])
            worker.validate_manifest(manifest, Path(receipt["manifest"]))
            self.assertEqual(manifest["total_workers"], 40)
            self.assertEqual(len({row["target_system_sha256"] for row in manifest["templates"].values()}), 1)
            adapter = (bundle / "tools/run_e5f_preference_share_candidate.py").read_text()
            self.assertIn(prepare.ADAPTER_GUARD_NEW, adapter)
            self.assertNotIn(prepare.ADAPTER_GUARD_OLD, adapter)
            selected = json.loads(SELECTED_B.read_text())["structural_parameters"]
            for cell, item in manifest["templates"].items():
                plan = worker.read(item["path"])
                self.assertEqual(plan["entry_specification"]["rule"], "fixed_reference_marginal")
                self.assertEqual(plan["income_specification"]["constructor_arguments"]["n_persistent"], 15 if cell[0] == "B" else 7)
                self.assertEqual(plan["income_specification"]["constructor_arguments"]["n_iid"], 1 if cell[0] == "B" else 3)
                self.assertEqual(len(plan["parameter_bounds"]), 9 if cell.endswith("floor") else 10)
                self.assertEqual({k: plan["structural_parameters"][k] for k in worker.SHARED},
                                 {k: selected[k] for k in worker.SHARED})
                result = subprocess.run([sys.executable, plan["adapter_path"], "--plan", item["path"],
                    "--output", str(bundle / "preflight" / cell), "--preflight"],
                    text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr[-1000:])
                self.assertEqual(worker.read(bundle / "preflight" / cell / "receipt.json")["solves"], 0)
                dynamic = bundle / "dynamic" / cell / "plan.json"
                dynamic.parent.mkdir(parents=True)
                worker.make_case(plan, plan["structural_parameters"], f"{cell}_dynamic",
                                 dynamic, 7200, {"design": "test_dynamic"})
                check = subprocess.run([sys.executable, plan["adapter_path"], "--plan", str(dynamic),
                    "--output", str(bundle / "dynamic_preflight" / cell), "--preflight"],
                    text=True, capture_output=True)
                self.assertEqual(check.returncode, 0, check.stderr[-1000:])
                self.assertEqual(worker.read(bundle / "dynamic_preflight" / cell / "receipt.json")["solves"], 0)

    def test_adaptive_proposal_and_bounds(self):
        center = {"H0": 7.5, "beta_annual": .986, "chi": 1., "first_birth_fixed_cost": .3,
                  "kappa_fert": .29, "kappa_fert_continuation": .45,
                  "theta0": .04, "theta1": .09, "h_P": 2.3}
        bounds = {name: [0., 10.] for name in center}
        bounds["beta_annual"] = [.94, .99]
        bounds["h_P"] = [.1, 2.3]
        a = worker.propose(center, bounds, 101, 2)
        changed = dict(center, H0=4.2)
        b = worker.propose(changed, bounds, 101, 2)
        self.assertNotEqual(a["H0"], b["H0"])
        self.assertTrue(all(bounds[k][0] <= v <= bounds[k][1] for k, v in a.items()))
        share_center = {k: v for k, v in center.items() if k != "h_P"}
        share_center.update(delta_alpha_jump=.02, delta_alpha=.017)
        share_bounds = {k: v for k, v in bounds.items() if k != "h_P"}
        share_bounds.update(delta_alpha_jump=[0., .25], delta_alpha=[0., .25])
        share = worker.propose(share_center, share_bounds, 101, 2)
        self.assertEqual({k: a[k] for k in worker.SHARED}, {k: share[k] for k in worker.SHARED})
        coordinate = worker.propose(center, bounds, 101, 5)
        self.assertEqual(sum(coordinate[k] != center[k] for k in center), 1)
        with tempfile.TemporaryDirectory() as scratch:
            plan = {"parameter_bounds": bounds, "structural_parameters": center,
                    "starting_structural_parameters": center}
            path = Path(scratch) / "case.json"
            worker.make_case(plan, a, "test_case", path, 7200)
            case = worker.read(path)["cases"][0]
            self.assertEqual((case["native_seconds"], case["wrapper_seconds"], case["seconds"]),
                             (6900, 7100, 7200))

    def test_checkpoint_provenance_exclusion_preserves_economic_gap(self):
        first = {"target_fit": [{"model_checkpoint_sha256": "aaa", "model_source_path": "/a",
                                 "target": .72, "model": 1.34, "gap": .62, "actual_weight": 5.}]}
        second = {"target_fit": [{"model_checkpoint_sha256": "bbb", "model_source_path": "/b",
                                  "target": .72, "model": 1.34, "gap": .62, "actual_weight": 5.}]}
        self.assertEqual(worker.verified_target_fit(first, "aaa", "B_floor"),
                         worker.verified_target_fit(second, "bbb", "B_floor"))
        second["target_fit"][0]["gap"] = .63
        self.assertNotEqual(worker.verified_target_fit(first, "aaa", "B_floor"),
                            worker.verified_target_fit(second, "bbb", "B_floor"))
        with self.assertRaises(worker.ContractError):
            worker.verified_target_fit(first, "wrong", "B_floor")


if __name__ == "__main__":
    unittest.main()
