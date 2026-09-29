"""Pure stdlib tests for the fail-closed preference-estimation controller."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch


CONTROLLER = Path(__file__).parents[2] / "cluster" / "run_e5f_preference_estimation_batch.py"
STAGED_SOURCE = Path(__file__).parent
SPEC = importlib.util.spec_from_file_location("preference_batch", CONTROLLER)
batch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch)


class BatchControllerTests(unittest.TestCase):
    def plan_args(self):
        cache_budget = 64 * 1024 ** 3
        return SimpleNamespace(
            cache_max_bytes=cache_budget, kind="four_successive", horizon=[6, 7],
            seed_horizon=6, perturbed_date=1, log_step=1e-5,
            reference_manifest_sha256=batch.REFERENCE_MANIFEST_SHA256,
            fit_evaluations=5, endpoint_evaluations=2, path_evaluations=2,
            source_dir=STAGED_SOURCE, readiness=Path("ready.json"), output=Path("plan.json"),
            jacobian_receipt=None, total_seconds=10., candidate_seconds=2.,
            endpoint_seconds=2., mapping_seconds=2., path_seconds=2., jacobian_seconds=2.,
            empirical_blocks=Path("blocks.csv"), annual=Path("annual.csv"),
            housing="fixed_stock", enable_execution=False,
        )

    def test_plan_uses_explicit_positive_cache_budget(self):
        args = self.plan_args()
        with patch.object(batch, "source_pins", return_value={}), \
             patch.object(batch, "check_readiness", return_value={}), \
             patch.object(batch, "target_contract", return_value={}), \
             patch.object(batch, "preflight_plan"):
            plan = batch.plan_from_args(args)
        self.assertEqual(plan["path"]["cache_max_bytes"], args.cache_max_bytes)

    def test_real_driver_structural_preflight_accepts_complete_64gib_plan(self):
        driver = batch.load_staged_driver(STAGED_SOURCE)
        args = self.plan_args()
        with patch.object(batch, "source_pins", return_value={"source": "pin"}), \
             patch.object(batch, "check_readiness", return_value={"path": "ready", "sha256": "pin"}), \
             patch.object(batch, "target_contract", return_value={"rows": [1]}), \
             patch.object(batch, "load_staged_driver", return_value=driver), \
             patch.object(driver, "validate_launch_inputs") as shared:
            plan = batch.plan_from_args(args)
        self.assertFalse(plan["execution_enabled"])
        self.assertEqual(plan["path"]["cache_max_bytes"], 64 * 1024 ** 3)
        shared.assert_called_once()
        self.assertTrue(shared.call_args.args[0]["execution_enabled"])

    def test_incoherent_plan_fails_real_driver_preflight(self):
        driver = batch.load_staged_driver(STAGED_SOURCE)
        args = self.plan_args()
        with patch.object(batch, "source_pins", return_value={"source": "pin"}), \
             patch.object(batch, "check_readiness", return_value={"path": "ready", "sha256": "pin"}), \
             patch.object(batch, "target_contract", return_value={"rows": [1]}), \
             patch.object(batch, "load_staged_driver", return_value=driver), \
             patch.object(driver, "validate_launch_inputs"):
            plan = batch.plan_from_args(args)
            plan["budget"]["maximum_policy_calls"] -= 1
            with self.assertRaisesRegex(ValueError, "solve count"):
                batch.preflight_plan(plan, args.source_dir)

    def test_preflight_failure_writes_no_plan(self):
        driver = batch.load_staged_driver(STAGED_SOURCE)
        args = self.plan_args()
        with tempfile.TemporaryDirectory() as directory:
            args.output = Path(directory) / "plan.json"
            with patch.object(batch, "source_pins", return_value={"source": "pin"}), \
                 patch.object(batch, "check_readiness", return_value={"path": "ready", "sha256": "pin"}), \
                 patch.object(batch, "target_contract", return_value={"rows": [1]}), \
                 patch.object(batch, "load_staged_driver", return_value=driver), \
                 patch.object(driver, "validate_launch_inputs", side_effect=ValueError("preflight rejected")), \
                 self.assertRaisesRegex(ValueError, "preflight rejected"):
                batch.create_plan(args)
            self.assertFalse(args.output.exists())

    def test_overlay_relative_path_contract_rejects_escape_forms(self):
        # Mirrors the launcher contract: only paths mounted at the writable
        # fixed-reference overlay may be staged or created.
        prefix = "output/model/fixed_reference_transition_20260928/"
        def mapped(value):
            if value.startswith("/") or "//" in value:
                raise ValueError("outside overlay")
            parts = value.split("/")
            if not value.startswith(prefix) or any(x in ("", ".", "..") for x in parts):
                raise ValueError("outside overlay")
            return "/scratch/overlay/" + value[len(prefix):]
        self.assertEqual(mapped(prefix + "four_shock_v1/plan.json"),
                         "/scratch/overlay/four_shock_v1/plan.json")
        for bad in ("/tmp/plan.json", "code/cluster/x.py", prefix + "../plan.json", prefix + "a//b"):
            with self.assertRaises(ValueError): mapped(bad)

    def test_policy_count_is_explicit_and_includes_seed_measurement(self):
        # 10H is the driver's conservative seed allowance, exceeding the five mappings.
        expected = 10 * 10 + 4 * 12 * (24 + 2 + (2 * 10 * 16 + 2 * 14 * 16)) + 8
        self.assertEqual(batch.maximum_policy_calls("four_successive", 10, 12, 24, 16, [10, 14]), expected)
        self.assertEqual(batch.maximum_policy_calls("one_permanent", 10, 12, 24, 16, [10, 14]),
                         10 * 10 + 12 * (24 + 2 + (2 * 10 * 16 + 2 * 14 * 16)) + 8)

    def test_readiness_fails_closed_on_source_mismatch_or_disabled_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "readiness.json"; pins = {"a": "1"}
            path.write_text(json.dumps({"status": "PASS", "tests_passed": True,
                "native_endpoint_and_fertility_smoke_passed": True, "preference_changes": False,
                "historical_fit": False, "estimator_sources": pins}))
            self.assertEqual(batch.check_readiness(path, pins)["path"], str(path))
            with self.assertRaisesRegex(ValueError, "pins"):
                batch.check_readiness(path, {"a": "2"})
            path.write_text(json.dumps({"status": "PASS", "tests_passed": True,
                "native_endpoint_and_fertility_smoke_passed": True, "preference_changes": True,
                "historical_fit": False, "estimator_sources": pins}))
            with self.assertRaisesRegex(ValueError, "unchanged-preference"):
                batch.check_readiness(path, pins)

    def test_target_windows_recompute_and_reject_duplicate_year(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); blocks = root / "blocks.csv"; annual = root / "annual.csv"
            with blocks.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["decision_year", "birth_year_start", "birth_year_end", "period_tfr_arithmetic_mean"]); writer.writeheader()
                writer.writerows({"decision_year": y, "birth_year_start": y + 1, "birth_year_end": y + 4, "period_tfr_arithmetic_mean": 2 - i / 10} for i, y in enumerate((2007, 2011, 2015, 2019)))
            rows = []
            for i, year in enumerate((2007, 2011, 2015, 2019)):
                rows.extend({"year": y, "period_tfr_births_per_woman": 2 - i / 10, "status": "verified_published_final", "source_url": "https://example.test"} for y in range(year + 1, year + 5))
            with annual.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            self.assertEqual(len(batch.target_contract(blocks, annual)["rows"]), 4)
            rows[-1]["year"] = rows[-2]["year"]
            with annual.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "unique annual"):
                batch.target_contract(blocks, annual)

    def test_execute_rejects_disabled_plan_before_any_driver_call(self):
        with tempfile.TemporaryDirectory() as directory:
            plan = Path(directory) / "plan.json"
            plan.write_text(json.dumps({"execution_enabled": False}))
            args = SimpleNamespace(plan=plan, plan_sha256=batch.sha(plan), source_dir=Path(directory),
                                   output=Path(directory) / "out", python="python")
            with patch.object(batch, "source_pins", side_effect=AssertionError("must not inspect source")), \
                 self.assertRaisesRegex(ValueError, "disabled"):
                batch.execute(args)

    def test_nonfinite_budget_is_rejected(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaisesRegex(ValueError, "Finite positive"):
                batch.require_positive_budgets({"total_seconds": value})


if __name__ == "__main__":
    unittest.main()
