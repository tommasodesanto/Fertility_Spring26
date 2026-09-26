"""Fail-closed runner checks; execute on Torch with no native solve."""
import json
import os
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest import mock

import e5f_utility_comparison_design as design
import run_e5f_utility_comparison as runner


class RunnerGuards(unittest.TestCase):
    def test_only_prespecified_native_rejections_are_inadmissible(self):
        for text in runner.INADMISSIBLE_PREFIXES:
            self.assertEqual(runner.failure_status(RuntimeError(text)),"inadmissible_parameter_proposal")
            self.assertEqual(runner.failure_status(ValueError(text)),"unclassified_failure_stop")
        for error in (RuntimeError("occupied value failure"),TimeoutError("deadline"),AssertionError("fiscal")):
            self.assertEqual(runner.failure_status(error),"unclassified_failure_stop")

    def fixture(self, root):
        budget=design.size_search_budget(objective_solve_p90_seconds=1025.7431141,
            objective_overhead_seconds=120,workers_per_arm=10,total_limit_seconds=28800,
            repeat_reserve_seconds=4500,export_reserve_seconds=900,initial_population=40,
            de_generations=3,objective_timeout_seconds=3100)
        budget_path=root/"budget.json"
        budget_path.write_text(json.dumps(budget))
        contract=dict(status="author_approved_frozen_design",launch_permitted=True,
                      files={"budget":{"path":str(budget_path)}},
                      approved_budget=budget,approved_assumptions=["explicit fixture only"])
        case_path=root/"case_plan.json"
        case_path.write_text(json.dumps(dict(contract_sha256="a"*64,arm="floor_linear",
            deadline_epoch=4e9,objective_cap_seconds=3100,point={})))
        args=SimpleNamespace(contract=root/"contract.json",case_plan=case_path,
                             arm="floor_linear",output=root/"output")
        return contract,args

    def test_prepared_contract_cannot_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            contract,args=self.fixture(root)
            contract.update(status="prepared_awaiting_author_review",launch_permitted=False)
            with mock.patch.object(runner,"verified_contract",return_value=contract), mock.patch.object(runner,"setup") as setup:
                with self.assertRaisesRegex(RuntimeError,"no objective launch"):
                    runner.evaluate(args)
                setup.assert_not_called()
                self.assertFalse(args.output.exists())

    def test_duplicate_case_preserves_existing_failure_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            contract,args=self.fixture(root)
            args.output.mkdir()
            failure=args.output/"failure.json"
            failure.write_text("original result must survive")
            with mock.patch.object(runner,"verified_contract",return_value=contract), \
                 mock.patch.object(runner.adapter,"file_hash",return_value="a"*64), \
                 mock.patch.object(runner,"setup") as setup:
                with self.assertRaises(FileExistsError): runner.evaluate(args)
                setup.assert_not_called()
                self.assertEqual(failure.read_text(),"original result must survive")

    def test_bool_approval_does_not_approve_a_different_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            contract,args=self.fixture(root)
            contract["approved_budget"]=True
            with mock.patch.object(runner,"verified_contract",return_value=contract):
                with self.assertRaisesRegex(RuntimeError,"no objective launch"):
                    runner.evaluate(args)

    def test_executing_source_must_be_the_pinned_path(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"contract.json"
            path.write_text(json.dumps({"files":{"run_e5f_utility_comparison.py":{"path":"/wrong/runner.py"}}}))
            with mock.patch.dict(os.environ,{"EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256":runner.adapter.file_hash(path)}):
                with self.assertRaisesRegex(RuntimeError,"executing source"):
                    runner.verified_contract(path)

    def recovery_fixture(self, root):
        contract, args = self.fixture(root)
        contract.update(candidate_failure_policy=runner.RECOVERY_POLICY,
            arms={args.arm: {"objective": {"sha256": "b" * 64}}})
        contract["files"]["run_e5f_utility_comparison.py"] = {"sha256": "c" * 64}
        request = runner.read(args.case_plan)
        request.update(candidate_id=args.output.name, controller_stage="initial",
            runner_source_sha256="c" * 64, deadline_owner="controller", controller_pid=os.getppid())
        args.case_plan.write_text(json.dumps(request))
        return contract, args

    def test_unknown_policy_cannot_fall_back_to_legacy(self):
        with self.assertRaisesRegex(RuntimeError, "unreviewed"):
            runner.recovery_enabled({"candidate_failure_policy": "anything"})

    def test_recovery_context_requires_exact_case_source_stage_and_owner(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, args = self.recovery_fixture(Path(directory))
            request = runner.read(args.case_plan)
            with mock.patch.object(runner.adapter, "file_hash", return_value="a" * 64):
                context = runner.recovery_context(contract, args.contract, args.arm, request, args.output)
                self.assertEqual(context.target_sha256, "b" * 64)
                self.assertEqual(context.point_sha256, design.canonical_fingerprint({}))
                for key, value in (("candidate_id", "other"), ("runner_source_sha256", "d" * 64),
                                   ("controller_stage", "other"), ("controller_pid", -1),
                                   ("deadline_owner", "worker")):
                    with self.subTest(key=key), self.assertRaises((RuntimeError, ValueError)):
                        runner.recovery_context(contract, args.contract, args.arm,
                                                dict(request, **{key: value}), args.output)

    def native_fixture(self, root):
        relative = "code/model/intergen_eqscale_seq_optimized/solver.py"
        path = root / "source" / relative
        path.parent.mkdir(parents=True)
        source = ("DEAD_MASS_TOL=1e-12\nDEAD_VALUE_CUTOFF=-1e9\n"
                  "class InfeasibleThetaError(RuntimeError):\n"
                  "    def __init__(self):\n        super().__init__('synthetic')\n")
        path.write_text(source)
        model = ModuleType("synthetic_native_identity")
        model.__file__ = str(path)
        exec(compile(source, str(path), "exec"), model.__dict__)
        manifest = root / "manifest.json"
        manifest.write_text(json.dumps({"files": {relative: runner.adapter.file_hash(path)}}))
        contract = {"reference_root": str(root), "parent_source_inventory": runner.pin(manifest)}
        return contract, model

    def test_native_capture_authenticates_file_manifest_module_and_gates(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, model = self.native_fixture(Path(directory))
            with mock.patch.dict(sys.modules, {model.__name__: model}):
                native_type, source = runner.authenticated_native_failure_type(contract, {"model": model})
                self.assertIs(native_type, model.InfeasibleThetaError)
                self.assertEqual(source["native_gate_tolerance"], 1e-12)
                for key, value in (("DEAD_MASS_TOL", 1e-10), ("DEAD_VALUE_CUTOFF", -1e8),
                                   ("__file__", "/wrong/solver.py")):
                    with self.subTest(key=key), mock.patch.object(model, key, value):
                        with self.assertRaisesRegex(RuntimeError, "pinned solver"):
                            runner.authenticated_native_failure_type(contract, {"model": model})
            with self.assertRaisesRegex(RuntimeError, "pinned solver"):
                runner.authenticated_native_failure_type(contract, {"model": model})
            Path(contract["parent_source_inventory"]["path"]).write_text("{}")
            with self.assertRaisesRegex(RuntimeError, "manifest changed"):
                runner.authenticated_native_failure_type(contract, {"model": model})

    def test_repair_uses_external_watchdog_and_preserves_raw_exception(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, args = self.recovery_fixture(Path(directory))
            old = SimpleNamespace(evaluate_point=mock.Mock(side_effect=SystemError("unknown kernel defect")))
            setup_result = (None, old, None, None, None, None, {"model": None}, None)
            source = {"native_gate_tolerance": 1e-12}
            with mock.patch.object(runner, "verified_contract", return_value=contract), \
                 mock.patch.object(runner.adapter, "file_hash", return_value="a" * 64), \
                 mock.patch.object(runner, "setup", return_value=setup_result), \
                 mock.patch.object(runner, "authenticated_native_failure_type", return_value=(None, source)), \
                 mock.patch.object(runner.signal, "setitimer") as timer:
                with self.assertRaisesRegex(SystemError, "unknown kernel defect"):
                    runner.evaluate(args)
                timer.assert_not_called()
            failure = runner.read(args.output / "failure.json")
            self.assertEqual(failure["status"], "unclassified_failure_stop")
            self.assertEqual(failure["recovery_evidence"]["raw_type"], "SystemError")
            self.assertFalse(failure["recovery_evidence"]["narrow_infeasibility_verified"])
            self.assertEqual(failure["recovery_evidence"]["context"]["candidate_id"], args.output.name)
            provenance = runner.read(args.output / "attempt_provenance.json")
            self.assertEqual(provenance["arm"], args.arm)
            self.assertEqual(provenance["context"]["stage"], "initial")
            old.evaluate_point.assert_called_once()

    def test_setup_failure_has_no_authenticated_native_source(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, args = self.recovery_fixture(Path(directory))
            with mock.patch.object(runner, "verified_contract", return_value=contract), \
                 mock.patch.object(runner.adapter, "file_hash", return_value="a" * 64), \
                 mock.patch.object(runner, "setup", side_effect=RuntimeError("source mismatch")), \
                 mock.patch.object(runner.signal, "setitimer") as timer:
                with self.assertRaisesRegex(RuntimeError, "source mismatch"):
                    runner.evaluate(args)
                timer.assert_not_called()
            failure = runner.read(args.output / "failure.json")
            self.assertIsNone(failure["native_source"])
            self.assertFalse(failure["recovery_evidence"]["narrow_infeasibility_verified"])

    def test_unserializable_diagnostic_does_not_mask_the_original_error(self):
        with tempfile.TemporaryDirectory() as directory:
            contract, args = self.recovery_fixture(Path(directory))
            context = runner.recovery.Context(args.output.name, "initial", "a" * 64,
                "c" * 64, "b" * 64, design.canonical_fingerprint({}))
            bad = runner.recovery.FailureEvidence(context, "SystemError", "original failure",
                {"dead_mass": float("nan")}, False, "nonfinite diagnostic")
            with mock.patch.object(runner, "verified_contract", return_value=contract), \
                 mock.patch.object(runner.adapter, "file_hash", return_value="a" * 64), \
                 mock.patch.object(runner, "setup", side_effect=SystemError("original failure")), \
                 mock.patch.object(runner.recovery, "capture_native_failure", return_value=bad):
                with self.assertRaisesRegex(SystemError, "original failure"):
                    runner.evaluate(args)
            failure = runner.read(args.output / "failure.json")
            self.assertEqual(failure["error"], "original failure")
            self.assertIn("recovery_capture_error", failure)
            self.assertNotIn("recovery_evidence", failure)


if __name__ == "__main__": unittest.main()
