"""Fail-closed runner checks; execute on Torch with no native solve."""
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
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


if __name__ == "__main__": unittest.main()
