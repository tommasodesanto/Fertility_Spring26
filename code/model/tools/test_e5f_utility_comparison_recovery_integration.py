"""Actual controller methods with synthetic subprocesses; no native imports."""
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import Mock

import run_e5f_utility_comparison_search as module


RUNNER = r'''
import argparse,hashlib,json,os,pathlib,time
p=argparse.ArgumentParser()
for key in ('stage','contract','arm','case-plan','output'): p.add_argument('--'+key)
a=p.parse_args(); root=pathlib.Path(a.output); root.mkdir()
read=lambda path:json.loads(pathlib.Path(path).read_text())
digest=lambda path:hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
save=lambda path,value:path.write_text(json.dumps(value))
c=read(a.contract); plan=read(a.case_plan)
if plan['point']['kind']==10: time.sleep(20)
context=dict(candidate_id=plan['candidate_id'],stage=plan['controller_stage'],
contract_sha256=digest(a.contract),source_sha256=plan['runner_source_sha256'],
target_sha256=c['arms'][a.arm]['objective']['sha256'],
point_sha256=hashlib.sha256(json.dumps(plan['point'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest())
save(root/'attempt_provenance.json',dict(context=context,arm=a.arm,candidate_failure_policy='reviewed_failure_v1',case_plan_sha256=digest(a.case_plan)))
kind=plan['point']['kind']
if kind==4: time.sleep(20)
if kind==5: raise SystemError('unknown error without recovery evidence')
if kind==6: time.sleep(.15)
if kind in (2,3,7,8,9):
 msg='Old-steady-state fertility normalization is not bracketed: fixture' if kind==3 else 'unknown scientific error'
 ev=dict(context=context,raw_type='RuntimeError',raw_error=msg,native_payload={},narrow_infeasibility_verified=False,validation='exception_type_not_the_authenticated_native_class')
 if kind==7: ev['context']['point_sha256']='0'*64
 if kind==8:
  row=dict(age=34.,b=0.,child_state=1,income=.22,location=0,mass=1.2e-13,parity=1,slack=-.024,tenure=0,transfer=0.,unsecured_position=0.,z=.054)
  rows=[dict(row,b=-.116,unsecured_position=-.116),dict(row,mass=8.4e-13),dict(row,parity=2),dict(row,parity=2,child_state=2)]
  ev.update(raw_type='InfeasibleThetaError',narrow_infeasibility_verified=True,validation='exact_native_type_and_structured_gate_failure',native_payload=dict(stage='forward_age_34',dead_mass=1.2e-12,native_gate_tolerance=1e-12,census=rows))
 source=pathlib.Path(c['reference_root'])/'source/code/model/intergen_eqscale_seq_optimized/solver.py'
 native=dict(native_source_path=str(source),native_source_sha256=digest(source),native_manifest_sha256=c['parent_source_inventory']['sha256'],native_gate_tolerance=1e-12,native_value_cutoff=-1e9)
 if kind==9: native['native_source_sha256']='0'*64
 save(root/'failure.json',dict(error_type=ev['raw_type'],error=msg,status='inadmissible_parameter_proposal',arm=a.arm,candidate_failure_policy='reviewed_failure_v1',recovery_evidence=ev,native_source=native))
 raise SystemExit(1)
case=root/'case'; case.mkdir()
save(case/'receipt.json',dict(comparison_contract_sha256=digest(a.contract),utility_comparison_arm=a.arm,target_system_sha256=c['arms'][a.arm]['objective']['sha256'],point=plan['point'],loss=plan['point']['kind']))
for name in ('initial_state.pkl.gz','target_fit.csv','parameters.csv'): (case/name).write_text('synthetic artifact only')
'''

COLLECTOR = r'''
import argparse,json,pathlib
p=argparse.ArgumentParser()
for key in ('stage','contract','arm','run-root','output'): p.add_argument('--'+key)
a=p.parse_args(); selection=json.loads((pathlib.Path(a.run_root)/a.arm/'selected.json').read_text())
assert len(selection['repeat_case_outputs'])==2
for path in [selection['original_case_output']]+selection['repeat_case_outputs']:
 assert (pathlib.Path(path)/'case/receipt.json').is_file()
output=pathlib.Path(a.output); output.mkdir()
(output/'collection_receipt.json').write_text(json.dumps(dict(status='synthetic_export_only',native_solves=0)))
'''


class ActualControllerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        source = self.base / "reference/source/code/model/intergen_eqscale_seq_optimized/solver.py"
        source.parent.mkdir(parents=True)
        source.write_text("synthetic source identity; never imported")
        manifest = self.base / "manifest.json"
        module.create(manifest, {"files": {"code/model/intergen_eqscale_seq_optimized/solver.py": module.file_hash(source)}})
        self.runner = self.base / "runner.py"
        self.runner.write_text(RUNNER)
        self.collector = self.base / "collector.py"
        self.collector.write_text(COLLECTOR)
        self.contract = dict(candidate_failure_policy=module.RECOVERY_MODE,
            reference_root=str(self.base / "reference"),
            parent_source_inventory=dict(path=str(manifest), sha256=module.file_hash(manifest)),
            files={"run_e5f_utility_comparison.py": dict(path=str(self.runner), sha256=module.file_hash(self.runner))},
            arms={"floor_linear": {"objective": {"sha256": "a" * 64}}})
        self.contract_path = self.base / "contract.json"
        module.create(self.contract_path, self.contract)
        self.fingerprint = module.file_hash(self.contract_path)
        self.root = self.base / "run"
        self.output = self.root / "floor_linear"
        self.output.mkdir(parents=True)
        now = time.time()
        self.clock = dict(start_epoch=now - 1, search_cutoff_epoch=now + 20,
                          repeat_cutoff_epoch=now + 25, end_epoch=now + 30,
                          contract_sha256=self.fingerprint)
        module.create(self.root / "clock.json", self.clock)

    def tearDown(self):
        self.tmp.cleanup()

    def controller(self, *, resume=False, cap=3.):
        c = module.Controller.__new__(module.Controller)
        c.contract, c.contract_path, c.fingerprint = self.contract, self.contract_path, self.fingerprint
        c.arm, c.root, c.output = "floor_linear", self.root, self.output
        c.runner, c.collector = str(self.runner), str(self.collector)
        c.recovery_mode, c.resuming = True, resume
        c.clock = self.clock
        c.budget = dict(objective_timeout_seconds=cap, workers_per_arm=2)
        c.env = os.environ.copy()
        c.phase, c.last_heartbeat = "initial", 0.
        c.interrupted, c.stop_reason, c.collection = False, None, None
        c.records, c.attempts, c.reused, c.completed_by_id = [], [], [], {}
        c.best, c.smoke_verified = None, False
        c.planned = dict(smoke=2, initial=40, de=120, repeat=2)
        if resume:
            c.restore()
        return c

    def candidates(self, *kinds, prefix="initial"):
        return [dict(id=f"{prefix}_{i}", parameters={"kind": kind}, arm="floor_linear") for i, kind in enumerate(kinds)]

    def fast_batch(self, c, candidates, phase="initial", workers=1):
        # Exercise Controller.batch/launch/finish, changing only poll latency.
        from unittest.mock import patch
        original = module.run_batch
        with patch.object(module, "run_batch", side_effect=lambda *a, **kw: original(*a, **dict(kw, poll_seconds=.01))):
            return c.batch(candidates, phase, c.clock["search_cutoff_epoch"], workers)

    def test_actual_rejection_continues_next_distinct_candidate(self):
        c = self.controller()
        result = self.fast_batch(c, self.candidates(3, 1))
        self.assertEqual([r["status"] for r in result["results"]], ["inadmissible", "success"])
        self.assertEqual(result["results"][0]["inadmissible_kind"], "normalization_or_equilibrium_gate")
        self.assertEqual(len(c.attempts), 2)
        self.assertEqual(c.best["point"], {"kind": 1})

    def test_actual_native_evidence_rejection_continues_distinct_slot(self):
        c = self.controller()
        result = self.fast_batch(c, self.candidates(8, 1))
        self.assertEqual(result["results"][0]["status"], "inadmissible")
        self.assertEqual(result["results"][0]["inadmissible_kind"], "reviewed_native_age34_evaluation_rejection")
        self.assertEqual(result["results"][1]["status"], "success")

    def test_native_source_mismatch_remains_fatal(self):
        c = self.controller()
        result = self.fast_batch(c, self.candidates(9, 1))
        self.assertEqual(result["results"][0]["status"], "failed")
        self.assertEqual(len(c.attempts), 1)

    def test_actual_fatal_halts_new_dispatch_preserving_sibling(self):
        c = self.controller()
        result = self.fast_batch(c, self.candidates(2, 6, 1), workers=2)
        self.assertEqual({r["candidate_id"]: r["status"] for r in result["results"]},
                         {"initial_0": "failed", "initial_1": "success"})
        self.assertEqual(result["unrun_ids"], ["initial_2"])
        self.assertTrue((self.root / "barrier/floor_linear.fatal.json").exists())

    def test_finish_validation_exception_preserves_running_sibling(self):
        c = self.controller()
        original_finish = c.finish
        def corrupted_plan(candidate, process, code):
            if candidate["id"] == "initial_0":
                plan_path = self.output / "plans/initial_0.json"
                plan = module.read(plan_path)
                plan["point"] = {"kind": 99}
                module.write(plan_path, plan)
            return original_finish(candidate, process, code)
        c.finish = corrupted_plan
        result = self.fast_batch(c, self.candidates(1, 6, 1), workers=2)
        self.assertEqual({r["candidate_id"]: r["status"] for r in result["results"]},
                         {"initial_0": "failed", "initial_1": "success"})
        self.assertEqual(result["unrun_ids"], ["initial_2"])
        self.assertEqual(c.summary()["counts"]["initial"]["incomplete"], 1)

    def test_forged_status_or_wrong_context_stays_fatal(self):
        c = self.controller()
        result = self.fast_batch(c, self.candidates(7, 1))
        self.assertEqual(result["results"][0]["status"], "failed")
        self.assertEqual(len(c.attempts), 1)

    def test_owned_timeout_is_local_unscored_and_next_slot_runs(self):
        c = self.controller(cap=1.)
        result = self.fast_batch(c, self.candidates(4, 1))
        first = result["results"][0]
        self.assertEqual(first["status"], "censored_timeout")
        self.assertTrue(first["deadline_evidence"]["owned_sigkill_reaped"])
        self.assertEqual(result["results"][1]["status"], "success")
        self.assertEqual(c.best["point"], {"kind": 1})

    def test_owned_startup_timeout_before_sidecar_write_is_censored(self):
        c = self.controller(cap=.1)
        result = self.fast_batch(c, self.candidates(10))
        record = result["results"][0]
        self.assertEqual(record["status"], "censored_timeout")
        self.assertEqual(record["attempt_provenance_status"], "not_written_before_owned_deadline")
        self.assertFalse((self.output / "initial_0/attempt_provenance.json").exists())
        self.assertIsNone(c.best)

    def test_owned_timeout_does_not_excuse_conflicting_or_malformed_evidence(self):
        for index, (name, content) in enumerate((
                ("attempt_provenance.json", '{"wrong":"context"}'),
                ("attempt_provenance.json", '{"partial":'),
                ("failure.json", '{"error_type":"SystemError","error":"unknown"}'))):
            with self.subTest(name=name, content=content):
                c = self.controller(cap=.1)
                candidate = dict(id=f"bad_sidecar_{index}", parameters={"kind": 10}, controller_stage="initial")
                process = c.launch(candidate, self.clock["search_cutoff_epoch"])
                try:
                    while process.poll() is None:
                        time.sleep(.01)
                    output = self.output / candidate["id"]
                    output.mkdir(exist_ok=True)
                    (output / name).write_text(content)
                    record = c.finish(candidate, process, process.process.returncode)
                    self.assertTrue(record["deadline_evidence"]["owned_sigkill_reaped"])
                    self.assertEqual(record["status"], "failed")
                finally:
                    process.close()

    def test_late_unknown_exit_is_fatal_not_timeout(self):
        c = self.controller(cap=.1)
        candidate = self.candidates(5)[0]
        candidate["controller_stage"] = "initial"
        process = c.launch(candidate, self.clock["search_cutoff_epoch"])
        try:
            process.process.wait(timeout=2)
            time.sleep(max(0., process.deadline - time.time()) + .02)
            code = process.poll()
            record = c.finish(candidate, process, code)
            self.assertEqual(record["status"], "failed")
            self.assertFalse(record["deadline_evidence"]["owned_sigkill_reaped"])
        finally:
            process.close()

    def test_resume_reuses_completed_case_and_original_clock_no_replay(self):
        c = self.controller()
        points = self.candidates(1, 6)
        self.fast_batch(c, points[:1])
        clock_hash = module.file_hash(self.root / "clock.json")
        resumed = self.controller(resume=True)
        result = self.fast_batch(resumed, points)
        self.assertEqual(len(resumed.attempts), 2)
        self.assertEqual(len(result["results"]), 2)
        self.assertEqual(module.file_hash(self.root / "clock.json"), clock_hash)
        original = next(r for r in result["results"] if r["candidate_id"] == "initial_0")
        self.assertEqual(original, c.records[0])

    def test_resume_rejects_orphan_dispatch_without_retry(self):
        c = self.controller()
        self.fast_batch(c, self.candidates(1))
        module.create(self.output / "dispatch/orphan.json", {"candidate_id": "orphan"})
        with self.assertRaisesRegex(RuntimeError, "orphan"):
            self.controller(resume=True)

    def test_resume_rejects_changed_scientific_artifact(self):
        c = self.controller()
        self.fast_batch(c, self.candidates(1))
        (self.output / "initial_0/case/target_fit.csv").write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "artifact changed"):
            self.controller(resume=True)

    def test_retained_selection_and_repeats_reach_export_without_replay(self):
        c = self.controller()
        self.fast_batch(c, self.candidates(1))
        c.smoke_verified = True
        c.stop_reason = "finite_search_bank_complete"
        c.repeat_and_collect()
        self.assertEqual(len(c.attempts), 3)
        self.assertEqual(c.collection["status"], "success")
        export_hash = module.file_hash(self.output / "export/collection_receipt.json")
        c.checkpoint()
        resumed = self.controller(resume=True)
        self.assertEqual(resumed.stop_reason, "finite_search_bank_complete")
        resumed.smoke_verified = True
        resumed.external = Mock(wraps=resumed.external)
        resumed.repeat_and_collect()
        self.assertEqual(len(resumed.attempts), 3)
        resumed.external.assert_called_once()
        self.assertEqual(resumed.external.call_args.args[1], "collection")
        self.assertEqual(module.file_hash(self.output / "export/collection_receipt.json"), export_hash)

    def test_native_payload_validation_rejects_malformed_or_owner_census(self):
        self.assertFalse(module.narrow_payload_valid({"stage": "forward_age_34"}, module.recovery_policy.CENSUS_FIELDS))

    def test_actual_loop_interruption_is_fatal_not_local_timeout(self):
        c = self.controller(cap=3.)
        candidate = self.candidates(4)[0]
        candidate["controller_stage"] = "initial"
        started = time.monotonic()
        result = module.run_batch([candidate], workers=1, deadline=self.clock["search_cutoff_epoch"],
            launch=c.launch, finish=c.finish, heartbeat=c.heartbeat,
            interrupted=lambda: time.monotonic() - started > .1, poll_seconds=.01,
            allowed_statuses=module.COMPLETED_STATUSES | module.CENSORED_STATUSES)
        record = result["results"][0]
        self.assertEqual(record["status"], "failed")
        self.assertTrue(record["deadline_evidence"]["cancelled"])
        self.assertFalse(record["deadline_evidence"]["owned_sigkill_reaped"])
        plan = module.read(self.output / "plans/initial_0.json")
        self.assertEqual(record["deadline_evidence"]["deadline_epoch"], plan["deadline_epoch"])

    def test_late_zero_required_repeat_cannot_trigger_verified_export(self):
        c = self.controller(cap=3.)
        initial = self.fast_batch(c, self.candidates(1))
        self.assertTrue(initial["complete"])
        self.assertEqual(initial["results"][0]["status"], "success")
        self.assertIsNotNone(c.best)
        c.smoke_verified = True
        candidate = dict(id="repeat_0", parameters={"kind": 1}, controller_stage="repeat", arm="floor_linear")
        process = c.launch(candidate, self.clock["repeat_cutoff_epoch"])
        try:
            self.assertEqual(process.process.wait(timeout=5), 0)
            time.sleep(max(0., process.deadline - time.time()) + .02)
            code = process.poll()
            record = c.finish(candidate, process, code)
            self.assertEqual(record["status"], "censored_late_completion")
            self.assertTrue((self.output / "repeat_0/case/receipt.json").is_file())
        finally:
            process.close()
        c.external = Mock()
        c.repeat_and_collect()
        c.external.assert_not_called()
        self.assertEqual(c.collection["status"], "incomplete_repetitions")
        self.assertFalse(c.collection["exact_repeat_claim"])


if __name__ == "__main__":
    unittest.main()
