"""Lightweight fault injection for an unintegrated recovery-policy proposal."""
import copy
from dataclasses import replace
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

import e5f_utility_recovery_policy_v1 as policy


class InfeasibleThetaError(RuntimeError):
    def __init__(self, stage="forward_age_34", mass=1.2e-12, rows=None):
        self.stage, self.dead_mass = stage, mass
        self.census = rows if rows is not None else [dict(age=34., b=0., child_state=1,
            income=.22, location=0, mass=mass, parity=1, slack=-.024, tenure=0,
            transfer=0., unsecured_position=0., z=.054)]
        if rows is None:
            base = self.census[0]
            self.census = [dict(base, b=-.116, unsecured_position=-.116, mass=mass * .1),
                           dict(base, mass=mass * .7), dict(base, parity=2, mass=mass * .1),
                           dict(base, parity=2, child_state=2, mass=mass * .1)]
        super().__init__("raw native diagnostic retained without text parsing")


def context(name="case", stage="initial"):
    return policy.Context(name, stage, "a" * 64, "b" * 64, "c" * 64, "d" * 64)


def evidence(exc=None, ctx=None):
    return policy.capture_native_failure(exc or InfeasibleThetaError(),
        expected_native_type=InfeasibleThetaError, native_gate_tolerance=1e-12, context=ctx or context())


def execution(ctx=None, expired=False, code=1):
    return policy.Execution(ctx or context(), 123, 456, 100., 90., 101. if expired else 95.,
                            expired, -9 if expired else code, "/synthetic/raw.log", expired, expired)


class PolicyTests(unittest.TestCase):
    def test_manifest_cannot_authorize_launch(self):
        self.assertFalse(policy.MANIFEST["launch_permitted"])
        self.assertFalse(policy.MANIFEST["integrated"])

    def test_valid_infeasibility_rejected_with_raw_provenance(self):
        result = policy.classify_completion(execution(), failure=evidence())
        self.assertEqual(result["status"], "rejected_infeasible")
        self.assertFalse(result["halt_new_dispatch"])
        self.assertFalse(result["score_eligible"])
        self.assertIsNone(result["score"])
        self.assertEqual(result["failure"]["context"], policy.asdict(context()))
        self.assertIn("raw native", result["failure"]["raw_error"])

    def test_both_observed_complete_censuses_reject_only_the_evaluation(self):
        # Structured transcription of run_001's two independently reviewed
        # floor initial_0005 receipts; no production parser trusts their prose.
        fixtures = (
            (1.19443474042e-12, (2.291315006454442e-52, 1.1944347404170972e-12,
                                2.0310845754475758e-24, 9.023225940780418e-25)),
            (1.19378693038e-12, (2.280040844535949e-52, 1.1937869303791421e-12,
                                2.0278724402322508e-24, 9.0089558200943e-25)),
        )
        for total, masses in fixtures:
            exc = InfeasibleThetaError(mass=total)
            for index, row in enumerate(exc.census):
                row.update(mass=masses[index], income=.219771479538767, z=.05388932006414422,
                           b=-.11627906976744207 if index == 0 else 0.,
                           unsecured_position=-.11627906976744207 if index == 0 else 0.,
                           slack=-.03350248616398435 if index == 0 else -.0239173512802634)
            result = policy.classify_completion(execution(), failure=evidence(exc))
            self.assertEqual(result["status"], "rejected_infeasible")
            self.assertFalse(result["score_eligible"])

    def test_arbitrary_exception_text_and_lookalike_class_are_not_trusted(self):
        fake = RuntimeError("forward_age_34: dead-node mass 2e-12 exceeds 1e-12")
        self.assertFalse(evidence(fake).narrow_infeasibility_verified)
        Lookalike = type("InfeasibleThetaError", (RuntimeError,), {})
        self.assertFalse(evidence(Lookalike(str(fake))).narrow_infeasibility_verified)

    def test_other_stage_gate_or_bad_census_stays_fatal(self):
        failures = [InfeasibleThetaError(stage="forward_age_38"), InfeasibleThetaError(mass=1e-12),
                    InfeasibleThetaError(mass=float("inf")), InfeasibleThetaError(rows=[])]
        for key, value in (("slack", 0.), ("slack", float("nan")), ("mass", -1.),
                           ("age", 38.), ("child_state", 3.), ("tenure", 1.)):
            exc = InfeasibleThetaError()
            exc.census[0][key] = value
            failures.append(exc)
        for exc in failures:
            with self.subTest(payload=repr(exc.census)):
                self.assertEqual(policy.classify_completion(execution(), failure=evidence(exc))["status"], "fatal")

    def test_census_schema_and_mass_accounting_are_strict(self):
        exc = InfeasibleThetaError()
        exc.census[0]["extra"] = 1.
        self.assertFalse(evidence(exc).narrow_infeasibility_verified)
        exc = InfeasibleThetaError()
        exc.census *= 2
        self.assertFalse(evidence(exc).narrow_infeasibility_verified)
        exc = InfeasibleThetaError()
        exc.census[0]["mass"] *= .5
        self.assertFalse(evidence(exc).narrow_infeasibility_verified)
        self.assertFalse(policy.capture_native_failure(InfeasibleThetaError(),
            expected_native_type=InfeasibleThetaError, native_gate_tolerance=1e-10,
            context=context()).narrow_infeasibility_verified)

    def test_timeout_distinct_even_when_wrapped_systemerror(self):
        raw = evidence(SystemError("CPUDispatcher returned exception set"))
        result = policy.classify_completion(execution(expired=True), failure=raw)
        self.assertEqual(result["status"], "censored_timeout")
        self.assertFalse(result["score_eligible"])
        self.assertEqual(result["failure"]["raw_type"], "SystemError")
        self.assertEqual(policy.classify_completion(execution(), failure=raw)["status"], "fatal")

    def test_smoke_and_repeat_timeouts_and_rejections_halt(self):
        for stage in ("smoke", "repeat"):
            ctx = context(stage=stage)
            self.assertTrue(policy.classify_completion(execution(ctx, True))["halt_new_dispatch"])
            self.assertTrue(policy.classify_completion(execution(ctx), failure=evidence(ctx=ctx))["halt_new_dispatch"])

    def test_provenance_and_integrity_errors_override_local_recovery(self):
        self.assertEqual(policy.classify_completion(execution(), failure=evidence(ctx=context("other")))["status"], "fatal")
        for error in ("source_hash", "target_hash", "accounting_gate"):
            self.assertEqual(policy.classify_completion(execution(expired=True), integrity_error=error)["status"], "fatal")

    def test_no_success_at_or_after_deadline_or_without_verification(self):
        self.assertEqual(policy.classify_completion(execution(code=0), result_verified=True)["status"], "success")
        self.assertEqual(policy.classify_completion(execution(expired=True, code=0), result_verified=True)["status"], "censored_timeout")
        self.assertEqual(policy.classify_completion(execution(code=0))["status"], "fatal")


class ProcessTests(unittest.TestCase):
    def case(self, name, code="pass", cap=2.):
        return dict(id=name, context=context(name), command=[sys.executable, "-c", code], cap_seconds=cap)

    def run_cases(self, cases, verifier=None, workers=1, seconds=5.):
        with tempfile.TemporaryDirectory() as directory:
            return policy.run_synthetic_schedule(cases, workers=workers, cutoff_epoch=time.time() + seconds,
                directory=directory, verifier=verifier or (lambda c, e: dict(result_verified=e.returncode == 0)))

    def test_local_rejection_and_timeout_continue_distinct_slots_without_retry(self):
        def verifier(case, result):
            if case["id"] == "rejected":
                return dict(failure=evidence(ctx=case["context"]))
            return dict(result_verified=result.returncode == 0)
        result = self.run_cases([self.case("rejected", "raise SystemExit(1)"),
            self.case("timeout", "import time; time.sleep(2)", .08), self.case("next")], verifier)
        self.assertEqual(result["launched"], ["rejected", "timeout", "next"])
        self.assertEqual([row["status"] for row in result["completed"]], ["rejected_infeasible", "censored_timeout", "success"])
        self.assertEqual((result["retry_count"], result["replacement_slots"]), (0, 0))

    def test_fatal_stop_preserves_already_running_sibling(self):
        result = self.run_cases([self.case("fatal", "raise SystemExit(1)"),
            self.case("sibling", "import time; time.sleep(.15)"), self.case("unrun")], workers=2)
        self.assertEqual(result["launched"], ["fatal", "sibling"])
        self.assertEqual(result["unrun"], ["unrun"])
        self.assertEqual({row["candidate_id"]: row["status"] for row in result["completed"]},
                         {"fatal": "fatal", "sibling": "success"})

    def test_absolute_cutoff_stops_dispatch_and_bounds_running_child(self):
        result = self.run_cases([self.case("running", "import time; time.sleep(2)"), self.case("unrun")], seconds=.08)
        self.assertEqual(result["unrun"], ["unrun"])
        self.assertEqual(result["completed"][0]["status"], "censored_timeout")

    def test_verifier_exception_halts_new_work_preserving_sibling(self):
        def verifier(case, result):
            if case["id"] == "bad_provenance":
                raise ValueError("source fingerprint mismatch")
            return dict(result_verified=True)
        result = self.run_cases([self.case("bad_provenance"),
            self.case("sibling", "import time; time.sleep(.15)"), self.case("unrun")], verifier, workers=2)
        self.assertEqual(result["unrun"], ["unrun"])
        self.assertEqual({row["candidate_id"]: row["status"] for row in result["completed"]},
                         {"bad_provenance": "fatal", "sibling": "success"})

    def test_already_exited_zero_is_not_success_when_observed_after_deadline(self):
        with tempfile.TemporaryDirectory() as directory:
            process = policy.OwnedProcess([sys.executable, "-c", "pass"], context=context(),
                log_path=Path(directory) / "raw.log", cutoff_epoch=time.time() + 5, cap_seconds=.1)
            try:
                time.sleep(.16)
                outcome = process.poll()
                self.assertTrue(outcome.deadline_expired)
                self.assertEqual(policy.classify_completion(outcome, result_verified=True)["status"], "censored_late_completion")
            finally:
                process.close()

    def test_unknown_systemerror_exit_before_deadline_observed_late_stays_fatal(self):
        with tempfile.TemporaryDirectory() as directory:
            process = policy.OwnedProcess([sys.executable, "-c", "raise SystemError('unknown kernel bug')"],
                context=context(), log_path=Path(directory) / "raw.log",
                cutoff_epoch=time.time() + 5, cap_seconds=.1)
            try:
                process.process.wait(timeout=2)
                time.sleep(max(0., process.deadline_epoch - time.time()) + .02)
                outcome = process.poll()
                self.assertFalse(outcome.deadline_kill_reaped)
                self.assertEqual(policy.classify_completion(outcome)["status"], "fatal")
            finally:
                process.close()

    def test_duplicate_slots_rejected_before_execution(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.run_cases([self.case("same"), self.case("same")])


if __name__ == "__main__":
    unittest.main()
