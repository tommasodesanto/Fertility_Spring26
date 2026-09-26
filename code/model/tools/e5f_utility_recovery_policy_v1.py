"""Narrow native capture plus a review-only process-policy harness.

The capture predicate is used by the opt-in reviewed_failure_v1 runner and
controller. The standalone OwnedProcess/synthetic schedule below remain test
harnesses, not production entry points. New reviewed source/contract pins are
required; this module grants no launch permission. The caller authenticates the
native class's defining source and captures its attributes directly. Historical
error strings cannot authorize a new native feasibility rejection. Supervised
runner mode enforces wall time outside compiled kernels without an inner alarm.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import time


MANIFEST = dict(schema="utility_recovery_policy_v1", launch_permitted=False,
                integrated=False, requires_new_reviewed_contract=True,
                native_capture_integration="opt_in_reviewed_failure_v1",
                unchanged_dead_mass_gate=1e-12,
                existing_inadmissibility_barrier="unchanged; integration decision pending")
SEARCH_STAGES = {"initial", "de"}
CENSUS_FIELDS = {"age", "b", "child_state", "income", "location", "mass", "parity",
                 "slack", "tenure", "transfer", "unsecured_position", "z"}


def finite(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


@dataclass(frozen=True)
class Context:
    candidate_id: str
    stage: str
    contract_sha256: str
    source_sha256: str
    target_sha256: str
    point_sha256: str

    def __post_init__(self):
        if not self.candidate_id or self.stage not in {"initial", "de", "smoke", "repeat"}:
            raise ValueError("candidate identity/stage required")
        for name in ("contract_sha256", "source_sha256", "target_sha256", "point_sha256"):
            if re.fullmatch(r"[0-9a-f]{64}", getattr(self, name)) is None:
                raise ValueError("full provenance fingerprints required")


@dataclass(frozen=True)
class FailureEvidence:
    context: Context
    raw_type: str
    raw_error: str
    native_payload: dict
    narrow_infeasibility_verified: bool
    validation: str


def capture_native_failure(exc, *, expected_native_type, native_gate_tolerance, context):
    """Use exact native type identity and structured attributes, never regex prose.

    expected_native_type must come from the caller's source-pinned runtime, not
    from an error receipt. This prototype intentionally imports no native class.
    """
    payload = {}
    valid, reason = False, "exception_type_not_the_authenticated_native_class"
    if type(exc) is expected_native_type and expected_native_type.__name__ == "InfeasibleThetaError":
        try:
            payload = dict(stage=exc.stage, dead_mass=exc.dead_mass,
                           native_gate_tolerance=native_gate_tolerance,
                           census=[dict(row) for row in exc.census])
            mass, rows = payload["dead_mass"], payload["census"]
            if native_gate_tolerance != 1e-12:
                raise ValueError("native_gate_changed")
            if payload["stage"] != "forward_age_34":
                raise ValueError("other_forward_stage")
            if not finite(mass) or not 1e-12 < mass <= 1.:
                raise ValueError("mass_does_not_exceed_unchanged_gate")
            if len(rows) != 4:
                raise ValueError("census_not_the_reviewed_complete_four_rows")
            for row in rows:
                if set(row) != CENSUS_FIELDS or not all(finite(x) for x in row.values()):
                    raise ValueError("malformed_or_nonfinite_census")
                if row["age"] != 34. or row["slack"] >= 0. or not 0. < row["mass"] <= mass:
                    raise ValueError("census_not_age34_positive_mass_negative_slack")
                for key in ("child_state", "location", "parity", "tenure"):
                    if row[key] < 0 or int(row[key]) != row[key]:
                        raise ValueError("invalid_discrete_state")
                if row["child_state"] > row["parity"] or row["income"] < 0. or row["z"] <= 0.:
                    raise ValueError("invalid_family_or_income_state")
                if (row["location"] != 0 or row["tenure"] != 0 or row["transfer"] != 0.
                        or row["b"] > 0. or row["unsecured_position"] != row["b"]):
                    raise ValueError("unreviewed_location_tenure_or_budget_state")
            states = sorted((row["b"] < 0., row["parity"], row["child_state"]) for row in rows)
            if states != [(False, 1, 1), (False, 2, 1), (False, 2, 2), (True, 1, 1)]:
                raise ValueError("unreviewed_or_duplicate_census_state")
            # Require the complete census to agree at the native .12g precision;
            # a truncated example of negative-slack nodes is insufficient.
            if format(math.fsum(row["mass"] for row in rows), ".12g") != format(mass, ".12g"):
                raise ValueError("census_mass_disagrees_with_native_dead_mass")
            valid, reason = True, "exact_native_type_and_structured_gate_failure"
        except (AttributeError, TypeError, ValueError) as error:
            reason = str(error)
    return FailureEvidence(context, type(exc).__name__, str(exc), payload, valid, reason)


@dataclass(frozen=True)
class Execution:
    context: Context
    supervisor_pid: int
    child_pid: int
    deadline_epoch: float
    started_epoch: float
    observed_epoch: float
    deadline_expired: bool
    returncode: int
    raw_log_path: str
    observed_running_at_expiry: bool = False
    deadline_kill_reaped: bool = False


def classify_completion(execution, *, failure=None, result_verified=False, integrity_error=None):
    """Timeout evidence comes from the process owner, not error text or class.

    A wrapped SystemError alone remains fatal. Explicit source/target/accounting
    failures take precedence even when the deadline also expired. This API must
    receive Execution directly from OwnedProcess, not deserialize arbitrary JSON.
    """
    context = execution.context
    status, reason = "fatal", "unverified_or_unknown_failure"
    if integrity_error is not None:
        reason = "integrity_failure:" + str(integrity_error)
    elif failure is not None and failure.context != context:
        reason = "failure_provenance_mismatch"
    elif (execution.deadline_expired and execution.observed_running_at_expiry
          and execution.deadline_kill_reaped and execution.returncode == -signal.SIGKILL):
        status, reason = "censored_timeout", "supervisor_owned_wall_deadline"
    elif execution.deadline_expired or execution.observed_epoch >= execution.deadline_epoch:
        if execution.returncode == 0 and failure is None:
            status, reason = "censored_late_completion", "completion_first_observed_after_deadline"
        else:
            reason = "late_nonzero_exit_without_positive_timeout_provenance"
    elif execution.returncode == 0 and result_verified and failure is None:
        status, reason = "success", "verified_completion_before_deadline"
    elif (execution.returncode != 0 and failure is not None
          and failure.narrow_infeasibility_verified):
        status, reason = "rejected_infeasible", "unchanged_dead_mass_gate"
    local = context.stage in SEARCH_STAGES
    halt = status == "fatal" or (status in {"censored_timeout", "censored_late_completion", "rejected_infeasible"} and not local)
    return dict(status=status, reason=reason, halt_new_dispatch=halt,
                score_eligible=status == "success", score=None, automatic_retry=False,
                execution=asdict(execution), failure=asdict(failure) if failure else None,
                policy=MANIFEST)


class OwnedProcess:
    """External watchdog; existing exit is distinguished from a verified kill."""
    def __init__(self, command, *, context, log_path, cutoff_epoch, cap_seconds):
        if not finite(cutoff_epoch) or not finite(cap_seconds) or cap_seconds <= 0:
            raise ValueError("finite positive deadline inputs required")
        self.context, self.log_path = context, str(log_path)
        self.started_epoch = time.time()
        self.deadline_epoch = min(cutoff_epoch, self.started_epoch + cap_seconds)
        remaining = self.deadline_epoch - self.started_epoch
        if remaining <= 0:
            raise TimeoutError("no dispatch after the absolute cutoff")
        self.monotonic_deadline = time.monotonic() + remaining
        self.log = Path(log_path).open("xb")
        self.execution = None
        try:
            self.process = subprocess.Popen(command, stdout=self.log, stderr=subprocess.STDOUT,
                                            start_new_session=True)
        except BaseException:
            self.log.close()
            raise

    def close(self):
        sent = False
        try:
            os.killpg(self.process.pid, signal.SIGKILL)
            sent = True
        except ProcessLookupError:
            pass
        self.process.wait(timeout=5)
        self.log.close()
        return sent and self.process.returncode == -signal.SIGKILL

    def poll(self):
        if self.execution is not None:
            return self.execution
        observed = time.time()
        code = self.process.poll()
        expired = time.monotonic() >= self.monotonic_deadline or observed >= self.deadline_epoch
        running_at_expiry, killed_reaped = expired and code is None, False
        if running_at_expiry:
            killed_reaped = self.close()
            code = self.process.returncode
        if code is None:
            return None
        self.execution = Execution(self.context, os.getpid(), self.process.pid, self.deadline_epoch,
                                   self.started_epoch, observed, expired, code, self.log_path,
                                   running_at_expiry, killed_reaped)
        return self.execution


def run_synthetic_schedule(cases, *, workers, cutoff_epoch, directory, verifier,
                           poll_seconds=.01):
    """Fault-injection harness only: finite unique slots, no optimizer or model.

    Each case supplies id, context, command and cap_seconds. verifier(case,
    execution) returns classification kwargs. Rejection/timeout never supplies a
    score; the existing production >50% barrier is deliberately not implemented.
    """
    if workers < 1 or not finite(cutoff_epoch) or poll_seconds <= 0:
        raise ValueError("invalid scheduler bounds")
    identifiers = [row["id"] for row in cases]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("duplicate candidate IDs would retry a slot")
    if any(row["id"] != row["context"].candidate_id for row in cases):
        raise ValueError("case identity differs from its provenance")
    pending, active, completed, launched = list(cases), {}, [], []
    halt = None
    try:
        while pending or active:
            for key, (case, process) in list(active.items()):
                execution = process.poll()
                if execution is None:
                    continue
                del active[key]
                try:
                    try:
                        result = classify_completion(execution, **verifier(case, execution))
                    except Exception as error:
                        result = classify_completion(execution,
                            integrity_error=type(error).__name__ + ":" + str(error))
                finally:
                    process.close()
                completed.append(dict(candidate_id=key, **result))
                if result["halt_new_dispatch"]:
                    halt = halt or "fatal_case_no_new_dispatch"
            if time.time() >= cutoff_epoch:
                halt = halt or "absolute_cutoff"
            while pending and len(active) < workers and halt is None:
                if time.time() >= cutoff_epoch:
                    halt = "absolute_cutoff"
                    break
                case = pending.pop(0)
                process = OwnedProcess(case["command"], context=case["context"],
                    log_path=Path(directory) / (case["id"] + ".log"),
                    cutoff_epoch=cutoff_epoch, cap_seconds=case["cap_seconds"])
                active[case["id"]] = (case, process)
                launched.append(case["id"])
            if not active and (halt or not pending):
                break
            time.sleep(poll_seconds)
    finally:
        for _, process in active.values():
            process.close()
    return dict(launched=launched, completed=completed, unrun=[row["id"] for row in pending],
                halt=halt, retry_count=0, replacement_slots=0, policy=MANIFEST)
