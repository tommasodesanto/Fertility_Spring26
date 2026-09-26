#!/usr/bin/env python3
"""Torch-only bounded four-arm controller; never grants launch approval.

One Slurm array member owns one arm. Four ready members establish one absolute
clock immediately before the first smoke. The same subprocess scheduler runs
both exact smokes, the finite search, and both selected repetitions. No failed
case is retried. All numerical work, collection, and deadline enforcement remain
on Torch when the submitting laptop sleeps.
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import e5f_utility_comparison_design as design
import e5f_utility_recovery_policy_v1 as recovery_policy


ARMS = design.ARM_NAMES
THREAD_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
HEARTBEAT_SECONDS = 30.
READINESS_SECONDS = 1800.
COLLECTOR_NAME = "collect_e5f_utility_comparison.py"
COMPLETED_STATUSES = {"success", "inadmissible"}
CENSORED_STATUSES = {"censored_timeout", "censored_late_completion"}
RECOVERY_NAME = "e5f_utility_recovery_policy_v1.py"
RECOVERY_MODE = "reviewed_failure_v1"
OLD_INADMISSIBLE_PREFIXES = (
    "Old-steady-state fertility normalization is not bracketed:",
    "Old-steady-state fertility normalization missed tolerance:",
    "Initial housing equilibrium failed its unchanged strict gate",
)


class BarrierStop(RuntimeError):
    """A completed finite generation crossed a prespecified stop threshold."""


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def create(path, value):
    """Reserve a new artifact without overwriting a previous run or approval."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".new.{os.getpid()}")
    try:
        with temporary.open("x") as stream:
            json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write("\n")
        # Linking publishes complete content atomically and refuses overwrite.
        # This matters when the other Slurm array members inspect the barrier.
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_contract(path):
    """Check all pinned preparation files before any native model subprocess."""
    path = Path(path).resolve()
    contract = read(path)
    fingerprint = file_hash(path)
    if os.environ.get("EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256") != fingerprint:
        raise RuntimeError("explicit approved contract SHA256 is required")
    if (contract.get("status") != "author_approved_frozen_design"
            or contract.get("launch_permitted") is not True
            or not contract.get("approved_assumptions")):
        raise RuntimeError("frozen design and common assumptions are not approved")
    required = (Path(__file__).name, Path(design.__file__).name,
                "run_e5f_utility_comparison.py", COLLECTOR_NAME,
                RECOVERY_NAME, "budget", "arm_catalog", "proposal_bank")
    for key in required:
        if key not in contract["files"]:
            raise RuntimeError("missing required source/design pin: " + key)
    for pin in contract["files"].values():
        if file_hash(pin["path"]) != pin["sha256"]:
            raise RuntimeError("pinned source/input changed: " + pin["path"])
    for actual in (__file__, design.__file__, recovery_policy.__file__):
        if Path(actual).resolve() != Path(contract["files"][Path(actual).name]["path"]).resolve():
            raise RuntimeError("executing source is not the pinned file: " + actual)
    mode = contract.get("candidate_failure_policy")
    if mode not in (None, RECOVERY_MODE):
        raise RuntimeError("unknown candidate failure policy")
    if mode == RECOVERY_MODE:
        if RECOVERY_NAME not in contract["files"]:
            raise RuntimeError("reviewed failure policy requires its own source pin and new contract")
        module = importlib.import_module(RECOVERY_NAME[:-3])
        if Path(module.__file__).resolve() != Path(contract["files"][RECOVERY_NAME]["path"]).resolve():
            raise RuntimeError("executing recovery policy differs from the new contract")
    budget = read(contract["files"]["budget"]["path"])
    design.validate_budget(budget)
    if contract.get("approved_budget") != budget:
        raise RuntimeError("approved_budget must contain the complete pinned budget")
    expected = {"arm_count": 4, "workers_per_arm": 10, "total_limit_seconds": 28800.,
                "repeat_reserve_seconds": 4500., "export_reserve_seconds": 900.,
                "objective_timeout_seconds": 3100., "initial_population": 40,
                "de_generations": 3, "identical_smokes_per_arm": 2,
                "selected_repeats_per_arm": 2}
    for key, value in expected.items():
        if budget[key] != value:
            raise RuntimeError("controller frozen schedule differs: " + key)
    catalog = read(contract["files"]["arm_catalog"]["path"])
    bank = read(contract["files"]["proposal_bank"]["path"])
    if bank["budget_sha256"] != design.canonical_fingerprint(budget):
        raise RuntimeError("proposal bank budget fingerprint changed")
    if bank["arm_catalog_sha256"] != design.canonical_fingerprint(catalog):
        raise RuntimeError("proposal bank arm catalog fingerprint changed")
    if set(contract["arms"]) != set(ARMS) or set(bank["arms"]) != set(ARMS):
        raise RuntimeError("exactly the four approved arms are required")
    de = contract["differential_evolution"]
    if de != {"generations": 3, "mutation_factor": .7, "crossover_probability": 1.}:
        raise RuntimeError("DE schedule differs from the frozen controller design")
    for arm in ARMS:
        if len(bank["arms"][arm]) != 40:
            raise RuntimeError("each arm requires the complete forty-member bank")
        if catalog[arm]["free_parameter_count"] != (8 if arm.startswith("floor") else 9):
            raise RuntimeError("arm coordinate count changed")
    return contract, fingerprint, budget, catalog, bank


class ManagedProcess:
    """One owned process group, with a hard absolute subprocess deadline."""
    def __init__(self, command, log_path, deadline, env):
        if deadline <= time.time():
            raise TimeoutError("subprocess deadline already exhausted")
        self.command, self.deadline = list(command), float(deadline)
        self.scheduled_deadline = self.deadline
        self.cancelled = False
        self.started_epoch = time.time()
        self.log = Path(log_path).open("xb")
        self.timed_out = False
        self.monotonic_deadline = time.monotonic() + self.deadline - self.started_epoch
        self.observed_running_at_expiry = False
        self.deadline_kill_reaped = False
        self.observed_epoch = None
        self.deadline_expired = False
        try:
            self.process = subprocess.Popen(self.command, stdout=self.log, stderr=subprocess.STDOUT,
                                            env=env, start_new_session=True)
        except BaseException:
            self.log.close()
            raise

    def kill(self):
        # Only a group created here with start_new_session; never a Slurm job or
        # a discovered PID. SIGKILL prevents an unbounded graceful-shutdown tail.
        sent = False
        try:
            os.killpg(self.process.pid, signal.SIGKILL)
            sent = True
        except ProcessLookupError:
            pass
        self.process.wait(timeout=10)
        return sent and self.process.returncode == -signal.SIGKILL

    def poll(self):
        observed = time.time()
        code = self.process.poll()
        expired = observed >= self.deadline or time.monotonic() >= self.monotonic_deadline
        if code is None and expired:
            self.timed_out = True
            self.observed_running_at_expiry = True
            self.deadline_kill_reaped = self.kill()
            code = self.process.returncode
        if code is not None:
            if self.observed_epoch is None:
                self.observed_epoch, self.deadline_expired = observed, expired
            self.log.close()
        return code

    def close(self):
        # Also reap a background descendant if its immediate parent exited.
        self.kill()
        self.log.close()

    def cancel(self):
        """Cancellation is fatal, distinct from exhausting the approved cap."""
        self.cancelled = True
        self.close()


def run_batch(candidates, *, workers, deadline, launch, finish, heartbeat,
              interrupted=lambda: False, poll_seconds=1., allowed_statuses=None,
              guard=lambda: None, on_finish_error=lambda candidate, error: None):
    """Run a finite barrier. Any failed result stops new dispatch, not siblings.

    Missing/unrun trials are never submitted to DE selection as rejected scores.
    Dependency injection permits synthetic process tests without model imports.
    """
    pending, active, results = list(candidates), {}, []
    reason = None
    allowed_statuses = COMPLETED_STATUSES if allowed_statuses is None else set(allowed_statuses)
    try:
        while pending or active:
            now = time.time()
            if interrupted() and reason is None:
                reason = "controller_interrupted"
            if now >= deadline and reason is None:
                reason = "stage_deadline_exhausted"
            external_stop = guard()
            if external_stop and reason is None:
                reason = external_stop
            # Reap every completed child before filling vacancies, so an already
            # observed failure cannot be followed by an unnecessary new dispatch.
            for key, (candidate, process) in list(active.items()):
                code = process.poll()
                if code is None:
                    continue
                del active[key]
                try:
                    try:
                        result = finish(candidate, process, code)
                    except Exception as exc:
                        # Preserve already-running siblings even when artifact
                        # validation itself fails. The missing completed receipt
                        # remains an explicit orphan and cannot be resumed.
                        on_finish_error(candidate, exc)
                        result = dict(candidate_id=candidate["id"], status="failed",
                            halt_new_dispatch=True, error_type=type(exc).__name__, error=str(exc))
                finally:
                    process.close()
                results.append(result)
                if (result.get("halt_new_dispatch") or result["status"] not in allowed_statuses) and reason is None:
                    reason = "case_failed_no_retry"
            if interrupted():
                for _, process in active.values():
                    process.cancel()
            while pending and len(active) < workers and reason is None:
                external_stop = guard()
                if external_stop:
                    reason = external_stop
                    break
                if time.time() >= deadline:
                    reason = "stage_deadline_exhausted"
                    break
                candidate = pending.pop(0)
                process = launch(candidate, deadline)
                active[candidate["id"]] = (candidate, process)
            heartbeat(active=len(active), pending=len(pending), stop_reason=reason)
            if not active and (reason is not None or not pending):
                break
            time.sleep(min(poll_seconds, max(.001, deadline - time.time())) if reason is None else poll_seconds)
    finally:
        for _, process in active.values():
            process.close()
    return {"results": results, "unrun_ids": [row["id"] for row in pending],
            "stop_reason": reason, "complete": len(results) == len(candidates)}


def complete_scores(candidates, results):
    """Require complete proposals and the frozen inadmissibility threshold."""
    if len(results) != len(candidates) or any(row["status"] not in COMPLETED_STATUSES for row in results):
        raise RuntimeError("incomplete or failed DE barrier; no replacement proposals")
    scores = {row["candidate_id"]: row["loss"] for row in results}
    if set(scores) != {row["id"] for row in candidates}:
        raise RuntimeError("DE result IDs differ from the finite proposal bank")
    rejected = sum(row["status"] == "inadmissible" for row in results)
    if rejected > len(results) / 2:
        raise BarrierStop("more_than_half_of_completed_barrier_inadmissible")
    if not any(row["status"] == "success" for row in results):
        raise BarrierStop("no_valid_parent_in_completed_barrier")
    return scores


def classified_inadmissible(failure, stage):
    """Accept only the pinned runner's exact whitelist classification.

    The runner classifies three explicitly frozen normalization/equilibrium gate
    messages. This consumes a distinct planned proposal, never a retry or a gate
    relaxation. Smoke and selected-repeat failures cannot use this route.
    """
    return (stage in {"initial", "de"}
            and failure.get("status") == "inadmissible_parameter_proposal"
            and failure.get("error_type") == "RuntimeError")


def narrow_payload_valid(payload, fields):
    """Independently recheck the serialized native census; trust no status flag."""
    try:
        mass, rows = payload["dead_mass"], payload["census"]
        if (set(payload) != {"stage", "dead_mass", "native_gate_tolerance", "census"}
                or payload["stage"] != "forward_age_34" or payload["native_gate_tolerance"] != 1e-12
                or isinstance(mass, bool) or not math.isfinite(mass) or not 1e-12 < mass <= 1
                or len(rows) != 4):
            return False
        for row in rows:
            if set(row) != set(fields) or any(isinstance(x, bool) or not isinstance(x, (int, float))
                                             or not math.isfinite(x) for x in row.values()):
                return False
            if (row["age"] != 34 or row["slack"] >= 0 or not 0 < row["mass"] <= mass
                    or row["location"] != 0 or row["tenure"] != 0 or row["transfer"] != 0
                    or row["b"] > 0 or row["unsecured_position"] != row["b"]
                    or row["income"] < 0 or row["z"] <= 0):
                return False
        states = sorted((r["b"] < 0, r["parity"], r["child_state"]) for r in rows)
        return (states == [(False, 1, 1), (False, 2, 1), (False, 2, 2), (True, 1, 1)]
                and format(math.fsum(r["mass"] for r in rows), ".12g") == format(mass, ".12g"))
    except (KeyError, TypeError, ValueError, OverflowError):
        return False


class Controller:
    def __init__(self, args):
        self.args = args
        self.contract_path = args.contract.resolve()
        (self.contract, self.fingerprint, self.budget,
         self.catalog, self.bank) = load_contract(self.contract_path)
        self.arm = args.arm
        self.recovery_mode = self.contract.get("candidate_failure_policy") == RECOVERY_MODE
        self.resuming = bool(getattr(args, "resume", False))
        if self.resuming and not self.recovery_mode:
            raise RuntimeError("resume requires the new reviewed failure contract")
        self.root = args.run_root.resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.output = self.root / self.arm
        if self.resuming and not self.output.is_dir():
            raise RuntimeError("resume requires the original arm directory")
        self.output.mkdir(exist_ok=self.resuming)
        self.ownership = (self.output / "controller.lock").open("a")
        fcntl.flock(self.ownership, fcntl.LOCK_EX | fcntl.LOCK_NB)
        self.env = os.environ.copy()
        self.env.update({name: "1" for name in THREAD_ENV})
        self.env.update(PYTHONUNBUFFERED="1", PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg",
                        NUMBA_CACHE_DIR=str(self.output / "numba_cache"),
                        UTILITY_COMPARISON_NUMBA_CACHE=str(self.output / "numba_cache"))
        (self.output / "numba_cache").mkdir(exist_ok=self.resuming)
        self.runner = self.contract["files"]["run_e5f_utility_comparison.py"]["path"]
        self.collector = self.contract["files"][COLLECTOR_NAME]["path"]
        self.phase = "readiness"
        self.last_heartbeat = 0.
        self.interrupted = False
        self.best = None
        self.clock = None
        self.stop_reason = None
        self.records = []
        self.attempts = []
        self.reused = []
        self.planned = {"smoke": 2, "initial": 40, "de": 120, "repeat": 2}
        self.collection = None
        self.smoke_verified = False
        self.completed_by_id = {}
        for name in ("latest_completed.json", "best_so_far.json"):
            if not self.resuming:
                create(self.output / name, dict(status="no_completed_objective", arm=self.arm,
                                               contract_sha256=self.fingerprint))
        if self.resuming:
            try:
                self.restore()
            except Exception:
                self.mark_fatal("invalid_resume_checkpoint")
                raise
        for signum in (signal.SIGTERM, signal.SIGINT):
            signal.signal(signum, self.interrupt)

    def expected_context(self, candidate):
        return dict(candidate_id=candidate["id"], stage=candidate["controller_stage"],
            contract_sha256=self.fingerprint,
            source_sha256=self.contract["files"]["run_e5f_utility_comparison.py"]["sha256"],
            target_sha256=self.contract["arms"][self.arm]["objective"]["sha256"],
            point_sha256=design.canonical_fingerprint(candidate["parameters"]))

    def validate_plan(self, candidate, plan):
        expected = dict(contract_sha256=self.fingerprint, arm=self.arm, point=candidate["parameters"],
            candidate_id=candidate["id"], controller_stage=candidate["controller_stage"],
            runner_source_sha256=self.contract["files"]["run_e5f_utility_comparison.py"]["sha256"],
            graphs=candidate["controller_stage"] == "smoke",
            objective_cap_seconds=self.budget["objective_timeout_seconds"], deadline_owner="controller")
        if any(plan.get(key) != value for key, value in expected.items()):
            raise RuntimeError("saved case plan differs from exact candidate/contract")
        if (set(plan) != set(expected) | {"deadline_epoch", "controller_pid"}
                or not math.isfinite(plan["deadline_epoch"])
                or type(plan["controller_pid"]) is not int or plan["controller_pid"] <= 0):
            raise RuntimeError("unexpected case plan fields/deadline")
        stage_cutoff = self.clock["repeat_cutoff_epoch" if candidate["controller_stage"] == "repeat" else "search_cutoff_epoch"]
        if not self.clock["start_epoch"] <= plan["deadline_epoch"] <= stage_cutoff:
            raise RuntimeError("saved case deadline lies outside the original shared clock")

    def checkpoint(self):
        if not self.recovery_mode or self.clock is None:
            return
        files = {}
        for pattern in ("records/*.json", "dispatch/*.json", "plans/*.json", "population_*.json",
                        "proposals_*.json", "external_receipts/*.json", "selected.json",
                        "initial_seed_reuse.json", "smoke_comparison.json", "complete.json"):
            for path in self.output.glob(pattern):
                files[str(path.relative_to(self.output))] = file_hash(path)
        write(self.output / "resume_checkpoint.json", dict(contract_sha256=self.fingerprint,
            clock_sha256=file_hash(self.root / "clock.json"), files=files))

    def retain(self, path, value):
        if self.resuming and Path(path).exists():
            if read(path) != value:
                raise RuntimeError("resume checkpoint differs from reconstructed state: " + str(path))
        else:
            create(path, value)
        self.checkpoint()

    def restore(self):
        checkpoint = read(self.output / "resume_checkpoint.json")
        if (checkpoint["contract_sha256"] != self.fingerprint
                or checkpoint["clock_sha256"] != file_hash(self.root / "clock.json")):
            raise RuntimeError("resume requires the identical contract and unchanged shared clock")
        self.clock = read(self.root / "clock.json")
        if self.clock["contract_sha256"] != self.fingerprint:
            raise RuntimeError("resume clock belongs to another contract")
        actual_files = {str(p.relative_to(self.output)) for pattern in (
            "records/*.json", "dispatch/*.json", "plans/*.json", "population_*.json",
            "proposals_*.json", "external_receipts/*.json", "selected.json", "initial_seed_reuse.json",
            "smoke_comparison.json", "complete.json") for p in self.output.glob(pattern)}
        if actual_files != set(checkpoint["files"]):
            raise RuntimeError("uncheckpointed/orphan durable state cannot be resumed")
        for relative, digest in checkpoint["files"].items():
            path = self.output / relative
            if path.resolve().parent != self.output.resolve() and self.output.resolve() not in path.resolve().parents:
                raise RuntimeError("checkpoint path escapes arm output")
            if file_hash(path) != digest:
                raise RuntimeError("resume checkpoint artifact changed: " + relative)
        dispatch = {p.stem: p for p in self.output.glob("dispatch/*.json")}
        records = {p.stem: p for p in self.output.glob("records/*.json")}
        plans = {p.stem: p for p in self.output.glob("plans/*.json")}
        if set(dispatch) != set(records) or set(plans) != set(records):
            raise RuntimeError("orphan dispatched/planned case: no automatic retry or adoption")
        for identifier in records:
            for path in (dispatch[identifier], records[identifier], plans[identifier]):
                if str(path.relative_to(self.output)) not in checkpoint["files"]:
                    raise RuntimeError("uncheckpointed case cannot be resumed")
            record, attempt, plan = read(records[identifier]), read(dispatch[identifier]), read(plans[identifier])
            candidate = dict(id=identifier, parameters=record["point"], controller_stage=record["stage"])
            self.validate_plan(candidate, plan)
            if (record["contract_sha256"] != self.fingerprint or record["candidate_id"] != identifier
                    or attempt["candidate_id"] != identifier or attempt["stage"] != record["stage"]
                    or record["deadline_evidence"]["supervisor_pid"] != plan["controller_pid"]
                    or record["plan_sha256"] != file_hash(plans[identifier])
                    or record["dispatch_sha256"] != file_hash(dispatch[identifier])):
                raise RuntimeError("record/dispatch/plan provenance differs")
            for relative, digest in record["artifact_sha256"].items():
                path = self.output / relative
                if self.output.resolve() not in path.resolve().parents or file_hash(path) != digest:
                    raise RuntimeError("saved case artifact changed")
            self.records.append(record)
            self.attempts.append(attempt)
            self.completed_by_id[identifier] = record
        valid = [r for r in self.records if r["status"] == "success" and r["stage"] != "repeat"]
        self.best = min(valid, key=lambda r: (r["loss"], r["finished_epoch"])) if valid else None
        reuse = self.output / "initial_seed_reuse.json"
        if reuse.exists():
            self.reused = [read(reuse)]
        for name in ("smoke_comparison", "collection"):
            receipt = self.output / "external_receipts" / (name + ".json")
            if (self.output / (name + ".log")).exists() and not receipt.exists():
                raise RuntimeError("orphan external operation cannot be repeated")
            if receipt.exists():
                if str(receipt.relative_to(self.output)) not in checkpoint["files"]:
                    raise RuntimeError("external receipt was not checkpointed")
                for relative, digest in read(receipt).get("artifact_sha256", {}).items():
                    if file_hash(self.output / relative) != digest:
                        raise RuntimeError("external output changed")
        smoke = self.output / "external_receipts/smoke_comparison.json"
        self.smoke_verified = smoke.exists() and read(smoke)["status"] == "success"
        selection = self.output / "selected.json"
        if selection.exists():
            self.stop_reason = read(selection)["search_stop_reason"]

    def interrupt(self, signum, _frame):
        self.interrupted = True
        self.stop_reason = "signal_" + str(signum)

    def heartbeat(self, *, force=False, **state):
        now = time.time()
        if force or now - self.last_heartbeat >= HEARTBEAT_SECONDS:
            value = dict(arm=self.arm, phase=self.phase, epoch=now,
                         contract_sha256=self.fingerprint, pid=os.getpid(),
                         completed_cases=len(self.records), stop_reason=self.stop_reason,
                         **{key: value for key, value in state.items() if key != "stop_reason"})
            if state.get("stop_reason"):
                value["stop_reason"] = state["stop_reason"]
            write(self.output / "heartbeat.json", value)
            self.last_heartbeat = now

    def readiness(self):
        """Bound startup queue skew before the shared eight-hour clock begins."""
        if self.resuming:
            self.heartbeat(force=True)
            return  # restore authenticated the original clock; never reset it.
        barrier = self.root / "barrier"
        barrier.mkdir(exist_ok=True)
        create(barrier / (self.arm + ".ready.json"),
               dict(arm=self.arm, contract_sha256=self.fingerprint, ready_epoch=time.time(),
                    array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID")))
        with (barrier / "clock.lock").open("a") as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            first = barrier / "readiness.json"
            if not first.exists():
                create(first, dict(start_epoch=time.time(), maximum_wait_seconds=READINESS_SECONDS))
        cutoff = read(first)["start_epoch"] + READINESS_SECONDS
        while True:
            if self.interrupted:
                raise RuntimeError("controller interrupted before the shared clock")
            if time.time() >= cutoff and not (self.root / "clock.json").exists():
                raise TimeoutError("four-controller readiness exceeded thirty minutes")
            paths = [barrier / (arm + ".ready.json") for arm in ARMS]
            if all(path.exists() for path in paths):
                rows = [read(path) for path in paths]
                if any(row["contract_sha256"] != self.fingerprint for row in rows):
                    raise RuntimeError("mixed contracts at readiness barrier")
                if len({row["array_job_id"] for row in rows}) != 1:
                    raise RuntimeError("readiness files belong to different Slurm arrays")
                with (barrier / "clock.lock").open("a") as stream:
                    fcntl.flock(stream, fcntl.LOCK_EX)
                    clock_path = self.root / "clock.json"
                    if not clock_path.exists():
                        start = time.time()
                        total = self.budget["total_limit_seconds"]
                        create(clock_path, dict(start_epoch=start, end_epoch=start + total,
                            search_cutoff_epoch=start + total - self.budget["repeat_reserve_seconds"]
                                                - self.budget["export_reserve_seconds"],
                            repeat_cutoff_epoch=start + total - self.budget["export_reserve_seconds"],
                            contract_sha256=self.fingerprint,
                            rule="clock begins immediately before first exact-loop smoke dispatch"))
                    self.clock = read(clock_path)
                if self.clock["contract_sha256"] != self.fingerprint:
                    raise RuntimeError("shared clock contract mismatch")
                self.checkpoint()
                return
            if time.time() >= cutoff:
                raise TimeoutError("four-controller readiness exceeded thirty minutes")
            self.heartbeat()
            time.sleep(min(2., max(.001, cutoff - time.time())))

    def launch(self, candidate, deadline):
        case_id = candidate["id"]
        output = self.output / case_id
        plan = self.output / "plans" / (case_id + ".json")
        cap = self.budget["objective_timeout_seconds"]
        cutoff = min(deadline, time.time() + cap)
        request = dict(contract_sha256=self.fingerprint, arm=self.arm,
                          point=candidate["parameters"], deadline_epoch=cutoff,
                          objective_cap_seconds=cap, graphs=candidate["controller_stage"] == "smoke")
        if self.recovery_mode:
            request.update(candidate_id=case_id, controller_stage=candidate["controller_stage"],
                runner_source_sha256=self.contract["files"]["run_e5f_utility_comparison.py"]["sha256"],
                deadline_owner="controller", controller_pid=os.getpid())
            self.validate_plan(candidate, request)
        create(plan, request)
        command = [sys.executable, self.runner, "--stage", "evaluate", "--contract",
                   str(self.contract_path), "--arm", self.arm, "--case-plan", str(plan),
                   "--output", str(output)]
        attempt = dict(candidate_id=case_id, stage=candidate["controller_stage"],
                       dispatch_epoch=time.time(), deadline_epoch=cutoff)
        self.attempts.append(attempt)
        create(self.output / "dispatch" / (case_id + ".json"), attempt)
        return ManagedProcess(command, self.output / (case_id + ".log"), cutoff, self.env)

    def recovery_failure_kind(self, candidate, failure):
        evidence = failure["recovery_evidence"]
        if (evidence["context"] != self.expected_context(candidate) or failure["arm"] != self.arm
                or failure["candidate_failure_policy"] != RECOVERY_MODE
                or evidence["raw_type"] != failure["error_type"] or evidence["raw_error"] != failure["error"]):
            raise RuntimeError("failure evidence provenance mismatch")
        pin = self.contract["parent_source_inventory"]
        if file_hash(pin["path"]) != pin["sha256"]:
            raise RuntimeError("native source inventory changed")
        relative = "code/model/intergen_eqscale_seq_optimized/solver.py"
        entry = read(pin["path"])["files"][relative]
        digest = entry if isinstance(entry, str) else entry["sha256"]
        native = failure["native_source"]
        expected_path = Path(self.contract["reference_root"]) / "source" / relative
        if (native is None or Path(native["native_source_path"]).resolve() != expected_path.resolve()
                or native["native_source_sha256"] != digest
                or file_hash(expected_path) != digest
                or native["native_manifest_sha256"] != pin["sha256"]
                or native["native_gate_tolerance"] != 1e-12 or native["native_value_cutoff"] != -1e9):
            raise RuntimeError("native failure producer authentication failed")
        if (evidence["raw_type"] == "RuntimeError"
                and evidence["raw_error"].startswith(OLD_INADMISSIBLE_PREFIXES)
                and evidence["narrow_infeasibility_verified"] is False):
            return "normalization_or_equilibrium_gate"
        policy = importlib.import_module(RECOVERY_NAME[:-3])
        if (evidence["raw_type"] == "InfeasibleThetaError"
                and evidence["narrow_infeasibility_verified"] is True
                and evidence["validation"] == "exact_native_type_and_structured_gate_failure"
                and narrow_payload_valid(evidence["native_payload"], policy.CENSUS_FIELDS)):
            return "reviewed_native_age34_evaluation_rejection"
        return None

    def mark_fatal(self, reason, candidate_id=None):
        path = self.root / "barrier" / (self.arm + ".fatal.json")
        if not path.exists():
            create(path, dict(contract_sha256=self.fingerprint, arm=self.arm,
                              reason=reason, candidate_id=candidate_id, epoch=time.time()))

    def global_guard(self, phase):
        if not self.recovery_mode or phase not in {"smoke", "initial", "de"}:
            return None
        for arm in ARMS:
            path = self.root / "barrier" / (arm + ".fatal.json")
            if path.exists():
                try:
                    if read(path)["contract_sha256"] != self.fingerprint:
                        return "global_guard_contract_mismatch"
                except (ValueError, KeyError):
                    return "global_guard_malformed"
                return "another_arm_fatal_stop"
        return None

    def phase_barrier(self, name, success):
        if not self.recovery_mode:
            return success
        path = self.root / "barrier" / (name + "." + self.arm + ".json")
        value = dict(contract_sha256=self.fingerprint, phase=name, arm=self.arm, success=success)
        if path.exists():
            if read(path) != value:
                raise RuntimeError("phase checkpoint differs on resume")
        else:
            create(path, value)
        if not success:
            self.mark_fatal("failed_or_incomplete_phase:" + name)
            return False
        while time.time() < self.clock["search_cutoff_epoch"]:
            if self.global_guard("initial") or self.interrupted:
                return False
            paths = [self.root / "barrier" / (name + "." + arm + ".json") for arm in ARMS]
            if all(p.exists() for p in paths):
                rows = [read(p) for p in paths]
                if any(r["contract_sha256"] != self.fingerprint or r["phase"] != name
                       or not r["success"] for r in rows):
                    self.mark_fatal("invalid_phase_barrier:" + name)
                    return False
                return True
            self.heartbeat()
            time.sleep(1.)
        self.mark_fatal("missing_phase_at_absolute_cutoff:" + name)
        return False

    def finish(self, candidate, process, code):
        output = self.output / candidate["id"]
        record = dict(candidate_id=candidate["id"], stage=candidate["controller_stage"],
                      case_output=str(output), point=candidate["parameters"],
                      started_epoch=process.started_epoch, finished_epoch=time.time(),
                      returncode=code, status="timed_out" if process.timed_out else "failed")
        failure_path = output / "failure.json"
        if self.recovery_mode:
            plan_path = self.output / "plans" / (candidate["id"] + ".json")
            saved_plan = read(plan_path)
            self.validate_plan(candidate, saved_plan)
            dispatch_path = self.output / "dispatch" / (candidate["id"] + ".json")
            dispatched = [a for a in self.attempts if a["candidate_id"] == candidate["id"]]
            if len(dispatched) != 1 or read(dispatch_path) != dispatched[0]:
                raise RuntimeError("dispatch receipt differs from the owned attempt")
            record.update(contract_sha256=self.fingerprint, plan_sha256=file_hash(plan_path),
                dispatch_sha256=file_hash(self.output / "dispatch" / (candidate["id"] + ".json")),
                deadline_evidence=dict(supervisor_pid=os.getpid(), child_pid=process.process.pid,
                    deadline_epoch=process.deadline, observed_epoch=process.observed_epoch,
                    expired=process.deadline_expired, observed_running_at_expiry=process.observed_running_at_expiry,
                    owned_sigkill_reaped=process.deadline_kill_reaped, cancelled=process.cancelled), status="failed")
            known_kind = None
            invalid_evidence = process.cancelled or process.deadline != saved_plan["deadline_epoch"]
            if invalid_evidence:
                record.update(error_type="CancelledOrChangedDeadline", error="approved objective deadline was not exhausted normally")
            owned_timeout = (not invalid_evidence and process.deadline_expired
                and process.observed_running_at_expiry and process.deadline_kill_reaped
                and code == -signal.SIGKILL and process.scheduled_deadline == saved_plan["deadline_epoch"])
            provenance_path = output / "attempt_provenance.json"
            try:
                provenance = read(provenance_path)
                expected_provenance = dict(context=self.expected_context(candidate), arm=self.arm,
                    candidate_failure_policy=RECOVERY_MODE, case_plan_sha256=file_hash(plan_path))
                if provenance != expected_provenance:
                    raise RuntimeError("attempt provenance differs from pinned request")
                record["attempt_provenance_status"] = "verified"
            except FileNotFoundError as exc:
                # Startup can exhaust its cap before the pinned runner writes
                # its sidecar. Only absence plus positive controller-owned kill
                # evidence is censored; conflicting/partial evidence stays fatal.
                if owned_timeout and not provenance_path.is_symlink() and not output.is_symlink():
                    record["attempt_provenance_status"] = "not_written_before_owned_deadline"
                else:
                    invalid_evidence = True
                    record.update(error_type=type(exc).__name__, error=str(exc))
            except (ValueError, OSError, RuntimeError) as exc:
                invalid_evidence = True
                record.update(error_type=type(exc).__name__, error=str(exc))
            if failure_path.is_file():
                try:
                    failure = read(failure_path)
                    record["failure"] = failure
                    known_kind = self.recovery_failure_kind(candidate, failure)
                    invalid_evidence = invalid_evidence or known_kind is None
                except (KeyError, ValueError, TypeError, OSError, RuntimeError) as exc:
                    invalid_evidence = True
                    record.update(error_type=type(exc).__name__, error=str(exc))
            local = candidate["controller_stage"] in {"initial", "de"}
            if not invalid_evidence:
                if owned_timeout:
                    record.update(status="censored_timeout", loss=None)
                elif process.deadline_expired and code == 0:
                    record.update(status="censored_late_completion", loss=None)
                elif not process.deadline_expired and code != 0 and known_kind and local:
                    record.update(status="inadmissible", loss=None, inadmissible_kind=known_kind)
            record["halt_new_dispatch"] = record["status"] == "failed" or (not local and record["status"] in CENSORED_STATUSES)
        elif code != 0 and not process.timed_out and failure_path.is_file():
            try:
                failure = read(failure_path)
                if classified_inadmissible(failure, candidate["controller_stage"]):
                    record.update(status="inadmissible", loss=None, failure=failure)
            except (ValueError, OSError):
                pass  # Malformed failure receipts remain unknown failures.
        if code == 0 and not process.timed_out and (not self.recovery_mode or (
                not process.deadline_expired and not failure_path.exists() and not invalid_evidence)):
            try:
                case = output / "case"
                for name in ("receipt.json", "initial_state.pkl.gz", "target_fit.csv", "parameters.csv"):
                    if not (case / name).is_file():
                        raise RuntimeError("missing successful-case artifact: " + name)
                receipt = read(case / "receipt.json")
                if (receipt["comparison_contract_sha256"] != self.fingerprint
                        or receipt["utility_comparison_arm"] != self.arm
                        or receipt["target_system_sha256"] != self.contract["arms"][self.arm]["objective"]["sha256"]
                        or receipt["point"] != candidate["parameters"]):
                    raise RuntimeError("case result differs from its contract, arm, target or point")
                loss = float(receipt["loss"])
                if not math.isfinite(loss) or loss < 0.:
                    raise RuntimeError("invalid scored loss")
                record.update(status="success", loss=loss,
                              halt_new_dispatch=False,
                              receipt_sha256=file_hash(case / "receipt.json"),
                              target_system_sha256=receipt["target_system_sha256"])
            except Exception as exc:
                record.update(error_type=type(exc).__name__, error=str(exc))
        if self.recovery_mode:
            paths = [p for p in (output / "case").glob("*") if p.is_file() and p.name in {
                "receipt.json", "initial_state.pkl.gz", "target_fit.csv", "parameters.csv"}]
            paths += [p for p in (failure_path, output / "attempt_provenance.json",
                                 self.output / (candidate["id"] + ".log")) if p.is_file()]
            record["artifact_sha256"] = {str(p.relative_to(self.output)): file_hash(p) for p in paths}
            if record.get("halt_new_dispatch"):
                self.mark_fatal("case_failed_no_retry", candidate["id"])
        self.records.append(record)
        self.completed_by_id[candidate["id"]] = record
        create(self.output / "records" / (candidate["id"] + ".json"), record)
        write(self.output / "latest_completed.json", record)
        if record["status"] == "success" and candidate["controller_stage"] != "repeat":
            if self.best is None or record["loss"] < self.best["loss"]:
                self.best = copy.deepcopy(record)
                write(self.output / "best_so_far.json", self.best)
        self.summary()
        self.checkpoint()
        self.heartbeat(force=True)
        return record

    def batch(self, candidates, phase, deadline, workers=None):
        self.phase = phase
        for candidate in candidates:
            candidate["controller_stage"] = phase
        retained, pending = [], []
        for candidate in candidates:
            previous = self.completed_by_id.get(candidate["id"])
            if previous is None:
                pending.append(candidate)
            else:
                if previous["stage"] != phase or previous["point"] != candidate["parameters"]:
                    raise RuntimeError("saved candidate differs from the exact planned slot")
                retained.append(previous)
        allowed = COMPLETED_STATUSES | (CENSORED_STATUSES if self.recovery_mode and phase in {"initial", "de"} else set())
        if any(r["status"] not in allowed or r.get("halt_new_dispatch") for r in retained):
            return dict(results=retained, unrun_ids=[r["id"] for r in pending],
                        stop_reason="retained_fatal_case_no_retry", complete=not pending)
        result = run_batch(pending, workers=workers or self.budget["workers_per_arm"],
                         deadline=deadline, launch=self.launch, finish=self.finish,
                         heartbeat=self.heartbeat, interrupted=lambda: self.interrupted,
                         allowed_statuses=allowed, guard=lambda: self.global_guard(phase),
                         on_finish_error=lambda candidate, exc: self.mark_fatal(
                             "finish_validation_error:" + type(exc).__name__ + ":" + str(exc), candidate["id"]))
        result["results"] = retained + result["results"]
        result["complete"] = len(result["results"]) == len(candidates)
        return result

    def external(self, command, name, deadline):
        self.phase = name
        saved = self.output / "external_receipts" / (name + ".json")
        if self.resuming and saved.exists():
            return read(saved)
        if time.time() >= deadline or self.interrupted:
            return dict(status="not_run", reason="deadline_or_interruption")
        process = ManagedProcess(command, self.output / (name + ".log"), deadline, self.env)
        try:
            while process.poll() is None:
                if self.interrupted:
                    process.cancel()
                self.heartbeat()
                time.sleep(.5)
            result = dict(status="success" if process.process.returncode == 0 and not process.cancelled and (
                not self.recovery_mode or not process.deadline_expired) else "failed",
                returncode=process.process.returncode, timed_out=process.timed_out, cancelled=process.cancelled)
        finally:
            process.close()
        if self.recovery_mode:
            outputs = [self.output / "smoke_comparison.json"] if name == "smoke_comparison" else list((self.output / "export").rglob("*"))
            marker = self.output / "smoke_comparison.json" if name == "smoke_comparison" else self.output / "export/collection_receipt.json"
            if result["status"] == "success" and not marker.is_file():
                result.update(status="failed", error="external process omitted its required completion receipt")
            result["artifact_sha256"] = {str(p.relative_to(self.output)): file_hash(p) for p in outputs if p.is_file()}
            create(saved, result)
            self.checkpoint()
        return result

    def smoke(self):
        seed = copy.deepcopy(self.bank["arms"][self.arm][0])
        completed = []
        for index in range(2):
            candidate = copy.deepcopy(seed)
            candidate["id"] = f"smoke_{index}"
            result = self.batch([candidate], "smoke", self.clock["search_cutoff_epoch"], workers=1)
            completed.extend(result["results"])
            if not result["complete"] or result["stop_reason"]:
                self.stop_reason = result["stop_reason"] or "incomplete_smoke"
                break
        passed = len(completed) == 2 and all(row["status"] == "success" for row in completed)
        if passed:
            command = [sys.executable, self.collector, "--stage", "smoke", "--contract",
                       str(self.contract_path), "--arm", self.arm, "--original",
                       completed[0]["case_output"], "--repeat", completed[1]["case_output"],
                       "--required-count", "1", "--output", str(self.output / "smoke_comparison.json")]
            comparison = self.external(command, "smoke_comparison", min(self.clock["search_cutoff_epoch"],
                time.time() + self.budget["export_reserve_seconds"]))
            passed = comparison["status"] == "success"
            if not passed:
                self.stop_reason = "exact_smoke_comparison_failed"
        marker = self.root / "barrier" / (self.arm + ".smoke.json")
        value = dict(arm=self.arm, status="success" if passed else "failed", contract_sha256=self.fingerprint)
        if self.resuming and marker.exists():
            if read(marker) != value:
                raise RuntimeError("saved smoke barrier differs")
        else:
            create(marker, value)
        if not passed:
            return None
        self.phase = "common_smoke_barrier"
        while time.time() < self.clock["search_cutoff_epoch"] and not self.interrupted:
            stop = self.global_guard("smoke")
            if stop:
                self.stop_reason = stop
                return None
            paths = [self.root / "barrier" / (arm + ".smoke.json") for arm in ARMS]
            rows = [read(path) for path in paths if path.exists()]
            if any(row["contract_sha256"] != self.fingerprint or row["status"] != "success" for row in rows):
                self.stop_reason = "another_arm_smoke_failed_no_search"
                return None
            if len(rows) == len(ARMS):
                return completed[0]
            self.heartbeat()
            time.sleep(2.)
        self.stop_reason = "common_smoke_barrier_budget_exhausted"
        return None

    def search(self, smoke_result):
        population = copy.deepcopy(self.bank["arms"][self.arm])
        # The exact smoke already scored this identical initial seed. Reuse that
        # objective; the pinned budget separately records this reused slot.
        reused = dict(smoke_result, candidate_id=population[0]["id"], stage="initial",
                      reused_from=smoke_result["candidate_id"])
        if not self.reused:
            self.reused.append(reused)
        self.retain(self.output / "initial_seed_reuse.json", reused)
        result = self.batch(population[1:], "initial", self.clock["search_cutoff_epoch"])
        records = [reused] + result["results"]
        if not result["complete"] or result["stop_reason"]:
            self.stop_reason = result["stop_reason"] or "incomplete_initial_population"
            self.phase_barrier("initial", False)
            return
        if self.recovery_mode and any(r["status"] in CENSORED_STATUSES for r in records):
            self.stop_reason = "typed_censored_DE_selection_not_integrated"
            self.phase_barrier("initial", False)
            return
        try:
            scores = complete_scores(population, records)
        except BarrierStop as exc:
            self.stop_reason = str(exc)
            self.phase_barrier("initial", False)
            return
        self.retain(self.output / "population_000.json", dict(population=population, scores=scores))
        if not self.phase_barrier("initial", True):
            self.stop_reason = "global_initial_phase_stop"
            return
        for generation in range(1, self.budget["de_generations"] + 1):
            trials = design.make_de_generation(self.catalog, self.arm, population, scores,
                generation=generation, budget=self.budget, rng_seed=self.bank["rng_seed"],
                mutation_factor=self.contract["differential_evolution"]["mutation_factor"])
            self.retain(self.output / f"proposals_{generation:03d}.json", trials)
            result = self.batch(trials, "de", self.clock["search_cutoff_epoch"])
            if not result["complete"] or result["stop_reason"]:
                self.stop_reason = result["stop_reason"] or "incomplete_de_generation"
                self.phase_barrier(f"de_{generation:03d}", False)
                return
            if self.recovery_mode and any(r["status"] in CENSORED_STATUSES for r in result["results"]):
                self.stop_reason = "typed_censored_DE_selection_not_integrated"
                self.phase_barrier(f"de_{generation:03d}", False)
                return
            try:
                trial_scores = complete_scores(trials, result["results"])
            except BarrierStop as exc:
                self.stop_reason = str(exc)
                self.phase_barrier(f"de_{generation:03d}", False)
                return
            population, scores = design.select_de_generation(population, scores, trials, trial_scores)
            self.retain(self.output / f"population_{generation:03d}.json", dict(population=population, scores=scores))
            if not self.phase_barrier(f"de_{generation:03d}", True):
                self.stop_reason = "global_de_phase_stop"
                return
        self.stop_reason = "finite_search_bank_complete"

    def summary(self):
        counts = {}
        for stage, planned in self.planned.items():
            rows = [row for row in self.records if row["stage"] == stage]
            attempted = sum(row["stage"] == stage for row in self.attempts)
            reused = sum(row["stage"] == stage for row in self.reused)
            counts[stage] = dict(planned=planned, attempted=attempted, reused=reused,
                succeeded=sum(row["status"] == "success" for row in rows),
                inadmissible=sum(row["status"] == "inadmissible" for row in rows),
                failed=sum(row["status"] == "failed" for row in rows),
                timed_out=sum(row["status"] == "timed_out" for row in rows),
                censored_timeout=sum(row["status"] == "censored_timeout" for row in rows),
                censored_late_completion=sum(row["status"] == "censored_late_completion" for row in rows),
                incomplete=attempted - len(rows), unrun=planned - attempted - reused)
        value = dict(schema="e5f_utility_arm_controller_v1", arm=self.arm,
            contract_sha256=self.fingerprint, phase=self.phase, stop_reason=self.stop_reason,
            clock=self.clock, counts=counts, best=self.best, collection=self.collection,
            elapsed_since_shared_clock_seconds=(time.time() - self.clock["start_epoch"]) if self.clock else None,
            automatic_retry=False, gate_relaxation=False, local_numerical_workers=0,
            identification="moment/parameter counts are not a numerical rank or identification result")
        write(self.output / "controller_summary.json", value)
        return value

    def repeat_and_collect(self):
        if self.best is None or self.clock is None or not self.smoke_verified:
            self.phase = "no_valid_selection"
            self.summary()
            return
        repeats = [dict(id=f"repeat_{index}", arm=self.arm, parameters=self.best["point"]) for index in range(2)]
        selection = dict(contract_sha256=self.fingerprint, arm=self.arm,
            candidate_id=self.best["candidate_id"], point=self.best["point"], loss=self.best["loss"],
            original_case_output=self.best["case_output"],
            repeat_case_outputs=[str(self.output / row["id"]) for row in repeats],
            selected_epoch=time.time(), search_stop_reason=self.stop_reason)
        if self.resuming and (self.output / "selected.json").exists():
            old = read(self.output / "selected.json")
            for key in set(selection) - {"selected_epoch", "search_stop_reason"}:
                if old[key] != selection[key]:
                    raise RuntimeError("retained frozen selection differs from verified best")
        else:
            create(self.output / "selected.json", selection)
            self.checkpoint()
        result = self.batch(repeats, "repeat", self.clock["repeat_cutoff_epoch"], workers=2)
        self.summary()
        if self.recovery_mode and (not result["complete"] or len(result["results"]) != 2
                or any(row["status"] != "success" for row in result["results"])):
            self.collection = dict(status="incomplete_repetitions", exact_repeat_claim=False,
                reason="both_required_repeats_must_succeed_before_their_approved_deadlines",
                statuses={row["candidate_id"]: row["status"] for row in result["results"]},
                unrun_ids=result["unrun_ids"], automatic_retry=False)
            self.retain(self.output / "external_receipts/collection.json", self.collection)
            self.summary()
            return
        command = [sys.executable, self.collector, "--stage", "collect", "--contract",
                   str(self.contract_path), "--arm", self.arm, "--run-root", str(self.root),
                   "--output", str(self.output / "export")]
        self.collection = self.external(command, "collection", self.clock["end_epoch"])

    def run(self):
        error = None
        if self.resuming and (self.output / "complete.json").exists():
            return read(self.output / "complete.json")["controller_exit_code"]
        try:
            self.readiness()
            if not (self.resuming and (self.output / "selected.json").exists()):
                smoke = self.smoke()
                if smoke is not None:
                    self.smoke_verified = True
                    self.search(smoke)
        except Exception as exc:
            error = dict(error_type=type(exc).__name__, error=str(exc))
            self.stop_reason = "controller_failure_no_retry"
            if not (self.output / "controller_failure.json").exists():
                create(self.output / "controller_failure.json", error)
            if self.recovery_mode:
                self.mark_fatal("controller_failure_no_retry")
        finally:
            try:
                self.repeat_and_collect()
            except Exception as exc:
                self.collection = dict(status="failed", error_type=type(exc).__name__, error=str(exc))
            self.phase = "finished"
            summary = self.summary()
            self.heartbeat(force=True)
        fully_scored = all(row["failed"] == row["timed_out"] == row["unrun"] == row["incomplete"]
                           == row["censored_timeout"] == row["censored_late_completion"] == 0
                           for row in summary["counts"].values())
        ok = (error is None and fully_scored and self.stop_reason == "finite_search_bank_complete"
              and (self.collection or {}).get("status") == "success")
        exit_code = 0 if ok else 2
        summary["controller_exit_code"] = exit_code
        create(self.output / "complete.json", summary)
        self.checkpoint()
        return exit_code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--arm", choices=ARMS, required=True)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--resume", action="store_true", help="same-contract checkpoint recovery; never retries a dispatched case")
    args = parser.parse_args()
    if args.check_only:
        _, fingerprint, budget, _, _ = load_contract(args.contract)
        print(json.dumps(dict(status="approved_frozen_controller_contract_checked", contract_sha256=fingerprint,
                              total_limit_seconds=budget["total_limit_seconds"], native_solves=0)))
        return 0
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("native controller requires a Torch Slurm allocation")
    if int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) < 10:
        raise RuntimeError("ten allocated CPUs are required for this frozen arm schedule")
    return Controller(args).run()


if __name__ == "__main__":
    raise SystemExit(main())
