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
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import e5f_utility_comparison_design as design


ARMS = design.ARM_NAMES
THREAD_ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
HEARTBEAT_SECONDS = 30.
READINESS_SECONDS = 1800.
COLLECTOR_NAME = "collect_e5f_utility_comparison.py"
COMPLETED_STATUSES = {"success", "inadmissible"}


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
                "budget", "arm_catalog", "proposal_bank")
    for key in required:
        if key not in contract["files"]:
            raise RuntimeError("missing required source/design pin: " + key)
    for pin in contract["files"].values():
        if file_hash(pin["path"]) != pin["sha256"]:
            raise RuntimeError("pinned source/input changed: " + pin["path"])
    for actual in (__file__, design.__file__):
        if Path(actual).resolve() != Path(contract["files"][Path(actual).name]["path"]).resolve():
            raise RuntimeError("executing source is not the pinned file: " + actual)
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
        self.started_epoch = time.time()
        self.log = Path(log_path).open("xb")
        self.timed_out = False
        try:
            self.process = subprocess.Popen(self.command, stdout=self.log, stderr=subprocess.STDOUT,
                                            env=env, start_new_session=True)
        except BaseException:
            self.log.close()
            raise

    def kill(self):
        # Only a group created here with start_new_session; never a Slurm job or
        # a discovered PID. SIGKILL prevents an unbounded graceful-shutdown tail.
        try:
            os.killpg(self.process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        self.process.wait(timeout=10)

    def poll(self):
        code = self.process.poll()
        if code is None and time.time() >= self.deadline:
            self.timed_out = True
            self.kill()
            code = self.process.returncode
        if code is not None:
            self.log.close()
        return code

    def close(self):
        # Also reap a background descendant if its immediate parent exited.
        self.kill()
        self.log.close()


def run_batch(candidates, *, workers, deadline, launch, finish, heartbeat,
              interrupted=lambda: False, poll_seconds=1.):
    """Run a finite barrier. Any failed result stops new dispatch, not siblings.

    Missing/unrun trials are never submitted to DE selection as rejected scores.
    Dependency injection permits synthetic process tests without model imports.
    """
    pending, active, results = list(candidates), {}, []
    reason = None
    try:
        while pending or active:
            now = time.time()
            if interrupted() and reason is None:
                reason = "controller_interrupted"
            if now >= deadline and reason is None:
                reason = "stage_deadline_exhausted"
            # Reap every completed child before filling vacancies, so an already
            # observed failure cannot be followed by an unnecessary new dispatch.
            for key, (candidate, process) in list(active.items()):
                code = process.poll()
                if code is None:
                    continue
                del active[key]
                try:
                    result = finish(candidate, process, code)
                finally:
                    process.close()
                results.append(result)
                if result["status"] not in COMPLETED_STATUSES and reason is None:
                    reason = "case_failed_no_retry"
            if interrupted():
                for _, process in active.values():
                    process.deadline = min(process.deadline, time.time())
            while pending and len(active) < workers and reason is None:
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


class Controller:
    def __init__(self, args):
        self.args = args
        self.contract_path = args.contract.resolve()
        (self.contract, self.fingerprint, self.budget,
         self.catalog, self.bank) = load_contract(self.contract_path)
        self.arm = args.arm
        self.root = args.run_root.resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.output = self.root / self.arm
        self.output.mkdir(exist_ok=False)
        self.env = os.environ.copy()
        self.env.update({name: "1" for name in THREAD_ENV})
        self.env.update(PYTHONUNBUFFERED="1", PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg",
                        NUMBA_CACHE_DIR=str(self.output / "numba_cache"),
                        UTILITY_COMPARISON_NUMBA_CACHE=str(self.output / "numba_cache"))
        (self.output / "numba_cache").mkdir()
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
        for name in ("latest_completed.json", "best_so_far.json"):
            create(self.output / name, dict(status="no_completed_objective", arm=self.arm,
                                           contract_sha256=self.fingerprint))
        for signum in (signal.SIGTERM, signal.SIGINT):
            signal.signal(signum, self.interrupt)

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
        create(plan, dict(contract_sha256=self.fingerprint, arm=self.arm,
                          point=candidate["parameters"], deadline_epoch=cutoff,
                          objective_cap_seconds=cap, graphs=candidate["controller_stage"] == "smoke"))
        command = [sys.executable, self.runner, "--stage", "evaluate", "--contract",
                   str(self.contract_path), "--arm", self.arm, "--case-plan", str(plan),
                   "--output", str(output)]
        attempt = dict(candidate_id=case_id, stage=candidate["controller_stage"],
                       dispatch_epoch=time.time(), deadline_epoch=cutoff)
        self.attempts.append(attempt)
        create(self.output / "dispatch" / (case_id + ".json"), attempt)
        return ManagedProcess(command, self.output / (case_id + ".log"), cutoff, self.env)

    def finish(self, candidate, process, code):
        output = self.output / candidate["id"]
        record = dict(candidate_id=candidate["id"], stage=candidate["controller_stage"],
                      case_output=str(output), point=candidate["parameters"],
                      started_epoch=process.started_epoch, finished_epoch=time.time(),
                      returncode=code, status="timed_out" if process.timed_out else "failed")
        failure_path = output / "failure.json"
        if code != 0 and not process.timed_out and failure_path.is_file():
            try:
                failure = read(failure_path)
                if classified_inadmissible(failure, candidate["controller_stage"]):
                    record.update(status="inadmissible", loss=None, failure=failure)
            except (ValueError, OSError):
                pass  # Malformed failure receipts remain unknown failures.
        if code == 0 and not process.timed_out:
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
                              receipt_sha256=file_hash(case / "receipt.json"),
                              target_system_sha256=receipt["target_system_sha256"])
            except Exception as exc:
                record.update(error_type=type(exc).__name__, error=str(exc))
        self.records.append(record)
        create(self.output / "records" / (candidate["id"] + ".json"), record)
        write(self.output / "latest_completed.json", record)
        if record["status"] == "success" and candidate["controller_stage"] != "repeat":
            if self.best is None or record["loss"] < self.best["loss"]:
                self.best = copy.deepcopy(record)
                write(self.output / "best_so_far.json", self.best)
        self.summary()
        self.heartbeat(force=True)
        return record

    def batch(self, candidates, phase, deadline, workers=None):
        self.phase = phase
        for candidate in candidates:
            candidate["controller_stage"] = phase
        return run_batch(candidates, workers=workers or self.budget["workers_per_arm"],
                         deadline=deadline, launch=self.launch, finish=self.finish,
                         heartbeat=self.heartbeat, interrupted=lambda: self.interrupted)

    def external(self, command, name, deadline):
        self.phase = name
        if time.time() >= deadline or self.interrupted:
            return dict(status="not_run", reason="deadline_or_interruption")
        process = ManagedProcess(command, self.output / (name + ".log"), deadline, self.env)
        try:
            while process.poll() is None:
                if self.interrupted:
                    process.deadline = time.time()
                self.heartbeat()
                time.sleep(.5)
            return dict(status="success" if process.process.returncode == 0 else "failed",
                        returncode=process.process.returncode, timed_out=process.timed_out)
        finally:
            process.close()

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
        create(self.root / "barrier" / (self.arm + ".smoke.json"),
               dict(arm=self.arm, status="success" if passed else "failed", contract_sha256=self.fingerprint))
        if not passed:
            return None
        self.phase = "common_smoke_barrier"
        while time.time() < self.clock["search_cutoff_epoch"] and not self.interrupted:
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
        self.reused.append(reused)
        create(self.output / "initial_seed_reuse.json", reused)
        result = self.batch(population[1:], "initial", self.clock["search_cutoff_epoch"])
        records = [reused] + result["results"]
        if not result["complete"] or result["stop_reason"]:
            self.stop_reason = result["stop_reason"] or "incomplete_initial_population"
            return
        try:
            scores = complete_scores(population, records)
        except BarrierStop as exc:
            self.stop_reason = str(exc)
            return
        write(self.output / "population_000.json", dict(population=population, scores=scores))
        for generation in range(1, self.budget["de_generations"] + 1):
            trials = design.make_de_generation(self.catalog, self.arm, population, scores,
                generation=generation, budget=self.budget, rng_seed=self.bank["rng_seed"],
                mutation_factor=self.contract["differential_evolution"]["mutation_factor"])
            create(self.output / f"proposals_{generation:03d}.json", trials)
            result = self.batch(trials, "de", self.clock["search_cutoff_epoch"])
            if not result["complete"] or result["stop_reason"]:
                self.stop_reason = result["stop_reason"] or "incomplete_de_generation"
                return
            try:
                trial_scores = complete_scores(trials, result["results"])
            except BarrierStop as exc:
                self.stop_reason = str(exc)
                return
            population, scores = design.select_de_generation(population, scores, trials, trial_scores)
            create(self.output / f"population_{generation:03d}.json", dict(population=population, scores=scores))
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
        create(self.output / "selected.json", dict(contract_sha256=self.fingerprint, arm=self.arm,
            candidate_id=self.best["candidate_id"], point=self.best["point"], loss=self.best["loss"],
            original_case_output=self.best["case_output"],
            repeat_case_outputs=[str(self.output / row["id"]) for row in repeats],
            selected_epoch=time.time(), search_stop_reason=self.stop_reason))
        self.batch(repeats, "repeat", self.clock["repeat_cutoff_epoch"], workers=2)
        self.summary()
        command = [sys.executable, self.collector, "--stage", "collect", "--contract",
                   str(self.contract_path), "--arm", self.arm, "--run-root", str(self.root),
                   "--output", str(self.output / "export")]
        self.collection = self.external(command, "collection", self.clock["end_epoch"])

    def run(self):
        error = None
        try:
            self.readiness()
            smoke = self.smoke()
            if smoke is not None:
                self.smoke_verified = True
                self.search(smoke)
        except Exception as exc:
            error = dict(error_type=type(exc).__name__, error=str(exc))
            self.stop_reason = "controller_failure_no_retry"
            create(self.output / "controller_failure.json", error)
        finally:
            try:
                self.repeat_and_collect()
            except Exception as exc:
                self.collection = dict(status="failed", error_type=type(exc).__name__, error=str(exc))
            self.phase = "finished"
            summary = self.summary()
            self.heartbeat(force=True)
            create(self.output / "complete.json", summary)
        fully_scored = all(row["failed"] == row["timed_out"] == row["unrun"] == row["incomplete"] == 0
                           for row in summary["counts"].values())
        ok = (error is None and fully_scored and self.stop_reason == "finite_search_bank_complete"
              and (self.collection or {}).get("status") == "success")
        return 0 if ok else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--arm", choices=ARMS, required=True)
    parser.add_argument("--check-only", action="store_true")
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
