"""Bounded sequential controller for the earnings and wealth smoke."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

SCHEMA = "e5f_earnings_wealth_smoke_v1"
SCORE_SCHEMA = "e5f_initial_minimum_distance_result_v1"
MAX_SECONDS = 7200.0
POLL_SECONDS = 0.5
HEARTBEAT_SECONDS = 30.0
_CASE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class ContractError(ValueError):
    pass


def read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def validate_plan(plan: dict[str, Any], plan_path: Path) -> None:
    if plan.get("schema") != SCHEMA:
        raise ContractError("wrong smoke plan schema")
    total_seconds = plan.get("total_seconds")
    if not finite(total_seconds) or float(total_seconds) <= 0 or float(total_seconds) > MAX_SECONDS:
        raise ContractError("total_seconds must be positive and at most 7200")
    cases = plan.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ContractError("cases must be a nonempty list")
    adapter = Path(plan.get("adapter_path", ""))
    if not adapter.is_file() or plan.get("adapter_sha256") != sha256(adapter):
        raise ContractError("adapter path/hash invalid")
    if adapter.resolve() == plan_path.resolve():
        raise ContractError("plan and adapter must differ")
    seen: set[str] = set()
    budget = 0.0
    for case in cases:
        case_id = case.get("id") if isinstance(case, dict) else None
        if not isinstance(case_id, str) or not _CASE_ID.fullmatch(case_id) or case_id in seen:
            raise ContractError("invalid or duplicate case id")
        seen.add(case_id)
        if not isinstance(case.get("arm"), str) or not case["arm"]:
            raise ContractError(f"case {case_id} needs an arm")
        repetitions = case.get("repetitions")
        if not isinstance(repetitions, int) or isinstance(repetitions, bool) or repetitions not in (1, 2):
            raise ContractError(f"case {case_id} repetitions must be 1 or 2")
        seconds = case.get("seconds")
        if not finite(seconds) or float(seconds) <= 0:
            raise ContractError(f"case {case_id} needs a positive seconds budget")
        budget += float(seconds)
    if budget > float(total_seconds):
        raise ContractError("case budgets exceed total_seconds")


def _kill_group(process: subprocess.Popen[Any]) -> None:
    try:
        import psutil
        parent = psutil.Process(process.pid)
        descendants = parent.children(recursive=True)
        for child in reversed(descendants):
            try:
                child.terminate()
            except psutil.Error:
                pass
        _, alive = psutil.wait_procs(descendants, timeout=2)
        for child in alive:
            try:
                child.kill()
            except psutil.Error:
                pass
    except Exception:
        # The adapter may exit while descendants are being enumerated; still
        # continue to process-group cleanup below.
        pass
    try:
        os.killpg(process.pid, signal.SIGTERM)
        deadline = time.monotonic() + 2.0
        while process.poll() is None and time.monotonic() < deadline:
            time.sleep(0.1)
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
    except (ProcessLookupError, OSError):
        pass


def _stationary_counts(case_dir: Path) -> int | None:
    counts: list[int] = []
    for path in case_dir.glob("evaluation/raw/repetition_*/stationary_solves.json"):
        try:
            payload = json.loads(path.read_text())
            if isinstance(payload, list):
                items = payload
            elif isinstance(payload, dict):
                items = [payload]
            else:
                continue
            counts.append(len(items))
        except (OSError, json.JSONDecodeError):
            continue
    return sum(counts) if counts else None


def _score_path(case_dir: Path) -> Path:
    paths = (
        case_dir / "evaluation/scored_repetition_01/score.json",
        case_dir / "evaluation/score.json",
    )
    for path in paths:
        if path.is_file():
            return path
    raise ContractError("missing scored repetition score.json")


def validate_case_receipt(case_dir: Path, plan: dict[str, Any], case: dict[str, Any]) -> dict[str, Any]:
    evaluation = case_dir / "evaluation"
    summary_path = evaluation / "summary.json"
    if not summary_path.is_file():
        raise ContractError("missing evaluation/summary.json")
    summary = read(summary_path)
    score = read(_score_path(case_dir))
    if score.get("schema") != SCORE_SCHEMA:
        raise ContractError("native score schema mismatch")
    expected_contract = plan.get("objective_canonical_sha256", plan.get("contract_sha256"))
    if not expected_contract or score.get("contract_sha256") != expected_contract:
        raise ContractError("native objective contract mismatch")
    expected_sources = plan.get("source_fingerprints")
    if expected_sources is not None and score.get("source_fingerprints") != expected_sources:
        raise ContractError("native source fingerprint mismatch")
    if not finite(score.get("loss")):
        raise ContractError("native score loss is not finite")
    if len(score.get("target_fit", [])) != 13 or len(score.get("parameters", [])) != 17:
        raise ContractError("native receipt must contain 13 targets and 17 parameters")
    if summary.get("status") != "verified_scored_candidate":
        raise ContractError("native summary is not a verified success")
    if summary.get("case_id") != case["id"]:
        raise ContractError("native case id mismatch")
    if summary.get("repetitions") != case["repetitions"]:
        raise ContractError("native repetitions mismatch")
    if case["repetitions"] == 2 and summary.get("exact_loss_equality") is not True:
        raise ContractError("native exact repetition verification is absent")
    return {
        "loss": float(score["loss"]),
        "stationary_solves": _stationary_counts(case_dir),
    }


def run_plan(plan_path: Path, output: Path, *, preflight: bool = False) -> dict[str, Any]:
    try:
        import psutil  # noqa: F401
    except ImportError as error:
        raise ContractError("psutil is required for descendant cleanup") from error
    plan = read(plan_path)
    validate_plan(plan, plan_path)
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    last_heartbeat = started - HEARTBEAT_SECONDS
    results: list[dict[str, Any]] = []
    best_loss = math.inf
    heartbeat_path = output / "heartbeat.json"

    def beat(status: str, force: bool = False, **extra: Any) -> None:
        nonlocal last_heartbeat
        now = time.monotonic()
        if force or now - last_heartbeat >= HEARTBEAT_SECONDS:
            write(heartbeat_path, {
                "schema": SCHEMA,
                "status": status,
                "elapsed_seconds": now - started,
                "completed_cases": len(results),
                **extra,
            })
            last_heartbeat = now

    env = {**os.environ, "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    current_process: subprocess.Popen[Any] | None = None
    current_case: dict[str, Any] | None = None
    beat("started", force=True)
    if preflight:
        command = [
            plan.get("python", sys.executable), str(plan["adapter_path"]), "--preflight",
            "--plan", str(plan_path), "--output", str(output / "preflight"),
        ]
        completed = subprocess.run(command, env=env, check=False)
        if completed.returncode:
            beat("preflight_failed", force=True, returncode=completed.returncode)
            raise RuntimeError("adapter preflight failed")
        result = {
            "schema": SCHEMA,
            "status": "preflight_passed",
            "elapsed_seconds": time.monotonic() - started,
            "cases": [],
            "evaluations": 0,
        }
        write(output / "receipt.json", result)
        beat("preflight_passed", force=True)
        return result

    try:
        for index, case in enumerate(plan["cases"], 1):
            current_case = case
            if time.monotonic() - started > float(plan["total_seconds"]):
                raise TimeoutError("total smoke budget exhausted")
            case_dir = output / case["id"]
            command = [
                plan.get("python", sys.executable), str(plan["adapter_path"]), "--mode", "pilot",
                "--plan", str(plan_path), "--arm", case["arm"], "--output", str(case_dir),
                "--repetitions", str(case["repetitions"]),
            ]
            process = subprocess.Popen(command, env=env, start_new_session=True)
            current_process = process
            case_started = time.monotonic()
            beat("running", force=True, case_id=case["id"], case_index=index)
            while process.poll() is None:
                if time.monotonic() - case_started > float(case["seconds"]):
                    _kill_group(process)
                    current_process = None
                    raise TimeoutError(f"case timed out: {case['id']}")
                if time.monotonic() - started > float(plan["total_seconds"]):
                    _kill_group(process)
                    current_process = None
                    raise TimeoutError("total smoke budget exhausted")
                beat("running", case_id=case["id"], case_index=index)
                time.sleep(POLL_SECONDS)
            if process.returncode:
                current_process = None
                raise RuntimeError(f"case failed: {case['id']} (returncode {process.returncode})")
            current_process = None
            receipt = validate_case_receipt(case_dir, plan, case)
            item = {
                "id": case["id"], "arm": case["arm"], "repetitions": case["repetitions"],
                "status": "succeeded", **receipt,
            }
            results.append(item)
            current_case = None
            write(output / "latest.json", item)
            if receipt["loss"] < best_loss:
                best_loss = receipt["loss"]
                write(output / "best.json", {
                    "status": "lowest_loss_diagnostic_arm", **item,
                })
            beat("case_completed", force=True, case_id=case["id"])
    except (TimeoutError, RuntimeError, ContractError) as error:
        if current_case is not None:
            results.append({
                "id": current_case["id"], "arm": current_case["arm"],
                "repetitions": current_case["repetitions"], "status": "incomplete",
                "incomplete": True, "unknown_in_flight_solves": True,
                "stationary_solves": _stationary_counts(output / current_case["id"]),
            })
        failure = {
            "schema": SCHEMA,
            "status": "failed",
            "error": type(error).__name__,
            "message": str(error),
            "elapsed_seconds": time.monotonic() - started,
            "cases": results,
        }
        write(output / "latest.json", failure)
        write(output / "receipt.json", failure)
        beat("timeout" if isinstance(error, TimeoutError) else "failed", force=True,
             message=str(error))
        raise
    except BaseException as error:
        if current_process is not None and current_process.poll() is None:
            _kill_group(current_process)
        if current_case is not None:
            failure = {"schema": SCHEMA, "status": "interrupted",
                       "error": type(error).__name__, "message": str(error),
                       "elapsed_seconds": time.monotonic() - started,
                       "cases": results + [{"id": current_case["id"],
                                            "arm": current_case["arm"],
                                            "repetitions": current_case["repetitions"],
                                            "status": "incomplete", "incomplete": True,
                                            "unknown_in_flight_solves": True,
                                            "stationary_solves": _stationary_counts(output / current_case["id"])}]}
            write(output / "latest.json", failure)
            write(output / "receipt.json", failure)
        raise
    result = {
        "schema": SCHEMA,
        "status": "completed",
        "elapsed_seconds": time.monotonic() - started,
        "cases": results,
        "best_loss": best_loss,
        "selection_label": "lowest_loss_diagnostic_arm",
    }
    write(output / "receipt.json", result)
    beat("completed", force=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run_plan(args.plan, args.output, preflight=args.preflight), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
