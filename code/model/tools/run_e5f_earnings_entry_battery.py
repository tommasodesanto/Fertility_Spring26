#!/usr/bin/env python3
"""Run one pinned task in the four-cell earnings-by-entry diagnostic battery."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

SCHEMA = "e5f_earnings_entry_battery_manifest_v1"
PARAMETERS = ("beta_annual", "kappa_fert", "kappa_fert_continuation", "chi", "H0",
              "theta0", "theta1", "first_birth_fixed_cost", "h_P")
OBJECTIVE_SCHEMA = "e5f_initial_minimum_distance_result_v1"
CASE_SECONDS = 3200
TOTAL_SECONDS = 3600


class ContractError(ValueError):
    pass


def read(path: Path) -> Any:
    return json.loads(path.read_text())


def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    tmp.replace(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def cases_for(manifest: dict[str, Any], stage: str, cell_id: str) -> list[dict[str, Any]]:
    if cell_id not in manifest.get("cells", {}) or stage not in ("smoke", "production"):
        raise ContractError("unknown stage or battery cell")
    raw = manifest.get(stage, {}).get("cases", [])
    rows = [row for row in raw if row.get("cell_id") == cell_id]
    expected = 1 if stage == "smoke" else 10
    if len(rows) != expected:
        raise ContractError(f"{stage}/{cell_id} must contain exactly {expected} worker(s)")
    ids = [row.get("task_id") for row in rows]
    if set(ids) != set(range(1, expected + 1)) or len(ids) != len(set(ids)):
        raise ContractError(f"{stage}/{cell_id} task IDs must be 1..{expected}")
    for row in rows:
        proposals = row.get("proposals")
        if not isinstance(proposals, list) or not 1 <= len(proposals) <= 6:
            raise ContractError("each task requires one to six predeclared proposals")
        if stage == "smoke" and len(proposals) != 1:
            raise ContractError("each smoke task must run exactly one proposal")
        seen_ids: set[str] = set()
        for p in proposals:
            for key in ("proposal_id", "case_id", "plan_path", "plan_sha256"):
                if not isinstance(p.get(key), str) or not p[key]:
                    raise ContractError(f"proposal requires {key}")
            if p["proposal_id"] in seen_ids:
                raise ContractError("proposal IDs must be unique per worker")
            seen_ids.add(p["proposal_id"])
    return sorted(rows, key=lambda x: x["task_id"])


def validate_manifest(manifest: dict[str, Any]) -> None:
    if manifest.get("schema") != SCHEMA or manifest.get("status") != "ready":
        raise ContractError("manifest must be ready and use the battery schema")
    cells = manifest.get("cells", {})
    if not isinstance(cells, dict) or len(cells) != 4:
        raise ContractError("manifest must define the four named earnings-by-entry cells")
    combinations = {(c.get("earnings_specification"), c.get("entry_rule")) for c in cells.values() if isinstance(c, dict)}
    income_levels = {x[0] for x in combinations}
    entry_levels = {x[1] for x in combinations}
    if (len(combinations) != 4 or len(income_levels) != 2 or len(entry_levels) != 2
            or combinations != {(income, entry) for income in income_levels for entry in entry_levels}):
        raise ContractError("cells must cover every combination of two income mappings and two entry rules")
    if any(c.get("income_mapping") != "direct_period" for c in cells.values()):
        raise ContractError("both earnings candidates must use the direct-period constructor")
    state_pairs = {(int(c.get("persistent_states", -1)), int(c.get("iid_states", -1)))
                   for c in cells.values()}
    labels_by_grid = {(c.get("earnings_specification"), int(c.get("persistent_states", -1)), int(c.get("iid_states", -1)))
                      for c in cells.values()}
    if state_pairs != {(7, 1), (7, 3)} or len(labels_by_grid) != 2:
        raise ContractError("earnings cells must distinguish the 7x1 and 7x3 period-income grids")
    if manifest.get("workers_per_cell") != 10 or manifest.get("case_seconds") != CASE_SECONDS or manifest.get("total_seconds") != TOTAL_SECONDS:
        raise ContractError("worker counts or time budgets differ from the authorized design")
    minimum_next = float(manifest.get("minimum_next_proposal_seconds", 0))
    if not math.isfinite(minimum_next) or not 1 <= minimum_next <= CASE_SECONDS:
        raise ContractError("minimum_next_proposal_seconds must lie in [1,3200]")
    if not manifest.get("objective_canonical_sha256") or not manifest.get("target_system_sha256"):
        raise ContractError("objective and target-system fingerprints are required")
    if manifest.get("scored_moment_count") != 13 or manifest.get("parameter_row_count") != 17:
        raise ContractError("the unchanged objective requires 13 moments and 17 parameter rows")
    for stage in ("smoke", "production"):
        rows = manifest.get(stage, {}).get("cases", [])
        expected = 4 if stage == "smoke" else 40
        if len(rows) != expected:
            raise ContractError(f"{stage} must contain exactly {expected} tasks")
        for cell_id in cells:
            cases_for(manifest, stage, cell_id)


def task_case(manifest: dict[str, Any], stage: str, cell_id: str, task_id: int) -> dict[str, Any]:
    rows = cases_for(manifest, stage, cell_id)
    if task_id < 1 or task_id > len(rows):
        raise ContractError("task ID is outside this cell array")
    row = rows[task_id - 1]
    if row["task_id"] != task_id:
        raise ContractError("task ID mismatch")
    return row


def verify_source_inventory(plan: dict[str, Any], cell: dict[str, Any],
                            verified: set[tuple[str, str]] | None = None) -> dict[str, str]:
    files = plan.get("files", {})
    pin = files.get("initial_contract", {})
    path = Path(pin.get("path", ""))
    if not path.is_file() or sha256(path) != pin.get("sha256"):
        raise ContractError("initial contract path/hash mismatch")
    initial = read(path)
    inventory = initial.get("source_sha256")
    if not isinstance(inventory, dict) or len(inventory) != int(cell.get("source_file_count", 641)):
        raise ContractError("full native source inventory missing or wrong size")
    if canonical_sha256(inventory) != cell.get("source_inventory_sha256"):
        raise ContractError("source inventory fingerprint differs from cell contract")
    root = Path(plan.get("source_root", ""))
    key = (str(root.resolve()), canonical_sha256(inventory))
    if verified is None or key not in verified:
        for relative, expected in inventory.items():
            source = root / relative
            if not source.is_file() or sha256(source) != expected:
                raise ContractError(f"native source pin mismatch: {relative}")
        if verified is not None:
            verified.add(key)
    return inventory


def validate_plan(manifest: dict[str, Any], cell_id: str, proposal: dict[str, Any],
                  verified_sources: set[tuple[str, str]] | None = None) -> tuple[Path, dict[str, Any]]:
    cell = manifest["cells"][cell_id]
    path = Path(proposal["plan_path"]).resolve()
    if not path.is_file() or sha256(path) != proposal["plan_sha256"]:
        raise ContractError("candidate plan missing or hash mismatch")
    plan = read(path)
    if plan.get("objective_canonical_sha256") != manifest["objective_canonical_sha256"]:
        raise ContractError("candidate objective fingerprint mismatch")
    if plan.get("target_system_sha256") != manifest["target_system_sha256"]:
        raise ContractError("candidate target-system fingerprint mismatch")
    params = plan.get("structural_parameters", plan.get("starting_structural_parameters", {}))
    if set(params) != set(PARAMETERS):
        raise ContractError("candidate must pin the same nine structural coordinates")
    if len(plan.get("cases", [])) != 1:
        raise ContractError("each pinned candidate plan must contain exactly one case")
    native_case = plan["cases"][0]
    if native_case.get("id") != proposal["case_id"] or native_case.get("repetitions") != 1:
        raise ContractError("candidate case ID/repetition contract mismatch")
    if plan.get("income_specification", {}).get("mapping") != cell.get("income_mapping"):
        raise ContractError("candidate income mapping differs from battery cell")
    constructor = plan.get("income_specification", {}).get("constructor_arguments", {})
    if (int(constructor.get("n_persistent", -1)) != int(cell.get("persistent_states", -2))
            or int(constructor.get("n_iid", -1)) != int(cell.get("iid_states", -2))):
        raise ContractError("candidate income-grid dimensions differ from earnings cell")
    if plan.get("entry_specification", {}).get("rule") != cell.get("entry_rule"):
        raise ContractError("candidate entry rule differs from battery cell")
    verify_source_inventory(plan, cell, verified_sources)
    adapter = Path(plan.get("adapter_path", ""))
    if not adapter.is_file() or plan.get("adapter_sha256") != sha256(adapter):
        raise ContractError("candidate adapter path/hash mismatch")
    for label, item in plan.get("files", {}).items():
        if not isinstance(item, dict) or not item.get("path") or not item.get("sha256"):
            raise ContractError(f"plan file pin malformed: {label}")
        pinned = Path(item["path"])
        if not pinned.is_file() or sha256(pinned) != item["sha256"]:
            raise ContractError(f"plan runtime/target file pin mismatch: {label}")
    # The native wrapper remains responsible for all frozen accounting, value,
    # market-clearing, feasibility, and objective-specific gates.
    return path, plan


def kill_group(proc: subprocess.Popen[Any]) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        deadline = time.monotonic() + 5
        while proc.poll() is None and time.monotonic() < deadline:
            time.sleep(.1)
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, OSError):
        pass


def solve_count(path: Path) -> int | None:
    files = list(path.glob("evaluation/raw/repetition_*/stationary_solves.json"))
    if not files:
        return None
    count = 0
    for file in files:
        try:
            raw = read(file)
            count += len(raw) if isinstance(raw, list) else 1
        except (OSError, ValueError, TypeError):
            continue
    return count


def update_rollups(stage_root: Path, result: dict[str, Any]) -> None:
    stage_root.mkdir(parents=True, exist_ok=True)
    with (stage_root / ".rollup.lock").open("w") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        write(stage_root / "latest.json", result)
        completed = []
        for path in stage_root.glob("**/status.json"):
            try:
                item = read(path)
            except (OSError, ValueError):
                continue
            if item.get("status") == "completed" and math.isfinite(float(item.get("loss", "nan"))):
                completed.append(item)
        if completed:
            write(stage_root / "best_so_far.json", min(completed, key=lambda x: float(x["loss"])))


def execute_task(manifest: dict[str, Any], cell_id: str, row: dict[str, Any], stage: str,
                 output_root: Path, *, popen=subprocess.Popen, monotonic=time.monotonic,
                 sleep=time.sleep, killer=kill_group) -> dict[str, Any]:
    worker_dir = output_root / stage / cell_id / f"worker_{row['task_id']:02d}"
    worker_dir.mkdir(parents=True, exist_ok=False)
    started, begin = time.time(), monotonic()
    write(worker_dir / "status.json", {"status": "started", "cell_id": cell_id, "task_id": row["task_id"], "started": started})
    update_rollups(output_root / stage / cell_id, read(worker_dir / "status.json"))
    results = []
    verified_sources: set[tuple[str, str]] = set()
    previous_seconds = 0.0
    try:
        minimum_next = float(manifest.get("minimum_next_proposal_seconds", 1.0))
        if not math.isfinite(minimum_next) or not 1.0 <= minimum_next <= CASE_SECONDS:
            raise ContractError("minimum_next_proposal_seconds must lie in [1,3200]")
        if stage == "production":
            smoke_path = output_root / "smoke" / cell_id / "worker_01" / "status.json"
            if not smoke_path.is_file():
                raise ContractError("cell smoke receipt is missing")
            smoke = read(smoke_path)
            smoke_seconds = float(smoke.get("elapsed_seconds", 0.0))
            if smoke.get("status") != "completed" or not math.isfinite(smoke_seconds) or smoke_seconds <= 0:
                raise ContractError("cell smoke receipt did not verify successfully")
            minimum_next = max(minimum_next, smoke_seconds * 1.1)
        for idx, proposal in enumerate(row["proposals"]):
            elapsed = monotonic() - begin
            remaining = TOTAL_SECONDS - elapsed
            reserve = CASE_SECONDS if idx == 0 else max(minimum_next, previous_seconds * 1.1)
            if remaining < reserve:
                break
            case_dir = worker_dir / f"proposal_{idx+1:02d}_{proposal['proposal_id']}"
            case_dir.mkdir(parents=True, exist_ok=False)
            plan_path, plan = validate_plan(manifest, cell_id, proposal, verified_sources)
            ready_elapsed = monotonic() - begin
            ready_remaining = TOTAL_SECONDS - ready_elapsed
            ready_reserve = CASE_SECONDS if idx == 0 else max(minimum_next, previous_seconds * 1.1)
            if ready_remaining < ready_reserve:
                write(case_dir / "status.json", {"status": "not_started_budget_stop", "cell_id": cell_id,
                    "task_id": row["task_id"], "proposal_id": proposal["proposal_id"],
                    "remaining_seconds": ready_remaining, "required_seconds": ready_reserve})
                break
            adapter = Path(plan["adapter_path"]).resolve()
            case_started = time.time()
            write(case_dir / "status.json", {"status": "started", "cell_id": cell_id, "task_id": row["task_id"],
                "proposal_id": proposal["proposal_id"], "plan_sha256": proposal["plan_sha256"], "started": case_started})
            cmd = [sys.executable, str(adapter), "--plan", str(plan_path), "--output", str(case_dir / "result"),
                   "--arm", plan["cases"][0]["arm"], "--repetitions", "1"]
            log = (case_dir / "worker.log").open("w")
            proc = popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            solve_start = monotonic()
            deadline = solve_start + min(CASE_SECONDS, ready_remaining)
            next_beat = 0.0
            try:
                while proc.poll() is None:
                    now = monotonic()
                    if now >= deadline:
                        raise TimeoutError("candidate exceeded its bounded solve time")
                    if now >= next_beat:
                        heartbeat = {"status": "running", "cell_id": cell_id, "task_id": row["task_id"],
                            "proposal_id": proposal["proposal_id"], "updated": time.time(), "pid": proc.pid,
                            "elapsed_seconds": now - solve_start}
                        write(case_dir / "heartbeat.json", heartbeat)
                        write(worker_dir / "heartbeat.json", heartbeat)
                        next_beat = now + 30
                    sleep(min(2.0, deadline-now))
            except BaseException:
                killer(proc)
                raise
            finally:
                log.close()
            if proc.returncode != 0:
                raise RuntimeError(f"candidate adapter exited {proc.returncode}")
            score_path = case_dir / "result/evaluation/scored_repetition_01/score.json"
            summary_path = case_dir / "result/evaluation/summary.json"
            if not score_path.is_file() or not summary_path.is_file():
                raise ContractError("candidate omitted score or wrapper summary")
            score, summary = read(score_path), read(summary_path)
            if score.get("schema") != OBJECTIVE_SCHEMA or len(score.get("target_fit", [])) != 13 or len(score.get("parameters", [])) != 17:
                raise ContractError("score schema or full target/parameter table dimensions differ")
            structural = {item.get("parameter") for item in score["parameters"] if item.get("structural_coordinate")}
            if structural != set(PARAMETERS):
                raise ContractError("score structural coordinates differ from the nine-parameter contract")
            if summary.get("status") != "verified_scored_candidate" or summary.get("objective_canonical_sha256") != manifest["objective_canonical_sha256"] or summary.get("repetitions") != 1:
                raise ContractError("unchanged native score/quality gate failed")
            loss = float(score.get("loss"))
            if not math.isfinite(loss):
                raise ContractError("candidate loss is not finite")
            duration = monotonic() - solve_start
            receipt = {"status": "completed", "cell_id": cell_id, "task_id": row["task_id"],
                "proposal_id": proposal["proposal_id"], "objective": loss, "target_count": 13,
                "parameter_count": 17, "stationary_solve_count": solve_count(case_dir / "result"),
                "elapsed_seconds": duration, "score_sha256": sha256(score_path),
                "wrapper_summary_sha256": sha256(summary_path), "completed": time.time()}
            write(case_dir / "status.json", {**receipt, "loss": loss})
            results.append(receipt)
            previous_seconds = duration
            write(worker_dir / "latest.json", receipt)
            best_path = worker_dir / "best_so_far.json"
            best = read(best_path) if best_path.is_file() else None
            if best is None or loss < float(best["objective"]):
                write(best_path, receipt)
            update_rollups(output_root / stage / cell_id, {**receipt, "loss": loss})
        elapsed = monotonic() - begin
        if not results:
            raise ContractError("worker completed zero scored proposals")
        final = {"status": "completed", "cell_id": cell_id, "task_id": row["task_id"],
            "started": started, "completed": time.time(), "elapsed_seconds": elapsed,
            "proposal_count": len(results), "proposals": results,
            "stop_reason": "proposal_pool_exhausted" if len(results) == len(row["proposals"]) else "insufficient_budget_for_next_proposal"}
        write(worker_dir / "status.json", final)
        update_rollups(output_root / stage / cell_id, final)
        return final
    except BaseException as exc:
        failure = {"status": "incomplete", "cell_id": cell_id, "task_id": row["task_id"],
            "started": started, "completed": time.time(), "elapsed_seconds": monotonic()-begin,
            "proposal_count": len(results), "proposals": results, "error_type": type(exc).__name__, "error": str(exc)}
        write(worker_dir / "status.json", failure)
        update_rollups(output_root / stage / cell_id, failure)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--stage", choices=("smoke", "production"), required=True)
    parser.add_argument("--cell-id", required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    manifest = read(args.manifest.resolve())
    validate_manifest(manifest)
    row = task_case(manifest, args.stage, args.cell_id, args.task_id)
    print(json.dumps(execute_task(manifest, args.cell_id, row, args.stage, args.output_root.resolve()), sort_keys=True))


if __name__ == "__main__":
    main()
