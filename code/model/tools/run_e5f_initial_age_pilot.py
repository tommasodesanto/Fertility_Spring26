#!/usr/bin/env python3
"""Run the bounded six-case age-profile pilot under the joint rebated closure."""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
from typing import Any
from unittest.mock import patch

import e5f_initial_age_profile as age


REBATED_HELPER_SHA256 = "d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790"
JOINT_ADAPTER_SHA256 = "9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb"
SCORED_EVALUATOR_SHA256 = "9da35b55466d74dc10a6a85ff91dabe887a807aab68417eabf63a9d7a42ca967"
AGE_HELPER_SHA256 = "5ef6be93fbf28f929d0e09f3eb4ffbe932d61aee5cd66eb707efdea6acd0d801"
CPS_AGE_DATA_SHA256 = "b415b75c916a61113f1ae054ba0baadecd35a26a00c88a21451262937271f210"
MAX_SECONDS = 10_800
MAX_WORKERS = 2
PRIMARY_CASES = 6
FINAL_RESERVE_SECONDS = 2_100


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path: Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def load_module(name: str, path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(name, Path(path).resolve())
    if specification is None or specification.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def require_source_pins(helper: Path, joint: Path, evaluator: Path, data: Path) -> None:
    expected = (
        (helper, REBATED_HELPER_SHA256, "rebated helper"),
        (joint, JOINT_ADAPTER_SHA256, "joint adapter"),
        (evaluator, SCORED_EVALUATOR_SHA256, "scored evaluator"),
        (Path(age.__file__).resolve(), AGE_HELPER_SHA256, "age-profile helper"),
        (data, CPS_AGE_DATA_SHA256, "CPS age-profile packet"),
    )
    for path, expected_sha, label in expected:
        if digest(Path(path)) != expected_sha:
            raise ValueError(f"fixed {label} source pin changed")


def proposals(initial: dict[str, Any], objective: dict[str, Any]) -> list[dict[str, Any]]:
    """Return baseline and five clipped half-percent-of-bound-span proposals."""
    center = {name: float(value) for name, value in initial["structural_candidate"].items()}
    restrictions = {row["parameter"]: row for row in objective["parameter_restrictions"]}
    required = ("first_birth_fixed_cost", "kappa_fert", "kappa_fert_continuation")
    if any(name not in center or name not in restrictions for name in required):
        raise ValueError("candidate or objective omits an age-pilot coordinate")
    result = [{"case_id": "baseline", "parameters": center.copy()}]
    changes = (("first_birth_fixed_cost", -1.0), ("first_birth_fixed_cost", 1.0),
               ("kappa_fert", -1.0), ("kappa_fert", 1.0),
               ("kappa_fert_continuation", 1.0))
    for index, (name, direction) in enumerate(changes, 1):
        restriction = restrictions[name]
        lower, upper = float(restriction["lower"]), float(restriction["upper"])
        requested = center[name] + direction * 0.005 * (upper - lower)
        candidate = center.copy()
        candidate[name] = min(upper, max(lower, requested))
        result.append({
            "case_id": f"proposal_{index:02d}_{name}_{'plus' if direction > 0 else 'minus'}",
            "parameters": candidate, "changed_parameter": name,
            "requested_span_fraction": direction * 0.005,
            "actual_change": candidate[name] - center[name],
            "clipped_to_existing_bounds": abs(candidate[name] - requested) > 1e-14,
        })
    if len(result) != PRIMARY_CASES:
        raise RuntimeError("pilot must contain baseline plus five proposals")
    return result


def extra_score(packet: dict[str, Any], data_path: Path) -> dict[str, Any]:
    original = age._empirical_profiles
    with patch.object(age, "_empirical_profiles", lambda: original(Path(data_path))):
        score = age.score_extra(packet)
    score["metadata"]["empirical_source"] = str(Path(data_path).resolve())
    score["metadata"]["empirical_source_sha256"] = digest(Path(data_path))
    return score


def evaluator_command(evaluator: Path, helper: Path, joint: Path, template: Path,
                      proposal: Path, output: Path) -> list[str]:
    return [sys.executable, "-B", str(Path(evaluator).resolve()),
            "--helper", str(Path(helper).resolve()), "--helper-sha256", REBATED_HELPER_SHA256,
            "--joint", str(Path(joint).resolve()), "--joint-sha256", JOINT_ADAPTER_SHA256,
            "--template", str(Path(template).resolve()), "--output", str(Path(output).resolve()),
            "--proposal", str(Path(proposal).resolve())]


def validate_smoke(summary_path: Path, packet: dict[str, Any],
                   restrictions: dict[str, Any]) -> dict[str, Any]:
    """Verify that the supplied smoke used this closure and full original objective."""
    summary_path = Path(summary_path).resolve()
    root = summary_path.parent
    summary, launch = read(summary_path), read(root / "launch_contract.json")
    if (summary.get("status") != "verified_rebated_initial_smoke"
            or summary.get("method") != "joint_price_fertility_rebate"
            or launch.get("helper_sha256") != REBATED_HELPER_SHA256
            or launch.get("joint_sha256") != JOINT_ADAPTER_SHA256
            or launch.get("objective_sha256") != packet["run"]["working_objective"]["canonical_sha256"]
            or Path(summary.get("source_root", "")).resolve() != packet["source_root"]):
        raise ValueError("--smoke is not the pinned verified joint rebated initial smoke")
    accounting = read(root / "case/evaluation/rebate_accounting.json")
    score = read(root / "case/evaluation/scored_repetition_01/score.json")
    if accounting.get("status") != "verified" or len(accounting.get("receipts", [])) != 1:
        raise ValueError("--smoke has no complete independent rebate-accounting receipt")
    if score.get("loss") != summary.get("loss") or summary.get("checkpoint") != accounting["receipts"][-1]:
        raise ValueError("--smoke summary, original score, and checkpoint disagree")
    controller_path = Path(packet["plan_path"]).parent / "run_capped_beta.py"
    controller = load_module("age_pilot_saved_controller", controller_path)
    controller.validate_score(score, launch["proposal"], restrictions)
    return {"summary": summary, "launch": launch, "score": score,
            "accounting": accounting["receipts"], "summary_sha256": digest(summary_path)}


def score_signature(score: dict[str, Any]) -> str:
    retained = {"loss": score["loss"], "target_fit": score["target_fit"],
                "parameters": score["parameters"]}
    return hashlib.sha256(json.dumps(retained, sort_keys=True, separators=(",", ":"),
                                 allow_nan=False).encode()).hexdigest()


def augmented_target_rows(score: dict[str, Any], extra: dict[str, Any]) -> list[dict[str, Any]]:
    """Combine the twelve scored original rows with the six experimental rows."""
    rows = [dict(system="original_12", **row) for row in score["target_fit"]
            if row.get("scored") is True]
    rows.extend(dict(system="experimental_extra_6", **row) for row in extra["rows"])
    if len(rows) != 18:
        raise RuntimeError("selected target table is not 12 scored plus 6 experimental rows")
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--helper", type=Path, required=True)
    parser.add_argument("--joint", type=Path, required=True)
    parser.add_argument("--evaluator", type=Path,
                        default=Path(__file__).with_name("run_e5f_joint_rebated_initial_scored.py"))
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=MAX_WORKERS)
    parser.add_argument("--seconds", type=int, default=MAX_SECONDS)
    arguments = parser.parse_args()
    if arguments.workers != MAX_WORKERS:
        raise ValueError("this bounded pilot requires exactly two workers")
    if not FINAL_RESERVE_SECONDS < arguments.seconds <= MAX_SECONDS:
        raise ValueError(f"--seconds must lie in [{FINAL_RESERVE_SECONDS + 1},{MAX_SECONDS}]")

    # Every source, target, objective, and smoke gate precedes output creation.
    helper_path, joint_path = arguments.helper.resolve(), arguments.joint.resolve()
    evaluator_path, template_path = arguments.evaluator.resolve(), arguments.template.resolve()
    data_path = arguments.data.resolve()
    require_source_pins(helper_path, joint_path, evaluator_path, data_path)
    helper = load_module("age_pilot_rebated_helper", helper_path)
    packet = helper.saved_packet(template_path)
    restrictions = helper.validate_scientific_contract(packet)
    smoke = validate_smoke(arguments.smoke, packet, restrictions)

    seed = copy.deepcopy(packet["plan"]["resume_proposal"])
    cases = proposals({"structural_candidate": seed["parameters"]}, packet["objective"])
    for item in cases:
        item.update(initial_psi=float(seed["initial_psi"]), repetitions=1)
    objective_hash = packet["run"]["working_objective"]["canonical_sha256"]
    components = {
        "original_objective_canonical_sha256": objective_hash,
        "original_objective_file_sha256": digest(packet["objective_path"]),
        "cps_age_profile_sha256": digest(data_path),
        "extra_age_helper_sha256": digest(Path(age.__file__).resolve()),
        "rebated_helper_sha256": digest(helper_path), "joint_adapter_sha256": digest(joint_path),
        "scored_evaluator_sha256": digest(evaluator_path),
        "verified_joint_rebated_smoke_summary_sha256": smoke["summary_sha256"],
        "extra_weight_rule": "diagonal inverse variance; synthetic scale=0.05*abs(target); not SE; overlap covariance ignored",
    }
    fingerprint = hashlib.sha256(json.dumps(components, sort_keys=True,
                                             separators=(",", ":")).encode()).hexdigest()
    output = arguments.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); deadline = started + arguments.seconds
    primary_deadline = deadline - FINAL_RESERVE_SECONDS
    state = {"status": "running", "phase": "primary_cases", "completed": 0, "failed": 0,
             "primary_deadline_seconds": arguments.seconds - FINAL_RESERVE_SECONDS,
             "exact_repetition_reserve_seconds": FINAL_RESERVE_SECONDS}
    records: list[dict[str, Any]] = []
    lock, stop = threading.Lock(), threading.Event()

    def compact(row):
        return {key: value for key, value in row.items() if key not in {"score", "extra_score"}}

    def persist():
        compact_records = [compact(row) for row in records]
        valid = [row for row in compact_records if row["status"] == "verified"]
        best = min(valid, key=lambda row: row["augmented_loss"]) if valid else None
        write(output / "latest_completed.json", {**state, "latest": compact_records[-1] if records else None})
        write(output / "best_so_far.json", {"status": "available" if best else "none", "best": best})
        write(output / "cases.json", compact_records)

    def heartbeat():
        while not stop.wait(30.0):
            write(output / "heartbeat.json", {**state, "elapsed_seconds": time.monotonic() - started,
                                                "deadline_seconds": arguments.seconds})

    write(output / "contract.json", {"status": "running", "fingerprint": fingerprint, **components,
          "template": str(template_path), "smoke": str(arguments.smoke.resolve()),
          "workers": arguments.workers, "seconds": arguments.seconds, "primary_cases": cases,
          "original_scored_rows": 12, "experimental_extra_rows": 6, "production_eligible": False})
    persist(); threading.Thread(target=heartbeat, daemon=True).start()

    def run_case(item: dict[str, Any], case_deadline: float) -> dict[str, Any]:
        case_id = item["case_id"]
        case_dir = output / "cases" / case_id; case_dir.mkdir(parents=True, exist_ok=False)
        write(case_dir / "phase.json", {"status": "running", "phase": "joint_rebated_original_score"})
        proposal_path = case_dir / "proposal.json"; write(proposal_path, item)
        evaluation = case_dir / "evaluation"
        remaining = case_deadline - time.monotonic()
        if remaining <= 0:
            result = {"case_id": case_id, "status": "skipped_deadline", "proposal": item}
            write(case_dir / "case_receipt.json", result); return result
        command = evaluator_command(evaluator_path, helper_path, joint_path, template_path,
                                    proposal_path, evaluation)
        environment = dict(os.environ, PYTHONOPTIMIZE="0", NUMBA_DISABLE_JIT="0",
                           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                           NUMBA_NUM_THREADS="1", MPLCONFIGDIR=str(case_dir / "mpl"))
        timed_out = False
        with (case_dir / "evaluator.log").open("w", encoding="utf-8") as log:
            child = subprocess.Popen(command, cwd=packet["source_root"], env=environment,
                                     stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                returncode = child.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                timed_out = True; os.killpg(child.pid, signal.SIGKILL); returncode = child.wait()
        receipt_path = evaluation / "candidate_result.json"
        if timed_out:
            result = {"case_id": case_id, "status": "timeout", "proposal": item,
                      "failure_detail": "candidate exhausted its stage wall-clock budget"}
            write(case_dir / "case_receipt.json", result); return result
        if not receipt_path.exists():
            raise RuntimeError(
                f"source/target preflight or evaluator failed for {case_id} (exit {returncode}); "
                f"no candidate_result.json; see {case_dir / 'evaluator.log'}"
            )
        candidate = read(receipt_path)
        if candidate.get("status") != "verified":
            result = {"case_id": case_id, "status": "failed_candidate", "proposal": item,
                      "returncode": returncode,
                      "failure_detail": candidate.get("failure_detail", "unspecified candidate failure"),
                      "output": str(evaluation)}
            write(case_dir / "case_receipt.json", result); return result
        summary, launch = read(evaluation / "summary.json"), read(evaluation / "launch_contract.json")
        if (summary.get("status") != "verified_rebated_initial_smoke"
                or summary.get("method") != "joint_price_fertility_rebate"
                or launch.get("helper_sha256") != REBATED_HELPER_SHA256
                or launch.get("joint_sha256") != JOINT_ADAPTER_SHA256
                or launch.get("objective_sha256") != objective_hash
                or candidate["loss"] != candidate["score"]["loss"]
                or candidate["loss"] != summary.get("loss")):
            raise RuntimeError(f"verified receipt gates failed for {case_id}")
        early_root = Path(candidate["output"]) / "raw"
        extras = []
        for repetition in range(1, int(item["repetitions"]) + 1):
            early = read(early_root / f"repetition_{repetition:02d}/early_measurement.json")
            extras.append(extra_score(early["fertility"]["uniform_birth_time"], data_path))
        if len(extras) == 2 and extras[0] != extras[1]:
            raise RuntimeError("selected two-repetition age-profile score did not reproduce exactly")
        write(case_dir / "extra_score.json", extras[0])
        if len(extras) == 2: write(case_dir / "extra_score_repetition_02.json", extras[1])
        result = {"case_id": case_id, "status": "verified", "original_loss": float(candidate["loss"]),
                  "extra_loss": float(extras[0]["extra_loss"]),
                  "augmented_loss": float(candidate["loss"]) + float(extras[0]["extra_loss"]),
                  "checkpoint": candidate["accounting"][-1], "repetitions": int(item["repetitions"]),
                  "second_signature_equal": candidate.get("second_signature_equal"), "proposal": item,
                  "output": str(evaluation), "score": candidate["score"], "extra_score": extras[0]}
        write(case_dir / "case_receipt.json", compact(result))
        write(case_dir / "phase.json", {"status": "verified", "phase": "complete",
                                         "checkpoint": result["checkpoint"]})
        return result

    def record(result):
        with lock:
            records.append(result); state["completed"] += 1
            state["failed"] += result["status"] != "verified"; persist()

    try:
        with futures.ThreadPoolExecutor(max_workers=arguments.workers) as pool:
            tasks = [pool.submit(run_case, item, primary_deadline) for item in cases]
            for task in futures.as_completed(tasks): record(task.result())
        primary = [row for row in records if row["case_id"] != "selected_exact_repetitions"]
        valid = [row for row in primary if row["status"] == "verified"]
        baseline = next((row for row in valid if row["case_id"] == "baseline"), None)
        if baseline is None: raise RuntimeError("the joint rebated baseline did not complete successfully")
        selected = min(valid, key=lambda row: row["augmented_loss"])
        improved = selected["case_id"] != "baseline" and selected["augmented_loss"] < baseline["augmented_loss"]
        repeated = None
        state["phase"] = "selected_exact_repetitions" if improved else "no_improved_candidate"; persist()
        if improved and time.monotonic() < deadline:
            repeated_item = copy.deepcopy(selected["proposal"])
            repeated_item.update(case_id="selected_exact_repetitions", repetitions=2)
            repeated = run_case(repeated_item, deadline); record(repeated)
            if (repeated["status"] != "verified" or repeated.get("second_signature_equal") is not True
                    or score_signature(repeated["score"]) != score_signature(selected["score"])
                    or repeated["extra_score"] != selected["extra_score"]):
                raise RuntimeError("improved candidate failed exact two-repetition reproduction")

        selected_record = repeated if repeated is not None else selected
        selected_evaluation = Path(selected_record["output"])
        shutil.copy2(selected_evaluation / "selected_target_fit.csv", output / "selected_original_target_fit.csv")
        shutil.copy2(selected_evaluation / "selected_parameters.csv", output / "selected_candidate_parameters.csv")
        _write_csv(output / "selected_extra_target_fit.csv", selected["extra_score"]["rows"])
        write(output / "selected_extra_age_profiles.json", selected["extra_score"]["age_profiles"])
        combined = augmented_target_rows(selected["score"], selected["extra_score"])
        _write_csv(output / "selected_augmented_target_fit.csv", combined)
        write(output / "selected_checkpoint.json", selected_record["checkpoint"])
        all_original, all_extra, all_parameters = [], [], []
        for row in valid:
            all_original.extend(dict(case_id=row["case_id"], **fit) for fit in row["score"]["target_fit"])
            all_extra.extend(dict(case_id=row["case_id"], **fit) for fit in row["extra_score"]["rows"])
            all_parameters.extend(dict(case_id=row["case_id"], **value) for value in row["score"]["parameters"])
        _write_csv(output / "all_original_target_fits.csv", all_original)
        _write_csv(output / "all_extra_target_fits.csv", all_extra)
        _write_csv(output / "all_candidate_parameters.csv", all_parameters)
        exact_verified = repeated is not None and repeated["status"] == "verified"
        final_status = ("completed_augmented_pilot" if exact_verified else
                        "stopped_before_exact_repetitions" if improved else
                        "completed_no_improved_candidate")
        final = {"status": final_status,
                 "phase": "finished", "elapsed_seconds": time.monotonic() - started,
                 "deadline_seconds": arguments.seconds, "attempted_primary_cases": len(primary),
                 "verified_primary_cases": len(valid), "failed_primary_cases": len(primary) - len(valid),
                 "selected_case": selected["case_id"], "selected_improves_baseline": improved,
                 "selected_augmented_loss": selected["augmented_loss"],
                 "selected_original_loss": selected["original_loss"], "selected_extra_loss": selected["extra_loss"],
                 "selected_checkpoint": selected_record["checkpoint"],
                 "selected_exact_repetitions_verified": exact_verified, "pilot_fingerprint": fingerprint,
                 "original_scored_rows": 12, "experimental_extra_rows": 6, "production_eligible": False}
        write(output / "summary.json", final); print(json.dumps(final), flush=True)
        if improved and not exact_verified:
            raise SystemExit(2)
    finally:
        stop.set()


if __name__ == "__main__":
    main()
