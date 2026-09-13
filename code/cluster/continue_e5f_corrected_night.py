#!/usr/bin/env python3
"""Validate corrected initial results, gate native smokes, and hand off arrays."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path: Path):
    return json.loads(path.read_text())


def save(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def require_pin(path: Path, expected: str, label: str) -> None:
    if not path.is_file() or sha(path) != expected:
        raise ValueError(f"{label} SHA-256 mismatch: {path}")


def checkpoint_identity(item: dict) -> tuple[str, str]:
    path = item.get("checkpoint", item.get("path"))
    digest = item.get("checkpoint_sha256", item.get("sha256"))
    if not path or not digest:
        raise ValueError("Checkpoint path and SHA-256 are required")
    return str(Path(path).resolve()), str(digest)


def validate_initial(initial_dir: Path) -> tuple[Path, str, dict]:
    summary_path = initial_dir / "summary.json"
    candidate_path = initial_dir / "candidate_result.json"
    if not summary_path.exists() or not candidate_path.exists():
        raise FileNotFoundError("Corrected initial summary and candidate_result are not complete")
    summary, candidate = load(summary_path), load(candidate_path)
    if summary.get("status") != "verified_rebated_initial_smoke":
        raise ValueError("Corrected initial summary is not verified")
    if (candidate.get("status") != "verified"
            or candidate.get("second_signature_equal") is not True
            or candidate.get("proposal", {}).get("repetitions") != 2):
        raise ValueError("Corrected initial candidate lacks its two-repetition verification")
    accounting = candidate.get("accounting")
    if not isinstance(accounting, list) or len(accounting) != 2:
        raise ValueError("Corrected initial candidate must have two accounting receipts")
    identities = []
    for item in accounting:
        identity = checkpoint_identity(item)
        require_pin(Path(identity[0]), identity[1], "initial repetition checkpoint")
        identities.append(identity)
    repetition_ids = [item.get("repetition") for item in accounting]
    if len(set(repetition_ids)) != 2 or any(not value for value in repetition_ids):
        raise ValueError("Initial repetition identifiers must be distinct")
    if checkpoint_identity(summary.get("checkpoint", {})) != identities[-1]:
        raise ValueError("Initial summary checkpoint is not the last verified repetition")
    return summary_path.resolve(), sha(summary_path), {
        "candidate_result": str(candidate_path.resolve()),
        "candidate_result_sha256": sha(candidate_path),
        "repetition_ids": repetition_ids,
        "checkpoint_sha256": identities[-1][1],
    }


def verify_source_pins(array_path: Path) -> None:
    array = load(array_path)
    for name, digest in array.get("source_pins", {}).items():
        require_pin(Path(name), digest, "array source")
    if not array.get("source_pins"):
        raise ValueError(f"Array source pins missing: {array_path}")


def run_preparer(args, state_dir: Path, summary_path: Path, summary_sha: str) -> dict:
    prepared_path = args.batch / "history_prepare_corrected_auto_v2_receipt.json"
    if prepared_path.exists():
        receipt = load(prepared_path)
    else:
        partials = list(args.batch.glob("history_*_corrected_auto_v2*.json"))
        if partials:
            raise ValueError("Partial corrected_auto_v2 preparation exists without its receipt")
        source_receipt = args.batch / "corrected_history_source_v2_receipt.json"
        command = [sys.executable, "-B", str(args.preparer), "--batch", str(args.batch),
                   "--initial-summary", str(summary_path), "--initial-summary-sha256", summary_sha,
                   "--source-receipt", str(source_receipt),
                   "--source-receipt-sha256", sha(source_receipt),
                   "--label", "corrected_auto_v2", "--deadline-epoch", str(args.deadline_epoch)]
        save(state_dir / "preparer_command.json", {"command": command})
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                timeout=max(1.0, min(600.0, args.seconds)))
        (state_dir / "preparer_stdout.json").write_text(result.stdout)
        (state_dir / "preparer_stderr.log").write_text(result.stderr)
        if result.returncode:
            raise RuntimeError(f"Preparer failed with exit status {result.returncode}")
        if not prepared_path.exists():
            raise ValueError("Preparer returned without its receipt")
        receipt = load(prepared_path)
    if receipt.get("status") != "prepared_not_submitted":
        raise ValueError("Preparation receipt is not complete")
    for key in ("plan", "manifest"):
        require_pin(Path(receipt[key]), receipt[f"{key}_sha256"], f"prepared {key}")
    manifest = load(Path(receipt["manifest"]))
    for name, digest in manifest.get("file_sha256", {}).items():
        require_pin(Path(name), digest, "scientific manifest input")
    source_receipt = args.batch / "corrected_history_source_v2_receipt.json"
    source_item = manifest.get("corrected_source_receipt", {})
    if (Path(manifest.get("initial_summary", "")).resolve() != summary_path
            or manifest.get("file_sha256", {}).get(str(summary_path)) != summary_sha
            or Path(source_item.get("path", "")).resolve() != source_receipt.resolve()
            or source_item.get("sha256") != sha(source_receipt)
            or manifest.get("seed_step") != -0.02
            or manifest.get("reuse_forecast_jacobian") is not True
            or manifest.get("root_controls", {}).get("automatic_fiscal_polish") is not True
            or "resume_history" in manifest):
        raise ValueError("Prepared manifest does not match the corrected no-resume contract")
    arrays = receipt.get("array_manifests", {})
    if set(arrays) != {"6", "24", "100"}:
        raise ValueError("Preparation receipt must contain separate 6/24/100 arrays")
    for count, item in arrays.items():
        path = Path(item["path"])
        require_pin(path, item["sha256"], f"prepared {count}-date array")
        verify_source_pins(path)
        manifest = load(path)
        if manifest.get("absolute_deadline_epoch") != args.deadline_epoch:
            raise ValueError("Prepared array has the wrong absolute deadline")
    return receipt


def parse_submitter_output(text: str) -> dict:
    value = json.loads(text)
    if not isinstance(value, dict) or not isinstance(value.get("command"), list):
        raise ValueError("Submitter returned an invalid receipt")
    return value


def submit_array(args, state_dir: Path, count: int, array_item: dict) -> dict:
    record_path = state_dir / f"submission_{count}.json"
    if record_path.exists():
        return load(record_path)
    array_path = Path(array_item["path"])
    array = load(array_path)
    native_path = Path(array["runroot"]) / "submission.json"
    if native_path.exists():
        native = load(native_path)
        if native.get("manifest_sha256") != sha(array_path) or not native.get("job"):
            raise ValueError(f"Existing native submission receipt is invalid for count {count}")
        record = {"status": "submitted_existing", "count": count, **native}
        save(record_path, record)
        return record
    command = [sys.executable, "-B", str(args.submitter), str(array_path)]
    if args.submit:
        command.append("--submit")
    save(state_dir / f"submitter_command_{count}.json", {"command": command})
    try:
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                timeout=max(1.0, min(120.0, args.seconds)))
    except subprocess.TimeoutExpired as exc:
        record = {"status": "submission_uncertain" if args.submit else "dry_run_failed",
                  "count": count, "command": command, "error": "submitter timeout",
                  "stdout": str(exc.stdout or "")[-4000:], "stderr": str(exc.stderr or "")[-4000:],
                  "retry_permitted": False if args.submit else None}
        save(record_path if args.submit else state_dir / f"dry_run_{count}.json", record)
        return record
    if not args.submit:
        receipt = parse_submitter_output(result.stdout) if result.returncode == 0 else {}
        record = {"status": "dry_run_ready" if result.returncode == 0 else "dry_run_failed",
                  "count": count, "command": command, "submitter": receipt,
                  "returncode": result.returncode, "stderr": result.stderr[-4000:]}
        save(state_dir / f"dry_run_{count}.json", record)
        return record
    if result.returncode:
        record = {"status": "submission_uncertain", "count": count, "command": command,
                  "returncode": result.returncode, "stdout": result.stdout[-4000:],
                  "stderr": result.stderr[-4000:], "retry_permitted": False}
        save(record_path, record)
        return record
    receipt = parse_submitter_output(result.stdout)
    if not receipt.get("job"):
        raise ValueError("Successful submitter call did not return a job id")
    record = {"status": "submitted", "count": count, "command": command,
              "job": str(receipt["job"]), "submitter": receipt}
    save(record_path, record)
    return record


def all_true(value) -> bool:
    return isinstance(value, dict) and bool(value) and all(item is True for item in value.values())


def smoke_evidence(array_item: dict, case: str) -> dict:
    label = "A0" if case == "A0" else "Aplus"
    root = Path(load(Path(array_item["path"]))["runroot"]) / f"{label}_6"
    failure = root / "failure.json"
    if failure.exists():
        return {"status": "failed", "case": case, "failure": str(failure)}
    required = [root / "native_smoke.json", root / "realized_fit.json"]
    if not all(path.exists() for path in required):
        return {"status": "pending", "case": case, "root": str(root)}
    native, realized = load(required[0]), load(required[1])
    if native.get("passed") is not True or native.get("count") != 6:
        return {"status": "failed", "case": case, "reason": "native smoke did not pass"}
    if not isinstance(realized, list) or not realized:
        return {"status": "failed", "case": case, "reason": "first realized fit is absent"}
    fit = realized[0]
    try:
        gap = float(fit["gap"])
    except (KeyError, TypeError, ValueError):
        return {"status": "failed", "case": case, "reason": "invalid realized gap"}
    if fit.get("year") != 2007 or not math.isfinite(gap) or abs(gap) > 0.005:
        return {"status": "failed", "case": case, "reason": "first-window fertility gate failed",
                "gap": gap}
    folder = fit.get("folder")
    if not isinstance(folder, str) or not Path(folder).is_absolute():
        return {"status": "failed", "case": case, "reason": "realized fit folder is invalid"}
    selected = Path(folder).resolve()
    candidates, incomplete = [], []
    for candidate in (selected, selected / "alternative"):
        root_receipt_path = candidate / "root_receipt.json"
        accepted_path = candidate / "accepted_forecast.pkl.gz"
        if not root_receipt_path.exists():
            incomplete.append(str(candidate))
            continue
        try:
            receipt = load(root_receipt_path)
            final = receipt.get("final") or {}
            payload = final.get("payload") or {}
            dated_gates = [item.get("gates") for item in payload.get("dated_audits", [])]
            receipt_psi = float(receipt.get("psi", math.nan))
            fit_psi = float(fit.get("psi", math.nan))
            passed = (
                receipt.get("converged") is True and receipt.get("status") == "converged"
                and receipt.get("finite_horizon_market_fiscal_converged") is True
                and receipt.get("final_reproduction_max_abs") == 0
                and receipt.get("start_year") == fit.get("year") == 2007
                and receipt.get("case") == case and receipt.get("count") == 6
                and math.isfinite(receipt_psi) and receipt_psi == fit_psi
                and final.get("mapping_valid") is True
                and math.isfinite(float(final.get("score", math.nan)))
                and all_true(payload.get("boundary_household_gates"))
                and len(dated_gates) == 6 and all(all_true(gates) for gates in dated_gates)
            )
        except json.JSONDecodeError:
            incomplete.append(str(candidate))
            continue
        except (AttributeError, TypeError, ValueError, KeyError):
            passed = False
        if passed:
            if accepted_path.exists():
                candidates.append((root_receipt_path, accepted_path))
            else:
                incomplete.append(str(candidate))
    if len(candidates) > 1:
        return {"status": "failed", "case": case,
                "reason": "ambiguous passing parent and alternative roots",
                "root_receipts": [str(item[0]) for item in candidates]}
    if not candidates:
        if incomplete:
            return {"status": "pending", "case": case,
                    "reason": "accepted root artifacts still incomplete", "candidates": incomplete}
        return {"status": "failed", "case": case,
                "reason": "all completed accepted-root candidates failed validation"}
    root_receipt_path, accepted_path = candidates[0]
    return {"status": "passed", "case": case, "gap": gap,
            "root_receipt": str(root_receipt_path), "root_receipt_sha256": sha(root_receipt_path),
            "accepted_forecast": str(accepted_path),
            "native_smoke": str(required[0]), "native_smoke_sha256": sha(required[0])}


def update_progress(state_dir: Path, started: float, phase: str, smoke: dict, submissions: dict) -> None:
    state = {"phase": phase, "elapsed_seconds": time.monotonic() - started,
             "smoke": smoke, "submissions": submissions, "updated_unix": time.time()}
    save(state_dir / "latest.json", state)
    save(state_dir / "heartbeat.json", state)
    passed = {case: item for case, item in smoke.items() if item.get("status") == "passed"}
    save(state_dir / "best_so_far.json", {"passed_smoke_cases": passed,
                                           "submissions": submissions})


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--initial-dir", type=Path, required=True)
    parser.add_argument("--preparer", type=Path, required=True)
    parser.add_argument("--preparer-sha256", required=True)
    parser.add_argument("--submitter", type=Path, required=True)
    parser.add_argument("--submitter-sha256", required=True)
    parser.add_argument("--deadline-epoch", type=float, default=1789322400.0)
    parser.add_argument("--seconds", type=float, default=5400)
    parser.add_argument("--poll-seconds", type=float, default=60.0, help=argparse.SUPPRESS)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)
    args.batch, args.initial_dir = args.batch.resolve(), args.initial_dir.resolve()
    args.preparer, args.submitter = args.preparer.resolve(), args.submitter.resolve()
    if not (0 < args.seconds <= 5400) or not (0 < args.poll_seconds <= 60):
        parser.error("Controller and polling seconds are outside their bounded ranges")
    require_pin(args.preparer, args.preparer_sha256, "preparer")
    require_pin(args.submitter, args.submitter_sha256, "submitter")
    state_dir = args.batch / "corrected_night_handoff_auto_v2"
    state_dir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    end = min(started + args.seconds, started + max(0.0, args.deadline_epoch - time.time()))
    submissions, smoke = {}, {case: {"status": "pending", "case": case} for case in ("A0", "A+")}
    update_progress(state_dir, started, "validating_initial", smoke, submissions)
    try:
        summary_path, summary_sha, initial = validate_initial(args.initial_dir)
        save(state_dir / "initial_receipt.json", initial)
        prepared = run_preparer(args, state_dir, summary_path, summary_sha)
        if not args.submit:
            dry_runs = {count: submit_array(args, state_dir, int(count), prepared["array_manifests"][count])
                        for count in ("6", "24", "100")}
            final = {"status": "dry_run_ready", "submitted": False, "prepared_receipt": prepared,
                     "dry_runs": dry_runs, "initial": initial}
            save(state_dir / "final_receipt.json", final)
            print(json.dumps({"status": final["status"], "receipt": str(state_dir / "final_receipt.json")}))
            return
        if args.deadline_epoch - time.time() <= 10860:
            final = {"status": "bounded_partial_deadline_before_smoke", "submitted": False,
                     "absolute_deadline_epoch": args.deadline_epoch, "submissions": submissions,
                     "smoke": smoke, "unmet_cases": ["A0", "A+"]}
            save(state_dir / "final_receipt.json", final)
            update_progress(state_dir, started, final["status"], smoke, submissions)
            print(json.dumps({"status": final["status"], "receipt": str(state_dir / "final_receipt.json")}))
            return
        submissions["6"] = submit_array(args, state_dir, 6, prepared["array_manifests"]["6"])
        if submissions["6"]["status"] == "submission_uncertain":
            final = {"status": "count6_submission_uncertain", "submitted": True,
                     "submissions": submissions, "smoke": smoke, "retry_permitted": False}
            save(state_dir / "final_receipt.json", final)
            print(json.dumps({"status": final["status"], "receipt": str(state_dir / "final_receipt.json")}))
            return
        while time.monotonic() < end:
            for case in ("A0", "A+"):
                if smoke[case]["status"] == "pending":
                    smoke[case] = smoke_evidence(prepared["array_manifests"]["6"], case)
            update_progress(state_dir, started, "gating_count6", smoke, submissions)
            if all(item["status"] == "passed" for item in smoke.values()):
                break
            if all(item["status"] != "pending" for item in smoke.values()):
                break
            time.sleep(min(args.poll_seconds, max(0.0, end - time.monotonic())))
        if all(item["status"] == "passed" for item in smoke.values()):
            if args.deadline_epoch - time.time() <= 10860:
                status = "bounded_partial_deadline_before_long_arrays"
            else:
                for count in (24, 100):
                    submissions[str(count)] = submit_array(
                        args, state_dir, count, prepared["array_manifests"][str(count)])
                uncertain = any(item["status"] == "submission_uncertain" for item in submissions.values())
                status = "long_submission_uncertain" if uncertain else "long_arrays_submitted"
        else:
            status = "bounded_partial_unmet_smoke_prerequisite"
        final = {"status": status, "submitted": True, "elapsed_seconds": time.monotonic() - started,
                 "absolute_deadline_epoch": args.deadline_epoch, "submissions": submissions,
                 "smoke": smoke, "unmet_cases": [k for k, v in smoke.items() if v["status"] != "passed"]}
        save(state_dir / "final_receipt.json", final)
        update_progress(state_dir, started, status, smoke, submissions)
        print(json.dumps({"status": status, "receipt": str(state_dir / "final_receipt.json")}))
    except Exception as exc:
        failure = {"status": "required_prerequisite_failed", "error": f"{type(exc).__name__}: {exc}",
                   "elapsed_seconds": time.monotonic() - started, "submissions": submissions,
                   "smoke": smoke}
        save(state_dir / "final_receipt.json", failure)
        update_progress(state_dir, started, failure["status"], smoke, submissions)
        raise


if __name__ == "__main__":
    main()
