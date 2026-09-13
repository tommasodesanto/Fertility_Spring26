#!/usr/bin/env python3
"""Prepare and collect corrected-history final-window readouts.

This utility only observes explicitly supplied case directories.  It never
solves, submits, or batches model work; the strict observer remains the sole
reader of native pickle artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import math
import time

WINDOW = [2007, 2011, 2015, 2019]

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def load(path: Path):
    return json.loads(path.read_text())

def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)

def resolve(value, base: Path) -> Path:
    p = Path(value)
    return (base / p).resolve() if not p.is_absolute() else p.resolve()

def explicit_pins(manifest: dict, path: Path) -> dict[Path, str]:
    raw = manifest.get("file_sha256") or manifest.get("source_pins")
    if not isinstance(raw, dict) or not raw:
        raise ValueError("scientific manifest has no explicit file_sha256 pins")
    out = {}
    for name, digest in raw.items():
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError("invalid scientific SHA-256 pin")
        p = resolve(name, path.parent)
        if p in out and out[p] != digest.lower():
            raise ValueError("duplicate scientific pin with conflicting digest")
        out[p] = digest.lower()
    return out

def verify_pins(pins: dict[Path, str]) -> None:
    for p, expected in pins.items():
        if not p.is_file() or sha(p) != expected:
            raise ValueError(f"pinned source changed or missing: {p}")

def completed_history(case: Path):
    rows = load(case / "realized_fit.json")
    complete = load(case / "finite_history_complete.json")
    if not isinstance(rows, list) or [r.get("year") for r in rows] != WINDOW:
        raise ValueError("Exactly four ordered realized windows are required")
    if complete.get("realized") != rows:
        raise ValueError("Completed history receipt differs from realized fit")
    for r in rows:
        gap = float(r["gap"])
        if not math.isfinite(gap) or abs(gap) > .005:
            raise ValueError("A historical fit gap fails the retained tolerance")
        if abs(float(r["model"]) - float(r["target"]) - gap) > 1e-12:
            raise ValueError("Historical fit arithmetic differs")
    folder = resolve(rows[-1]["folder"], case)
    if not folder.is_relative_to(case):
        raise ValueError("Accepted final folder escapes its case")
    return rows, folder

def build_manifest(case: Path, out: Path, scientific: Path, helper: Path,
                   reader: Path, profile: Path) -> tuple[Path, str]:
    sr = case / "contract_receipt.json"
    complete = case / "finite_history_complete.json"
    if not sr.is_file() or not complete.is_file():
        raise FileNotFoundError("case contract and finite_history_complete are required")
    contract = load(sr)
    sm = load(scientific)
    if contract.get("manifest_sha256") != sha(scientific):
        raise ValueError("Case scientific manifest differs")
    source_root = resolve(contract["source_root"], sr.parent)
    pins = explicit_pins(sm, scientific)
    prior_path = resolve(sm["prior_plan"], scientific.parent)
    if prior_path not in pins or sha(prior_path) != pins[prior_path]:
        raise ValueError("Scientific prior plan is not pinned")
    prior = load(prior_path)
    if source_root != resolve(prior["source_root"], prior_path.parent):
        raise ValueError("Case model source differs from the scientific plan")
    for path, digest in explicit_pins(prior, prior_path).items():
        if path in pins and pins[path] != digest:
            raise ValueError("Scientific manifest and prior plan pins conflict")
        pins[path] = digest
    verify_pins(pins)
    rows, fit = completed_history(case)
    row = rows[-1]
    selected = []
    for candidate in [fit, fit / "alternative"]:
        if not (candidate / "accepted_forecast.pkl.gz").is_file():
            continue
        root = load(candidate / "root_receipt.json")
        observations = load(candidate / "fertility.json")
        if (root.get("start_year") == 2019 and root.get("count") == contract.get("count")
                and root.get("case") == contract.get("case") and root.get("converged")
                and root.get("finite_horizon_market_fiscal_converged")
                and (root.get("final") or {}).get("mapping_valid")
                and math.isfinite(float(root.get("final_reproduction_max_abs", math.inf)))
                and float(root.get("final_reproduction_max_abs", math.inf)) <= 2e-10
                and abs(float(root.get("psi", math.inf)) - float(row["psi"])) <= 1e-12
                and observations and observations[0].get("calendar_year") == 2019
                and abs(float(observations[0]["period_tfr_topcode_adjusted"]) - float(row["model"])) <= 2e-10):
            selected.append(candidate)
    if len(selected) != 1:
        raise ValueError("Final fit requires exactly one matching accepted primary/alternative forecast")
    fit = selected[0]
    artifacts = {"accepted_forecast": fit / "accepted_forecast.pkl.gz",
                 "first_period_diagnostics": fit / "first_period_diagnostics.pkl.gz",
                 "native_2023_snapshot": fit / "native_2023_snapshot.pkl.gz",
                 "root_receipt": fit / "root_receipt.json"}
    for p in artifacts.values():
        if not p.is_file(): raise FileNotFoundError(f"native artifact missing: {p}")
        pins[p] = sha(p)
    helper = helper.resolve(); reader = reader.resolve(); profile = profile.resolve()
    # Reader/profile digests have already been checked against explicit CLI pins.
    # Model and helper pins come from the actual frozen scientific manifest.
    for p in [reader, profile]:
        if p in pins and pins[p] != sha(p):
            raise ValueError("Reader/profile contradicts the scientific source pins")
        pins[p] = sha(p)
    if not any(p.is_relative_to(helper) for p in pins):
        raise ValueError("No scientific kernel pins cover the requested helper directory")
    kernels = sorted(str(p) for p in pins if p.suffix == ".py")
    payload = dict(sm)
    payload.update({k: str(v) for k, v in artifacts.items()})
    payload.update(source_root=str(source_root), helper_root=str(helper),
                   profile_kernel=str(profile), kernel_files=kernels,
                   runtime_paths=[str(helper), str(profile.parent),
                                  str(source_root / "code/model/tools"), str(source_root / "code/model")])
    for path in [sr, complete, case / "realized_fit.json", scientific, fit / "fertility.json"]:
        pins[path.resolve()] = sha(path)
    payload["file_sha256"] = {str(p): d for p, d in pins.items()}
    payload["historical_fit_status"] = dict(complete=True, realized=rows,
        source=str(case / "realized_fit.json"),
        verification_method="completed native driver receipt and four fit-gap checks; final root independently checked")
    payload["preparation"] = {"case_dir": str(case), "fit_2019": row,
                              "window": WINDOW, "prepared_unix": time.time()}
    out.mkdir(parents=True, exist_ok=True)
    manifest = out / "readout_manifest.json"
    atomic_json(manifest, payload)
    return manifest, sha(manifest)

def inspect_case(case: Path, args, output: Path) -> dict:
    state = {"case_dir": str(case), "checked_unix": time.time()}
    try:
        if not (case / "finite_history_complete.json").is_file():
            state.update(status="PENDING", reason="Four-window history is not yet complete")
            return state
        work = output / case.name
        prior_path = work / "status.json"
        if prior_path.exists():
            prior = load(prior_path)
            if prior.get("status") in {"COMPLETE", "ERROR", "RUNNING"}:
                # A completed or attempted immutable input is never blindly rerun.
                return prior
        manifest, digest = build_manifest(case, work, args.scientific_manifest,
                                           args.helper_root, args.reader, args.profile)
        verify_pins({Path(k): v for k, v in load(manifest)["file_sha256"].items()})
        readout = work / "readout"
        if readout.exists() or (work / "attempt_started.json").exists():
            state.update(status="ERROR", reason="existing readout; refusing blind rerun")
            return state
        cmd = [sys.executable, "-B", str(args.reader), "--manifest", str(manifest),
               "--manifest-sha256", digest, "--out", str(readout)]
        atomic_json(work / "command_receipt.json", {"command": cmd, "manifest_sha256": digest})
        if args.once:
            state.update(status="READY", manifest=str(manifest), manifest_sha256=digest)
            return state
        atomic_json(work / "attempt_started.json", dict(command=cmd, started_unix=time.time()))
        proc = subprocess.run(cmd, text=True, capture_output=True,
                              timeout=max(1, int(args.command_timeout)))
        verification = load(readout / "verification.json") if (readout / "verification.json").exists() else {}
        passed = (proc.returncode == 0 and verification.get("status") == "PASS"
                  and verification.get("manifest_sha256") == digest)
        state.update(status="COMPLETE" if passed else "ERROR",
                     manifest=str(manifest), manifest_sha256=digest,
                     returncode=proc.returncode, stdout=proc.stdout[-4000:], stderr=proc.stderr[-4000:])
    except FileNotFoundError as exc:
        state.update(status="PENDING", reason=str(exc))
    except OSError as exc:
        state.update(status="WAIT_FILESYSTEM", error_type=type(exc).__name__, reason=str(exc))
    except Exception as exc:
        state.update(status="ERROR", error_type=type(exc).__name__, reason=str(exc))
    return state

def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--case-dir", action="append", required=True)
    p.add_argument("--deadline-epoch", type=float, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--scientific-manifest", type=Path, required=True)
    p.add_argument("--scientific-manifest-sha256", required=True)
    p.add_argument("--helper-root", type=Path, required=True)
    p.add_argument("--reader", type=Path, required=True)
    p.add_argument("--reader-sha256", required=True)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--profile-sha256", required=True)
    p.add_argument("--interval", type=int, default=300)
    p.add_argument("--command-timeout", type=int, default=900)
    p.add_argument("--once", action="store_true")
    a = p.parse_args(argv); cases = [Path(x).resolve() for x in a.case_dir]
    if len(set(c.name for c in cases)) != len(cases):
        raise ValueError("Case output names must be unique within each collector group")
    for path, expected in [(a.scientific_manifest, a.scientific_manifest_sha256),
                           (a.reader, a.reader_sha256), (a.profile, a.profile_sha256)]:
        if sha(path) != expected:
            raise ValueError("Explicit collection source pin differs: " + str(path))
    if a.interval < 60: raise ValueError("interval must be at least 60 seconds")
    while True:
        for path, expected in [(a.scientific_manifest, a.scientific_manifest_sha256),
                               (a.reader, a.reader_sha256), (a.profile, a.profile_sha256)]:
            if sha(path) != expected:
                raise ValueError("Collection source changed while monitoring")
        rows = [inspect_case(c, a, a.output) for c in cases]
        for row in rows:
            case_state = a.output / Path(row["case_dir"]).name
            atomic_json(case_state / "status.json", row)
            atomic_json(case_state / "latest.json", row)
            if row["status"] == "COMPLETE" and not (case_state / "best_so_far.json").exists():
                atomic_json(case_state / "best_so_far.json", row)
        atomic_json(a.output / "latest.json", {"deadline_epoch": a.deadline_epoch,
                    "cases": rows, "updated_unix": time.time()})
        if a.once or all(x["status"] in ("COMPLETE", "ERROR") for x in rows) or time.time() >= a.deadline_epoch:
            break
        time.sleep(min(a.interval, max(1, a.deadline_epoch - time.time())))

if __name__ == "__main__": main()
