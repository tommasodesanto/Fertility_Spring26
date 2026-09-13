#!/usr/bin/env python3
"""Prepare, but never submit, corrected-source finite-history manifests."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import time


EXPECTED_DRIVER_SHA256 = "312d17bbf499388d114b3370e4cde194022913026b997e3bde6abf9338d41ca7"
PYTHON = "/share/apps/anaconda3/2025.06/bin/python"


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load(path: Path) -> dict:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def require_pin(path: Path, expected: str, label: str) -> None:
    if not path.is_file() or sha(path) != expected:
        raise ValueError(f"{label} SHA-256 mismatch: {path}")


def write_new(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def relative_source_pins(receipt: dict) -> tuple[Path, dict[str, str]]:
    source_root = Path(receipt["source_root"]).resolve()
    pins = receipt["source_sha256"]
    if receipt.get("source_files") != len(pins) or not pins:
        raise ValueError("Corrected source receipt has an inconsistent file count")
    relative = {}
    for name, digest in pins.items():
        path = Path(name).resolve()
        try:
            key = str(path.relative_to(source_root))
        except ValueError as exc:
            raise ValueError(f"Corrected source pin escapes source_root: {path}") from exc
        require_pin(path, digest, "corrected source")
        relative[key] = digest
    if len(relative) != len(pins):
        raise ValueError("Corrected source receipt has duplicate relative paths")
    return source_root, dict(sorted(relative.items()))


def validate_kernel(receipt: dict, corrected_root: Path) -> tuple[Path, Path, Path]:
    kernel_path = Path(receipt["kernel_equivalence"]).resolve()
    require_pin(kernel_path, receipt["kernel_equivalence_sha256"], "kernel receipt")
    kernel = load(kernel_path)
    pairs = kernel.get("pairs")
    if kernel.get("status") != "verified" or not isinstance(pairs, list) or len(pairs) != 77:
        raise ValueError("Kernel receipt must be verified and contain exactly 77 pairs")
    initial_roots, history_roots = set(), set()
    marker = Path("code/model/intergen_eqscale_seq_optimized")
    for pair in pairs:
        first, second = Path(pair["initial"]).resolve(), Path(pair["history"]).resolve()
        require_pin(first, pair["sha256"], "initial kernel")
        require_pin(second, pair["sha256"], "history kernel")
        first_parts, second_parts = first.parts, second.parts
        mark = marker.parts
        try:
            i = next(i for i in range(len(first_parts)) if first_parts[i:i + len(mark)] == mark)
            j = next(i for i in range(len(second_parts)) if second_parts[i:i + len(mark)] == mark)
        except StopIteration as exc:
            raise ValueError("Kernel pair does not lie below the model kernel directory") from exc
        initial_roots.add(Path(*first_parts[:i]))
        history_roots.add(Path(*second_parts[:j]))
    if len(initial_roots) != 1 or history_roots != {corrected_root}:
        raise ValueError("Kernel receipt roots do not match the corrected history source")
    return kernel_path, next(iter(initial_roots)), corrected_root


def validate_initial(path: Path, expected: str, initial_root: Path) -> dict:
    require_pin(path, expected, "corrected initial summary")
    summary = load(path)
    if (summary.get("status") != "verified_rebated_initial_smoke"
            or Path(summary.get("source_root", "")).resolve() != initial_root):
        raise ValueError("Initial summary is not a verified rebated result from corrected initial source")
    checkpoint = summary.get("checkpoint")
    if not isinstance(checkpoint, dict):
        raise ValueError("Initial summary checkpoint receipt missing")
    checkpoint_path = Path(checkpoint.get("checkpoint", "")).resolve()
    require_pin(checkpoint_path, checkpoint.get("checkpoint_sha256", ""), "initial checkpoint")
    return summary


def validate_helper_snapshot(source: Path, cached: Path) -> dict[str, str]:
    source_files = {p.name: p for p in source.glob("*.py")}
    cached_files = {p.name: p for p in cached.glob("*.py")}
    if source_files.keys() != cached_files.keys() or not source_files:
        raise ValueError("Auto source must have exactly the cached helper filenames")
    changed = [name for name in sorted(source_files)
               if sha(source_files[name]) != sha(cached_files[name])]
    if changed != ["run_e5f_final_rebated_history.py"]:
        raise ValueError(f"Only the history driver may differ from cached source: {changed}")
    require_pin(source_files[changed[0]], EXPECTED_DRIVER_SHA256, "auto history driver")
    return {str(path.resolve()): sha(path) for path in source_files.values()}


def assert_plan_diff(base: dict, new: dict) -> None:
    normalized = copy.deepcopy(new)
    normalized["source_root"] = base["source_root"]
    normalized["terminal_template"]["source_sha256"] = base["terminal_template"]["source_sha256"]
    if normalized != base:
        raise AssertionError("Corrected history plan changed outside the two source relocation fields")


def assert_manifest_diff(base: dict, new: dict) -> None:
    normalized = copy.deepcopy(new)
    for key in ("prior_plan", "initial_summary", "initial_source_root", "kernel_equivalence", "file_sha256"):
        normalized[key] = copy.deepcopy(base[key])
    normalized["seed_step"] = base["seed_step"]
    normalized["root_controls"].pop("automatic_fiscal_polish")
    normalized.pop("reuse_forecast_jacobian")
    normalized.pop("corrected_source_receipt")
    normalized.pop("kernel_equivalence_pairs")
    if normalized != base:
        raise AssertionError("Corrected history manifest changed an unapproved field")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--initial-summary", type=Path, required=True)
    parser.add_argument("--initial-summary-sha256", required=True)
    parser.add_argument("--source-receipt", type=Path, required=True)
    parser.add_argument("--source-receipt-sha256", required=True)
    parser.add_argument("--label", default="corrected_auto_v1")
    parser.add_argument("--deadline-epoch", type=float, default=1789322400.0)
    args = parser.parse_args(argv)
    if not args.label.replace("_", "").isalnum():
        raise ValueError("Simple label required")

    batch = args.batch.resolve()
    source = batch / "history_source_auto_v1"
    cached = batch / "history_source_cached"
    base_plan_path = batch / "history_plan.json"
    base_manifest_path = batch / "history_manifest_refit.json"
    require_pin(args.source_receipt.resolve(), args.source_receipt_sha256, "corrected source receipt")
    source_receipt = load(args.source_receipt.resolve())
    corrected_root, corrected_pins = relative_source_pins(source_receipt)
    kernel_path, initial_root, _ = validate_kernel(source_receipt, corrected_root)
    validate_initial(args.initial_summary.resolve(), args.initial_summary_sha256, initial_root)
    helper_pins = validate_helper_snapshot(source, cached)
    cache_source = batch / "cache_source_v2"
    cache_runner = cache_source / "run_e5f_cached_history.py"
    cache_module = cache_source / "e5f_exact_policy_cache.py"
    cache_proof_path = batch / "policy_cache_probe_v2" / "summary.json"
    cache_proof = load(cache_proof_path)
    cache_sha = sha(cache_module)
    if (cache_proof.get("status") != "verified"
            or cache_proof.get("cache_sha256") != cache_sha
            or cache_proof.get("exact_mapping_equal") is not True
            or cache_proof.get("household_and_accounting_mapping_valid") is not True):
        raise ValueError("Exact policy cache lacks its matching native forecast certificate")

    base_plan = load(base_plan_path)
    for name, digest in base_plan["file_sha256"].items():
        require_pin(Path(name), digest, "base plan input")
    plan = copy.deepcopy(base_plan)
    plan["source_root"] = str(corrected_root)
    plan["terminal_template"]["source_sha256"] = corrected_pins
    assert_plan_diff(base_plan, plan)
    plan_path = batch / f"history_plan_{args.label}.json"
    write_new(plan_path, plan)

    base_manifest = load(base_manifest_path)
    manifest = copy.deepcopy(base_manifest)
    manifest["prior_plan"] = str(plan_path)
    manifest["initial_summary"] = str(args.initial_summary.resolve())
    manifest["initial_source_root"] = str(initial_root)
    manifest["kernel_equivalence"] = str(kernel_path)
    pins = {**helper_pins,
            **{str(corrected_root / name): digest for name, digest in corrected_pins.items()}}
    pins[str(plan_path)] = sha(plan_path)
    pins[str(kernel_path)] = sha(kernel_path)
    pins[str(args.source_receipt.resolve())] = sha(args.source_receipt.resolve())
    pins[str(args.initial_summary.resolve())] = sha(args.initial_summary.resolve())
    for path in base_plan["empirical_blocks"],:
        pins[path] = sha(Path(path))
    manifest["file_sha256"] = dict(sorted(pins.items()))
    manifest["seed_step"] = -0.02
    manifest["root_controls"]["automatic_fiscal_polish"] = True
    manifest["reuse_forecast_jacobian"] = True
    manifest["corrected_source_receipt"] = {
        "path": str(args.source_receipt.resolve()),
        "sha256": args.source_receipt_sha256,
    }
    manifest["kernel_equivalence_pairs"] = 77
    manifest.pop("resume_history", None)
    assert_manifest_diff(base_manifest, manifest)
    manifest_path = batch / f"history_manifest_{args.label}.json"
    write_new(manifest_path, manifest)

    remaining = math.floor(args.deadline_epoch - time.time())
    stage_seconds = min(43200, remaining)
    if stage_seconds <= 10860:
        raise ValueError("Insufficient absolute-deadline budget for policy reserve and a forecast")
    driver_seconds = stage_seconds - 60
    runroot = batch / f"histories_{args.label}"
    cache_pins = {
        str(cache_runner): sha(cache_runner), str(cache_module): cache_sha,
        str(cache_proof_path): sha(cache_proof_path),
    }
    array_paths = {}
    smoke_commands = []
    deadline_wrapper = (
        "import json,math,pathlib,subprocess,sys,time; deadline=float(sys.argv[1]); "
        "argv=sys.argv[2:]; seconds=min(int(argv[argv.index('--seconds')+1]),"
        "math.floor(deadline-time.time())-60); "
        "out=pathlib.Path(argv[argv.index('--output')+1]); "
        "out.mkdir(parents=True,exist_ok=True); "
        "(out/'queue_deadline_receipt.json').write_text(json.dumps("
        "{'remaining_driver_seconds':seconds,'deadline_unix':deadline})); "
        "argv[argv.index('--seconds')+1]=str(seconds); "
        "sys.exit(subprocess.call(argv) if seconds>10800 else 0)"
    )
    for count, memory in ((6, "16G"), (24, "32G"), (100, "48G")):
        runroot = batch / f"histories_{args.label}_{count}"
        stages = []
        for case, case_label in (("A0", "A0"), ("A+", "Aplus")):
            name = f"{case_label}_{count}"
            driver_argv = [
                PYTHON, "-B", str(cache_runner), "--cache-proof", str(cache_proof_path),
                "--cache-sha256", cache_sha, "--", "--manifest", str(manifest_path),
                "--case", case, "--count", str(count), "--output", str(runroot / name),
                "--seconds", str(driver_seconds),
            ]
            command = [PYTHON, "-c", deadline_wrapper, str(args.deadline_epoch), *driver_argv]
            stages.append({
                "name": name, "seconds": stage_seconds, "cpus": 1, "command": command,
            })
            if count == 6:
                smoke_commands.append(command)
        array = {
            "runroot": str(runroot), "horizon_hours": stage_seconds / 3600,
            "absolute_deadline_epoch": args.deadline_epoch, "max_workers": 2, "memory": memory,
            "source_pins": dict(sorted({**pins, **cache_pins,
                                        str(manifest_path): sha(manifest_path)}.items())),
            "stages": stages,
            "environment": {
                "PYTHONPATH": f"{cache_source}:{source}:{corrected_root / 'code/model/tools'}:{corrected_root / 'code/model'}",
                "NUMBA_CACHE_DIR": str(corrected_root / "output/cache/numba"),
            },
        }
        array_path = batch / f"history_array_manifest_{args.label}_{count}.json"
        write_new(array_path, array)
        array_paths[str(count)] = {"path": str(array_path), "sha256": sha(array_path), "memory": memory}
    receipt = {
        "status": "prepared_not_submitted", "plan": str(plan_path), "plan_sha256": sha(plan_path),
        "manifest": str(manifest_path), "manifest_sha256": sha(manifest_path),
        "array_manifests": array_paths, "absolute_deadline_epoch": args.deadline_epoch,
        "driver_sha256": EXPECTED_DRIVER_SHA256, "kernel_equivalence_pairs": 77,
        "cache_sha256": cache_sha, "cache_proof_sha256": sha(cache_proof_path),
        "economic_diff_assertions": {"plan": "passed", "manifest": "passed"},
        "smoke_commands": smoke_commands,
    }
    receipt_path = batch / f"history_prepare_{args.label}_receipt.json"
    write_new(receipt_path, receipt)
    print(json.dumps({**receipt, "receipt": str(receipt_path)}, indent=2))


if __name__ == "__main__":
    main()
