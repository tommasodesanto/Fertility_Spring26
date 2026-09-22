"""Build a self-contained, hash-pinned earnings-wealth run bundle.

This is packaging only.  It never imports the model, runs a preflight, stages
to Torch, or submits a scheduler job.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

EXPECTED_SOURCE_FILES = 641
CODE_KEYS = ("adapter", "accounting", "income", "period_income", "wrapper",
             "search_controller", "controller")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canonical(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def read(path: Path) -> Any:
    return json.loads(path.read_text())


def write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def require_hash(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise ValueError(f"missing {label}: {path}")
    actual = digest(path)
    if actual != expected:
        raise ValueError(f"{label} hash mismatch: expected {expected}, got {actual}")


def copy_pinned(src: Path, dst: Path, expected: str, label: str, manifest: list[dict[str, Any]]) -> None:
    require_hash(src, expected, label)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    actual = digest(dst)
    if actual != expected:
        raise ValueError(f"copied {label} hash mismatch")
    manifest.append({"label": label, "path": str(dst), "sha256": actual, "source": str(src)})


def _basename_path(entry: dict[str, Any], root: Path, label: str) -> Path:
    path = Path(entry["path"])
    if not path.is_absolute():
        path = (root / path).resolve()
    return path


def prepare(plan_path: Path, bundle: Path, execution_root: Path, python_path: str) -> dict[str, Any]:
    plan_path = plan_path.resolve()
    plan = read(plan_path)
    if plan.get("status") != "validated_by_lead":
        raise ValueError("plan status must be validated_by_lead")
    if plan.get("config_fixed") is not True and plan.get("configuration_status") != "fixed":
        raise ValueError("plan must declare config_fixed=true or configuration_status=fixed")
    if bundle.exists():
        raise ValueError(f"bundle must be a new directory: {bundle}")
    files = plan.get("files") or {}
    missing = [key for key in CODE_KEYS if key not in files]
    if missing:
        raise ValueError(f"plan missing pinned runtime files: {', '.join(missing)}")
    bundle = bundle.resolve()
    execution_root = execution_root.resolve()
    bundle.mkdir(parents=True)
    manifest: list[dict[str, Any]] = []
    root = plan_path.parents[5] if len(plan_path.parents) > 5 else Path.cwd()
    tools = bundle / "tools"
    inputs = bundle / "inputs"

    rewritten = json.loads(json.dumps(plan))
    rewritten["python"] = python_path
    rewritten["bundle_execution_root"] = str(execution_root)
    rewritten["source_root"] = str(execution_root / "source")
    rewritten["bundle_status"] = "prepared_hash_verified"
    rewritten["bundle_parent_plan_sha256"] = digest(plan_path)

    for key in (*CODE_KEYS, *(["local_supervisor"] if "local_supervisor" in files else [])):
        src = _basename_path(files[key], root, key)
        dst = tools / src.name
        copy_pinned(src, dst, files[key]["sha256"], key, manifest)
        rewritten["files"][key] = {**files[key], "path": str(execution_root / "tools" / src.name)}

    initial_src = _basename_path(files["initial_contract"], root, "initial_contract")
    require_hash(initial_src, files["initial_contract"]["sha256"], "initial_contract")
    initial_dst = inputs / "initial_contract.json"
    initial = read(initial_src)
    normalized = Path(initial["normalized_checkpoint"])
    normalized_hash = initial["normalized_checkpoint_sha256"]
    copy_pinned(normalized, inputs / "normalized_checkpoint" / normalized.name, normalized_hash, "normalized_checkpoint", manifest)
    initial["normalized_checkpoint"] = str(execution_root / "inputs" / "normalized_checkpoint" / normalized.name)
    initial_dst.write_text(json.dumps(initial, indent=2, sort_keys=True) + "\n")
    if digest(initial_dst) == files["initial_contract"]["sha256"]:
        raise ValueError("initial contract unexpectedly unchanged after checkpoint relocation")
    manifest.append({"label": "initial_contract_relocated", "path": str(initial_dst), "sha256": digest(initial_dst), "source": str(initial_src), "original_sha256": files["initial_contract"]["sha256"]})
    rewritten["files"]["initial_contract"] = {**files["initial_contract"], "path": str(execution_root / "inputs" / initial_dst.name), "sha256": digest(initial_dst)}

    for key in ("period_estimates", "period_input_receipt"):
        src = _basename_path(files[key], root, key)
        dst = inputs / key / src.name
        copy_pinned(src, dst, files[key]["sha256"], key, manifest)
        rewritten["files"][key] = {**files[key], "path": str(execution_root / "inputs" / key / src.name)}

    run_src = _basename_path(files["run_contract"], root, "run_contract")
    require_hash(run_src, files["run_contract"]["sha256"], "run_contract")
    run = read(run_src)
    objective_src = _basename_path(run["working_objective"], root, "working_objective")
    require_hash(objective_src, run["working_objective"]["sha256"], "working_objective")
    objective = read(objective_src)
    if canonical(objective) != plan["objective_canonical_sha256"]:
        raise ValueError("canonical objective mismatch before packaging")
    run["source_root"] = str(execution_root / "source")
    run["initial_solve_contract"] = {"path": str(execution_root / "inputs" / "initial_contract.json"), "sha256": digest(initial_dst)}
    run["wrapper_sha256"] = files["wrapper"]["sha256"]
    obj = run.get("objective_source_files", {})
    for label, item in obj.items():
        src = _basename_path(item, root, label)
        dst = inputs / "objective_source_files" / label / src.name
        kind = item.get("hash_kind", "bytes")
        expected = objective["source_fingerprints"][label]
        if (digest(src) if kind == "bytes" else canonical(read(src))) != expected:
            raise ValueError(f"objective source pin mismatch before packaging: {label}")
        copy_pinned(src, dst, digest(src), f"objective_source:{label}", manifest)
        if (digest(dst) if kind == "bytes" else canonical(read(dst))) != expected:
            raise ValueError(f"objective source fingerprint changed: {label}")
        item["path"] = str(execution_root / "inputs" / "objective_source_files" / label / src.name)
    for label in ("scorer", "validator", "working_objective"):
        item = run[label]
        src = _basename_path(item, root, label)
        dst = inputs / label / src.name
        copy_pinned(src, dst, item["sha256"], label, manifest)
        item["path"] = str(execution_root / "inputs" / label / src.name)
    run_dst = inputs / "run_contract.json"
    write(run_dst, run)
    manifest.append({"label": "run_contract_relocated", "path": str(run_dst), "sha256": digest(run_dst), "source": str(run_src), "original_sha256": files["run_contract"]["sha256"]})
    rewritten["files"]["run_contract"] = {**files["run_contract"], "path": str(execution_root / "inputs" / run_dst.name), "sha256": digest(run_dst)}

    source_root = Path(plan["source_root"])
    listed = read(initial_src).get("source_sha256") or {}
    if isinstance(listed, dict):
        names = sorted(listed)
    else:
        names = sorted(listed)
    if len(names) != EXPECTED_SOURCE_FILES:
        raise ValueError(f"source manifest must contain exactly {EXPECTED_SOURCE_FILES} files, got {len(names)}")
    for name in names:
        src = source_root / name
        expected = listed[name] if isinstance(listed, dict) else None
        if expected:
            copy_pinned(src, bundle / "source" / name, expected, f"source:{name}", manifest)
        else:
            if not src.is_file():
                raise ValueError(f"missing source manifest file: {name}")
            copy_pinned(src, bundle / "source" / name, digest(src), f"source:{name}", manifest)
    rewritten["source_root"] = str(execution_root / "source")
    rewritten["adapter_path"] = rewritten["files"]["adapter"]["path"]
    rewritten["adapter_sha256"] = rewritten["files"]["adapter"]["sha256"]
    rewritten["source_manifest"] = {"file_count": len(names), "files": listed}
    write(bundle / "plan.json", rewritten)
    write(bundle / "hash_manifest.json", {"schema": "e5f_earnings_wealth_bundle_v1", "plan_sha256": digest(bundle / "plan.json"), "files": manifest})
    return {"status": "prepared_hash_verified", "bundle": str(bundle), "plan": str(bundle / "plan.json"), "file_count": len(manifest), "source_file_count": len(names)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--execution-root", type=Path, required=True)
    parser.add_argument("--python", dest="python_path", required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.plan, args.bundle, args.execution_root, args.python_path), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
