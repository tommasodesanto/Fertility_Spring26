"""Package the passed revised-timing source plus pinned continuation inputs."""
from __future__ import annotations

import ast
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PACKETS = Path("output/model/fixed_reference_economics_20260928")
DEPLOY = ROOT / PACKETS / "soft_timing_continuation_20261003_v1/deployment"
PARENT = ROOT / PACKETS / "soft_timing_calibration_20261002_v1/deployment/expanded/stage.tar.gz"
PARENT_SHA = "604096fe0bc9d8d7f7f3de7b52fbf3a7cad411dd4478f4405e4add68f54371a0"
DRIVER = Path("code/cluster/soft_timing_calibration/continuation/calibrate.py")
STARTS = PACKETS / "soft_timing_continuation_20261003_v1/start_plan.json"
REMOTE = "/scratch/td2248/projects/soft_timing_continuation_20261003_v1"
TARGET = "db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1"
WEIGHT = "2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0"


def sha(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def function_ast(blob: bytes, name: str) -> str:
    nodes = [n for n in ast.walk(ast.parse(blob)) if isinstance(n, ast.FunctionDef) and n.name == name]
    if len(nodes) != 1:
        raise SystemExit(f"Expected one {name}")
    return ast.dump(nodes[0], include_attributes=False)


def main() -> None:
    if sha(PARENT.read_bytes()) != PARENT_SHA:
        raise SystemExit("Passed v3 archive drift")
    with tarfile.open(PARENT) as archive:
        parent_inventory = json.load(archive.extractfile("inventory.json"))
        source = {name.removeprefix("source/"): archive.extractfile(name).read()
                  for name in archive.getnames() if name.startswith("source/")}
    if {key: sha(source[key]) for key in source} != parent_inventory["files"]:
        raise SystemExit("Passed v3 source inventory drift")
    original_driver = source["code/model/experiments/purchase_timing_sandbox/calibrate.py"]
    new_driver = (ROOT / DRIVER).read_bytes()
    for name in ("checked_inputs", "install_timing_observer", "completion_receipt", "selected_repeat_path", "objective"):
        if function_ast(original_driver, name) != function_ast(new_driver, name):
            raise SystemExit(f"Economic/receipt function changed: {name}")
    table = json.loads((ROOT / STARTS).read_text())
    plan = json.loads(source[str(PACKETS / "normalized_calibration_v2/plan.json")])
    if table["target_contract"] != plan["base_target_contract"]:
        raise SystemExit("Original target contract changed")
    if (not (table["target_fingerprint"] == parent_inventory["target_fingerprint"] == TARGET) or
            not (table["weight_fingerprint"] == parent_inventory["weight_fingerprint"] == WEIGHT)):
        raise SystemExit("Target or weight fingerprint changed")
    if table["source_checkpoint_sha256"] != parent_inventory["selected_source_sha256"]:
        raise SystemExit("Original soft source anchor changed")
    old_bounds = json.loads(source[str(PACKETS / "soft_timing_calibration_20261002_v1/expanded_start_plan.json")])["bounds"]
    if table["bounds"] != old_bounds or len(table["starts"]) != 10 or len(set(tuple(sorted(r.items())) for r in table["starts"])) != 10:
        raise SystemExit("Ten-start original-bound contract changed")
    added = [DRIVER, STARTS]
    for i, provenance in enumerate(table["start_provenance"]):
        path = Path(provenance["source_path"])
        blob = (ROOT / path).read_bytes()
        receipt = json.loads(blob)
        if (sha(blob) != provenance["source_sha256"] or
                receipt["status"] != "selected_numerically_verified" or receipt["arm"] != "alternative" or
                receipt["target_fingerprint"] != TARGET or receipt["weight_fingerprint"] != WEIGHT or
                receipt["selected"]["parameters"] != table["starts"][i] or
                receipt["native_loss"] != provenance["native_loss"]):
            raise SystemExit(f"Verified start {i} source drift")
        added.append(path)
    if table["selected_continuation_checkpoint_sha256"] != table["start_provenance"][0]["source_sha256"]:
        raise SystemExit("Best continuation source drift")
    for path in added:
        if str(path) in source:
            raise SystemExit(f"Unexpected source replacement: {path}")
        source[str(path)] = (ROOT / path).read_bytes()
    entrypoints = {p.name: sha(p.read_bytes()) for p in HERE.iterdir()
                   if p.is_file() and p.suffix in (".py", ".sh") and p.name != "build_stage.py"}
    inventory = dict(files={key: sha(blob) for key, blob in sorted(source.items())},
                     entrypoints=entrypoints, parent_archive_sha256=PARENT_SHA,
                     parent_inventory_sha256=sha(json.dumps(parent_inventory, sort_keys=True).encode()),
                     new_inputs=[str(p) for p in added], target_fingerprint=TARGET,
                     weight_fingerprint=WEIGHT, selected_source_sha256=parent_inventory["selected_source_sha256"],
                     selected_continuation_checkpoint_sha256=table["selected_continuation_checkpoint_sha256"],
                     start_plan_sha256=sha(source[str(STARTS)]), remote_root=REMOTE,
                     source_prefix="source/", no_cache_or_results=True)
    DEPLOY.mkdir(parents=True, exist_ok=True)
    (DEPLOY / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    entries = {"source/" + rel: blob for rel, blob in source.items()}
    entries["inventory.json"] = (DEPLOY / "inventory.json").read_bytes()
    entries.update({name: (HERE / name).read_bytes() for name in entrypoints})
    archive_path = DEPLOY / "stage.tar.gz"
    with archive_path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w") as archive:
            for name, blob in sorted(entries.items()):
                info = tarfile.TarInfo(name)
                info.size, info.mtime = len(blob), 0
                info.mode = 0o755 if name.endswith(".sh") else 0o644
                archive.addfile(info, io.BytesIO(blob))
    receipt = dict(archive=str(archive_path), sha256=sha(archive_path.read_bytes()),
                   bytes=archive_path.stat().st_size, source_files=len(source), new_inputs=inventory["new_inputs"],
                   target_fingerprint=TARGET, weight_fingerprint=WEIGHT, start_plan_sha256=inventory["start_plan_sha256"])
    (DEPLOY / "stage_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
