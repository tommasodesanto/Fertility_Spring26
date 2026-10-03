"""Zero-solve parent and continuation contract gate."""
from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

PACKETS = Path("output/model/fixed_reference_economics_20260928")
OLD_DRIVER = Path("code/model/experiments/purchase_timing_sandbox/calibrate.py")
NEW_DRIVER = Path("code/cluster/soft_timing_calibration/continuation/calibrate.py")
STARTS = PACKETS / "soft_timing_continuation_20261003_v1/start_plan.json"
PARENT = Path("/scratch/td2248/projects/soft_timing_calibration_20261002_v3")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def func(path: Path, name: str) -> str:
    nodes = [n for n in ast.walk(ast.parse(path.read_text())) if isinstance(n, ast.FunctionDef) and n.name == name]
    if len(nodes) != 1:
        raise AssertionError(f"Expected one {name}")
    return ast.dump(nodes[0], include_attributes=False)


def verify(stage: Path, parent: Path = PARENT) -> dict:
    inv = json.loads((stage / "inventory.json").read_text())
    p_inv = json.loads((parent / "inventory.json").read_text())
    if sha(parent / "stage.tar.gz") != inv["parent_archive_sha256"]:
        raise AssertionError("Passed parent archive drift")
    for rel, digest in p_inv["files"].items():
        if inv["files"].get(rel) != digest:
            raise AssertionError(f"Passed source changed: {rel}")
    new_inputs = set(inv["new_inputs"])
    if set(inv["files"]) - set(p_inv["files"]) != new_inputs:
        raise AssertionError("Unexpected source in continuation stage")
    for name in ("checked_inputs", "install_timing_observer", "completion_receipt", "selected_repeat_path", "objective"):
        if func(parent / "source" / OLD_DRIVER, name) != func(stage / "source" / NEW_DRIVER, name):
            raise AssertionError(f"Economic/receipt function changed: {name}")
    table_path = stage / "source" / STARTS
    if sha(table_path) != inv["start_plan_sha256"]:
        raise AssertionError("Start plan SHA drift")
    table = json.loads(table_path.read_text())
    parent_starts = json.loads((parent / "source" / PACKETS / "soft_timing_calibration_20261002_v1/expanded_start_plan.json").read_text())
    original = json.loads((parent / "source" / PACKETS / "normalized_calibration_v2/plan.json").read_text())
    if table["target_contract"] != original["base_target_contract"] or table["bounds"] != parent_starts["bounds"]:
        raise AssertionError("Original target or bound changed")
    if table["target_fingerprint"] != inv["target_fingerprint"] or table["weight_fingerprint"] != inv["weight_fingerprint"]:
        raise AssertionError("Objective fingerprint drift")
    if len(table["starts"]) != 10 or len(table["start_provenance"]) != 10:
        raise AssertionError("Ten-start table incomplete")
    for i, (point, provenance) in enumerate(zip(table["starts"], table["start_provenance"])):
        if set(point) != set(table["bounds"]) or any(not lo <= point[key] <= hi for key, (lo, hi) in table["bounds"].items()):
            raise AssertionError(f"Start {i} violates bounds")
        path = stage / "source" / provenance["source_path"]
        if sha(path) != provenance["source_sha256"]:
            raise AssertionError(f"Checkpoint {i} SHA drift")
        receipt = json.loads(path.read_text())
        if receipt["status"] != "selected_numerically_verified" or receipt["selected"]["parameters"] != point:
            raise AssertionError(f"Checkpoint {i} selected point drift")
    if len({tuple(sorted(row.items())) for row in table["starts"]}) != 10:
        raise AssertionError("Duplicate continuation starts")
    return dict(status="continuation_zero_solve_gate_passed", chains=10,
                target_fingerprint=inv["target_fingerprint"], weight_fingerprint=inv["weight_fingerprint"],
                selected_continuation_checkpoint_sha256=inv["selected_continuation_checkpoint_sha256"])


if __name__ == "__main__":
    stage = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    parent = Path(sys.argv[2]) if len(sys.argv) > 2 else PARENT
    print(json.dumps(verify(stage, parent), sort_keys=True))
