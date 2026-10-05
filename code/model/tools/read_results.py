#!/usr/bin/env python3
"""Constant-scope reader for explicitly indexed saved model results.

Reads the index and requested small JSON/CSV only. It never discovers results,
opens solution arrays, imports the model, or changes a selected slot implicitly.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[3]
INDEX = ROOT / "output/model/production/readout_index.json"
REQUIRED_IDENTITY = ("classification", "reference", "contract", "source", "status")


def read_index() -> dict:
    with INDEX.open(encoding="utf-8") as handle:
        index = json.load(handle)
    if index.get("schema") != "saved_result_readout_v1" or not isinstance(index.get("slots"), dict):
        raise ValueError(f"Invalid readout index: {INDEX}")
    return index


def path_for(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def status(slot: dict, kind: str) -> dict:
    value = slot.get("artifacts", {}).get(kind)
    if not value:
        return {"ready": False, "reason": "not indexed"}
    paths = value if isinstance(value, list) else [value]
    missing = [str(path_for(p)) for p in paths if not path_for(p).is_file()]
    return {"ready": not missing, "paths": [str(path_for(p)) for p in paths], "missing": missing}


def describe(name: str, slot: dict) -> dict:
    result = {"slot": name, **{k: slot.get(k) for k in REQUIRED_IDENTITY}}
    result["epoch"] = slot.get("epoch")
    result["summary"] = slot.get("summary", {})
    result["source_artifacts"] = slot.get("source_artifacts", {})
    source_paths = [slot.get("source")] + list(result["source_artifacts"].values())
    result["missing_source_artifacts"] = [str(path_for(p)) for p in source_paths if p and not path_for(p).exists()]
    result["staleness"] = "source identity must be reviewed before explicit reselection; no automatic freshness claim"
    result["readiness"] = {k: status(slot, k) for k in ("fit", "parameters", "plots")}
    result["warnings"] = slot.get("warnings", [])
    return result


def table(slot: dict, kind: str) -> list[dict]:
    detail = status(slot, kind)
    if not detail["ready"]:
        raise ValueError(f"{kind} unavailable: {detail}")
    paths = detail["paths"]
    if len(paths) != 1:
        raise ValueError(f"{kind} requires exactly one CSV")
    with Path(paths[0]).open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Empty {kind} table: {paths[0]}")
    return rows


def supplemental_table(slot: dict) -> list[dict]:
    value = slot.get("artifacts", {}).get("supplemental_parameters")
    if not value:
        return []
    path = path_for(value)
    if not path.is_file():
        raise ValueError(f"Supplemental parameter table unavailable: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def select_slot(index: dict, slot_name: str, manifest_path: Path) -> dict:
    """Explicitly replace one pointer from a reviewed, small JSON manifest."""
    if slot_name not in index["slots"]:
        raise ValueError(f"Unknown slot: {slot_name}")
    with manifest_path.open(encoding="utf-8") as handle:
        replacement = json.load(handle)
    if replacement.get("slot") != slot_name:
        raise ValueError("Manifest slot does not match requested slot")
    if any(not replacement.get(key) for key in REQUIRED_IDENTITY):
        raise ValueError(f"Manifest requires {', '.join(REQUIRED_IDENTITY)}")
    if not replacement.get("epoch") or not isinstance(replacement.get("artifacts"), dict):
        raise ValueError("Manifest requires epoch and artifacts")
    for kind in ("fit", "parameters", "plots"):
        detail = status(replacement, kind)
        if detail.get("paths") and not detail["ready"]:
            raise ValueError(f"Missing {kind} artifact: {detail['missing']}")
    before = index["slots"][slot_name]
    history = index.setdefault("selection_history", [])
    history.append({"slot": slot_name, "replaced_at_utc": datetime.now(timezone.utc).isoformat(), "previous": before})
    replacement.pop("slot")
    index["slots"][slot_name] = replacement
    fd, scratch = tempfile.mkstemp(prefix=".readout_index_", suffix=".json", dir=INDEX.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(index, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(scratch, INDEX)
    finally:
        if os.path.exists(scratch):
            os.unlink(scratch)
    return describe(slot_name, replacement)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    for command in ("show", "fit", "parameters", "plots"):
        sub.add_parser(command).add_argument("slot")
    update = sub.add_parser("select", help="explicitly replace a slot with a reviewed manifest")
    update.add_argument("slot")
    update.add_argument("manifest", type=Path)
    args = parser.parse_args()
    try:
        index = read_index()
        if args.command == "list":
            output = {"index": str(INDEX), "slots": [describe(k, v) for k, v in index["slots"].items()]}
        elif args.command == "select":
            output = select_slot(index, args.slot, args.manifest)
        else:
            if args.slot not in index["slots"]:
                raise ValueError(f"Unknown slot: {args.slot}")
            slot = index["slots"][args.slot]
            output = describe(args.slot, slot)
            if args.command in ("fit", "parameters"):
                output["rows"] = table(slot, args.command)
                if args.command == "parameters":
                    output["supplemental_rows"] = supplemental_table(slot)
            elif args.command == "plots":
                if not output["readiness"]["plots"]["ready"]:
                    raise ValueError(f"Plots unavailable: {output['readiness']['plots']}")
        print(json.dumps(output, separators=(",", ":"), ensure_ascii=False))
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(json.dumps({"error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
