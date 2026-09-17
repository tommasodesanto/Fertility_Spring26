"""Check the September 14 reference hashes without importing the model."""
from __future__ import annotations

import hashlib
import json
import argparse
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision", choices=("working-tree", "HEAD"), default="working-tree")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    packet = root / "output/model/paper_baseline_sep14"
    manifest_name = "main_expected_source_manifest.json" if (packet / "main_expected_source_manifest.json").exists() else "manifest.json"
    manifest = json.loads((packet / manifest_name).read_text())
    failures = []
    counts = {}
    head_hashes = {}
    if args.revision == "HEAD":
        tree = subprocess.check_output(["git", "ls-tree", "-r", "--full-tree", "HEAD"], cwd=root, text=True)
        head_hashes = {line.split("\t", 1)[1]: line.split()[2] for line in tree.splitlines() if "\t" in line}
    for section in ("source_files", "artifacts"):
        entries = manifest[section]
        counts[section] = len(entries)
        for relative, expected in entries.items():
            path = root / relative
            if args.revision == "HEAD":
                blob = head_hashes.get(relative)
                actual = hashlib.sha256(subprocess.check_output(["git", "cat-file", "blob", blob], cwd=root)).hexdigest() if blob else None
            else:
                actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
            if actual is None:
                failures.append({"path": relative, "reason": "missing"})
            elif actual != expected:
                failures.append({"path": relative, "reason": "hash differs"})
    result = dict(status="failed" if failures else "passed", **counts,
                  model_solves=0, failures=failures)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
