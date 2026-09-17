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
    packet_relative = "output/model/paper_baseline_sep14"
    if args.revision == "HEAD":
        # Keep the committed check independent of local edits to the packet.
        tree_names = set(subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", "HEAD", "--", packet_relative],
            cwd=root,
            text=True,
        ).splitlines())
        manifest_name = next(
            (name for name in ("main_expected_source_manifest.json", "manifest.json")
             if f"{packet_relative}/{name}" in tree_names),
            None,
        )
        if manifest_name is None:
            raise FileNotFoundError("no paper-baseline manifest is committed in HEAD")
        manifest = json.loads(subprocess.check_output(
            ["git", "show", f"HEAD:{packet_relative}/{manifest_name}"], cwd=root, text=True
        ))
    else:
        manifest_name = ("main_expected_source_manifest.json"
                         if (packet / "main_expected_source_manifest.json").exists()
                         else "manifest.json")
        manifest = json.loads((packet / manifest_name).read_text())
    failures = []
    counts = {}
    head_contents = {}
    if args.revision == "HEAD":
        tree = subprocess.check_output(["git", "ls-tree", "-r", "--full-tree", "HEAD"], cwd=root, text=True)
        head_blobs = {line.split("\t", 1)[1]: line.split()[2] for line in tree.splitlines() if "\t" in line}
        requested = [head_blobs[path] for section in ("source_files", "artifacts")
                     for path in manifest[section] if path in head_blobs]
        batch = subprocess.run(
            ["git", "cat-file", "--batch"],
            cwd=root,
            input="".join(f"{blob}\n" for blob in requested).encode(),
            stdout=subprocess.PIPE,
            check=True,
        ).stdout
        offset = 0
        for relative, blob in ((path, head_blobs[path]) for section in ("source_files", "artifacts")
                               for path in manifest[section] if path in head_blobs):
            header_end = batch.index(b"\n", offset)
            _object, kind, size = batch[offset:header_end].decode().split()
            if kind != "blob":
                raise RuntimeError(f"HEAD object for {relative} is not a blob")
            start = header_end + 1
            end = start + int(size)
            head_contents[relative] = batch[start:end]
            offset = end + 1
    for section in ("source_files", "artifacts"):
        entries = manifest[section]
        counts[section] = len(entries)
        for relative, expected in entries.items():
            path = root / relative
            if args.revision == "HEAD":
                content = head_contents.get(relative)
                actual = hashlib.sha256(content).hexdigest() if content is not None else None
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
