"""Check the September 14 reference hashes without importing the model."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[3]
    packet = root / "output/model/paper_baseline_sep14"
    manifest = json.loads((packet / "manifest.json").read_text())
    failures = []
    counts = {}
    for section in ("source_files", "artifacts"):
        entries = manifest[section]
        counts[section] = len(entries)
        for relative, expected in entries.items():
            path = root / relative
            if not path.is_file():
                failures.append({"path": relative, "reason": "missing"})
            elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                failures.append({"path": relative, "reason": "hash differs"})
    result = dict(status="failed" if failures else "passed", **counts,
                  model_solves=0, failures=failures)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
