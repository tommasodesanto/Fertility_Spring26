"""Hash every source and launcher input in the isolated continuation stage."""
import hashlib
import json
import sys
from pathlib import Path

REMOTE = Path("/scratch/td2248/projects/soft_timing_continuation_20261003_v1")
STAGE = Path("/work/deployment")
REPO = Path("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    mounted = sys.argv[1:] == ["--container"]
    if sys.argv[1:] not in ([], ["--host"], ["--container"]):
        raise SystemExit("Usage: verify_stage.py [--host|--container]")
    stage = STAGE if mounted else REMOTE
    inventory = json.loads((stage / "inventory.json").read_text())
    for rel, digest in inventory["files"].items():
        path = (REPO if mounted else stage / "source") / rel
        if not path.is_file() or sha(path) != digest:
            raise SystemExit(f"Source pin drift: {rel}")
    for rel, digest in inventory["entrypoints"].items():
        if sha(stage / rel) != digest:
            raise SystemExit(f"Entrypoint pin drift: {rel}")
    plan = json.loads(((REPO if mounted else stage / "source") /
                       "output/model/fixed_reference_economics_20260928/soft_timing_continuation_20261003_v1/start_plan.json").read_text())
    if plan["target_fingerprint"] != inventory["target_fingerprint"] or plan["weight_fingerprint"] != inventory["weight_fingerprint"]:
        raise SystemExit("Target/weight fingerprint mismatch")
    if len(plan["starts"]) != 10 or plan["bounds"]["beta_annual"] != [.94, .99]:
        raise SystemExit("Ten-chain or beta-bound contract mismatch")
    print(json.dumps(dict(status="passed", mode="container" if mounted else "host",
                          source_files=len(inventory["files"]), target_fingerprint=inventory["target_fingerprint"],
                          weight_fingerprint=inventory["weight_fingerprint"])))


if __name__ == "__main__":
    main()
