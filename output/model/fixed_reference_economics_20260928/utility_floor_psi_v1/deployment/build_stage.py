"""Build a deterministic, compact overlay; never solve or regenerate source pins."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = PACKET.parents[3]
BASE_REMOTE = "/scratch/td2248/projects/grid_resolution_credit053_v2"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extra-source", action="append", default=[],
                        help="Additional compact dependency, relative to repository root")
    args = parser.parse_args()
    for name in ("run_psi.py", "plan.json", "source_pins.json"):
        if not (PACKET / name).is_file():
            raise SystemExit(f"Missing final prepared source: {name}")
    # Pins are authored exactly once by the preparation owner, after driver finalization.
    pins = json.loads((PACKET / "source_pins.json").read_text())
    for rel, expected in pins.items():
        path = ROOT / rel
        if digest(path) != expected:
            raise SystemExit(f"Source pin drift: {rel}")
    sources = set(PACKET.glob("*.py")) | set(HERE.glob("*.py")) | set(HERE.glob("*.sh"))
    sources |= {PACKET / "plan.json", PACKET / "source_pins.json"}
    for rel in pins:
        path = ROOT / rel
        if path.suffix in {".py", ".sh", ".json", ".csv"} and path.stat().st_size <= 2_000_000:
            sources.add(path)
    for rel in args.extra_source:
        path = (ROOT / rel).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise SystemExit(f"Invalid extra source: {rel}")
        sources.add(path)
    # Large caches/checkpoints/results cannot accidentally enter this source-only overlay.
    for path in sources:
        if path.suffix not in {".py", ".sh", ".json", ".csv"} or path.stat().st_size > 2_000_000:
            raise SystemExit(f"Not a compact driver source: {path}")
    files = {str(path.relative_to(ROOT)): digest(path) for path in sorted(sources)}
    inventory = dict(files=files, includes_checkpoints=False, includes_caches=False,
                     reused_remote_stage=BASE_REMOTE, source_prefix="source/")
    inventory_path = HERE / "inventory.json"
    inventory_path.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    archive = HERE / "utility_floor_psi_v2_stage.tar.gz"
    with archive.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w") as tar:
                entries = [("source/" + rel, ROOT / rel) for rel in sorted(files)]
                entries += [("inventory.json", inventory_path), ("launch_torch.sh", HERE / "launch_torch.sh")]
                for name, path in entries:
                    blob = path.read_bytes()
                    info = tarfile.TarInfo(name)
                    info.size, info.mtime = len(blob), 0
                    info.mode = 0o755 if path.suffix == ".sh" else 0o644
                    tar.addfile(info, io.BytesIO(blob))
    receipt = dict(archive=str(archive), sha256=digest(archive), bytes=archive.stat().st_size,
                   source_files=len(files), reused_remote_stage=BASE_REMOTE)
    (HERE / "stage_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
