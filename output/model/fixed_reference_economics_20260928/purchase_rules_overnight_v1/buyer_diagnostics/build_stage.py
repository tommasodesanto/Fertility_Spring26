"""Build an immutable, source-only Torch addendum; no solve or submission."""
from __future__ import annotations

import hashlib
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = PACKET.parents[3]
SOURCES = ("summarize_saved_buyers.py", "financial_access.py", "run_selected.py",
           "run_dated.py", "preflight_native.py", "preflight_torch.sh", "stage_torch.sh",
           "launch_readout.sh", "README.md")
REL = HERE.relative_to(ROOT)


def main() -> None:
    files = {str(REL / name): hashlib.sha256((HERE / name).read_bytes()).hexdigest()
             for name in SOURCES}
    inventory = {"schema": "purchase_buyer_diagnostics_source_v1", "files": files,
                 "no_model_solve": True,
                 "separate_from_calibration_and_mechanism_source_pins": True}
    (HERE / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    archive = HERE / "buyer_diagnostics_stage_v4.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for name in SOURCES:
            tar.add(HERE / name, arcname=str(Path("source") / REL / name), recursive=False)
        tar.add(HERE / "inventory.json", arcname="inventory.json", recursive=False)
    print(json.dumps({"archive": str(archive),
                      "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                      "files": files}, indent=2))


if __name__ == "__main__":
    main()
