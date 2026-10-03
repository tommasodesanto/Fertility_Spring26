"""Package the passed soft-timing source snapshot plus two reviewed new inputs."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACKETS = Path("output/model/fixed_reference_economics_20260928")
DEPLOY = ROOT / PACKETS / "alternative_wealth_cluster_20261003_v1/deployment"
PARENT = ROOT / PACKETS / "soft_timing_calibration_20261002_v1/deployment/attempt2/stage.tar.gz"
PARENT_SHA = "f7a8fec4ff370fd3690c0d0068ca595b75a17dd8aaac3bd47f6009ef73ecd68b"
DRIVER = Path("code/model/experiments/alternative_wealth_calibration/cluster_calibrate.py")
STARTS = PACKETS / "alternative_wealth_cluster_20261003_v1/start_plan.json"
WINNER = PACKETS / "alternative_wealth_local_20261003_v1/overnight_two_wave/chain_02/completed.json"
REMOTE = "/scratch/td2248/projects/alternative_wealth_calibration_20261003_v1"
TARGET = "c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70"
WEIGHT = "f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4"


def sha(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def main() -> None:
    if sha(PARENT.read_bytes()) != PARENT_SHA:
        raise SystemExit("Passed parent archive drift")
    with tarfile.open(PARENT) as archive:
        parent_inventory = json.load(archive.extractfile("inventory.json"))
        source = {name.removeprefix("source/"): archive.extractfile(name).read()
                  for name in archive.getnames() if name.startswith("source/")}
    if {key: sha(source[key]) for key in source} != parent_inventory["files"]:
        raise SystemExit("Parent source inventory drift")
    for rel in (DRIVER, STARTS, WINNER):
        if str(rel) in source or not (ROOT / rel).is_file():
            raise SystemExit(f"Unexpected or absent new input: {rel}")
        source[str(rel)] = (ROOT / rel).read_bytes()
    table = json.loads(source[str(STARTS)])
    old = json.loads(source[str(PACKETS / "normalized_calibration_v2/plan.json")])
    expected = [dict(row) for row in old["base_target_contract"]]
    matches = [row for row in expected if row["moment"] == "wealth_earnings"]
    if len(matches) != 1 or matches[0]["target"] != "6.92658379107299" or matches[0]["weight"] != "7.595098472533724":
        raise SystemExit("Parent wealth target or weight drift")
    matches[0]["target"] = "4.45838713455674"
    if table["target_contract"] != expected or table["target_fingerprint"] != TARGET or table["weight_fingerprint"] != WEIGHT:
        raise SystemExit("Experimental target/weight contract drift")
    if table["source_checkpoint_sha256"] != sha(source[str(WINNER)]):
        raise SystemExit("New-wealth winner checkpoint drift")
    winner = json.loads(source[str(WINNER)])
    if winner.get("status") != "selected_numerically_verified" or winner.get("target_fingerprint") != TARGET or winner.get("weight_fingerprint") != WEIGHT:
        raise SystemExit("New-wealth winner contract drift")
    if table["starts"][0] != winner["selected"]["parameters"]:
        raise SystemExit("First start differs from verified new-wealth winner")
    if table["bounds"]["beta_annual"] != [0.93, 0.99] or len(table["starts"]) != 10:
        raise SystemExit("Search bounds or chain count drift")
    if len({tuple(sorted(row.items())) for row in table["starts"]}) != 10:
        raise SystemExit("Duplicate starts")
    for row in table["starts"]:
        if set(row) != set(table["bounds"]) or any(not lo <= row[key] <= hi for key, (lo, hi) in table["bounds"].items()):
            raise SystemExit("Start violates bounds")
    if not any(row["beta_annual"] < .94 for row in table["starts"]):
        raise SystemExit("No new low-beta start")
    entrypoints = {p.name: sha(p.read_bytes()) for p in HERE.iterdir()
                   if p.is_file() and p.suffix in (".py", ".sh") and p.name != "build_stage.py"}
    inventory = dict(files={key: sha(blob) for key, blob in sorted(source.items())},
                     entrypoints=entrypoints, parent_archive_sha256=PARENT_SHA,
                     parent_inventory_sha256=sha(json.dumps(parent_inventory, indent=2, sort_keys=True).encode()),
                     new_inputs=[str(DRIVER), str(STARTS), str(WINNER)], target_fingerprint=TARGET,
                     weight_fingerprint=WEIGHT, selected_source_sha256=sha(source[str(WINNER)]),
                     original_soft_selected_source_sha256=parent_inventory["selected_source_sha256"],
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
