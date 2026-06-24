#!/usr/bin/env python3
"""One-click runner for the current June 2026 intergen model.

Open this file and run it. It runs the active one-market intergenerational
housing-fertility model from the saved current diagnostic theta, then writes a
mechanics packet with plots and tables under `output/model/`.

This is the file to use when someone asks: "run the model and show me the
results." It is diagnostic, not a production calibration search.
"""

from __future__ import annotations

import datetime as dt
import csv
import json
import subprocess
import sys
import time
from pathlib import Path


# =============================================================================
# Current run settings
# =============================================================================

# Saved June 23 current/balanced diagnostic theta.
SOURCE_RECORD = "output/model/intergen_room_distribution_current_best_20260623/summary.json"

# Output folder. Change this if you want to preserve several named runs.
OUTPUT_FOLDER = "output/model/intergen_model_run_current"

# Active target set used for target/model/gap tables.
TARGET_SET = "candidate_replacement_young_old_roomgap_v1"

# Main model grid. The saved source record already contains the current owner
# ladder H_own=[2,4,6,8,10], hR_max=8, J=16, Nb=60, and Nz=5. These are listed
# here so the run is legible when you open this file.
J = 16
NB = 60
INCOME_STATES = 5
OWNER_LADDER = "2,4,6,8,10"
RENTER_MAX_ROOMS = 8.0
MAX_ITER_EQ = 10

# Turn this on only when you want the diagnostic policy proof-of-concept cases:
# baseline, parent LTV relief, property tax increase, and estate tax wedge.
RUN_POLICY_CASES = False


def main() -> None:
    run_start = time.perf_counter()
    start_stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    repo_root = Path(__file__).resolve().parents[2]
    model_dir = Path(__file__).resolve().parent
    venv_python = model_dir / ".venv/bin/python"
    packet_python = venv_python if venv_python.exists() else Path(sys.executable)

    packet_tool = model_dir / "tools/build_intergen_mechanics_packet.py"
    cmd = [
        str(packet_python),
        str(packet_tool),
        "--source",
        str(repo_root / SOURCE_RECORD),
        "--outdir",
        str(repo_root / OUTPUT_FOLDER),
        "--target-set",
        TARGET_SET,
        "--J",
        str(J),
        "--Nb",
        str(NB),
        "--income-states",
        str(INCOME_STATES),
        "--H-own",
        OWNER_LADDER,
        "--hR-max",
        str(RENTER_MAX_ROOMS),
        "--max-iter-eq",
        str(MAX_ITER_EQ),
    ]
    if RUN_POLICY_CASES:
        cmd.append("--run-policy-cases")

    print("Running current one-market intergen model", flush=True)
    print(f"Started: {start_stamp}", flush=True)
    print(f"Source theta: {repo_root / SOURCE_RECORD}", flush=True)
    print(f"Output folder: {repo_root / OUTPUT_FOLDER}", flush=True)
    print(f"Target set: {TARGET_SET}", flush=True)
    print(f"Python used for solve: {packet_python}", flush=True)
    print(
        f"Grid: J={J}, Nb={NB}, Nz={INCOME_STATES}, H_own=[{OWNER_LADDER}], hR_max={RENTER_MAX_ROOMS}",
        flush=True,
    )
    print(f"max_iter_eq={MAX_ITER_EQ}", flush=True)
    print(flush=True)

    subprocess.run(cmd, cwd=str(model_dir), check=True)

    load_outputs_for_spyder(repo_root / OUTPUT_FOLDER)
    elapsed = time.perf_counter() - run_start
    end_stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    readme = repo_root / OUTPUT_FOLDER / "README.md"
    contact_sheet = repo_root / OUTPUT_FOLDER / "contact_sheet.png"
    print()
    print("Done.")
    print(f"Finished: {end_stamp}")
    print(f"Total runtime: {format_elapsed(elapsed)}")
    print(f"Start with: {readme}")
    print(f"Visual summary: {contact_sheet}")
    print("Loaded Spyder variables: solution_summary, moments, target_fit, age_profiles, room_bin_fit")


def format_elapsed(seconds: float) -> str:
    total = int(round(seconds))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m {secs}s"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def load_outputs_for_spyder(outdir: Path) -> None:
    """Expose the main run artifacts as globals for Spyder's Variable Explorer."""
    global output_folder, readme_path, contact_sheet_path
    global solution_summary, moments, target_fit, age_profiles, room_bin_fit

    output_folder = outdir
    readme_path = outdir / "README.md"
    contact_sheet_path = outdir / "contact_sheet.png"
    solution_summary = read_json(outdir / "solution_summary.json")
    moments = read_json(outdir / "moments.json")
    target_fit = read_csv_table(outdir / "target_fit.csv")
    age_profiles = read_csv_table(outdir / "age_profiles.csv")
    room_bin_fit = read_csv_table(outdir / "room_bin_fit_prime30_55_childless.csv")


def read_json(path: Path):
    with path.open() as fh:
        return json.load(fh)


def read_csv_table(path: Path):
    try:
        import pandas as pd

        return pd.read_csv(path)
    except Exception:
        with path.open(newline="") as fh:
            return list(csv.DictReader(fh))


if __name__ == "__main__":
    main()
