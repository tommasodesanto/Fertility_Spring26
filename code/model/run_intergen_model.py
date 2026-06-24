#!/usr/bin/env python3
"""One-click runner for the current June 2026 intergen model.

Open this file and run it. It runs the active one-market intergenerational
housing-fertility model from the saved current diagnostic theta, then writes a
mechanics packet with plots and tables under `output/model/`.

This is the file to use when someone asks: "run the model and show me the
results." It is diagnostic, not a production calibration search.
"""

from __future__ import annotations

import os
import subprocess
import sys
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
    repo_root = Path(__file__).resolve().parents[2]
    model_dir = Path(__file__).resolve().parent
    venv_python = model_dir / ".venv/bin/python"

    # If this file is launched from VS Code, PyCharm, Finder, or another
    # generic Python, rerun it inside the project venv so imports are stable.
    if venv_python.exists() and Path(sys.executable).resolve() != venv_python.resolve():
        os.execv(str(venv_python), [str(venv_python), str(Path(__file__).resolve())])

    packet_tool = model_dir / "tools/build_intergen_mechanics_packet.py"
    cmd = [
        sys.executable,
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
    print(f"Source theta: {repo_root / SOURCE_RECORD}", flush=True)
    print(f"Output folder: {repo_root / OUTPUT_FOLDER}", flush=True)
    print(f"Target set: {TARGET_SET}", flush=True)
    print(
        f"Grid: J={J}, Nb={NB}, Nz={INCOME_STATES}, H_own=[{OWNER_LADDER}], hR_max={RENTER_MAX_ROOMS}",
        flush=True,
    )
    print(f"max_iter_eq={MAX_ITER_EQ}", flush=True)
    print(flush=True)

    subprocess.run(cmd, cwd=str(model_dir), check=True)

    readme = repo_root / OUTPUT_FOLDER / "README.md"
    contact_sheet = repo_root / OUTPUT_FOLDER / "contact_sheet.png"
    print()
    print("Done.")
    print(f"Start with: {readme}")
    print(f"Visual summary: {contact_sheet}")


if __name__ == "__main__":
    main()
