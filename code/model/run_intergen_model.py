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
import os
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

# Set this to True after at least one successful solve has written
# output/model/intergen_model_run_current/solution_cache.pkl. It rebuilds all
# plots and CSV diagnostics from that saved solution object without solving.
REFRESH_PLOTS_FROM_SAVED_SOLUTION = False

# Keep cache refreshes interactive: update first-look CSVs/plots, but do not
# spend time rebuilding the larger diagnostics folder or contact sheet.
FAST_REFRESH_FROM_SAVED_SOLUTION = True


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
    if REFRESH_PLOTS_FROM_SAVED_SOLUTION:
        cmd.append("--refresh-plots-from-cache")
        if FAST_REFRESH_FROM_SAVED_SOLUTION:
            cmd.extend(["--skip-standard-diagnostics", "--skip-contact-sheet", "--quick-first-look-only"])

    if REFRESH_PLOTS_FROM_SAVED_SOLUTION:
        print("Refreshing current one-market intergen graphs from saved solution", flush=True)
    else:
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
    print(f"refresh_plots_from_saved_solution={REFRESH_PLOTS_FROM_SAVED_SOLUTION}", flush=True)
    print(f"fast_refresh_from_saved_solution={FAST_REFRESH_FROM_SAVED_SOLUTION}", flush=True)
    if REFRESH_PLOTS_FROM_SAVED_SOLUTION:
        print(f"Saved solution cache: {repo_root / OUTPUT_FOLDER / 'solution_cache.pkl'}", flush=True)
    else:
        print_runtime_diagnostics("Before solve")
    print(flush=True)

    subprocess.run(cmd, cwd=str(model_dir), check=True)

    load_outputs_for_spyder(repo_root / OUTPUT_FOLDER)
    print_solver_timing_summary(solution_summary, cached=REFRESH_PLOTS_FROM_SAVED_SOLUTION)
    if not REFRESH_PLOTS_FROM_SAVED_SOLUTION:
        print_runtime_diagnostics("After solve")
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
    print(
        "Loaded Spyder variables: solution_summary, moments, target_fit, age_profiles, "
        "room_bin_fit, first_look_path, first_look_full_path, first_look_density_path, "
        "first_look_total_wealth_path, first_look_total_wealth_density_path, "
        "solution_cache_path, first_look_policy_lines, first_look_market_summary"
    )


def print_runtime_diagnostics(label: str) -> None:
    print()
    print(f"{label} runtime diagnostics:")
    print(f"  CPU power: {mac_cpu_power_status()}")
    print(f"  Thread env: {thread_environment_summary()}")


def mac_cpu_power_status() -> str:
    """Return compact macOS CPU power limits when available."""
    try:
        proc = subprocess.run(
            ["pmset", "-g", "therm"],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (FileNotFoundError, subprocess.SubprocessError):
        return "not available"
    text = (proc.stdout or proc.stderr or "").splitlines()
    values: dict[str, str] = {}
    for line in text:
        if "=" not in line or "CPU_" not in line:
            continue
        name, value = line.split("=", 1)
        values[name.strip()] = value.strip()
    if not values:
        return "not available"
    pieces = []
    for key in ["CPU_Speed_Limit", "CPU_Scheduler_Limit", "CPU_Available_CPUs"]:
        if key in values:
            pieces.append(f"{key}={values[key]}")
    return ", ".join(pieces) if pieces else "not available"


def thread_environment_summary() -> str:
    keys = [
        "NUMBA_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ]
    return ", ".join(f"{key}={os.environ.get(key, 'unset')}" for key in keys)


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
    global solution_cache_path
    global solution_summary, moments, target_fit, age_profiles, room_bin_fit
    global first_look_path, first_look_full_path, first_look_density_path
    global first_look_total_wealth_path, first_look_total_wealth_full_path
    global first_look_total_wealth_density_path, first_look_total_wealth_density
    global total_wealth_grid_coverage
    global first_look_policy_lines, first_look_market_summary

    output_folder = outdir
    readme_path = outdir / "README.md"
    contact_sheet_path = outdir / "contact_sheet.png"
    solution_cache_path = outdir / "solution_cache.pkl"
    first_look_path = outdir / "first_look_policies_markets.png"
    first_look_full_path = outdir / "first_look_policies_markets_full.png"
    first_look_density_path = outdir / "first_look_wealth_density.png"
    first_look_total_wealth_path = outdir / "first_look_policies_markets_total_wealth.png"
    first_look_total_wealth_full_path = outdir / "first_look_policies_markets_total_wealth_full.png"
    first_look_total_wealth_density_path = outdir / "first_look_total_wealth_density.png"
    solution_summary = read_json(outdir / "solution_summary.json")
    moments = read_json(outdir / "moments.json")
    target_fit = read_csv_table(outdir / "target_fit.csv")
    age_profiles = read_csv_table(outdir / "age_profiles.csv")
    room_bin_fit = read_csv_table(outdir / "room_bin_fit_prime30_55_childless.csv")
    first_look_policy_lines = read_csv_table(outdir / "first_look_policy_lines.csv")
    first_look_market_summary = read_csv_table(outdir / "first_look_market_summary.csv")
    first_look_total_wealth_density = read_csv_table(outdir / "first_look_total_wealth_density.csv")
    total_wealth_grid_coverage = read_csv_table(outdir / "total_wealth_grid_coverage.csv")


def print_solver_timing_summary(summary: dict, *, cached: bool = False) -> None:
    timings = summary.get("solver_timings") if isinstance(summary, dict) else None
    if not isinstance(timings, dict) or not timings:
        return
    refine = timings.get("scalar_market_refine")
    print()
    label = "Cached solve runtime" if cached else "Solve runtime"
    print(f"{label}: {format_elapsed(float(summary.get('elapsed_sec', 0.0)))}")
    print(f"Equilibrium iterations: {timings.get('iterations_completed')}")
    if isinstance(refine, dict) and refine.get("used"):
        print(
            "Scalar price refinement: "
            f"{refine.get('iterations')} iterations, "
            f"{refine.get('price_evaluations')} price evaluations, "
            f"{format_elapsed(float(refine.get('price_evaluation_time_sec', 0.0)))}"
        )
    if "markov_income_solve_time" in timings:
        print(f"Damped-loop solve time: {format_elapsed(float(timings['markov_income_solve_time']))}")
    if "distribution" in timings:
        print(f"Final distribution/statistics time: {format_elapsed(float(timings['distribution']))}")


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
