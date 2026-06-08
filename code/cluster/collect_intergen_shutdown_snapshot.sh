#!/bin/bash
# Snapshot intergen_housing_fertility panel outputs before Torch downtime.
# Run from:
#   $SCRATCH/projects/Fertility_Spring26/code/cluster

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"

SNAPSHOT_TAG="${1:-$(date +%Y%m%d_%H%M%S)}"
SNAPSHOT_DIR="${SCRIPT_DIR}/shutdown_snapshots/${SNAPSHOT_TAG}"
mkdir -p "${SNAPSHOT_DIR}"

JOB_IDS="${INTERGEN_SHUTDOWN_JOB_IDS:-10626531 10626676 10628379 10628390 10628391 10628618}"

if command -v module >/dev/null 2>&1; then
    module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || module load python/3.12 2>/dev/null || module load python/3.11 2>/dev/null || true
fi

PYTHON_BIN="${INTERGEN_PYTHON:-$(command -v python3 || command -v python)}"
export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"

echo "Snapshot tag: ${SNAPSHOT_TAG}" | tee "${SNAPSHOT_DIR}/snapshot.log"
echo "Started: $(date)" | tee -a "${SNAPSHOT_DIR}/snapshot.log"
echo "Python: ${PYTHON_BIN}" | tee -a "${SNAPSHOT_DIR}/snapshot.log"
echo "Job IDs: ${JOB_IDS}" | tee -a "${SNAPSHOT_DIR}/snapshot.log"

squeue -u "${USER:-td2248}" > "${SNAPSHOT_DIR}/squeue.txt" 2>&1 || true
sacct -j "$(echo "${JOB_IDS}" | tr ' ' ',')" --format=JobID,JobName%12,State,ExitCode,Elapsed,Submit,Start,End,ReqMem,MaxRSS%12 -P \
    > "${SNAPSHOT_DIR}/sacct.txt" 2>&1 || true

mapfile -t RUN_DIRS < <(
    find "${SCRIPT_DIR}" -maxdepth 1 -type d \
        \( -name 'results_intergen_housing_fertility_intergen_candidate_no_timing_v0_twohour_20260608' \
        -o -name 'results_intergen_housing_fertility_intergen_candidate_no_timing_v0_overnight*_3g_20260608' \
        -o -name 'results_intergen_housing_fertility_intergen_candidate_no_timing_v0_burst*_3g_20260608' \) \
        | sort
)

printf "%s\n" "${RUN_DIRS[@]}" > "${SNAPSHOT_DIR}/run_dirs.txt"

for run_dir in "${RUN_DIRS[@]}"; do
    echo "Collecting ${run_dir}" | tee -a "${SNAPSHOT_DIR}/snapshot.log"
    "${PYTHON_BIN}" "${MODEL_DIR}/tools/collect_intergen_panel_results.py" \
        --results-dir "${run_dir}" \
        --outdir "${run_dir}" \
        --top-n 200 \
        > "${SNAPSHOT_DIR}/collect_$(basename "${run_dir}").log" 2>&1 || true
done

"${PYTHON_BIN}" - "${SNAPSHOT_DIR}" "${RUN_DIRS[@]}" <<'PY'
import csv
import json
import math
import sys
from pathlib import Path

snapshot_dir = Path(sys.argv[1])
run_dirs = [Path(p) for p in sys.argv[2:]]
records = []
for run_dir in run_dirs:
    for case_path in sorted(run_dir.glob("task_*/cases.jsonl")):
        task_id = case_path.parent.name.removeprefix("task_")
        for line in case_path.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            try:
                loss = float(rec.get("rank_loss", math.inf))
            except (TypeError, ValueError):
                loss = math.inf
            if not math.isfinite(loss):
                continue
            rec["run_dir"] = run_dir.name
            rec["task_id"] = task_id
            records.append(rec)

records.sort(key=lambda r: float(r.get("rank_loss", math.inf)))
fields = [
    "run_dir", "task_id", "case", "label", "rank_loss",
    "full_old_nonlocation_loss", "market_residual", "elapsed_sec",
    "tfr", "childless_rate", "own_rate", "own_family_gap",
    "housing_increment_0to1", "housing_increment_1to2",
    "young_liquid_wealth_to_income", "old_age_own_rate",
    "old_age_parent_childless_gap", "liquid_wealth_to_income",
    "housing_user_cost_share", "prime_childless_renter_median_rooms",
    "prime_childless_owner_median_rooms", "mean_age_first_birth",
    "parity_share_0", "parity_share_1", "parity_share_2plus",
    "theta_json",
]
with (snapshot_dir / "combined_top200.csv").open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields)
    writer.writeheader()
    for rec in records[:200]:
        moments = dict(rec.get("moments", {}))
        row = {
            "run_dir": rec.get("run_dir"),
            "task_id": rec.get("task_id"),
            "case": rec.get("case"),
            "label": rec.get("label"),
            "rank_loss": rec.get("rank_loss"),
            "full_old_nonlocation_loss": rec.get("full_old_nonlocation_loss"),
            "market_residual": rec.get("market_residual"),
            "elapsed_sec": rec.get("elapsed_sec"),
            "theta_json": json.dumps(rec.get("theta", {}), sort_keys=True),
        }
        for name in fields:
            if name in moments:
                row[name] = moments.get(name)
        writer.writerow(row)

summary = {
    "n_records": len(records),
    "n_run_dirs": len(run_dirs),
    "run_dirs": [p.name for p in run_dirs],
    "best": None,
}
if records:
    best = records[0]
    summary["best"] = {
        "run_dir": best.get("run_dir"),
        "task_id": best.get("task_id"),
        "case": best.get("case"),
        "label": best.get("label"),
        "rank_loss": best.get("rank_loss"),
        "full_old_nonlocation_loss": best.get("full_old_nonlocation_loss"),
        "market_residual": best.get("market_residual"),
        "moments": best.get("moments", {}),
        "theta": best.get("theta", {}),
    }
(snapshot_dir / "combined_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
print(json.dumps(summary, indent=2, sort_keys=True))
PY

TAR_PATH="${SCRIPT_DIR}/shutdown_snapshots/${SNAPSHOT_TAG}.tar.gz"
tar -czf "${TAR_PATH}" \
    -C "${SCRIPT_DIR}" \
    "shutdown_snapshots/${SNAPSHOT_TAG}" \
    $(for d in "${RUN_DIRS[@]}"; do basename "$d"; done) \
    logs/slurm_ihf_2hr_10626531_*.out logs/slurm_ihf_2hr_10626531_*.err \
    logs/slurm_ihf_2hr_10626676_*.out logs/slurm_ihf_2hr_10626676_*.err \
    logs/slurm_ihf_2hr_10628379_*.out logs/slurm_ihf_2hr_10628379_*.err \
    logs/slurm_ihf_2hr_10628390_*.out logs/slurm_ihf_2hr_10628390_*.err \
    logs/slurm_ihf_2hr_10628391_*.out logs/slurm_ihf_2hr_10628391_*.err \
    logs/slurm_ihf_2hr_10628618_*.out logs/slurm_ihf_2hr_10628618_*.err \
    2> "${SNAPSHOT_DIR}/tar_warnings.log" || true

echo "Tarball: ${TAR_PATH}" | tee -a "${SNAPSHOT_DIR}/snapshot.log"

if [ "${INTERGEN_SHUTDOWN_CANCEL:-0}" = "1" ]; then
    echo "Cancelling jobs: ${JOB_IDS}" | tee -a "${SNAPSHOT_DIR}/snapshot.log"
    scancel ${JOB_IDS} 2>&1 | tee -a "${SNAPSHOT_DIR}/snapshot.log" || true
fi

echo "Finished: $(date)" | tee -a "${SNAPSHOT_DIR}/snapshot.log"
