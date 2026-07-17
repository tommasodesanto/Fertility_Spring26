#!/bin/bash
# Five-cell conditional profile for the internally estimated M4 luxury shifter.
# Each noncentral cell fixes theta1 and re-estimates the other 12 coordinates;
# the central cell reproduces the collected M4 winner exactly.
#SBATCH --job-name=ihfstdbqp
#SBATCH --output=logs/slurm_ihf_stdbqp_%A_%a.out
#SBATCH --error=logs/slurm_ihf_stdbqp_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:35:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-5

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RUN_ROOT="$PROJECT_ROOT/output/model/intergen_standard_bequest_recalibration_20260716"
WINNER_JSON="${STANDARD_BEQUEST_WINNER_JSON:-$RUN_ROOT/report/results.json}"
PROFILE_ROOT="${STANDARD_BEQUEST_PROFILE_ROOT:-$RUN_ROOT/theta1_profile}"
ACCEPTANCE_CSV="${STANDARD_BEQUEST_ACCEPTANCE_CSV:-$PROJECT_ROOT/code/data/mms_center_periphery/output_ownership_audit/acs_ownership_4year_acceptance_bins_6284.csv}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 5 ]; then
  echo "M4 theta1 profile requires SLURM_ARRAY_TASK_ID in 1..5" >&2
  exit 2
fi

CELL_VALUES="$($PYTHON_BIN - "$WINNER_JSON" "$TASK_ID" <<'PY'
import json
import math
import sys

winner_path, task_id = sys.argv[1], int(sys.argv[2])
winner = json.load(open(winner_path))["winners"]["M4"]
theta1 = float(winner["theta"]["theta1"])
lo, hi = 0.02, 16.0
unit = (math.log(theta1) - math.log(lo)) / (math.log(hi) - math.log(lo))
if unit <= 0.10:
    offsets = (0.0, 0.025, 0.05, 0.10, 0.20)
elif unit >= 0.90:
    offsets = (-0.20, -0.10, -0.05, -0.025, 0.0)
else:
    offsets = (-0.10, -0.05, 0.0, 0.05, 0.10)
offset = offsets[task_id - 1]
profile_unit = unit + offset
profile_theta1 = math.exp(math.log(lo) + profile_unit * (math.log(hi) - math.log(lo)))
print(f"{profile_theta1:.17g} {offset:.17g} {theta1:.17g}")
PY
)"
read -r THETA1_VALUE UNIT_OFFSET WINNER_THETA1 <<< "$CELL_VALUES"
OUTDIR="$PROFILE_ROOT/task_${TASK_ID}_offset_${UNIT_OFFSET}"
mkdir -p "$OUTDIR"
$PYTHON_BIN - "$OUTDIR/profile_cell.json" "$TASK_ID" "$UNIT_OFFSET" "$THETA1_VALUE" "$WINNER_THETA1" <<'PY'
import json
import sys

path, task_id, offset, theta1, winner_theta1 = sys.argv[1:]
with open(path, "w") as handle:
    json.dump(
        {
            "task_id": int(task_id),
            "theta1_unit_offset": float(offset),
            "theta1": float(theta1),
            "winner_theta1": float(winner_theta1),
        },
        handle,
        indent=2,
        sort_keys=True,
    )
PY

MAX_EVALS="${STANDARD_BEQUEST_PROFILE_MAX_EVALS:-400}"
MINUTES="${STANDARD_BEQUEST_PROFILE_MINUTES:-25}"
if [ "$UNIT_OFFSET" = "0" ]; then
  MAX_EVALS=1
  MINUTES=6
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" \
  --seed-record "$WINNER_JSON" \
  --seed-arm M4 \
  --acceptance-csv "$ACCEPTANCE_CSV" \
  --outdir "$OUTDIR" \
  --arm M4_PROFILE \
  --theta1 "$THETA1_VALUE" \
  --seed "$((2026071680 + TASK_ID))" \
  --start-mix 0 \
  --max-evals "$MAX_EVALS" \
  --minutes "$MINUTES" \
  --max-iter-eq 10 \
  --tol-eq 1e-4
