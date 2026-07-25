#!/bin/bash
# E5 wealth slice: one annual-beta row per task, all theta0 nodes per row.
#SBATCH --job-name=ihfe5w
#SBATCH --output=logs/slurm_ihfe5w_%A_%a.out
#SBATCH --error=logs/slurm_ihfe5w_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-8

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
BETA_VALUES=(0.94 0.9483 0.9566 0.9649 0.9732 0.9815 0.9898 0.9981)
if (( TASK_ID < 1 || TASK_ID > ${#BETA_VALUES[@]} )); then
  echo "Expected SLURM_ARRAY_TASK_ID in 1..${#BETA_VALUES[@]}" >&2
  exit 2
fi
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e5_wealth_slice.py" \
  --source "$PROJECT_ROOT/output/model/eqscale_seq_e5_recalibration_20260724/report/results.json" \
  --outdir "$PROJECT_ROOT/output/model/eqscale_seq_e5_wealth_slice_20260725/task_${TASK_ID}" \
  --beta-annual "${BETA_VALUES[$((TASK_ID - 1))]}" \
  --theta0 "0.05,0.1117,0.25,0.5,1.0,2.0,4.0,8.0"
