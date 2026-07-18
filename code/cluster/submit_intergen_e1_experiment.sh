#!/bin/bash
# E1 eqscale/sequential-fertility recalibration: eight independent NM chains.
#SBATCH --job-name=ihfe1
#SBATCH --output=logs/slurm_ihfe1_%A_%a.out
#SBATCH --error=logs/slurm_ihfe1_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
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
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 8 ]; then
  echo "E1 requires SLURM_ARRAY_TASK_ID in 1..8" >&2
  exit 2
fi

RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_recalibration_20260718/production"
OUTDIR="$RUN_ROOT/chain_${TASK_ID}"
ARGS=(
  --outdir "$OUTDIR"
  --seed "$((2026071800 + SLURM_ARRAY_TASK_ID))"
  --minutes 225
)
if [ "${E1_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke --minutes 10)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq/run_e1_chain.py" "${ARGS[@]}"
