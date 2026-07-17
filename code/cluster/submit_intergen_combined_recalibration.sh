#!/bin/bash
#SBATCH --job-name=ihfrecal
#SBATCH --output=logs/slurm_ihf_recal_%A_%a.out
#SBATCH --error=logs/slurm_ihf_recal_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RESULTS_ROOT="${COMBINED_RECAL_RESULTS_ROOT:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ -n "${COMBINED_PARENT_RESULTS_ROOT:-}" ]; then
  SEED_THETA="$COMBINED_PARENT_RESULTS_ROOT/task_$TASK_ID/best.json"
  if [ ! -s "$SEED_THETA" ]; then
    echo "Missing parent-wave seed: $SEED_THETA" >&2
    exit 2
  fi
else
  SEED_THETA="${COMBINED_SEED_THETA:?required when COMBINED_PARENT_RESULTS_ROOT is unset}"
fi
OUTDIR="$RESULTS_ROOT/task_$TASK_ID"
mkdir -p "$OUTDIR"
SEED_BASE="${COMBINED_SEED_BASE:-20260710}"
export COMBINED_SEED=$((SEED_BASE + 1000 * TASK_ID))
STEPS=(0.020 0.035 0.050 0.075)
export COMBINED_INITIAL_STEP="${COMBINED_INITIAL_STEP:-${STEPS[$(((TASK_ID - 1) % 4))]}}"
if [ "${COMBINED_DIVERSIFY_STARTS:-0}" = "1" ]; then
  MIXES=(0.000 0.030 0.060 0.100 0.150 0.220)
  export COMBINED_START_MIX="${COMBINED_START_MIX:-${MIXES[$(((TASK_ID - 1) / 4))]}}"
else
  export COMBINED_START_MIX="${COMBINED_START_MIX:-0.0}"
fi
export COMBINED_METHOD="${COMBINED_METHOD:-pattern}"
export COMBINED_MAX_EVALS="${COMBINED_MAX_EVALS:-1000}"
export COMBINED_MINUTES="${COMBINED_MINUTES:-230}"
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_combined_recalibration.py" \
  --seed-theta "$SEED_THETA" --outdir "$OUTDIR"
