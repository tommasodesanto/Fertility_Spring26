#!/bin/bash
#SBATCH --job-name=ihfcombo
#SBATCH --output=logs/slurm_ihf_combo_%A_%a.out
#SBATCH --error=logs/slurm_ihf_combo_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

MODE="${COMBINED_MODE:-fixed}"
RESULTS_ROOT="${COMBINED_RESULTS_ROOT:?required}"
SEED_THETA="${COMBINED_SEED_THETA:?required}"
TARGET_CSV="${COMBINED_TARGET_CSV:?required}"
mkdir -p "$RESULTS_ROOT"

if [ "$MODE" = smoke ]; then
  exec "$PYTHON_BIN" "$MODEL_DIR/tools/compare_intergen_combined_specification.py" \
    --mode smoke --theta-json "$SEED_THETA" --outdir "$RESULTS_ROOT" --quiet
fi

ARMS=(current bequest_normalization_only finance_supply_revisions_only all_revisions)
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
ARM="${ARMS[$((TASK_ID - 1))]}"
OUTDIR="$RESULTS_ROOT/$ARM"
mkdir -p "$OUTDIR"
exec "$PYTHON_BIN" "$MODEL_DIR/tools/compare_intergen_combined_specification.py" \
  --mode fixed --arm "$ARM" --theta-json "$SEED_THETA" \
  --housing-targets-csv "$TARGET_CSV" --outdir "$OUTDIR" \
  --smoke-path "$RESULTS_ROOT/inversion_smoke.json" --quiet
