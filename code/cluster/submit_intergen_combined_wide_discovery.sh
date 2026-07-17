#!/bin/bash
#SBATCH --job-name=ihfwide
#SBATCH --output=logs/slurm_ihf_wide_%A_%a.out
#SBATCH --error=logs/slurm_ihf_wide_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RESULTS_ROOT="${WIDE_RESULTS_ROOT:?required}"
SEED_THETA="${WIDE_SEED_THETA:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
ARMS=(positive_all deterministic_tenure deterministic_tenure_no_bequest)
ARM="${WIDE_ARM:-${ARMS[$(((TASK_ID - 1) % 3))]}}"
OUTDIR="$RESULTS_ROOT/task_$TASK_ID"
mkdir -p "$OUTDIR"
export WIDE_SEED=$((${WIDE_SEED_BASE:-2026071300} + 1000 * TASK_ID))
export WIDE_MUTATION="${WIDE_MUTATION:-0.80}"
export WIDE_CROSSOVER="${WIDE_CROSSOVER:-0.75}"

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_combined_wide_discovery.py" \
  --seed-theta "$SEED_THETA" --outdir "$OUTDIR" --arm "$ARM"
