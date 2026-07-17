#!/bin/bash
#SBATCH --job-name=ihfwpol
#SBATCH --output=logs/slurm_ihf_wpol_%A_%a.out
#SBATCH --error=logs/slurm_ihf_wpol_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=04:00:00
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

RESULTS_ROOT="${WIDE_POLISH_RESULTS_ROOT:?required}"
SEED_THETA="${WIDE_POLISH_SEED_THETA:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
OUTDIR="$RESULTS_ROOT/task_$TASK_ID"
mkdir -p "$OUTDIR"
export WIDE_POLISH_SEED=$((${WIDE_POLISH_SEED_BASE:-2026071300} + 1000 * TASK_ID))

ARGS=(
  --seed-theta "$SEED_THETA"
  --outdir "$OUTDIR"
  --arm "${WIDE_POLISH_ARM:-deterministic_tenure}"
)
if [ "${WIDE_POLISH_PROFILE:-0}" = "1" ]; then
  H0_FLOORS=(0.25 0.40 0.55 0.75 1.00)
  CHI_RUNGS=(0.90 1.00 1.10 1.25 1.50)
  H_IDX=$(((TASK_ID - 1) / 5))
  C_IDX=$(((TASK_ID - 1) % 5))
  ARGS+=(--fixed-h-bar-0 "${H0_FLOORS[$H_IDX]}" --fixed-chi "${CHI_RUNGS[$C_IDX]}")
else
  [ -z "${WIDE_POLISH_FIXED_H_BAR_0:-}" ] || ARGS+=(--fixed-h-bar-0 "$WIDE_POLISH_FIXED_H_BAR_0")
  [ -z "${WIDE_POLISH_FIXED_CHI:-}" ] || ARGS+=(--fixed-chi "$WIDE_POLISH_FIXED_CHI")
fi

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_combined_wide_polish.py" "${ARGS[@]}"
