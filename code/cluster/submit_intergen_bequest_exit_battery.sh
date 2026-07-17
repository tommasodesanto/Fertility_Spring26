#!/bin/bash
#SBATCH --job-name=ihfbqx
#SBATCH --output=logs/slurm_ihf_bqx_%A_%a.out
#SBATCH --error=logs/slurm_ihf_bqx_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RESULTS_ROOT="${BEQUEST_EXIT_RESULTS_ROOT:?required}"
SEED_RECORD="${BEQUEST_EXIT_SEED_RECORD:?required}"
ACCEPTANCE_CSV="${BEQUEST_EXIT_ACCEPTANCE_CSV:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"

# Eleven specifications, each with two explicit starts.  The LTV variants are
# external robustness cases and are never selected by scalar loss.
ARMS=(A0 A1 A2 A2 A2 A3 A3 A3 A4 A5 A3)
LBARS=(0.4 0.4 0.2 0.4 0.6 0.2 0.4 0.6 0.4 0.4 0.4)
THETA1S=(0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.50)
CASE_IDX=$(((TASK_ID - 1) / 2))
CHAIN_IDX=$(((TASK_ID - 1) % 2))
ARM="${ARMS[$CASE_IDX]}"
LBAR="${LBARS[$CASE_IDX]}"
THETA1="${THETA1S[$CASE_IDX]}"
START_MIX=0
[ "$CHAIN_IDX" -eq 0 ] || START_MIX=0.08

OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_${ARM}_lbar${LBAR}_theta1${THETA1}_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"
ARGS=(
  --seed-record "$SEED_RECORD"
  --acceptance-csv "$ACCEPTANCE_CSV"
  --outdir "$OUTDIR"
  --arm "$ARM"
  --ltv-terminal "$LBAR"
  --theta1 "$THETA1"
  --seed "$((${BEQUEST_EXIT_SEED_BASE:-2026071400} + TASK_ID))"
  --start-mix "$START_MIX"
  --max-evals "${BEQUEST_EXIT_MAX_EVALS:-2000}"
  --minutes "${BEQUEST_EXIT_MINUTES:-230}"
  --max-iter-eq "${BEQUEST_EXIT_MAX_ITER_EQ:-10}"
  --tol-eq "${BEQUEST_EXIT_TOL_EQ:-1e-4}"
)
[ "${BEQUEST_EXIT_SMOKE:-0}" = "1" ] || :
if [ "${BEQUEST_EXIT_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" "${ARGS[@]}"
