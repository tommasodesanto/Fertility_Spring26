#!/bin/bash
# Joint internal calibration of theta0, theta1, and theta_n after the M1
# mortality-only specification. Four chains use distinct transformed starts.
#SBATCH --job-name=ihfintbq
#SBATCH --output=logs/slurm_ihf_intbq_%A_%a.out
#SBATCH --error=logs/slurm_ihf_intbq_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:40:00
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

RESULTS_ROOT="${INTERNAL_BEQUEST_RESULTS_ROOT:?required}"
SEED_RECORD="${INTERNAL_BEQUEST_SEED_RECORD:?required}"
ACCEPTANCE_CSV="${INTERNAL_BEQUEST_ACCEPTANCE_CSV:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"

CHAIN_IDX=$((TASK_ID - 1))
START_MIXES=(0.00 0.04 0.08 0.12 0.02)
START_MIX="${INTERNAL_BEQUEST_START_MIX:-${START_MIXES[$CHAIN_IDX]}}"
OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_M3_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"

ARGS=(
  --seed-record "$SEED_RECORD"
  --seed-arm M1
  --acceptance-csv "$ACCEPTANCE_CSV"
  --outdir "$OUTDIR"
  --arm M3
  --theta1 "${INTERNAL_BEQUEST_SEED_THETA1:-0.55}"
  --seed-theta0 "${INTERNAL_BEQUEST_SEED_THETA0:-0.30}"
  --seed-theta-n "${INTERNAL_BEQUEST_SEED_THETA_N:-0.75}"
  --seed "$((${INTERNAL_BEQUEST_SEED_BASE:-2026071530} + TASK_ID))"
  --start-mix "$START_MIX"
  --max-evals "${INTERNAL_BEQUEST_MAX_EVALS:-1000}"
  --minutes "${INTERNAL_BEQUEST_MINUTES:-90}"
  --max-iter-eq "${INTERNAL_BEQUEST_MAX_ITER_EQ:-10}"
  --tol-eq "${INTERNAL_BEQUEST_TOL_EQ:-1e-4}"
)
if [ "${INTERNAL_BEQUEST_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" "${ARGS[@]}"
