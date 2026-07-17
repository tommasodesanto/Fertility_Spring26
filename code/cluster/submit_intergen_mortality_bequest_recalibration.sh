#!/bin/bash
#SBATCH --job-name=ihfmortbq
#SBATCH --output=logs/slurm_ihf_mortbq_%A_%a.out
#SBATCH --error=logs/slurm_ihf_mortbq_%A_%a.err
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

RESULTS_ROOT="${MORTALITY_BEQUEST_RESULTS_ROOT:?required}"
SEED_RECORD="${MORTALITY_BEQUEST_SEED_RECORD:?required}"
ACCEPTANCE_CSV="${MORTALITY_ACCEPTANCE_CSV:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"

CHAIN_IDX=$((TASK_ID - 1))
START_MIXES=(0.00 0.04 0.08 0.12)
START_MIX="${START_MIXES[$CHAIN_IDX]}"
OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_M2_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"

ARGS=(
  --seed-record "$SEED_RECORD"
  --seed-arm M1
  --acceptance-csv "$ACCEPTANCE_CSV"
  --outdir "$OUTDIR"
  --arm M2
  --theta1 "${MORTALITY_BEQUEST_THETA1:-0.25}"
  --seed "$((${MORTALITY_BEQUEST_SEED_BASE:-2026071510} + TASK_ID))"
  --start-mix "$START_MIX"
  --max-evals "${MORTALITY_BEQUEST_MAX_EVALS:-1000}"
  --minutes "${MORTALITY_BEQUEST_MINUTES:-90}"
  --max-iter-eq "${MORTALITY_BEQUEST_MAX_ITER_EQ:-10}"
  --tol-eq "${MORTALITY_BEQUEST_TOL_EQ:-1e-4}"
)
if [ "${MORTALITY_BEQUEST_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" "${ARGS[@]}"
