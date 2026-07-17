#!/bin/bash
#SBATCH --job-name=ihfbqprof
#SBATCH --output=logs/slurm_ihf_bqprof_%A_%a.out
#SBATCH --error=logs/slurm_ihf_bqprof_%A_%a.err
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

RESULTS_ROOT="${BEQUEST_EXIT_PROFILE_ROOT:?required}"
WINNER_JSON="${BEQUEST_EXIT_WINNER_JSON:?required}"
ACCEPTANCE_CSV="${BEQUEST_EXIT_ACCEPTANCE_CSV:?required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
THETA0S=(0 0.05 0.10 0.20 0.40)
PROFILE_IDX=$(((TASK_ID - 1) / 2))
CHAIN_IDX=$(((TASK_ID - 1) % 2))
THETA0="${THETA0S[$PROFILE_IDX]}"
START_MIX=0
[ "$CHAIN_IDX" -eq 0 ] || START_MIX=0.08
OUTDIR="$RESULTS_ROOT/task_${TASK_ID}_theta0${THETA0}_chain${CHAIN_IDX}"
mkdir -p "$OUTDIR"

exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_intergen_bequest_exit_chain.py" \
  --seed-record "$WINNER_JSON" \
  --acceptance-csv "$ACCEPTANCE_CSV" \
  --outdir "$OUTDIR" \
  --arm A3_PROFILE \
  --ltv-terminal 0.4 \
  --theta1 0.25 \
  --fixed-theta0 "$THETA0" \
  --seed "$((${BEQUEST_EXIT_PROFILE_SEED_BASE:-2026071500} + TASK_ID))" \
  --start-mix "$START_MIX" \
  --max-evals "${BEQUEST_EXIT_PROFILE_MAX_EVALS:-2000}" \
  --minutes "${BEQUEST_EXIT_PROFILE_MINUTES:-230}" \
  --max-iter-eq "${BEQUEST_EXIT_MAX_ITER_EQ:-10}" \
  --tol-eq "${BEQUEST_EXIT_TOL_EQ:-1e-4}"
