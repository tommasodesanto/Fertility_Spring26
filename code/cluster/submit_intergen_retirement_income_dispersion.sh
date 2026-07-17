#!/bin/bash
#SBATCH --job-name=ihfretinc
#SBATCH --output=logs/slurm_ihf_retinc_%A_%a.out
#SBATCH --error=logs/slurm_ihf_retinc_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:55:00
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

RESULTS_ROOT="${RETIREMENT_DISP_RESULTS_ROOT:?required}"
WINNER_JSON="${RETIREMENT_DISP_WINNER_JSON:?required}"
CELL="${SLURM_ARRAY_TASK_ID:?array task required}"
OUTDIR="$RESULTS_ROOT/cell_$(printf '%02d' "$CELL")"
mkdir -p "$OUTDIR"

ARGS=(
  --winner-json "$WINNER_JSON"
  --winner-arm M3
  --outdir "$OUTDIR"
  --cell "$CELL"
  --experiment retirement_income_dispersion
  --seed "$((${RETIREMENT_DISP_SEED_BASE:-2026071700} + CELL))"
  --max-evals "${RETIREMENT_DISP_MAX_EVALS:-120}"
  --minutes "${RETIREMENT_DISP_MINUTES:-85}"
  --initial-step "${RETIREMENT_DISP_INITIAL_STEP:-0.10}"
  --min-step "${RETIREMENT_DISP_MIN_STEP:-0.003}"
  --max-iter-eq "${RETIREMENT_DISP_MAX_ITER_EQ:-10}"
  --tol-eq "${RETIREMENT_DISP_TOL_EQ:-1e-4}"
)
if [ "${RETIREMENT_DISP_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/profile_intergen_bequest_reachability.py" "${ARGS[@]}"
