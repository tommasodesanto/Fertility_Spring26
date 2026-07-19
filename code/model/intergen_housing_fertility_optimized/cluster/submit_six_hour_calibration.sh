#!/bin/bash
# Eight-chain, six-hour continuation of the live M5 calibration in the isolated
# optimized package. Submit from this package's cluster directory.
#SBATCH --job-name=ihfopt6
#SBATCH --output=slurm_ihfopt6_%A_%a.out
#SBATCH --error=slurm_ihfopt6_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-8%8

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}/intergen_housing_fertility_optimized}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?array task required}"
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 8 ]; then
  echo "six-hour continuation requires array task 1..8" >&2
  exit 2
fi

METHODS=(nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern)
STEPS=(0.003 0.006 0.010 0.015 0.025 0.040 0.060 0.090)
MIXES=(0.0000 0.0000 0.0025 0.0050 0.0100 0.0200 0.0400 0.0800)
INDEX=$((TASK_ID - 1))
METHOD="${METHODS[$INDEX]}"
STEP="${STEPS[$INDEX]}"
MIX="${MIXES[$INDEX]}"
RUN_TAG="${OPTIMIZED_M5_RUN_TAG:-intergen_optimized_m5_continuation_20260719}"
RUN_ROOT="${PROJECT_ROOT}/output/model/${RUN_TAG}"
OUTDIR="${RUN_ROOT}/chains/chain_${TASK_ID}"
mkdir -p "$NUMBA_CACHE_DIR" "$OUTDIR"

ARGS=(
  --outdir "$OUTDIR"
  --seed "$((2026071900 + TASK_ID))"
  --method "$METHOD"
  --initial-step "$STEP"
  --start-mix "$MIX"
  --minimum-step "${OPTIMIZED_M5_MINIMUM_STEP:-0.00025}"
  --minutes "${OPTIMIZED_M5_MINUTES:-350}"
  --strict-reserve-minutes "${OPTIMIZED_M5_STRICT_RESERVE_MINUTES:-10}"
  --max-evals "${OPTIMIZED_M5_MAX_EVALS:-2000}"
)
if [ "${OPTIMIZED_M5_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi

echo "task=$TASK_ID method=$METHOD step=$STEP mix=$MIX outdir=$OUTDIR python=$PYTHON_BIN"
exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.calibration_search "${ARGS[@]}"
