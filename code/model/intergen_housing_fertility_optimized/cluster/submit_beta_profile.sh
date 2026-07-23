#!/bin/bash
# Historical seven-cell annual-beta profile, two nuisance chains per cell.
# New launches are blocked because the July 22 saving/bequest rows mix timing.
# Values above 0.9995 are diagnostic-only tests of the imposed calibration cap.
#SBATCH --job-name=ihfnewbp
#SBATCH --output=slurm_ihfnewbp_%A_%a.out
#SBATCH --error=slurm_ihfnewbp_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-14%14

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
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 14 ]; then
  echo "conditional beta profile requires array task 1..14" >&2
  exit 2
fi

BETAS=(0.9800 0.9800 0.9900 0.9900 0.9950 0.9950 0.9980 0.9980 0.9995 0.9995 0.9997 0.9997 0.9999 0.9999)
METHODS=(nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern)
STEPS=(0.020 0.030 0.020 0.030 0.020 0.030 0.020 0.030 0.020 0.030 0.020 0.030 0.020 0.030)
MIXES=(0.000 0.015 0.000 0.015 0.000 0.015 0.000 0.015 0.000 0.015 0.000 0.015 0.000 0.015)
INDEX=$((TASK_ID - 1))
BETA="${BETAS[$INDEX]}"
METHOD="${METHODS[$INDEX]}"
STEP="${STEPS[$INDEX]}"
MIX="${MIXES[$INDEX]}"

SOURCE_TAG="${NEW_MOMENT_DIAG_SOURCE_TAG:-intergen_new_moment_weight_tilt_20260722_canonical_cleanup}"
SEED_JSON="${NEW_MOMENT_DIAG_SEED_JSON:-${PROJECT_ROOT}/output/model/${SOURCE_TAG}/report/results.json}"
if [ ! -f "$SEED_JSON" ]; then
  echo "missing certified strict winner: $SEED_JSON" >&2
  exit 3
fi
RUN_TAG="${NEW_MOMENT_BETA_PROFILE_RUN_TAG:-intergen_new_moment_beta_profile_20260723}"
OUTDIR="${PROJECT_ROOT}/output/model/${RUN_TAG}/tasks/task_$(printf '%02d' "$TASK_ID")"
mkdir -p "$NUMBA_CACHE_DIR" "$OUTDIR"

ARGS=(
  --profile new-moments
  --outdir "$OUTDIR"
  --seed "$((2026072700 + TASK_ID))"
  --method "$METHOD"
  --start-theta-json "$SEED_JSON"
  --fixed-beta-annual "$BETA"
  --initial-step "$STEP"
  --start-mix "$MIX"
  --minimum-step "${NEW_MOMENT_BETA_PROFILE_MINIMUM_STEP:-0.00020}"
  --minutes "${NEW_MOMENT_BETA_PROFILE_MINUTES:-110}"
  --strict-reserve-minutes "${NEW_MOMENT_BETA_PROFILE_STRICT_RESERVE_MINUTES:-10}"
  --max-evals "${NEW_MOMENT_BETA_PROFILE_MAX_EVALS:-1200}"
)
if [ "${NEW_MOMENT_BETA_PROFILE_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi

echo "task=$TASK_ID beta_annual=$BETA method=$METHOD step=$STEP mix=$MIX seed=$SEED_JSON outdir=$OUTDIR"
exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.calibration_search "${ARGS[@]}"
