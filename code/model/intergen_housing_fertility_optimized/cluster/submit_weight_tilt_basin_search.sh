#!/bin/bash
# Diagnostic basin search for the corrected one-shot 14-moment system.
# Search navigation temporarily tilts the wealth and/or TFR rows, but every
# candidate retains its canonical loss and the two strict repeats are scored
# only on the unchanged canonical objective.
#SBATCH --job-name=ihfnewtilt
#SBATCH --output=slurm_ihfnewtilt_%A_%a.out
#SBATCH --error=slurm_ihfnewtilt_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:30:00
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
  echo "weight-tilt basin search requires array task 1..8" >&2
  exit 2
fi

# Four predeclared search tilts, each run with both algorithms. Broader start
# mixing than the local continuation is deliberate: this is a basin probe.
METHODS=(nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern)
STEPS=(0.030 0.030 0.040 0.040 0.060 0.060 0.080 0.080)
MIXES=(0.020 0.020 0.040 0.040 0.080 0.080 0.120 0.120)
WEALTH_FACTORS=(4 4 1 1 4 4 9 9)
TFR_FACTORS=(1 1 4 4 4 4 9 9)
INDEX=$((TASK_ID - 1))
METHOD="${METHODS[$INDEX]}"
STEP="${STEPS[$INDEX]}"
MIX="${MIXES[$INDEX]}"
WEALTH_FACTOR="${WEALTH_FACTORS[$INDEX]}"
TFR_FACTOR="${TFR_FACTORS[$INDEX]}"

SOURCE_TAG="${NEW_MOMENT_TILT_SOURCE_TAG:-intergen_new_moment_overnight_20260722_wave2}"
SEED_JSON="${NEW_MOMENT_TILT_SEED_JSON:-${PROJECT_ROOT}/output/model/${SOURCE_TAG}/report/results.json}"
if [ ! -f "$SEED_JSON" ]; then
  echo "missing collected canonical seed: $SEED_JSON" >&2
  exit 3
fi
RUN_TAG="${NEW_MOMENT_TILT_RUN_TAG:-intergen_new_moment_weight_tilt_20260722}"
RUN_ROOT="${PROJECT_ROOT}/output/model/${RUN_TAG}"
OUTDIR="${RUN_ROOT}/chains/chain_${TASK_ID}"
mkdir -p "$NUMBA_CACHE_DIR" "$OUTDIR"

ARGS=(
  --profile new-moments
  --outdir "$OUTDIR"
  --seed "$((${NEW_MOMENT_TILT_SEED_BASE:-2026072500} + TASK_ID))"
  --method "$METHOD"
  --start-theta-json "$SEED_JSON"
  --initial-step "$STEP"
  --start-mix "$MIX"
  --minimum-step "${NEW_MOMENT_TILT_MINIMUM_STEP:-0.00020}"
  --minutes "${NEW_MOMENT_TILT_MINUTES:-140}"
  --strict-reserve-minutes "${NEW_MOMENT_TILT_STRICT_RESERVE_MINUTES:-10}"
  --max-evals "${NEW_MOMENT_TILT_MAX_EVALS:-1800}"
)
if [ "$WEALTH_FACTOR" != "1" ]; then
  ARGS+=(
    --search-weight-multiplier
    "aggregate_wealth_to_annual_after_tax_earnings=${WEALTH_FACTOR}"
  )
fi
if [ "$TFR_FACTOR" != "1" ]; then
  ARGS+=(--search-weight-multiplier "tfr=${TFR_FACTOR}")
fi
if [ "${NEW_MOMENT_TILT_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi

echo "task=$TASK_ID method=$METHOD step=$STEP mix=$MIX wealth_factor=$WEALTH_FACTOR tfr_factor=$TFR_FACTOR seed=$SEED_JSON outdir=$OUTDIR"
exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.calibration_search "${ARGS[@]}"
