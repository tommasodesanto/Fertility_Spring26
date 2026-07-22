#!/bin/bash
# One eight-chain overnight continuation wave for the corrected one-shot
# 14-moment system. The start point must be a collected, twice-repeated strict
# winner from the preceding wave. Each chain checkpoints every evaluation and
# reserves ten minutes for two fresh strict winner repeats. Submit two waves
# with a collector between them to respect Torch's 3:55 cpu_short cap.
#SBATCH --job-name=ihfnewov
#SBATCH --output=slurm_ihfnewov_%A_%a.out
#SBATCH --error=slurm_ihfnewov_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
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
  echo "overnight continuation requires array task 1..8" >&2
  exit 2
fi

# Fine and medium transformed-coordinate neighborhoods around the collected
# winner. Alternating Nelder-Mead and pattern chains provide independent local
# search geometries rather than eight replicas of one optimizer.
METHODS=(nelder-mead pattern nelder-mead pattern nelder-mead pattern nelder-mead pattern)
STEPS=(0.0015 0.0030 0.0050 0.0080 0.0120 0.0200 0.0350 0.0600)
MIXES=(0.0000 0.0000 0.0010 0.0025 0.0050 0.0100 0.0200 0.0400)
INDEX=$((TASK_ID - 1))
METHOD="${METHODS[$INDEX]}"
STEP="${STEPS[$INDEX]}"
MIX="${MIXES[$INDEX]}"

SOURCE_TAG="${NEW_MOMENT_OVERNIGHT_SOURCE_TAG:-intergen_new_moment_calibration_20260722_corrected_3h}"
SEED_JSON="${NEW_MOMENT_OVERNIGHT_SEED_JSON:-${PROJECT_ROOT}/output/model/${SOURCE_TAG}/report/results.json}"
if [ ! -f "$SEED_JSON" ]; then
  echo "missing collected strict-winner seed: $SEED_JSON" >&2
  exit 3
fi
RUN_TAG="${NEW_MOMENT_OVERNIGHT_RUN_TAG:-intergen_new_moment_overnight_20260722}"
RUN_ROOT="${PROJECT_ROOT}/output/model/${RUN_TAG}"
OUTDIR="${RUN_ROOT}/chains/chain_${TASK_ID}"
mkdir -p "$NUMBA_CACHE_DIR" "$OUTDIR"

ARGS=(
  --profile new-moments
  --outdir "$OUTDIR"
  --seed "$((${NEW_MOMENT_OVERNIGHT_SEED_BASE:-2026072300} + TASK_ID))"
  --method "$METHOD"
  --start-theta-json "$SEED_JSON"
  --initial-step "$STEP"
  --start-mix "$MIX"
  --minimum-step "${NEW_MOMENT_OVERNIGHT_MINIMUM_STEP:-0.00010}"
  --minutes "${NEW_MOMENT_OVERNIGHT_MINUTES:-225}"
  --strict-reserve-minutes "${NEW_MOMENT_OVERNIGHT_STRICT_RESERVE_MINUTES:-10}"
  --max-evals "${NEW_MOMENT_OVERNIGHT_MAX_EVALS:-4000}"
)
if [ "${NEW_MOMENT_OVERNIGHT_SMOKE:-0}" = "1" ]; then
  ARGS+=(--smoke)
fi

echo "task=$TASK_ID method=$METHOD step=$STEP mix=$MIX seed=$SEED_JSON outdir=$OUTDIR python=$PYTHON_BIN"
exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.calibration_search "${ARGS[@]}"
