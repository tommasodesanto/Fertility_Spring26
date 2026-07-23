#!/bin/bash
# Boundary-aware tight Jacobian at the certified corrected one-shot winner.
#SBATCH --job-name=ihfnewjac
#SBATCH --output=slurm_ihfnewjac_%j.out
#SBATCH --error=slurm_ihfnewjac_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}/intergen_housing_fertility_optimized}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

SOURCE_TAG="${NEW_MOMENT_DIAG_SOURCE_TAG:-intergen_new_moment_weight_tilt_20260722_canonical_cleanup}"
SEED_JSON="${NEW_MOMENT_DIAG_SEED_JSON:-${PROJECT_ROOT}/output/model/${SOURCE_TAG}/report/results.json}"
if [ ! -f "$SEED_JSON" ]; then
  echo "missing certified strict winner: $SEED_JSON" >&2
  exit 3
fi
RUN_TAG="${NEW_MOMENT_JACOBIAN_RUN_TAG:-intergen_new_moment_final_jacobian_20260723}"
OUTDIR="${PROJECT_ROOT}/output/model/${RUN_TAG}/full"
ARGS=(
  --outdir "$OUTDIR"
  --start-theta-json "$SEED_JSON"
  --unit-step "${NEW_MOMENT_JACOBIAN_UNIT_STEP:-0.010}"
  --rank-relative-tol 1e-6
  --tight
)
if [ "${NEW_MOMENT_JACOBIAN_SMOKE:-0}" = "1" ]; then
  ARGS+=(--max-parameters 2)
fi

exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.audit_new_moment_jacobian "${ARGS[@]}"
