#!/bin/bash
# Collect strict winners after every continuation chain finishes.
#SBATCH --job-name=ihfopt6c
#SBATCH --output=slurm_ihfopt6_collect_%j.out
#SBATCH --error=slurm_ihfopt6_collect_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

RUN_TAG="${OPTIMIZED_M5_RUN_TAG:-intergen_optimized_m5_continuation_20260719}"
RUN_ROOT="${PROJECT_ROOT}/output/model/${RUN_TAG}"
exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.calibration_collect \
  --run-root "$RUN_ROOT/chains" \
  --output "$RUN_ROOT/report"
