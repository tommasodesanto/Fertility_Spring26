#!/bin/bash
# Collect the corrected recent-vintage gross/gross conditional-beta profile.
#SBATCH --job-name=ihfnewbpc
#SBATCH --output=slurm_ihfnewbp_collect_%j.out
#SBATCH --error=slurm_ihfnewbp_collect_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:10:00
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

RUN_TAG="${NEW_MOMENT_BETA_PROFILE_RUN_TAG:-intergen_new_moment_beta_recent_gross_20260723}"
RUN_ROOT="${PROJECT_ROOT}/output/model/${RUN_TAG}"
mkdir -p "${RUN_ROOT}/report"
exec "$PYTHON_BIN" -m intergen_housing_fertility_optimized.collect_beta_profile \
  --run-root "${RUN_ROOT}/tasks" \
  --output "${RUN_ROOT}/report" \
  --expected-betas 0.98 0.99 0.995 0.999 0.9995
