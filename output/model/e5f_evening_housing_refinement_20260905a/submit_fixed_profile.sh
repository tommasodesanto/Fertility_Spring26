#!/bin/bash
# Explicit fixed-jump diagnostic; retained calibration is unchanged.
#SBATCH --job-name=e5ffloorprofile
#SBATCH --output=logs/slurm_e5fmorning_%A_%a.out
#SBATCH --error=logs/slurm_e5fmorning_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-2%2
set -euo pipefail
CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:?}" && pwd)"
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
exec python3 "$MODEL_DIR/../../output/model/e5f_evening_housing_refinement_20260905a/fixed_jump_adapter.py" \
  --plan "${E5F_BOUNDED_PLAN:?}" --plan-sha256 "${E5F_BOUNDED_PLAN_SHA256:?}" \
  --case-id "${SLURM_ARRAY_TASK_ID:?}"
