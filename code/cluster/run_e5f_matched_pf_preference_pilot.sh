#!/bin/bash
#SBATCH --job-name=e5f_psi_pilot
#SBATCH --output=logs/psi_%A_%a.out
#SBATCH --error=logs/psi_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:05:00
#SBATCH --array=0-2%3
set -euo pipefail
: "${E5F_PILOT_ROOT:?Explicit isolated snapshot required}"
: "${E5F_PILOT_PHASE:?Explicit smoke or main required}"
: "${SLURM_ARRAY_TASK_ID:?Array case required}"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "$E5F_PILOT_ROOT"
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl_${E5F_PILOT_PHASE}_${SLURM_ARRAY_TASK_ID}"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba_${E5F_PILOT_PHASE}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
NUMBA_DISABLE_JIT=1 python -m unittest test_run_e5f_matched_pf_preference_pilot test_e5f_matched_pf_birth_path test_e5f_matched_pf_path_root test_run_e5f_matched_pf_historical_root
exec python code/model/tools/prepare_e5f_matched_pf_preference_case.py --phase "$E5F_PILOT_PHASE" --case "$SLURM_ARRAY_TASK_ID"
