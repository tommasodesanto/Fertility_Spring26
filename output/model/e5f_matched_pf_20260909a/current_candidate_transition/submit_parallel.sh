#!/bin/bash
#SBATCH --job-name=e5f_shock_bracket
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=04:15:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/bracket_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/bracket_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/batches/current_candidate/cache/bracket${SLURM_ARRAY_TASK_ID}/mpl" NUMBA_CACHE_DIR="$PWD/batches/current_candidate/cache/bracket${SLURM_ARRAY_TASK_ID}/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
if [[ "$SLURM_ARRAY_TASK_ID" == 0 ]]; then
 python batches/current_candidate/run_gated_shock_probe.py --delta -0.025 --label delta_m0025
else
 python batches/current_candidate/run_gated_shock_probe.py --delta -0.10 --label delta_m010
fi
