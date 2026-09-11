#!/bin/bash
#SBATCH --job-name=e5f_utility_fiscal
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:12:00
#SBATCH --array=0-3%4
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/output/utility_fiscal_decomposition_20260911a/slurm_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/output/utility_fiscal_decomposition_20260911a/slurm_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508
export PYTHONPATH=code/model/tools:code/model
# Reuse the verified compiled cache; no source or prior result is changed.
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
cases=(old_old old_balanced new_old new_balanced)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
python code/model/tools/run_e5f_initial_revision_probe.py --contract /scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/output/utility_fiscal_decomposition_20260911a/contract.json --contract-sha256 65c7b48a7d26a8c0ca67632acf499aeb53705d2b7436469833e61b6df8ece6a9 --case "$case_name" --output "/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/output/utility_fiscal_decomposition_20260911a/$case_name"
