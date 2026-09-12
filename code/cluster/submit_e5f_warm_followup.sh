#!/bin/bash
#SBATCH --job-name=e5f_warm_forecast
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:40:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_warm_followup_20260912/slurm_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_warm_followup_20260912/slurm_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/output/cache/numba
python -B /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_warm_followup_20260912/run_e5f_warm_followup.py --arm "$SLURM_ARRAY_TASK_ID"
