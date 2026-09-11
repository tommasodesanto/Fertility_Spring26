#!/bin/bash
#SBATCH --job-name=e5f_path_continue
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:40:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/resume_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/resume_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/batches/current_candidate/cache/resume/mpl" NUMBA_CACHE_DIR="$PWD/batches/current_candidate/cache/resume/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python -m unittest test_e5f_candidate_drivers -q
python batches/current_candidate/resume_history.py
