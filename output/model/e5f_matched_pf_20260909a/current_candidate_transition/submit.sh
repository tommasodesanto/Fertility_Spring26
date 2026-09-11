#!/bin/bash
#SBATCH --job-name=e5f_current_path
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=03:10:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/path_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/path_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/batches/current_candidate/cache/mpl" NUMBA_CACHE_DIR="$PWD/batches/current_candidate/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python -m unittest test_e5f_candidate_drivers test_e5f_balanced_history test_e5f_balanced_terminal test_e5f_approved_initial_state test_e5f_social_security_root -q
python batches/current_candidate/run_candidate_path.py --delta -0.05 --label delta_m005
