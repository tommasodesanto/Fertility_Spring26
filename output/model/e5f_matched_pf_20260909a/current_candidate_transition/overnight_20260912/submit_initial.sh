#!/bin/bash
#SBATCH --job-name=e5f_night_initial
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=18
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/night_20260912/initial_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/night_20260912/initial_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8
export PYTHONPATH="$PWD/batches/md_exact_loop/inputs:$PWD/code/model/tools:$PWD/code/model"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba"
python -B batches/night_20260912/run_e5f_overnight_initial_refinement.py --template batches/capped_beta_099_20260911 --output "$PWD/batches/night_20260912/search"
