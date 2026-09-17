#!/bin/bash
#SBATCH --job-name=e5f_beta099_cap
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=18
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/capped_beta_099_20260911/capped_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/capped_beta_099_20260911/capped_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONOPTIMIZE=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8
export PYTHONPATH="$PWD/batches/md_exact_loop/inputs:$PWD/code/model/tools:$PWD/code/model"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba"
python -B -m unittest discover -s batches/capped_beta_099_20260911 -p 'test_run_capped_beta.py' -q
python -B batches/capped_beta_099_20260911/run_capped_beta.py --plan batches/capped_beta_099_20260911/plan_capped_beta_099.json
