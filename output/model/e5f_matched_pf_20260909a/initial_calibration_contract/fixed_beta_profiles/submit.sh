#!/bin/bash
#SBATCH --job-name=e5f_beta_profiles
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/fixed_beta_profiles_20260911/profile_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/fixed_beta_profiles_20260911/profile_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONOPTIMIZE=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8
export PYTHONPATH="$PWD/batches/md_exact_loop/inputs:$PWD/code/model/tools:$PWD/code/model"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba"
python -B -m unittest discover -s batches/fixed_beta_profiles_20260911 -p 'test_run_profile.py' -q
python -B -m unittest discover -s batches/md_exact_loop/inputs -p 'test_*.py' -q
if [[ "$SLURM_ARRAY_TASK_ID" == 0 ]]; then
 python -B batches/fixed_beta_profiles_20260911/run_profile.py --plan batches/fixed_beta_profiles_20260911/plan_beta_098.json
else
 python -B batches/fixed_beta_profiles_20260911/run_profile.py --plan batches/fixed_beta_profiles_20260911/plan_beta_099.json
fi
