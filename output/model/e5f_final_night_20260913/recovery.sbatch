#!/bin/bash
#SBATCH --account=torch_pr_570_general
#SBATCH --job-name=e5f_initial_recovery
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/recovery_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/recovery_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/output/cache/numba
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_source:/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/initial_v2:/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batches/md_exact_loop/inputs:/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/code/model/tools:/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/code/model
python -B /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913/joint_source/recover_e5f_rebated_initial_search.py --batch /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/final_night_20260913 --seconds 2100 --original-search initial_search_joint
