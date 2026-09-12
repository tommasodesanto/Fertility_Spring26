#!/bin/bash
#SBATCH --job-name=e5f_surprise_fit
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-2
#SBATCH --cpus-per-task=3
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_surprises_20260912/slurm_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/night_surprises_20260912/slurm_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
cd /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a
batch="$PWD/batches/night_surprises_20260912"
export PYTHONPATH="$batch:$PWD/code/model/tools:$PWD/code/model"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba"
export MPLCONFIGDIR="$batch/mpl_${SLURM_ARRAY_TASK_ID}"
python -B -m unittest test_e5f_surprise_overnight test_e5f_successive_surprises -q
python -B "$batch/run_e5f_successive_surprises_overnight.py" --plan "$batch/plan.json" --arm "$SLURM_ARRAY_TASK_ID"
