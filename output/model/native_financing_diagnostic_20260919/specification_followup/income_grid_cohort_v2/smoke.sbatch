#!/usr/bin/env bash
#SBATCH --job-name=grid_smoke
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --chdir=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/logs/%x_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/logs/%x_%j.err
#SBATCH --kill-on-invalid-dep=yes
set -euo pipefail
module load anaconda3/2025.06
unset PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/numba_cache MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
cd /scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2
timeout --signal=TERM --kill-after=30s 1050s python -B code/model/tools/run_e5f_income_grid_cohort.py --mode smoke --checkpoint /scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/checkpoint.pkl.gz --source-root /scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/code/model --candidate-json /scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/candidate.json --baseline-summary /scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/baseline_summary.json --baseline-summary-sha256 a22c7cc58d0a937afa7414e2a90577d3fe86b922ddc601a6c71abb141008c37b --output /scratch/td2248/projects/Fertility_Spring26_specification_20260920/income_grid_cohort_v2/results/smoke
