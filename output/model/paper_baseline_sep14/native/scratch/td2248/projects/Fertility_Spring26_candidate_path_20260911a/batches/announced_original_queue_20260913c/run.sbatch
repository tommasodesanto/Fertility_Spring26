#!/bin/bash
#SBATCH --job-name=e5f_announced_run
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=370
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/announced_original_queue_20260913c/logs/run_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
/share/apps/anaconda3/2025.06/bin/python /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/announced_original_queue_20260913c/source/run_e5f_announced_original_queue.py --manifest /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches/announced_original_queue_20260913c/manifest.json --mode run
