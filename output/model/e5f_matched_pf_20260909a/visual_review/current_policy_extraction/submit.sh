#!/bin/bash
#SBATCH --job-name=e5f_policy_read
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_history_root_450bce1c/batches/current_policy_extraction/read_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_history_root_450bce1c/batches/current_policy_extraction/read_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONUNBUFFERED=1 NUMBA_DISABLE_JIT=1
python /scratch/td2248/projects/Fertility_Spring26_history_root_450bce1c/batches/current_policy_extraction/extract.py --output /scratch/td2248/projects/Fertility_Spring26_history_root_450bce1c/batches/current_policy_extraction/results
