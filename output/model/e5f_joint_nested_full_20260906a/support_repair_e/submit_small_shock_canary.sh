#!/bin/bash
#SBATCH --job-name=e5fnest_corner
#SBATCH --partition=cpu_short,cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/joint_support_${SLURM_JOB_ID}"
cd "${SLURM_SUBMIT_DIR:?}"
timeout --signal=TERM --kill-after=30 3600 python3 -u run_small_shock_canary.py
