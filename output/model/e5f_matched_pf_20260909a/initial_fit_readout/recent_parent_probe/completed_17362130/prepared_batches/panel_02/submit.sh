#!/bin/bash
# Prepared only; lead review and staging precede submission.
#SBATCH --job-name=e5f_recent_panel_02
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batch_panel_02_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batch_panel_02_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8
export PYTHONPATH="$PWD/code/model/tools:$PWD/code/model"
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
timeout --signal=KILL 240s python batches/panel_02/run.py --contract batches/panel_02/contract.json --contract-sha256 8691abb304897f4c23aebcef63858f8754a6227c89decbc314c0bbb5ef1e5d16 --source-root "$PWD" --output "output/batch_panel_02_${SLURM_JOB_ID}"
