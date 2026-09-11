#!/bin/bash
# Prepared only. Submit explicitly after lead review and archive staging.
#SBATCH --job-name=e5f_recent_parent_probe
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/probe_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/probe_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8
export PYTHONPATH="$PWD/code/model/tools:$PWD/code/model"
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
timeout --signal=KILL 240s python probe/run.py \
  --contract probe/contract.json \
  --contract-sha256 48a18bb61d5c83006a76c3517aaafc666e203e65965008c2763a6484c2ab9e81 \
  --source-root "$PWD" --output "output/probe_${SLURM_JOB_ID}"
