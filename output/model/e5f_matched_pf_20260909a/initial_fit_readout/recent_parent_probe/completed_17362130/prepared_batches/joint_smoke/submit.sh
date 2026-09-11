#!/bin/bash
# Prepared only; lead review and staging precede submission.
#SBATCH --job-name=e5f_recent_joint_smoke
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batch_joint_smoke_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8/batch_joint_smoke_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_recent_parent_probe_70abd4a8
export PYTHONPATH="$PWD/code/model/tools:$PWD/code/model"
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
timeout --signal=KILL 240s python batches/joint_smoke/run.py --contract batches/joint_smoke/contract.json --contract-sha256 a710834c8f9a74d5aed5689540c9569cdb12bf6d5efa002b17591c2bf93cf829 --source-root "$PWD" --output "output/batch_joint_smoke_${SLURM_JOB_ID}"
