#!/bin/bash
#SBATCH --job-name=e5f_dead_trace
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_feasibility_20260907
mkdir -p tmp
export TMPDIR="$PWD/tmp" NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 -u diagnose.py
