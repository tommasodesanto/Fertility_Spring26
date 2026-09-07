#!/bin/bash
#SBATCH --job-name=e5f_overlap_search
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=352G
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907n
export TMPDIR="$PWD/tmp" NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$TMPDIR/mpl"
python3 -u verify_and_prepare.py --run-search-after-verification
