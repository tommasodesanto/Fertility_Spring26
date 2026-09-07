#!/bin/bash
#SBATCH --job-name=e5f_mass_repair
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907h
export TMPDIR="$PWD/tmp"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$TMPDIR/mpl"
timeout --signal=TERM --kill-after=20s 3600s python3 -u verify_and_replay_case18.py
