#!/bin/bash
#SBATCH --job-name=e5f_joined_pf
#SBATCH --output=logs/joined_%j.out
#SBATCH --error=logs/joined_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_history_20260909a
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_history_20260909a/output/cache/mpl
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_history_20260909a/output/cache/numba
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python code/model/tools/test_e5f_matched_pf_history.py
python code/model/tools/run_e5f_matched_pf_history.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_history_20260909a/history_contract.json --contract-sha256 ab9058de1bd67b442c8ee274cd1f74648a0d7a7927b3d4db75c21741c73d0999 --arm sequential --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_history_20260909a/output/history_01/sequential
