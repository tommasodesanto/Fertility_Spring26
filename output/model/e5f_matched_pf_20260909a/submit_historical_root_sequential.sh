#!/bin/bash
#SBATCH --job-name=e5f_historical_root
#SBATCH --output=logs/historical_root_%j.out
#SBATCH --error=logs/historical_root_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a/output/cache/mpl_root_sequential
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a/output/cache/numba_root_sequential
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python code/model/tools/run_e5f_matched_pf_historical_root.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a/historical_root_sequential_contract.json --contract-sha256 68dd56ef25ab0c57b04c1503ced07fb3fef9b8e50657342da360dbeb2a340570 --arm sequential --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909a/output/historical_root_01/sequential
