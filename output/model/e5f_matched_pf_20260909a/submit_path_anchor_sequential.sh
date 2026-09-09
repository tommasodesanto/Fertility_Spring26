#!/bin/bash
#SBATCH --job-name=e5f_path_anchor
#SBATCH --output=logs/path_anchor_%j.out
#SBATCH --error=logs/path_anchor_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/cache/mpl_sequential
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/cache/numba_sequential
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
PYTHONPATH=code/model/tools:code/model NUMBA_DISABLE_JIT=1 python -m unittest test_e5f_matched_pf_price_root test_e5f_matched_pf_path_root test_e5f_matched_pf_endpoint test_e5f_matched_pf_initial_state test_e5f_matched_pf_moments test_e5f_matched_pf_history test_e5f_matched_pf_smoke
python code/model/tools/run_e5f_matched_pf_baseline.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/path_anchor_sequential_contract.json --contract-sha256 0cd547e37b3ca9d4d2e2ae91f5e27ae22fafa50efd74142d5794d63aa82a1dc9 --arm sequential --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/path_anchor_02/sequential
