#!/bin/bash
#SBATCH --job-name=e5f_horizon100
#SBATCH --output=logs/horizon100_%j.out
#SBATCH --error=logs/horizon100_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a/output/cache/mpl
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a/output/cache/numba
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
NUMBA_DISABLE_JIT=1 PYTHONPATH=code/model/tools:code/model python -m unittest test_run_e5f_matched_pf_historical_root test_e5f_matched_pf_history test_e5f_matched_pf_path_root test_collect_e5f_matched_pf_price_jacobian
python code/model/tools/run_e5f_matched_pf_baseline.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a/horizon100_sequential_contract.json --contract-sha256 cd1fbfd7264f9354f3f4596fb2144d9f40a8c962eb9b77f0e0cf378ad4224ff3 --arm sequential --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_horizon_20260910a/output/horizon100_anchor_01/sequential
