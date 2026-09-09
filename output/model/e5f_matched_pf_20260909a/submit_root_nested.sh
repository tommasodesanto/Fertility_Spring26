#!/bin/bash
#SBATCH --job-name=e5f_terminal_root
#SBATCH --output=logs/terminal_root_%j.out
#SBATCH --error=logs/terminal_root_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_roots_20260909a
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_roots_20260909a/output/cache/mpl_nested
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_roots_20260909a/output/cache/numba_nested
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
PYTHONPATH=code/model/tools:code/model NUMBA_DISABLE_JIT=1 python -m unittest test_e5f_matched_pf_price_root test_e5f_matched_pf_path_root test_e5f_matched_pf_endpoint test_e5f_matched_pf_initial_state test_e5f_matched_pf_moments test_e5f_matched_pf_history test_e5f_matched_pf_smoke
python code/model/tools/run_e5f_matched_pf_baseline.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_roots_20260909a/root_nested_contract.json --contract-sha256 931a3f4c2dfa0583c0c4416ad53bdaa4b6965eb9172d6f93ef8fb661ca09434e --arm nested --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_roots_20260909a/output/root_01/nested
