#!/bin/bash
#SBATCH --job-name=e5f_tail_sizing
#SBATCH --output=logs/tail_sizing_%j.out
#SBATCH --error=logs/tail_sizing_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:08:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/mpl_tail_sizing
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/cache/numba_tail_sizing
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python forecast_cached_terminal_tail_v2.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/historical_root_long_sequential_contract.json --contract-sha256 442aa4823d10173b970759d4763d54ae56c28d8c0ae56fcf1568f262fd6e068d --anchor-summary /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/path_long_anchor_01/sequential/summary.json --anchor-summary-sha256 f82c148fc54f4ef92ad6dfdb660eceeda5c3a6e936bf61f62dc546478c53cc49 --periods 2 --seconds 60 --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/tail_sizing_02/stationary_smoke --initial-state stationary
python forecast_cached_terminal_tail_v2.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/historical_root_long_sequential_contract.json --contract-sha256 442aa4823d10173b970759d4763d54ae56c28d8c0ae56fcf1568f262fd6e068d --anchor-summary /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/path_long_anchor_01/sequential/summary.json --anchor-summary-sha256 f82c148fc54f4ef92ad6dfdb660eceeda5c3a6e936bf61f62dc546478c53cc49 --periods 2 --seconds 60 --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/tail_sizing_02/smoke
python forecast_cached_terminal_tail_v2.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/historical_root_long_sequential_contract.json --contract-sha256 442aa4823d10173b970759d4763d54ae56c28d8c0ae56fcf1568f262fd6e068d --anchor-summary /scratch/td2248/projects/Fertility_Spring26_matched_pf_paths_20260909b/output/path_long_anchor_01/sequential/summary.json --anchor-summary-sha256 f82c148fc54f4ef92ad6dfdb660eceeda5c3a6e936bf61f62dc546478c53cc49 --periods 100 --seconds 240 --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_path_roots_20260909b/output/tail_sizing_02/full
