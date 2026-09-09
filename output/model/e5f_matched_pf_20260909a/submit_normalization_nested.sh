#!/bin/bash
#SBATCH --job-name=e5f_old_init
#SBATCH --output=logs/old_init_%j.out
#SBATCH --error=logs/old_init_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_normalized_20260909b
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_normalized_20260909b/output/cache/mpl_nested
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_normalized_20260909b/output/cache/numba_nested
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
PYTHONPATH=code/model/tools:code/model NUMBA_DISABLE_JIT=1 python -m unittest test_e5f_matched_pf_endpoint test_e5f_matched_pf_initial_state test_e5f_matched_pf_moments test_e5f_matched_pf_history test_e5f_matched_pf_smoke
python code/model/tools/run_e5f_matched_pf_baseline.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_normalized_20260909b/normalization_nested_contract.json --contract-sha256 27778a1811c560f3c54b566ba439f4ede57d1362e3a63a15a63ac358cccfcee8 --arm nested --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_normalized_20260909b/output/initial_01/nested
