#!/bin/bash
#SBATCH --job-name=e5f_endpoint_smoke
#SBATCH --output=logs/endpoint_%j.out
#SBATCH --error=logs/endpoint_%j.err
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_matched_pf_baseline_20260909a
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_baseline_20260909a/output/cache/mpl_nested
export NUMBA_CACHE_DIR=/scratch/td2248/projects/Fertility_Spring26_matched_pf_baseline_20260909a/output/cache/numba_nested
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python code/model/tools/run_e5f_matched_pf_baseline.py --contract /scratch/td2248/projects/Fertility_Spring26_matched_pf_baseline_20260909a/endpoint_nested_contract.json --contract-sha256 3d94c0c8790c7643e69730963beb15e42012666d076c54c9eea4007813023628 --arm nested --output /scratch/td2248/projects/Fertility_Spring26_matched_pf_baseline_20260909a/output/endpoint_01/nested
