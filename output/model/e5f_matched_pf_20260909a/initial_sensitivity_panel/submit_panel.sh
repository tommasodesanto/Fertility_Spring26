#!/bin/bash
#SBATCH --job-name=e5f_initial_panel
#SBATCH --partition=cpu_short
#SBATCH --array=0-18%19
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:32:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/panel_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/panel_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python contracts/run_array_case.py
