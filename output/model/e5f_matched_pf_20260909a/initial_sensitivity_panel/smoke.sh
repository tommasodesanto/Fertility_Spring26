#!/bin/bash
#SBATCH --job-name=e5f_candidate_smoke
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:32:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/smoke_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/smoke_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python code/model/tools/run_e5f_initial_revision_probe.py --contract contracts/smoke_contract.json --contract-sha256 00ef31f4a34af5bc118ef2b6719611ccbfb81a03e9d1332d7730aa9633e0a198 --case new_balanced --output output/candidate_smoke
