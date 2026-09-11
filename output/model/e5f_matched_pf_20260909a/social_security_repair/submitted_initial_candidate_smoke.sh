#!/bin/bash
#SBATCH --job-name=e5f_initial_repeat
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:32:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/initial_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/initial_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python -m unittest test_e5f_parenthood_utility test_e5f_stationary_paygo -v
python code/model/tools/run_e5f_initial_revision_probe.py --contract initial_contract.json --contract-sha256 ba25eab4449fd8fba4b4d4bbd0e5ea3b5fffc70b29649572feb54d80a34706d4 --case new_balanced --output output/new_balanced_smoke
