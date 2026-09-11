#!/bin/bash
#SBATCH --job-name=e5f_joint_round01
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:32:00
#SBATCH --array=0-20,22-23%24
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/output/joint_round_01/array_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/output/joint_round_01/array_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
python /scratch/td2248/projects/Fertility_Spring26_initial_panel_7e872053/contracts/joint_round_01/launch_gate.py run-array --plan-sha256 "${ROUND_PLAN_SHA256:?Set the reviewed run-plan SHA256}"
