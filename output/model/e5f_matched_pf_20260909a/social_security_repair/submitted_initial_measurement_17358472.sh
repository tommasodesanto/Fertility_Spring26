#!/bin/bash
#SBATCH --job-name=e5f_measure_fiscal
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_measurement_9e2d0a5b/measurement_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_measurement_9e2d0a5b/measurement_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_measurement_9e2d0a5b
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python code/model/tools/run_e5f_initial_measurement_probe.py --contract contract.json --contract-sha256 3ada0b4ddf7b1c10a9ab4025607378b1411593eff3dd8799118f01ccacb6fcb1 --output output/measurement
python -m unittest test_e5f_parenthood_utility test_e5f_stationary_paygo test_e5f_social_security_compiled -v
