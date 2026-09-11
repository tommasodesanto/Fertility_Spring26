#!/bin/bash
#SBATCH --job-name=e5f_psi_rates
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=03:15:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/psi_rates_%A_%a.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/psi_rates_%A_%a.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a
export PYTHONPATH=code/model/tools:code/model:batches/current_candidate
export MPLCONFIGDIR="$PWD/batches/current_candidate/cache/rates${SLURM_ARRAY_TASK_ID}/mpl" NUMBA_CACHE_DIR="$PWD/batches/current_candidate/cache/rates${SLURM_ARRAY_TASK_ID}/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python -m unittest test_e5f_candidate_drivers test_e5f_initial_fertility_observer -q
python -c 'import sys, unittest, run_observed_bracket_v2 as v2; sys.modules["run_observed_bracket"]=v2; result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromName("test_observed_bracket_v2")); sys.exit(not result.wasSuccessful())'
python -c 'from test_e5f_transition_accounting import test_period_tfr_reconstructs_parity_flows; test_period_tfr_reconstructs_parity_flows(); print("Period fertility flow accounting passed")'
python batches/current_candidate/run_observed_bracket_v2.py --case "$SLURM_ARRAY_TASK_ID"
