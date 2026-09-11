#!/bin/bash
#SBATCH --job-name=e5f_history_paygo
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:32:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_history_probe_54c52386/history_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_history_probe_54c52386/history_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_history_probe_54c52386
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python -m unittest test_run_e5f_balanced_history_probe test_e5f_balanced_history test_e5f_balanced_terminal test_e5f_approved_initial_state test_e5f_social_security_root test_e5f_matched_pf_endpoint -q
python code/model/tools/run_e5f_balanced_history_probe.py --contract contracts/contract.json --contract-sha256 b8e63815abe002867a0d3a607f82a40791c4510647dcb17884ca92e1433d8882 --output output/history_smoke
