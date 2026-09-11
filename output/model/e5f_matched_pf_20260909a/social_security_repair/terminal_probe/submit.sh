#!/bin/bash
#SBATCH --job-name=e5f_terminal_paygo
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:32:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_terminal_paygo_b373b142/terminal_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_terminal_paygo_b373b142/terminal_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_terminal_paygo_b373b142
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl" NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python -m unittest test_e5f_approved_initial_state test_e5f_balanced_terminal test_e5f_social_security_root test_e5f_matched_pf_endpoint -q
python code/model/tools/run_e5f_balanced_terminal_probe.py --contract contracts/contract.json --contract-sha256 f6baa5832ea12d82bd9564bf43faea7b52ace370db82c31c231bab35a61663e1 --output output/terminal_smoke
