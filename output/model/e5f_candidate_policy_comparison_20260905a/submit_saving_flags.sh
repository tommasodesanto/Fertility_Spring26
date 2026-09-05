#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/code/model:/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/code/model/tools
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 /scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/output/model/e5f_candidate_policy_comparison_20260905a/inspect_saving_flags.py
