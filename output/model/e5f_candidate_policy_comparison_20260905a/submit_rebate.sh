#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/code/model:/scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/code/model/tools
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 /scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/code/model/tools/observe_e5f_candidate_rebate.py --plan /scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a/output/model/e5f_candidate_policy_comparison_20260905a/rebate_plan.json --plan-sha256 05b28547b2a442c4625d264e057e1411faeaf68b540d48248ca7d17439d03101
