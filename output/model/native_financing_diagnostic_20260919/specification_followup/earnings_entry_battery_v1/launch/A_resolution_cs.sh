#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
test "$(sha256sum /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/earnings_entry_resolution_20260922_v1/A_plan.json | cut -d " " -f 1)" = 99cb7b80a1f49f80e5f7661b172891dd9e0e25fde717888cb9832deccac6494a
mkdir -p /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/earnings_entry_resolution_20260922_v1/results/smoke/A/worker_01/proposal_01_resolution_smoke
exec /share/apps/anaconda3/2025.06/bin/python /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/earnings_entry_battery_20260922_v1/tools/run_e5f_earnings_wealth_candidate.py --plan /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/earnings_entry_resolution_20260922_v1/A_plan.json --output /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/earnings_entry_resolution_20260922_v1/results/smoke/A/worker_01/proposal_01_resolution_smoke/result --arm literature_income_purchase --repetitions 1
