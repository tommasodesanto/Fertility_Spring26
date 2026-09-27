#!/usr/bin/env bash
# Run inside a Torch Slurm allocation; source must be a separately staged copy.
set -euo pipefail
: "${SLURM_JOB_ID:?Run code checks in a Torch Slurm allocation}"
source_root="${1:?staged source root}"
output_root="${2:?checks output directory}"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$source_root/code/model/tools:$source_root/code/model${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$output_root"
cd "$source_root/code/model/tools"
python -m unittest -v \
  test_e5f_warm_price \
  test_e5f_native_child_preferences \
  test_e5f_normalization_warm_start \
  test_e5f_early_fertility_observer \
  test_e5f_ahs_rooms_observer \
  test_e5f_initial_fertility_observer \
  test_e5f_initial_housing_observer \
  test_e5f_stationary_paygo > "$output_root/unittest.log" 2>&1
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$output_root/completed.txt"
tail -n 5 "$output_root/unittest.log"
