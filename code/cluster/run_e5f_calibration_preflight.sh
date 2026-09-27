#!/usr/bin/env bash
# Configuration and controller checks only: zero household solves, no search.
set -euo pipefail
: "${SLURM_JOB_ID:?Run in a Torch Slurm allocation}"
work="${1:?staged working directory}"
contract_dir="${2:?new contract output directory}"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$work/tools:$work/source/code/model/tools:$work/source/code/model"
export REVIEWED_RECOVERY_SEARCH="$work/tools/run_e5f_utility_comparison_search.py"
python -m unittest -v test_e5f_calibration_runtime test_e5f_utility_overnight_controller
python "$work/tools/prepare_e5f_calibration_launch.py" --work-root "$work" \
  --output "$contract_dir" --ahs-target "$work/inputs/ahs_2007_room_target.json" \
  --initial-parameters "$work/inputs/parameters.csv"
export EXPECTED_UTILITY_OVERNIGHT_SHA256="$(sha256sum "$contract_dir/contract.json" | cut -d' ' -f1)"
python "$work/tools/run_e5f_utility_overnight_calibration.py" --stage preflight \
  --contract "$contract_dir/contract.json" --output "$contract_dir/preflight"
