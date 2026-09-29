#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:12:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fertility_identification_20260928/numerical_pair_v1
mode=${1:?Use prepare, smoke, or run}
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_CACHE_DIR="$stage/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
if [ "$mode" = prepare ]; then
  args=("$packet/prepare.py")
else
  : "${EXPECTED_NUMERICAL_PAIR_SHA256:?Explicit reviewed pair contract pin required}"
  args=("$packet/runner.py" --mode "$mode" --contract "$packet/contract.json")
  if [ "$mode" = run ]; then
    : "${NUMERICAL_PAIR_APPROVAL:?Explicit reviewed approval JSON required}"
    : "${EXPECTED_NUMERICAL_PAIR_APPROVAL_SHA256:?Explicit approval pin required}"
    args+=(--approval "$NUMERICAL_PAIR_APPROVAL")
  fi
fi
exec apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "${args[@]}"
