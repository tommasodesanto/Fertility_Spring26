#!/usr/bin/env bash
#SBATCH --job-name=e5f_saved_render
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:10:00
set -euo pipefail
if [ "$#" -ne 6 ]; then
  echo 'usage: STAGE CONTRACT SMOKE_COMPLETE SMOKE_SHA OUTPUT HELPER_SHA' >&2
  exit 2
fi
stage=$1
contract=$2
smoke=$3
smoke_sha=$4
output=$5
helper_sha=$6
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
helper=output/model/overnight_calibration_20260928/cluster/check_saved_render.py
printf '%s  %s\n' "$helper_sha" "$stage/project/$helper" | sha256sum -c -
module load anaconda3/2025.06
unset NUMBA_DISABLE_JIT APPTAINERENV_NUMBA_DISABLE_JIT SINGULARITYENV_NUMBA_DISABLE_JIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export NUMBA_CACHE_DIR="$stage/numba_cache"
exec apptainer exec --bind "$stage/project:$original" --pwd "$original" \
 /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python \
 "$original/$helper" --contract "$contract" --smoke "$smoke" \
 --smoke-sha256 "$smoke_sha" --output "$output"
