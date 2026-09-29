#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:10:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fertility_identification_20260928/zero_fertility_taste_v1
mode=${1:?prepare, tests or search}
lane=${2:-}
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_CACHE_DIR="$stage/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
python=(apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python)
case "$mode" in
  prepare) exec "${python[@]}" "$packet/prepare.py" ;;
  tests) exec "${python[@]}" "$packet/tests_search.py" ;;
  search)
    : "${EXPECTED_ZERO_TASTE_CONFIG_SHA256:?Pinned configuration required}"
    [[ "$lane" == first_scale_zero || "$lane" == later_scale_zero ]] || { echo "Select one zero-scale lane" >&2; exit 2; }
    exec "${python[@]}" "$packet/search.py" --config "$original/$packet/config.json" --lane "$lane" --output "$original/$packet/run_v1/$lane"
    ;;
  *) exit 2 ;;
esac
