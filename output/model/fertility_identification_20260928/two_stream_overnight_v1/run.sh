#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=07:10:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fertility_identification_20260928/two_stream_overnight_v1
mode=${1:?prepare, tests, integration or search}
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_CACHE_DIR="$stage/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
python=(apptainer exec --bind "$stage/project:$original" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python)
case "$mode" in
  prepare) exec "${python[@]}" "$packet/prepare.py" ;;
  tests) exec "${python[@]}" "$packet/tests_search.py" ;;
  integration|search)
    : "${EXPECTED_TWO_STREAM_CONFIG_SHA256:?Pinned new config required}"
    index=${SLURM_ARRAY_TASK_ID:?Two independent one-worker tasks required}
    if [ "$index" = 0 ]; then lane=one_birth; elif [ "$index" = 1 ]; then lane=two_birth; else exit 2; fi
    if [ "$mode" = integration ]; then
      exec "${python[@]}" "$packet/integration_smoke.py" --lane "$lane"
    fi
    : "${TWO_STREAM_LAUNCH_APPROVAL:?Passed source-pinned readiness required}"
    exec "${python[@]}" "$packet/search.py" --config "$original/$packet/config.json" --lane "$lane" --output "$original/$packet/run_v1/$lane"
    ;;
  *) exit 2 ;;
esac
