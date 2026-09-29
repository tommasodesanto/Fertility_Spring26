#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
set -euo pipefail

reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
stage=/scratch/td2248/projects/fixed_reference_credit_ge_20260929
source_dir="$stage/source_preflight_v1"
results_dir="$stage/ge_results"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf

test -d "$source_dir" && test ! -e "$results_dir/preflight_ge_v1"
cd "$source_dir"
sha256sum -c source.sha256
exec apptainer exec --bind "$reference:$original:ro,$source_dir:/work/ge_source:ro,$results_dir:/work/ge_results" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/ge_source/preflight_ge.py /work/ge_results/preflight_ge_v1
