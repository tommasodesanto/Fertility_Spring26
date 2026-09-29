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
stage=/scratch/td2248/projects/fixed_reference_elasticity_20260929
old_credit=/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
source="$stage/source_v1"
results="$stage/results_v1"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$stage/cache"
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools:/work/elasticity_source"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
test -d "$reference" && test -d "$source" && test -d "$old_credit"
test ! -e "$results/preflight_v1"
mkdir -p "$results" "$stage/cache"
exec apptainer exec --bind "$reference:$original:ro,$source:/work/elasticity_source:ro,$results:/work/elasticity_results:rw,$old_credit:/work/credit_seed:ro" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/elasticity_source/preflight.py --output /work/elasticity_results/preflight_v1
