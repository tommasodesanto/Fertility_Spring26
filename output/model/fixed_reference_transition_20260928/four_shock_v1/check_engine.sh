#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:16:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
overlay=/scratch/td2248/projects/fixed_reference_transition_20260928
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
relative=output/model/fixed_reference_transition_20260928
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONDONTWRITEBYTECODE=1
export E5F_TRANSITION_SOURCE_SUBDIR=source_v2
export NUMBA_CACHE_DIR="$overlay/four_shock_v1/numba_cache"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
exec apptainer exec --bind "$stage/project:$original:ro" --bind "$overlay:$original/$relative:rw" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$relative/four_shock_v1/check_engine.py"
