#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:03:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
overlay=/scratch/td2248/projects/fixed_reference_transition_20260928/preparation_v1
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
relative=output/model/fixed_reference_transition_20260928/preparation_v1
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools${PYTHONPATH:+:$PYTHONPATH}"
export NUMBA_CACHE_DIR="$overlay/numba_cache"
exec apptainer exec --bind "$stage/project:$original:ro" --bind "$overlay:$original/$relative:rw" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$relative/inspect_parameter_difference.py"
