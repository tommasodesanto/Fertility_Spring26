#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/smoke.log
#SBATCH --error=/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/smoke.err
set -euo pipefail
reference=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
old=/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929
old_credit=/scratch/td2248/projects/fixed_reference_credit_20260929/results/solve_v1
stage=/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$original/tmp/e5f_overnight_local_20260927/portable/tools_v4:$original/code/model/tools:/work/recovery_source:/work/elasticity_source"
export EXPECTED_E5F_IDENTIFICATION_SHA256=68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf
exec apptainer exec --bind "$reference:$original:ro,$old/source_v2:/work/elasticity_source:ro,$old/results_v2:/work/elasticity_results:ro,$stage/source:/work/recovery_source:ro,$stage:/work/recovery_results:rw,$old_credit:/work/credit_seed:ro" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python /work/recovery_source/smoke.py
