#!/usr/bin/env bash
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
stage=/scratch/td2248/projects/fertility_night_calibration_20260928_v1
folder=output/model/fixed_reference_theory_20260928
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/scratch/td2248/projects/fertility_evening_calibration_20260927_v1/report_deps
exec apptainer exec --bind "$stage/project:$stage/project" --pwd "$stage/project" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$stage/project/$folder/build_note.py"
