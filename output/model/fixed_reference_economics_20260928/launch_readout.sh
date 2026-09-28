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
original=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
case_dir=output/model/fixed_reference_economics_20260928
source_dir="$case_dir/sources/readout_v1"
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/scratch/td2248/projects/fertility_evening_calibration_20260927_v1/report_deps
cd "$stage/project/$source_dir"
sha256sum -c source.sha256
exec apptainer exec --bind "$stage/project:$original:ro" --bind "$stage/project/output/pdf:$original/output/pdf:rw" --pwd "$original" /share/apps/images/ubuntu-24.04.4.sif /share/apps/anaconda3/2025.06/bin/python "$source_dir/build_readout.py" --findings "$source_dir/findings.json" --output output/pdf/fixed_reference_economics_block0506.pdf
