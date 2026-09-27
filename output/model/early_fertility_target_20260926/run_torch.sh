#!/bin/bash
#SBATCH --job-name=cps_early_fertility
#SBATCH --account=torch_pr_570_general
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
cd /scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/early_fertility_target_20260926
python build_early_fertility_target.py --receipt fertility_availability.json --partitions input --compressed-source input/cps_00003.dat.gz --schema loader.do --output output
