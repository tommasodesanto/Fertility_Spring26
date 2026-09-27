#!/bin/bash
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/housing_fertility_cost_diagnostic_20260926/continuation_001/slurm-%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
export PYTHONPATH=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/first_child_loading_probe_20260926/tools:/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2/tools:/scratch/td2248/commute_pdf_qa_deps
r=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/housing_fertility_cost_diagnostic_20260926
python "$r/tools/e5f_housing_fertility_cost_resume.py" --stage smoke --source "$r/run_001" --output "$r/continuation_001" --combined "$r/combined_001"
python "$r/tools/e5f_housing_fertility_cost_resume.py" --stage run --source "$r/run_001" --output "$r/continuation_001" --combined "$r/combined_001"
python "$r/tools/e5f_housing_fertility_cost_audit.py" --output "$r/combined_001"
python "$r/tools/run_e5f_housing_fertility_cost_diagnostic.py" --output "$r/combined_001" --stage render
