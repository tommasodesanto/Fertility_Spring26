#!/bin/bash
#SBATCH --job-name=simple_fertility
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_20260907a/output/model/e5f_simple_fertility_nest_20260907a/simple_fertility_%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_20260907a
python code/model/tools/run_e5f_bounded_calibration_refinement.py --plan output/model/e5f_simple_fertility_nest_20260907a/simple_fertility_plan.json --plan-sha256 d163e3812ac7a6623e0aed2c58ee7b0af52648fcb4e241597c78b3e2b31d09dc --case-id 1
