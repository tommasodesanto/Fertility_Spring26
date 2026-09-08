#!/bin/bash
#SBATCH --job-name=fertility-best-repeat
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=1-2
#SBATCH --time=01:40:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_recalibration_20260907a/output/model/e5f_simple_fertility_recalibration_20260907a/reviewed_repeat_%A_%a.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_recalibration_20260907a
python code/model/tools/run_e5f_bounded_calibration_refinement.py --plan output/model/e5f_simple_fertility_recalibration_20260907a/search/reviewed_repeats/plan.json --plan-sha256 8ceee024ffba823327496d128d383fe7bed838cc7bd966e5ec6a7606a0b68079 --case-id "$SLURM_ARRAY_TASK_ID"
