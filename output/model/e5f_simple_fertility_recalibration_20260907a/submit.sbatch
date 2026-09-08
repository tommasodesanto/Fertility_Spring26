#!/bin/bash
#SBATCH --job-name=fertility-nest-calibrate
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=23
#SBATCH --mem=322G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_recalibration_20260907a/output/model/e5f_simple_fertility_recalibration_20260907a/search_%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_recalibration_20260907a
python -m unittest discover -s code/model/tools -p 'test_e5f_simple_fertility*.py' -v
python code/model/tools/run_e5f_simple_fertility_search.py --contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_recalibration_20260907a/output/model/e5f_simple_fertility_recalibration_20260907a/contract.json --contract-sha256 fac7694e8ccf20a9b3c7f243e31f57bc17e804b5469180597e7b23ebe7e5a174
