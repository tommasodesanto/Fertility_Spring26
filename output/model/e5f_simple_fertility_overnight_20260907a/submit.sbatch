#!/bin/bash
#SBATCH --job-name=e5f_nest_overnight
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH --mem=322G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_overnight_20260907a/slurm-%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
python -m unittest discover -s code/model/tools -p 'test_e5f_simple_fertility*.py' -q
python code/model/tools/run_e5f_simple_fertility_overnight.py --contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_overnight_20260907a/contract.json --contract-sha256 a483fb6bd65230409062f82403038f3218a35a03e5b4e13c1f73d6846ac10f0c
