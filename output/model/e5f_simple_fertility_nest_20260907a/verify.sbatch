#!/bin/bash
#SBATCH --job-name=fertility-nest-check
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:35:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_20260907a/output/model/e5f_simple_fertility_nest_20260907a/verify_%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_20260907a
python -m unittest discover -s code/model/tools -p 'test_e5f_simple_fertility*.py' -v
python code/model/tools/run_e5f_simple_fertility_probe.py --checkpoint output/model/e5f_simple_fertility_nest_20260907a/dated_state.pkl --output output/model/e5f_simple_fertility_nest_20260907a/cluster_probe --seconds 1800 --skip-reference --cases fertility_nest_fixed_price sequential_exhaustive_fixed_price
