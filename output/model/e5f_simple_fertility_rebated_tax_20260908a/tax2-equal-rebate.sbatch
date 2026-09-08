#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=00:35:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_rebated_tax_20260908a/slurm-%j.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_rebated_tax_20260908a/mpl
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
python code/model/tools/run_e5f_simple_fertility_tax_impact.py --contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_rebated_tax_20260908a/contract.json --contract-sha256 931def1c35d3a6b4575bcb9e8cd65e345f71804d1b25d74f3671a8719b680a7d --case tax2-equal-rebate
