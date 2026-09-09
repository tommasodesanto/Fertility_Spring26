#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:35:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_channels_20260908a/slurm-%A_%a.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_channels_20260908a/mpl
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
python code/model/tools/run_e5f_simple_fertility_tax_channels.py --contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_channels_20260908a/contract.json --contract-sha256 de767353697c706dec4a7d960d1f3314dfb5cf2debad4c17e8438ce021966288 --cell "${SLURM_ARRAY_TASK_ID}"
