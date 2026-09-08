#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:35:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_policy_impact_20260908a/slurm-%A-%a.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_policy_impact_20260908a/mpl
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
cases=(baseline dependent-child-ltv95 supply-plus-20)
python code/model/tools/run_e5f_simple_fertility_policy_impact.py --contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_policy_impact_20260908a/contract.json --contract-sha256 1f7beed6c5ededada14cce5b7b480aee7fcc4619a3963fb09b0b6b7620ccd3ca --case "${cases[$SLURM_ARRAY_TASK_ID]}"
