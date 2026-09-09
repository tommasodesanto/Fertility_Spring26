#!/bin/bash
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=04:05:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/slurm-%A_%a.out
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/mpl
export TMPDIR=/scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/tmp-${SLURM_JOB_ID}
mkdir -p "$TMPDIR" "$MPLCONFIGDIR"
cd /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a
cases=(tax1-equal-rebate tax2-equal-rebate)
python code/model/tools/run_e5f_simple_fertility_tax_transition.py --contract /scratch/td2248/projects/Fertility_Spring26_simple_fertility_overnight_20260907a/output/model/e5f_simple_fertility_tax_transition_20260908a/contract.json --contract-sha256 ce80b6ad241bec9556d2f3d3cdccb9f89e0ccfdd5b2e80ddcd167ed5156b5b68 --stage full --case "${cases[$SLURM_ARRAY_TASK_ID]}"
