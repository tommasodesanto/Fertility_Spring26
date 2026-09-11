#!/bin/bash
#SBATCH --job-name=e5f_two_node_audit
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:06:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/branch_audit_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508/branch_audit_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=1 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_initial_revision_c6dd3508
timeout --signal=KILL 300s python contracts/branch_audit/run_branch_audit.py --contract contracts/branch_audit/input_contract.json --contract-sha256 1c63d2549fbb0113b787bcbdaa0d4e06da919f02330c9de8ad09173782899394 --output "output/branch_audit_${SLURM_JOB_ID}"
