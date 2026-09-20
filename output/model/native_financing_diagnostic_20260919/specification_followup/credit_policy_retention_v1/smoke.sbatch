#!/usr/bin/env bash
#SBATCH --job-name=cpr_smoke
#SBATCH --account=torch_pr_570_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --chdir=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1/logs/%x_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1/logs/%x_%j.err
#SBATCH --kill-on-invalid-dep=yes
set -euo pipefail
module load anaconda3/2025.06
unset PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0
export NUMBA_CACHE_DIR="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1/numba_cache" MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
cd "/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1"
W="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1/code/model/tools/run_e5f_credit_policy_retention.py"; M="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1/launch_manifest.json"; RES="/scratch/td2248/projects/Fertility_Spring26_specification_20260920/credit_policy_retention_v1/results"
python -B "$W" --mode verify --manifest "$M" --results "$RES" --receipt "$RES/verify_smoke_${SLURM_JOB_ID}.json"
DEADLINE=$(( $(date +%s) + 1800 ))
timeout --signal=TERM --kill-after=30s 600s python -B "$W" --mode case --case 0 --manifest "$M" --results "$RES" --deadline-epoch "$DEADLINE"
