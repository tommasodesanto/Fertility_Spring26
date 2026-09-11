#!/bin/bash
#SBATCH --job-name=e5f_ss_audit
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:10:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/audit_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/audit_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=1
export MPLCONFIGDIR=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/mpl
cd /scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a
python audit_e5f_social_security.py \
  --source-root /scratch/td2248/projects/Fertility_Spring26_preference_shape_20260910b \
  --contract contract.json \
  --contract-sha256 13ac77b34c5b71094bf98461cb641260b7b32d130d1b88c599003132598d54a3 \
  --helper e5f_social_security.py \
  --helper-sha256 f7bb7db2774be7ef7d2ca2a635eafa3cb27ba79a503fe69b6f00c1eef9658455 \
  --output budget_audit.json
