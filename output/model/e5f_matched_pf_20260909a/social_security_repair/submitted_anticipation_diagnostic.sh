#!/bin/bash
#SBATCH --job-name=e5f_ss_diagnose
#SBATCH --partition=cpu_short
#SBATCH --account=torch_pr_570_general
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/diagnose_%j.out
#SBATCH --error=/scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/diagnose_%j.err
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 PYTHONUNBUFFERED=1
cd /scratch/td2248/projects/Fertility_Spring26_paygo_smoke_a654219c
export PYTHONPATH=code/model/tools:code/model
export MPLCONFIGDIR="$PWD/output/cache/mpl"
export NUMBA_CACHE_DIR="$PWD/output/cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
python /scratch/td2248/projects/Fertility_Spring26_paygo_audit_20260911a/diagnose_compiled_anticipation.py
