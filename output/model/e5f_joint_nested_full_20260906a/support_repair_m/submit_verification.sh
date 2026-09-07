#!/bin/bash
#SBATCH --job-name=e5f_support_verify
#SBATCH --partition=cs
#SBATCH --qos=cpu48
#SBATCH --account=torch_pr_570_general
#SBATCH --time=01:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
set -euo pipefail
module load anaconda3/2025.06
cd /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907m
export TMPDIR="$PWD/tmp" NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="$TMPDIR/mpl"
python3 code/model/tools/test_e5f_exhaustive_saving.py
python3 code/model/tools/test_e5f_exhaustive_saving.py
python3 -u verify_original_policy.py
export E5F_JOINT_MODE=smoke E5F_JOINT_CONTRACT="$PWD/output/model/joint_nested_overnight/contract.json"
export E5F_JOINT_CONTRACT_SHA256="$(sha256sum "$E5F_JOINT_CONTRACT" | cut -d ' ' -f1)"
exec bash code/cluster/submit_e5f_joint_nested_long.sh
