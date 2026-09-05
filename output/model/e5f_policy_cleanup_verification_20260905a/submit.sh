#!/bin/bash
#SBATCH --job-name=e5f-reuse-check
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cpu_short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
set -euo pipefail
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /scratch/td2248/projects/Fertility_Spring26_policy_cleanup_20260905a
export PYTHONPATH="$PWD/code/model:$PWD/code/model/tools"
python3 code/model/tools/test_calendar_policy_reuse.py -v
python3 code/model/tools/test_e5f_transition_accounting.py
python3 output/model/e5f_policy_cleanup_verification_20260905a/replay.py \
 --manifest output/model/e5f_policy_cleanup_verification_20260905a/manifest.json \
 --manifest-sha256 6ceed97805435419286668acad3f543b543b10b8d92e78ddeb09c59528720cac \
 --outdir output/model/e5f_policy_cleanup_verification_20260905a/replay
