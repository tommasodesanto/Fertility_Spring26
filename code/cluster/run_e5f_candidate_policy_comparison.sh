#!/usr/bin/env bash
# Bounded policy comparison using an already verified scientific snapshot.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cpu_short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:55:00
set -euo pipefail
task_root=${1:?frozen project root}
plan=${2:?immutable plan}
plan_sha=${3:?plan SHA-256}
stage=${4:?smoke, full or rebate}
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$task_root/code/model:$task_root/code/model/tools"
cases=(baseline supply-plus-20 dependent-child-ltv95 combined property-tax-2pct-no-rebate)
case_index=${SLURM_ARRAY_TASK_ID:-1}
[[ "$case_index" =~ ^[1-5]$ ]] || { echo 'case index must be 1..5' >&2; exit 2; }
python3 "$task_root/code/model/tools/run_e5f_candidate_policy_comparison.py" \
  --plan "$plan" --plan-sha256 "$plan_sha" --stage "$stage" \
  --case "${cases[$((case_index-1))]}"
