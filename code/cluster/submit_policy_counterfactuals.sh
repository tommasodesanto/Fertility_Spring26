#!/bin/bash
# Policy counterfactuals from a saved calibration record.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26/code/cluster

#SBATCH --job-name=dtpolcf
#SBATCH --output=logs/slurm_policy_cf_%A_%a.out
#SBATCH --error=logs/slurm_policy_cf_%A_%a.err
#SBATCH --array=1-2%2
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MODEL_DIR="${REPO_DIR}/code/model"
mkdir -p "${SCRIPT_DIR}/logs"

if command -v module >/dev/null 2>&1; then
    module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || module load python/3.11 2>/dev/null || module load python/3.10 2>/dev/null || true
fi

if [ -n "${DT_POLICY_PYTHON:-}" ]; then
    PYTHON_BIN="${DT_POLICY_PYTHON}"
else
    PYTHON_BIN="$(command -v python3 || command -v python)"
fi

export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}}"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mkdir -p "${NUMBA_CACHE_DIR}"

RECORD="${DT_POLICY_RECORD:-${REPO_DIR}/output/model/alpha_share_a24_w250_20260528_compare_job011_job001/records/alpha_a24_job011_best.json}"
OUTROOT="${DT_POLICY_OUTROOT:-${REPO_DIR}/output/model/policy_counterfactuals_alpha_a24_job011_cluster_dezoning_20260528}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"

case "${TASK_ID}" in
    1)
        CASE_ID="parent_phi_relief"
        EXTRA_ARGS=(--case parent_phi_relief --parent-phi "${DT_POLICY_PARENT_PHI:-1.0}")
        ;;
    2)
        CASE_ID="center_xi_relaxed"
        EXTRA_ARGS=(
            --case center_xi_relaxed
            --center-supply-elasticity "${DT_POLICY_CENTER_SUPPLY_ELASTICITY:-2.0}"
        )
        ;;
    *)
        echo "Unsupported task id ${TASK_ID}" >&2
        exit 2
        ;;
esac

OUTDIR="${OUTROOT}/${CASE_ID}"
mkdir -p "${OUTDIR}"

echo "============================================"
echo "Policy counterfactual"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array Task ID: ${TASK_ID}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Case: ${CASE_ID}"
echo "Record: ${RECORD}"
echo "Outdir: ${OUTDIR}"
echo "Python: ${PYTHON_BIN}"
echo "Started: $(date)"
echo "============================================"

"${PYTHON_BIN}" "${MODEL_DIR}/tools/run_policy_counterfactuals_from_record.py" \
    --record "${RECORD}" \
    --outdir "${OUTDIR}" \
    --baseline-mismatch fail \
    --baseline-tol 1e-3 \
    --fixed-baseline-max-iter-eq "${DT_POLICY_FIXED_BASELINE_MAX_ITER_EQ:-200}" \
    --policy-max-iter-eq "${DT_POLICY_POLICY_MAX_ITER_EQ:-200}" \
    "${EXTRA_ARGS[@]}"

echo "Finished: $(date)"
