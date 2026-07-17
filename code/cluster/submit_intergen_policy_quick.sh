#!/bin/bash
# Parallel fixed-theta legacy policy battery for the repaired one-market model.

#SBATCH --job-name=ihfpolicy
#SBATCH --output=logs/slurm_ihf_policy_%A_%a.out
#SBATCH --error=logs/slurm_ihf_policy_%A_%a.err
#SBATCH --array=1-5%5
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MODEL_DIR="${REPO_DIR}/code/model"
mkdir -p "${SCRIPT_DIR}/logs"

PYTHON_BIN="${INTERGEN_POLICY_PYTHON:-/share/apps/anaconda3/2025.06/bin/python}"
SOURCE="${INTERGEN_POLICY_SOURCE:?set INTERGEN_POLICY_SOURCE to the saved best.json}"
OUTROOT="${INTERGEN_POLICY_OUTROOT:?set INTERGEN_POLICY_OUTROOT}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}"

case "${TASK_ID}" in
    1) POLICY_CASE="parent_ltv100_birth" ;;
    2) POLICY_CASE="universal_ltv95" ;;
    3) POLICY_CASE="birth_grant_A0p4_Hge6" ;;
    4) POLICY_CASE="property_tax_2pct" ;;
    5) POLICY_CASE="tax2_grant_A0p4_Hge6" ;;
    *) echo "Unsupported task ${TASK_ID}" >&2; exit 2 ;;
esac

OUTDIR="${OUTROOT}/task_${TASK_ID}_${POLICY_CASE}"
export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}}"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mkdir -p "${NUMBA_CACHE_DIR}" "${OUTDIR}"

"${PYTHON_BIN}" "${MODEL_DIR}/tools/build_intergen_mechanics_packet.py" \
    --source "${SOURCE}" \
    --outdir "${OUTDIR}" \
    --target-set candidate_replacement_post_audit_v1 \
    --J 17 --Nb 120 --income-states 5 --n-house 5 --hR-max 6.0 \
    --max-iter-eq 10 --combined-corrected-spec \
    --run-policy-cases --policy-battery legacy --policy-case "${POLICY_CASE}" \
    --skip-standard-diagnostics --skip-contact-sheet --quick-first-look-only --no-csv
