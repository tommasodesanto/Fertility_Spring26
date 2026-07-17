#!/bin/bash
# Full parallel fixed-theta policy battery for the repaired one-market model.

#SBATCH --job-name=ihfpolall
#SBATCH --output=logs/slurm_ihf_polall_%A_%a.out
#SBATCH --error=logs/slurm_ihf_polall_%A_%a.err
#SBATCH --array=1-21%21
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
SOURCE="${INTERGEN_POLICY_SOURCE:?set INTERGEN_POLICY_SOURCE}"
OUTROOT="${INTERGEN_POLICY_OUTROOT:?set INTERGEN_POLICY_OUTROOT}"
MAX_ITER_EQ="${INTERGEN_POLICY_MAX_ITER_EQ:-10}"
SCALAR_REFINE_ITER="${INTERGEN_POLICY_SCALAR_REFINE_ITER:-}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?missing SLURM_ARRAY_TASK_ID}"
CASES=(
  parent_ltv95_birth parent_ltv100_birth universal_ltv95
  birth_grant_A0p1_Hge6 birth_grant_A0p2_Hge6 birth_grant_A0p4_Hge6
  birth_grant_A0p58_Hge6 birth_grant_A1p0_Hge6 birth_grant_A0p4_all_rungs
  property_tax_2pct tax2_grant_A0p4_Hge6
  rental_hR6p5 rental_hR7 rental_hR7p5 rental_hR8
  supply_H0_plus10 supply_H0_plus20 grant_A0p4_Hge6_supply_plus20
  debt_line_lambda0p2 debt_line_lambda0p4 debt_line0p4_supply_plus20
)
POLICY_CASE="${CASES[$((TASK_ID - 1))]}"
OUTDIR="${OUTROOT}/task_${TASK_ID}_${POLICY_CASE}"

export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "${NUMBA_CACHE_DIR}" "${OUTDIR}"

EXTRA_ARGS=()
if [ -n "${SCALAR_REFINE_ITER}" ]; then
  EXTRA_ARGS+=(--scalar-market-refine-iter "${SCALAR_REFINE_ITER}")
fi

"${PYTHON_BIN}" "${MODEL_DIR}/tools/build_intergen_mechanics_packet.py" \
  --source "${SOURCE}" --outdir "${OUTDIR}" \
  --target-set candidate_replacement_post_audit_v1 \
  --J 17 --Nb 120 --income-states 5 --n-house 5 --hR-max 6.0 --max-iter-eq "${MAX_ITER_EQ}" \
  --combined-corrected-spec --run-policy-cases --policy-battery full --policy-case "${POLICY_CASE}" \
  --skip-standard-diagnostics --skip-contact-sheet --quick-first-look-only --no-csv \
  "${EXTRA_ARGS[@]}"
