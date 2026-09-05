#!/usr/bin/env bash
# Frozen-code dated reconstruction and exact-loop independent numerical audit.
#SBATCH --job-name=e5findaudit
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general
#SBATCH --output=logs/slurm_e5findaudit_%j.out
#SBATCH --error=logs/slurm_e5findaudit_%j.err
set -euo pipefail
: "${E5F_AUDIT_PRODUCTION_ROOT:?required original production snapshot}"
: "${E5F_AUDIT_OUTPUT:?required unique absolute output directory}"
: "${E5F_AUDIT_EXPECTED_DRIVER_SHA256:?required driver hash}"
CLUSTER_DIR="${SLURM_SUBMIT_DIR:?submit from frozen code/cluster}"
AUDIT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
KIND="${E5F_AUDIT_KIND:-numerical}"
case "$KIND" in
  numerical) DRIVER="$AUDIT_ROOT/code/model/tools/run_e5f_independent_numerical_audit.py" ;;
  global-saving) DRIVER="$AUDIT_ROOT/code/model/tools/run_e5f_global_saving_quantification.py" ;;
  *) echo 'Unknown audit kind' >&2; exit 2 ;;
esac
ACTUAL_HASH="$(sha256sum "$DRIVER" | awk '{print $1}')"
[[ "$ACTUAL_HASH" = "$E5F_AUDIT_EXPECTED_DRIVER_SHA256" ]] || { echo 'Audit driver hash mismatch' >&2; exit 2; }
module load anaconda3/2025.06
export PYTHONPATH="$AUDIT_ROOT/code/model:$AUDIT_ROOT/code/model/tools"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/e5findaudit_mpl_${SLURM_JOB_ID}"
cd "$AUDIT_ROOT"
ARGS=(--outdir "$E5F_AUDIT_OUTPUT")
if [[ "$KIND" = numerical ]]; then
  ARGS+=(--production-root "$E5F_AUDIT_PRODUCTION_ROOT" --draws "${E5F_AUDIT_DRAWS:-240}")
else
  : "${E5F_AUDIT_CHECKPOINT:?required saved dated state}"
  ARGS+=(--mode "${E5F_AUDIT_MODE:-smoke}")
  ARGS+=(--market-tol "${E5F_AUDIT_MARKET_TOL:-0.0002}")
  if [[ "${E5F_AUDIT_MODE:-smoke}" = equilibrium ]]; then
    : "${E5F_AUDIT_SMOKE_SUMMARY:?required exact-loop smoke}"
    : "${E5F_AUDIT_SMOKE_SUMMARY_SHA256:?required smoke hash}"
    ARGS+=(--smoke-summary "$E5F_AUDIT_SMOKE_SUMMARY" --smoke-summary-sha256 "$E5F_AUDIT_SMOKE_SUMMARY_SHA256")
  fi
fi
if [[ -n "${E5F_AUDIT_CHECKPOINT:-}" ]]; then
  : "${E5F_AUDIT_CHECKPOINT_SHA256:?required checkpoint hash}"
  ARGS+=(--checkpoint "$E5F_AUDIT_CHECKPOINT" --checkpoint-sha256 "$E5F_AUDIT_CHECKPOINT_SHA256")
fi
python3 -u "$DRIVER" "${ARGS[@]}"
