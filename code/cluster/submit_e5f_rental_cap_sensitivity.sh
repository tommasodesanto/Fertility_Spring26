#!/usr/bin/env bash
# Explicit cap-eight economic sensitivity; production remains unchanged.
#SBATCH --job-name=e5fcap8
#SBATCH --partition=cpu_short
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general
#SBATCH --output=logs/slurm_e5fcap8_%j.out
#SBATCH --error=logs/slurm_e5fcap8_%j.err
set -euo pipefail
: "${E5F_CAP_OUTPUT:?unique output directory required}"
: "${E5F_CAP_STAGE:?smoke or policies required}"
: "${E5F_CAP_DRIVER_SHA256:?driver hash required}"
: "${E5F_CAP_CHECKPOINT:?certified 239-node baseline required}"
CLUSTER_DIR="${SLURM_SUBMIT_DIR:?submit from frozen code/cluster}"
AUDIT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$AUDIT_ROOT/code/model/tools/run_e5f_rental_cap_sensitivity.py"
ACTUAL_HASH="$(sha256sum "$DRIVER" | awk '{print $1}')"
[[ "$ACTUAL_HASH" = "$E5F_CAP_DRIVER_SHA256" ]] || { echo 'Cap driver hash mismatch' >&2; exit 2; }
module load anaconda3/2025.06
export PYTHONPATH="$AUDIT_ROOT/code/model:$AUDIT_ROOT/code/model/tools"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/e5fcap8_mpl_${SLURM_JOB_ID}"
cd "$AUDIT_ROOT"
ARGS=(--stage "$E5F_CAP_STAGE" --outdir "$E5F_CAP_OUTPUT" --checkpoint "$E5F_CAP_CHECKPOINT")
if [[ "$E5F_CAP_STAGE" = policies ]]; then
  : "${E5F_CAP_SMOKE_SUMMARY:?passed exact-loop smoke required}"
  : "${E5F_CAP_SMOKE_SHA256:?smoke hash required}"
  ARGS+=(--smoke-summary "$E5F_CAP_SMOKE_SUMMARY" --smoke-summary-sha256 "$E5F_CAP_SMOKE_SHA256")
fi
python3 -u "$DRIVER" "${ARGS[@]}"
