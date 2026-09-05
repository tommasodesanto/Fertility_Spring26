#!/usr/bin/env bash
# Isolated nested-grid impact sensitivity; no production/calibration changes.
#SBATCH --job-name=e5fgrid
#SBATCH --partition=cpu_short
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=torch_pr_570_general
#SBATCH --output=logs/slurm_e5fgrid_%j.out
#SBATCH --error=logs/slurm_e5fgrid_%j.err
set -euo pipefail
: "${E5F_GRID_OUTPUT:?unique output directory required}"
: "${E5F_GRID_STAGE:?smoke, fine-policies, or finer required}"
: "${E5F_GRID_DRIVER_SHA256:?driver hash required}"
: "${E5F_GRID_CHECKPOINT:?certified inherited checkpoint required}"
: "${E5F_GRID_REFERENCE_CSV:?completed global comparison required}"
: "${E5F_GRID_REFERENCE_SHA256:?reference hash required}"
CLUSTER_DIR="${SLURM_SUBMIT_DIR:?submit from frozen code/cluster}"
AUDIT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$AUDIT_ROOT/code/model/tools/run_e5f_nested_wealth_grid_audit.py"
ACTUAL_HASH="$(sha256sum "$DRIVER" | awk '{print $1}')"
[[ "$ACTUAL_HASH" = "$E5F_GRID_DRIVER_SHA256" ]] || { echo 'Grid driver hash mismatch' >&2; exit 2; }
module load anaconda3/2025.06
export PYTHONPATH="$AUDIT_ROOT/code/model:$AUDIT_ROOT/code/model/tools"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/e5fgrid_mpl_${SLURM_JOB_ID}"
cd "$AUDIT_ROOT"
ARGS=(--stage "$E5F_GRID_STAGE" --outdir "$E5F_GRID_OUTPUT" --checkpoint "$E5F_GRID_CHECKPOINT"
      --reference-csv "$E5F_GRID_REFERENCE_CSV" --reference-csv-sha256 "$E5F_GRID_REFERENCE_SHA256")
if [[ "$E5F_GRID_STAGE" != smoke ]]; then
  : "${E5F_GRID_PREVIOUS_SUMMARY:?completed preceding stage required}"
  : "${E5F_GRID_PREVIOUS_SHA256:?preceding stage hash required}"
  ARGS+=(--previous-summary "$E5F_GRID_PREVIOUS_SUMMARY" --previous-summary-sha256 "$E5F_GRID_PREVIOUS_SHA256")
fi
python3 -u "$DRIVER" "${ARGS[@]}"
