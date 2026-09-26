#!/usr/bin/env bash
# Prepare-only launcher.  It never submits unless --submit is explicit and both
# successful Torch receipts are supplied; it never polls disconnected jobs.
#SBATCH --job-name=e5f_estate_probe
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --time=03:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then
  : "${E5F_ESTATE_PROBE_ROOT:?submission must preserve the diagnostic project root}"
  PROJECT_ROOT="$E5F_ESTATE_PROBE_ROOT"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
fi
DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_estate_receiver_probe.py"
need_receipt() { [ -s "$1" ] || { echo "missing required receipt: $1" >&2; exit 2; }; }
if [ "${1:-}" = "--submit" ]; then
  shift
  : "${E5F_ESTATE_PROBE_STAGE:?stage required before submission}"
  : "${E5F_ESTATE_PROBE_OUTPUT:?new remote output directory required}"
  test ! -e "$E5F_ESTATE_PROBE_OUTPUT"
  case "$E5F_ESTATE_PROBE_STAGE" in
    preflight) stage_limit=00:20:00;;
    smoke) stage_limit=00:20:00;;
    run) stage_limit=03:10:00;;
    *) echo "invalid stage" >&2; exit 2;;
  esac
  case "$E5F_ESTATE_PROBE_STAGE" in
    preflight) ;;
    smoke) : "${E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT:?successful source-preflight receipt required}"; need_receipt "$E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT" ;;
    run) : "${E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT:?successful source-preflight receipt required}"; : "${E5F_ESTATE_PROBE_SMOKE_RECEIPT:?successful exact-loop smoke receipt required}"; need_receipt "$E5F_ESTATE_PROBE_PREFLIGHT_RECEIPT"; need_receipt "$E5F_ESTATE_PROBE_SMOKE_RECEIPT" ;;
    *) echo "invalid stage" >&2; exit 2 ;;
  esac
  export E5F_ESTATE_PROBE_ROOT="$PROJECT_ROOT"
  exec sbatch --account=torch_pr_570_general --partition=cs --time="$stage_limit" --cpus-per-task=1 --mem=16G --export=ALL "$0" "$@"
fi
: "${SLURM_JOB_ID:?prepare only: invoke with --submit after review and successful receipts}"
: "${E5F_ESTATE_PROBE_STAGE:?set preflight, smoke, or run}"
: "${E5F_ESTATE_PROBE_OUTPUT:?new remote output directory required}"
case "$E5F_ESTATE_PROBE_STAGE" in preflight|smoke|run) ;; *) echo "invalid stage" >&2; exit 2;; esac
for key in NUMBA_NUM_THREADS OMP_NUM_THREADS MKL_NUM_THREADS OPENBLAS_NUM_THREADS; do export "$key=1"; done
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/scratch/td2248/commute_pdf_qa_deps${PYTHONPATH:+:$PYTHONPATH}
export NUMBA_CACHE_DIR="$E5F_ESTATE_PROBE_OUTPUT-numba-cache"
mkdir -p "$NUMBA_CACHE_DIR"
export E5F_ESTATE_RECEIVER_TORCH=1
ARGS=(--stage "$E5F_ESTATE_PROBE_STAGE" --output-dir "$E5F_ESTATE_PROBE_OUTPUT")
if [[ "$E5F_ESTATE_PROBE_STAGE" == run ]]; then stage_seconds=10800; else stage_seconds=1100; fi
exec timeout --signal=TERM --kill-after=30s "${stage_seconds}s" python3 "$DRIVER" "${ARGS[@]}"
