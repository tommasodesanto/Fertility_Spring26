#!/usr/bin/env bash
# Three fixed-price household solves on Torch; no search or recalibration.
#SBATCH --job-name=e5f_estate_risk
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then
  : "${E5F_ESTATE_RISK_ROOT:?preserve the diagnostic root at submission}"
else
  risk_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  export E5F_ESTATE_RISK_ROOT="$(cd "$risk_script_dir/../.." && pwd)"
fi
: "${E5F_ESTATE_RISK_OUTPUT:?set a new remote output directory}"
if [ "${1:-}" = "--submit" ]; then
  shift
  test ! -e "$E5F_ESTATE_RISK_OUTPUT"
  mkdir -p "$(dirname "$E5F_ESTATE_RISK_OUTPUT")"
  exec sbatch --output="$E5F_ESTATE_RISK_OUTPUT-%j.log" --export=ALL "$@" "$0"
fi
: "${SLURM_JOB_ID:?Torch allocation required; use --submit on Torch}"
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/scratch/td2248/commute_pdf_qa_deps${PYTHONPATH:+:$PYTHONPATH}
export NUMBA_CACHE_DIR="$E5F_ESTATE_RISK_OUTPUT-numba-cache"
mkdir -p "$NUMBA_CACHE_DIR"
cd "$E5F_ESTATE_RISK_ROOT"
exec timeout --signal=TERM --kill-after=30s 1500s python3 code/model/tools/run_e5f_estate_receipt_risk.py --output-dir "$E5F_ESTATE_RISK_OUTPUT" --subdivide-wealth-grid "${E5F_ESTATE_RISK_GRID_SUBDIVISION:-1}"
