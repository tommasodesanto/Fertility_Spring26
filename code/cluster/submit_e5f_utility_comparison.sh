#!/usr/bin/env bash
# Submit this file on Torch only, after the author-approved immutable contract
# has been frozen and the exact controller/collector loop has been checked.
# Required: E5F_UTILITY_TOOLS, E5F_UTILITY_CONTRACT, E5F_UTILITY_RUN_ROOT,
# EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256. This script never writes approval.
#SBATCH --job-name=e5f_utility_compare
#SBATCH --partition=cs
#SBATCH --account=torch_pr_570_general
#SBATCH --array=0-3%4
#SBATCH --time=08:35:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --no-requeue
#SBATCH --output=slurm_e5f_utility_%A_%a.out

set -euo pipefail
utility_tools="${E5F_UTILITY_TOOLS:?absolute isolated bundle tools directory}"
utility_contract="${E5F_UTILITY_CONTRACT:?absolute author-approved contract path}"
utility_results="${E5F_UTILITY_RUN_ROOT:?new absolute shared result directory}"
: "${EXPECTED_UTILITY_COMPARISON_CONTRACT_SHA256:?approved immutable contract SHA256}"
: "${SLURM_JOB_ID:?Torch Slurm allocation required}"
case "$utility_tools:$utility_contract:$utility_results" in
  /*:/*:/*) ;;
  *) printf '%s\n' 'All utility comparison paths must be absolute.' >&2; exit 2 ;;
esac
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg NUMBA_DISABLE_JIT=0
utility_reference="$(python3 - "$utility_contract" <<'PY'
import json, pathlib, sys
path = pathlib.Path(json.loads(pathlib.Path(sys.argv[1]).read_text())["reference_root"])
if not path.is_absolute():
    raise SystemExit("reference_root must be an absolute frozen source path")
print(path)
PY
)"
export PYTHONPATH="$utility_tools:$utility_reference/source/code/model/tools:$utility_reference/source/code/model:/scratch/td2248/commute_pdf_qa_deps"
utility_arms=(floor_linear floor_concave shares_linear shares_concave)
utility_slot="${SLURM_ARRAY_TASK_ID:?four-arm Slurm array required}"
case "$utility_slot" in 0|1|2|3) ;; *) exit 2 ;; esac
exec python3 "$utility_tools/run_e5f_utility_comparison_search.py" \
  --contract "$utility_contract" --run-root "$utility_results" \
  --arm "${utility_arms[$utility_slot]}"
