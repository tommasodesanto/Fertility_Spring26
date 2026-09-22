#!/usr/bin/env bash
# Package, stage, or explicitly submit the frozen earnings-wealth run.
# Default is stage-only; submission requires E5F_EW_SUBMIT=1.
set -euo pipefail

project_root="$(cd "$(dirname "$0")/../.." && pwd)"
plan="${E5F_EW_PLAN:-$project_root/output/model/native_financing_diagnostic_20260919/specification_followup/earnings_wealth_v1/overnight_plan.local.json}"
helper="$project_root/code/model/tools/prepare_e5f_earnings_wealth_run.py"
controller="run_e5f_earnings_wealth_search.py"
remote_root="${E5F_EW_REMOTE_ROOT:-/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/earnings_wealth_direct_period_v1}"
bundle="${E5F_EW_BUNDLE:?set E5F_EW_BUNDLE to a new local bundle directory}"
python_path="${E5F_EW_PYTHON:-/share/apps/anaconda3/2025.06/bin/python}"
mode="${1:-stage}"

prepare_local() {
  python3 "$helper" --plan "$plan" --bundle "$bundle" \
    --execution-root "$remote_root" --python "$python_path"
}

stage_bundle() {
  ssh torch "test ! -e '$remote_root' && mkdir -p '$remote_root'"
  rsync -a "$bundle/" "torch:$remote_root/"
  printf 'staged bundle=%s remote_root=%s\n' "$bundle" "$remote_root"
}

submit_jobs() {
  [[ "${E5F_EW_SUBMIT:-}" == 1 ]] || { echo 'refusing submission: set E5F_EW_SUBMIT=1' >&2; exit 2; }
  local smoke production
  smoke=$(ssh torch sbatch --parsable --account=torch_pr_570_general \
    --cpus-per-task=2 --mem=64G --time=02:00:00 \
    --output="$remote_root/slurm_smoke_%j.out" --error="$remote_root/slurm_smoke_%j.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
root='$remote_root'
timeout --signal=TERM --kill-after=60s 7000s "$python_path" -B "\$root/tools/$controller" --mode smoke --plan "\$root/plan.json" --output "\$root/output/smoke"
SBATCH
  )
  production=$(ssh torch sbatch --parsable --dependency="afterok:$smoke" --account=torch_pr_570_general \
    --cpus-per-task=4 --mem=256G --time=06:00:00 \
    --output="$remote_root/slurm_production_%j.out" --error="$remote_root/slurm_production_%j.err" <<SBATCH
#!/usr/bin/env bash
set -euo pipefail
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1
root='$remote_root'
timeout --signal=TERM --kill-after=60s 21000s "$python_path" -B "\$root/tools/$controller" --mode search --plan "\$root/plan.json" --output "\$root/output/search" --verified-smoke "\$root/output/smoke/smoke_receipt.json"
SBATCH
  )
  printf '%s\n' "smoke_job=$smoke production_job=$production" | tee "$bundle/submission_receipt.txt"
}

case "$mode" in
  prepare) prepare_local ;;
  stage) prepare_local; stage_bundle ;;
  submit) prepare_local; stage_bundle; submit_jobs ;;
  *) echo "usage: $0 [prepare|stage|submit]" >&2; exit 2 ;;
esac
