#!/usr/bin/env bash
set -euo pipefail
work=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1
run_root="$work/results/run_001"
if [[ "${1:-}" != --submit || "${APPROVED_PAIR_SUBMISSION:-}" != YES ]]; then
  echo "Preparation only. A lead-reviewed lock and explicit --submit with APPROVED_PAIR_SUBMISSION=YES are required." >&2
  exit 2
fi
: "${EXPECTED_PAIR_LOCK_SHA256:?exact reviewed lock SHA required}"
test -f "$work/render_pair.py"
test ! -e "$run_root"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHONDONTWRITEBYTECODE=1 python3 "$work/run_pair.py" --stage preflight
if squeue -h -u td2248 -o '%j' | grep -q '^e5f-nightpair-'; then
  echo "Existing paired Slurm job found; refusing duplicate submission" >&2
  exit 2
fi
mkdir -p "$run_root"
smoke=$(sbatch --parsable --account=torch_pr_570_general --partition=cs --array=1-2%2 \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=06:10:00 \
  --export=ALL,EXPECTED_PAIR_LOCK_SHA256="$EXPECTED_PAIR_LOCK_SHA256" \
  --job-name=e5f-nightpair-smoke --output="$run_root/smoke_%A_%a.out" \
  --error="$run_root/smoke_%A_%a.err" "$work/run_task.sh" smoke "$run_root")
workers=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterok:$smoke" \
  --account=torch_pr_570_general --partition=cs --array=1-40%40 \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=06:10:00 \
  --export=ALL,EXPECTED_PAIR_LOCK_SHA256="$EXPECTED_PAIR_LOCK_SHA256" \
  --job-name=e5f-nightpair-worker --output="$run_root/worker_%A_%a.out" \
  --error="$run_root/worker_%A_%a.err" "$work/run_task.sh" worker "$run_root")
repeats=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$workers" \
  --account=torch_pr_570_general --partition=cs --array=1-4%4 \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=01:15:00 \
  --export=ALL,EXPECTED_PAIR_LOCK_SHA256="$EXPECTED_PAIR_LOCK_SHA256" \
  --job-name=e5f-nightpair-repeat --output="$run_root/repeat_%A_%a.out" \
  --error="$run_root/repeat_%A_%a.err" "$work/run_task.sh" repeat "$run_root")
export_job=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency="afterany:$repeats" \
  --account=torch_pr_570_general --partition=cs --array=1-2%2 \
  --ntasks=1 --cpus-per-task=1 --mem=16G --time=00:20:00 \
  --export=ALL,EXPECTED_PAIR_LOCK_SHA256="$EXPECTED_PAIR_LOCK_SHA256" \
  --job-name=e5f-nightpair-export --output="$run_root/export_%A_%a.out" \
  --error="$run_root/export_%A_%a.err" "$work/run_task.sh" export "$run_root")
printf 'smoke=%s\nworkers=%s\nrepeats=%s\nexport=%s\nrun_root=%s\n' "$smoke" "$workers" "$repeats" "$export_job" "$run_root"
