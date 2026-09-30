#!/usr/bin/env bash
#SBATCH --job-name=small_credit_v1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/small_credit_v1/logs/%x-%j.out
set -euo pipefail
mode="${1:-}"
[[ "$mode" == smoke || "$mode" == full ]] || { echo 'use smoke|full' >&2; exit 2; }
start_epoch="$(date +%s)"
if [[ "$mode" == smoke ]]; then
  deadline_epoch="$((start_epoch + 300))"
else
  deadline_epoch="$((start_epoch + 2400))"
fi
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/small_credit_v1
source="$remote/source"
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
results="$remote/results"
mkdir -p "$results" "$remote/cache" "$remote/logs"
[[ -f "$source/source.sha256" && -f "$inputs/bundle.json" && -f "$inputs/arrays.npz" ]] || {
  echo 'staged source or authenticated bundle missing' >&2; exit 2;
}
(cd "$source" && sha256sum -c source.sha256)
actual_bundle="$(sha256sum "$inputs/bundle.json" | cut -d' ' -f1)"
[[ "$actual_bundle" == 427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7 ]] || {
  echo 'bundle hash differs' >&2; exit 2;
}
# User explicitly authorized immediate production launch without the smoke run.
[[ ! -e "$results/$mode" ]] || { echo 'refusing duplicate output' >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR=/work/cache
export PYTHONPATH=/work/source:/work/source/source
remaining="$((deadline_epoch - $(date +%s)))"
[[ "$remaining" -gt 0 ]] || { echo 'deadline exhausted before container start' >&2; exit 2; }
set +e
timeout --signal=TERM --kill-after=5s "${remaining}s" apptainer exec \
  --bind "$frozen:$repo:ro,$source:/work/source:ro,$inputs:/work/inputs:ro,$results:/work/results:rw,$remote/cache:/work/cache:rw" \
  --pwd "$repo" /share/apps/images/ubuntu-24.04.4.sif \
  /share/apps/anaconda3/2025.06/bin/python /work/source/driver.py "$mode" \
  --reference-root "$repo" --bundle /work/inputs \
  --out "/work/results/$mode" --deadline-epoch "$deadline_epoch"
rc=$?
set -e
if [[ $rc -ne 0 && ! -f "$results/$mode/failure.json" ]]; then
  mkdir -p "$results/$mode"
  /usr/bin/python3 - "$results/$mode/failure.json" "$rc" <<'PY'
import json,sys,time
json.dump({'status':'failed','returncode':int(sys.argv[2]),'time_epoch':time.time(),
           'no_auto_retry':True},open(sys.argv[1],'w'),indent=2)
PY
fi
exit "$rc"
