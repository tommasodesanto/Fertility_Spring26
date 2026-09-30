#!/usr/bin/env bash
set -euo pipefail
mode="${1:-}"
[[ "$mode" == smoke || "$mode" == control ]] || { echo "usage: sbatch .../launch_runtime_validation.sh smoke|control" >&2; exit 2; }
export LAUNCH_STARTED_EPOCH="$(date +%s)"
export CASE_DEADLINE_EPOCH="$(( LAUNCH_STARTED_EPOCH + 300 ))"
export TOTAL_DEADLINE_EPOCH="$(( LAUNCH_STARTED_EPOCH + 900 ))"
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
remote=/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v3
source="/scratch/td2248/projects/fixed_reference_credit_runtime_20260930_v3/source"
results="$remote/results"
[[ -d "$frozen/output/model/fertility_identification_20260928" && -f "$source/source.sha256" ]] || { echo "physical host inputs absent" >&2; exit 2; }
(cd "$source" && sha256sum -c source.sha256)
if [[ "$mode" == control ]]; then
  smoke_receipt="$results/smoke/receipt.json"
  [[ -f "$smoke_receipt" ]] || { echo "lead-reviewable PASS smoke is required before control" >&2; exit 2; }
  /usr/bin/python3 - "$smoke_receipt" <<'PY'
import json,sys
r=json.load(open(sys.argv[1])); assert r["status"]=="passed" and r["lifecycle_solves"]==0
PY
  [[ -f "$results/SMOKE_LEAD_REVIEWED" ]] || { echo "create SMOKE_LEAD_REVIEWED only after lead review" >&2; exit 2; }
fi
[[ ! -e "$results/$mode" ]] || { echo "refusing duplicate result directory" >&2; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR=/work/cache
remaining="$(( CASE_DEADLINE_EPOCH - $(date +%s) ))"
[[ "$remaining" -gt 0 ]] || { echo "case deadline exhausted before container start" >&2; exit 2; }
set +e
timeout --signal=TERM --kill-after=5s "${remaining}s" apptainer exec --bind "$frozen:$repo:ro,$source:/work/source:ro,$results:/work/results:rw,$remote/cache:/work/cache:rw" --pwd "$repo" /share/apps/images/ubuntu-24.04.4.sif \
  /share/apps/anaconda3/2025.06/bin/python /work/source/run_runtime_validation.py "$mode" --plan /work/source/plan.json --output "/work/results/$mode"
rc=$?
set -e
if [[ $rc -ne 0 && ! -f "$results/$mode/failure.json" ]]; then
  mkdir -p "$results/$mode"
  /usr/bin/python3 - "$results/$mode/failure.json" "$mode" "$rc" <<'PY'
import json,sys,time
json.dump({"status":"failed","mode":sys.argv[2],"returncode":int(sys.argv[3]),"time_epoch":time.time(),"no_retry":True},open(sys.argv[1],"w"),indent=2)
PY
fi
exit "$rc"
