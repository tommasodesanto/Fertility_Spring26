#!/usr/bin/env bash
set -euo pipefail
stage="${1:?stage}"
run_root="${2:?run_root}"
work=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/nightpair_20260925_v1
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/scratch/td2248/commute_pdf_qa_deps${PYTHONPATH:+:$PYTHONPATH}
: "${EXPECTED_PAIR_LOCK_SHA256:?reviewed lock pin required}"
case "$stage" in
  smoke|worker|repeat|export) ;;
  *) echo "invalid stage" >&2; exit 2 ;;
esac
if [[ "$stage" == smoke ]]; then
  task_id="${SLURM_ARRAY_TASK_ID:?array task id}"
  if (( task_id == 1 )); then arm=greaney_179; else arm=oasi_087510; fi
elif [[ "$stage" == worker ]]; then
  task_id="${SLURM_ARRAY_TASK_ID:?array task id}"
  if (( task_id <= 20 )); then arm=greaney_179; slot="$task_id"
  else arm=oasi_087510; slot=$((task_id-20)); fi
elif [[ "$stage" == repeat ]]; then
  task_id="${SLURM_ARRAY_TASK_ID:?array task id}"
  if (( task_id <= 2 )); then arm=greaney_179; rep="$task_id"
  else arm=oasi_087510; rep=$((task_id-2)); fi
else
  task_id="${SLURM_ARRAY_TASK_ID:?array task id}"
  if (( task_id == 1 )); then arm=greaney_179; else arm=oasi_087510; fi
fi
if [[ "$stage" == smoke ]]; then
  remaining=17100
else
  remaining=$(python3 - "$run_root/deadline.json" "$stage" <<'PY'
import json,sys,time
d=json.load(open(sys.argv[1])); stage=sys.argv[2]
key="search_cutoff_epoch" if stage=="worker" else "export_cutoff_epoch" if stage=="repeat" else "deadline_epoch"
reserve=30 if stage=="export" else 0
print(max(0,int(float(d[key])-time.time()-reserve)))
PY
)
fi
if (( remaining <= 0 )); then echo "absolute paired budget exhausted" >&2; exit 124; fi
if [[ "$stage" == smoke ]]; then
  exec timeout --signal=TERM --kill-after=30s "${remaining}s" python3 "$work/run_pair.py" --stage smoke --arm "$arm" --run-root "$run_root"
elif [[ "$stage" == worker ]]; then
  exec timeout --signal=TERM --kill-after=30s "${remaining}s" python3 "$work/run_pair.py" --stage worker --arm "$arm" --slot "$slot" --run-root "$run_root"
elif [[ "$stage" == repeat ]]; then
  exec timeout --signal=TERM --kill-after=30s "${remaining}s" python3 "$work/run_pair.py" --stage repeat --arm "$arm" --repeat-id "$rep" --run-root "$run_root"
fi
exec timeout --signal=TERM --kill-after=30s "${remaining}s" python3 "$work/render_pair.py" --arm "$arm" --run-root "$run_root"
