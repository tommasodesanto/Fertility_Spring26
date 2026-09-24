#!/usr/bin/env bash
set -euo pipefail
stage="${1:?stage required}"
run_root="${2:?run root required}"
bundle=/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/commute_calibration_20260924_v1
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export NUMBA_DISABLE_JIT=0 MPLBACKEND=Agg PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
if [[ "$stage" == smoke ]]; then
  remaining=3600
else
  remaining=$(python3 - "$run_root/deadline.json" "$stage" <<'PY'
import json,sys,time
data=json.load(open(sys.argv[1]));reserve=240 if sys.argv[2]=='worker' else 30
print(max(0,int(data['deadline_epoch']-time.time()-reserve)))
PY
)
fi
if (( remaining <= 0 )); then
  echo "Shared one-hour deadline exhausted before $stage" >&2
  exit 124
fi
if [[ "$stage" == worker ]]; then
  exec timeout --signal=TERM --kill-after=30s "${remaining}s" \
    python3 "$bundle/run_commute.py" --stage worker \
      --worker-id "${SLURM_ARRAY_TASK_ID:?}" --output "$run_root"
fi
exec timeout --signal=TERM --kill-after=30s "${remaining}s" \
  python3 "$bundle/run_commute.py" --stage "$stage" --output "$run_root"
