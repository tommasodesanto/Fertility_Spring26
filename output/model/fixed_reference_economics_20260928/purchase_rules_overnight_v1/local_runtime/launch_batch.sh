#!/usr/bin/env bash
# Detached supervised local ten-chain ramp; do not invoke until lead review.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
: "${ALLOW_LOCAL_CALIBRATION:?Set ALLOW_LOCAL_CALIBRATION=1 only after lead review}"
[[ "$ALLOW_LOCAL_CALIBRATION" == 1 ]] || exit 2
out="$here/runs/local10_v1"
[[ ! -e "$out" ]] || { echo "Refusing existing batch directory: $out" >&2; exit 2; }
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
runner=("$here/../../../../../code/model/.venv/bin/python" "$here/supervisor.py" --out-root "$out")
if command -v caffeinate >/dev/null 2>&1; then
  runner=(caffeinate -i "${runner[@]}")
fi
nohup "${runner[@]}" > "$here/supervisor.log" 2>&1 < /dev/null &
echo "$!" > "$here/supervisor.pid"
echo "supervisor_pid=$! output=$out"
