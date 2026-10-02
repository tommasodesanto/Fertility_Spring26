#!/usr/bin/env bash
# One local calibration chain; requires explicit launch flag from coordinator.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
: "${ALLOW_LOCAL_CALIBRATION:?Set ALLOW_LOCAL_CALIBRATION=1 only after lead authorizes the single useful run}"
[[ "$ALLOW_LOCAL_CALIBRATION" == 1 ]] || exit 2
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
exec "$here/../../../../../code/model/.venv/bin/python" "$here/worker.py" "$@"
