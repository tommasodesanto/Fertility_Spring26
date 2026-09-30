#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
runner_dir=$(cd "$(dirname "$0")" && pwd)
# Launcher wall clock includes imports, authentication, solves, observers and serialization.
exec "${PYTHON:-python}" - "$runner_dir" "$@" <<'PY'
import os, signal, subprocess, sys, time
from pathlib import Path
folder=Path(sys.argv[1]); args=sys.argv[2:]
mode=args[0] if args else ''
cap=300 if mode=='preflight' else 2400
started=time.time()
command=[sys.executable,str(folder/'run_comparison.py'),*args,'--deadline-epoch',str(started+2400)]
process=subprocess.Popen(command,start_new_session=True)
try:
    code=process.wait(timeout=cap)
except subprocess.TimeoutExpired:
    os.killpg(process.pid,signal.SIGKILL)
    process.wait()
    raise
sys.exit(code)
PY
