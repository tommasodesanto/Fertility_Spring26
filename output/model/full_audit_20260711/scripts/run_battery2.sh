#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/../../../../code/model"
PY=".venv/bin/python"
SCRIPT="../../output/model/full_audit_20260711/scripts/reproduce_candidate.py"
REPORT="../../output/model/combined_recalibration/overnight_20260710_report"
OUT="../../output/model/full_audit_20260711/repro"
run() { env NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" "$PY" "$SCRIPT" "$@" ; }
# housing candidate: local determinism check (repeat of exact config)
run --record "$REPORT/housing_relaxed_best.json" --out "$OUT/housing_repro_repeat.json" --label housing_repro_repeat --nb 120
# housing candidate: tight equilibrium (honest loss)
run --record "$REPORT/housing_relaxed_best.json" --out "$OUT/housing_eq30.json" --label housing_eq30 --nb 120 --max-iter-eq 30
# current candidate: tight equilibrium at Nb=240 (honest verification-grid loss)
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_nb240_eq30.json" --label current_nb240_eq30 --nb 240 --max-iter-eq 30
echo "BATTERY2 DONE"
