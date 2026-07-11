#!/usr/bin/env bash
# Repeatability battery for the July 10 combined-spec candidates (audit-only).
set -uo pipefail
cd "$(dirname "$0")/../../../../code/model"
PY=".venv/bin/python"
SCRIPT="../../output/model/full_audit_20260711/scripts/reproduce_candidate.py"
REPORT="../../output/model/combined_recalibration/overnight_20260710_report"
OUT="../../output/model/full_audit_20260711/repro"
ENVV="NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=$PWD"

run() { env NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH="$PWD" "$PY" "$SCRIPT" "$@" ; }

# 1. current-bound: second fresh-process run (cross-process determinism)
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_repro_repeat.json" --label current_repro_repeat --nb 120
# 2. housing-relaxed: exact reproduction
run --record "$REPORT/housing_relaxed_best.json" --out "$OUT/housing_repro_exact.json" --label housing_repro_exact --nb 120
# 3. current-bound: more equilibrium iterations
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_eq30.json" --label current_eq30 --nb 120 --max-iter-eq 30
# 4. current-bound: tighter tolerance
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_tol25e6.json" --label current_tol25e6 --nb 120 --max-iter-eq 40 --tol-eq 2.5e-5
# 5/6. current-bound: perturbed initial prices around p_eq=0.68976
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_pinit_lo.json" --label current_pinit_lo --nb 120 --p-init 0.655
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_pinit_hi.json" --label current_pinit_hi --nb 120 --p-init 0.725
# 7. current-bound: Nb=240 verification grid
run --record "$REPORT/current_bound_best.json" --out "$OUT/current_nb240.json" --label current_nb240 --nb 240
# 8. housing-relaxed: Nb=240 verification grid
run --record "$REPORT/housing_relaxed_best.json" --out "$OUT/housing_nb240.json" --label housing_nb240 --nb 240
echo "BATTERY DONE"
