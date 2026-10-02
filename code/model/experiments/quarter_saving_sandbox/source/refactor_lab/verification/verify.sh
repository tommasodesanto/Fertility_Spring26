#!/usr/bin/env bash
# Refactor-lab phase driver. MODE=local (this Mac, explicit PY=<NumPy-2 interpreter>,
# one process at a time) or MODE=torch (sbatch + apptainer). Lives in refactor_lab/verification/; see README.
#
#   MODE=local PHASE=old-fixed-price             ORIGINAL engine baseline on this machine (1 solve)
#   MODE=local PHASE=fixed-price BUNDLE_SHA=.. [COMPARISON_REFERENCE=<baseline dir> COMPARISON_PIN=<receipt sha>]
#
#   PHASE=smoke                                  no solve: sequencing, fail-fast, timeout paths (DRY_RUN)
#   PHASE=export                                 input bundle (written once)
#   PHASE=acceptance  BUNDLE_SHA=..              component tests                     cap 600 s
#   PHASE=fixed-price BUNDLE_SHA=..              2 lab FP solves + compare + oracle   cap 900 s
#   PHASE=oracle-fp   FP_RESULT=<results dir>    oracle only on saved FP artifacts    cap 900 s
#   PHASE=ge          BUNDLE_SHA=.. PRICE_FACTOR=1.05   lab GE, old GE, certify both, compare   cap 2700 s
#   PHASE=credit-d0                              refused here (other chat owns it; see README)
#
# Every step is capped by the remaining phase budget (and GE solve steps by
# ENGINE_CAP). A failed or timed-out step stops the phase immediately; a
# timeout is failure, never a partial certificate. No automatic retries.
# progress.json is rewritten every 300 s and after each step;
# latest_completed.json and summary.json after each step.
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:50:00
#SBATCH --job-name=refactor_lab
set -euo pipefail
need_var() { local v; for v in "$@"; do [ -n "${!v:-}" ] || { echo "required variable $v is unset" >&2; exit 2; }; done; }
need_var PHASE
MODE="${MODE:-torch}"
CKPT_SHA=b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d
MAX_RSS_GIB=12
if [ "$MODE" = local ]; then
  ROOT="${ROOT:-$(cd "$(dirname "$0")/../../../.." && pwd)}"
  LAB="${LAB:-$ROOT/output/model/publication_refactor_20260929/local_runs}"
  REF="$ROOT"; ORIG="$ROOT"      # reference artifacts are read only by convention; nothing writes under them
  CKPT_REL=output/model/fertility_identification_20260928/resume_v1/selected_export/primary/initial_state.pkl.gz
  SRC="${LAB_SRC:-$ROOT/code/model}"   # directory containing refactor_lab/ (e.g. an indexed candidate stage)
  PY="${PY:-}"                          # REQUIRED: explicit, separately verified NumPy-2-capable interpreter
  TESTS="$SRC/refactor_lab/tests"
  OVERLAY_ROOT="$ROOT/output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/overlay"
else
  LAB="${LAB:-/scratch/td2248/projects/publication_refactor_20260929}"
  REF=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
  ORIG=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26; ROOT="$ORIG"
  CKPT_REL=output/model/fertility_identification_20260928/resume_v1/repeat_0212_primary/case/initial_state.pkl.gz
  SRC="${LAB_SRC:-$LAB/lab_src}"; PY=/share/apps/anaconda3/2025.06/bin/python; OVERLAY_ROOT=""; TESTS=/lab/refactor_lab/tests
fi
EXPORT="${EXPORT_DIR:-$LAB/export_v1}"
CACHE_DIR="$LAB/numba_cache"
DRY_RUN="${DRY_RUN:-0}"
JOB="${SLURM_JOB_ID:-local$$}"
OUT="${OUT_OVERRIDE:-$LAB/results/${PHASE}_${JOB}}"
[ "$MODE" = local ] && need_var PY
now_s() { perl -MTime::HiRes=time -e 'printf "%.3f\n", time'; }   # portable (BSD date has no %N)
# Destinations may never lie inside frozen reference/source subtrees (reads only there).
for dest in "$OUT" "$EXPORT" "$LAB"; do
  python3 - "$ROOT" "$dest" <<'PY_GUARD' || { echo "refused destination inside frozen reference subtree: $dest" >&2; exit 2; }
import os, sys
root, dest = sys.argv[1], os.path.realpath(sys.argv[2])
frozen = ["output/model/fertility_identification_20260928", "output/model/fixed_reference_economics_20260928",
          "output/model/overnight_calibration_20260928", "code/model/intergen_eqscale_seq_optimized", "code/model/tools", "tmp"]
sys.exit(1 if any(dest == os.path.realpath(os.path.join(root, f)) or dest.startswith(os.path.realpath(os.path.join(root, f)) + os.sep) for f in frozen) else 0)
PY_GUARD
done
if [ -e "$OUT" ]; then echo "output directory already exists (fresh output required): $OUT" >&2; exit 2; fi
SOLVE_DEADLINE=900       # enforced IN-PROCESS on the engine solve stage (CallCounter), separate from import/output
IMPORT_OUTPUT_RESERVE=300
ENGINE_CAP=$((SOLVE_DEADLINE + IMPORT_OUTPUT_RESERVE))   # external backstop for one engine step
MAX_LIFECYCLE=18         # enforced BEFORE at-price call 19 starts (exit 4 = failed_budget)
case "$PHASE" in
  acceptance) CAP=600 ;; fixed-price|oracle-fp|old-fixed-price) CAP=900 ;; ge) CAP=2700 ;; export) CAP=1200 ;; *) CAP=600 ;;
esac
[ "$DRY_RUN" = 1 ] && CAP="${SMOKE_CAP:-$CAP}"
mkdir -p "$OUT"
T_START=$(date +%s)
if [ "$MODE" = local ] && [ -x "${PY:-}" ]; then
  "$PY" -c "import os,sys,json;print(json.dumps(dict(py_resolved=os.path.realpath(sys.executable))))" > "$OUT/interpreter.json" 2>&1 || true
fi
json_str() { printf '"%s"' "$(printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g')"; }
write_plan() {
  cat > "$OUT/plan.json" <<JSON
{"phase": $(json_str "$PHASE"), "cap_seconds": $CAP, "engine_cap_seconds": $ENGINE_CAP,
 "max_lifecycle_evaluations_per_engine": $MAX_LIFECYCLE, "planned_steps": $(json_str "$1"),
 "estimate": $(json_str "$2"), "stop": "first failed or timed-out step stops the phase; no retries; timeout is failure",
 "source_stage": $(json_str "$SRC"), "dry_run": $DRY_RUN}
JSON
}
CURRENT="start"; COMPLETED=""
progress() {
  local now; now=$(date +%s)
  printf '{"phase": %s, "current_step": %s, "elapsed_seconds": %s, "cap_seconds": %s, "completed": [%s], "status": %s}\n' \
    "$(json_str "$PHASE")" "$(json_str "$CURRENT")" "$((now - T_START))" "$CAP" "$COMPLETED" "$(json_str "$1")" > "$OUT/progress.json"
}
( while sleep "${HEARTBEAT_SECONDS:-300}"; do [ -f "$OUT/.done" ] && exit 0
    printf '{"heartbeat_epoch": %s, "elapsed_seconds": %s, "current_step": "%s"}\n' "$(date +%s)" \
      "$(( $(date +%s) - T_START ))" "$(cat "$OUT/.current" 2>/dev/null)" > "$OUT/progress_heartbeat.json"; done ) </dev/null >/dev/null 2>&1 &
HEARTBEAT=$!
trap 'rc=$?; touch "$OUT/.done"; kill $HEARTBEAT 2>/dev/null; exit $rc' EXIT   # preserve exit status
sha_check() { if command -v sha256sum >/dev/null; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
if [ "$DRY_RUN" != 1 ]; then
  if [ "$MODE" = torch ]; then (cd "$SRC" && sha256sum -c --quiet SHA256SUMS); module load anaconda3/2025.06
  else (cd "$SRC" && find refactor_lab \( -name "*.py" -o -name "*.json" -o -name "*.sh" \) | sort | xargs shasum -a 256) > "$OUT/source_sha256.txt"; fi
  echo "$SRC" > "$OUT/source_stage.txt"
  echo "$CKPT_SHA  $REF/$CKPT_REL" | sha_check -c --quiet -
  env PYTHONPATH="$SRC" "$PY" -c "import json;from refactor_lab.verification import baseline_identity as b;print(json.dumps(b.runtime_identity(),indent=1))" > "$OUT/runtime_identity.json"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
py() {
  if [ "$DRY_RUN" = 1 ]; then   # smoke: no container, no import, no solve
    if [ "$MODE" = local ]; then echo "DRY: PYTHONPATH=$SRC NUMBA_CACHE_DIR=$CACHE_DIR REFACTOR_ROOT=$ROOT $PY $*"; else echo "DRY: $*"; fi
    [ "${DRY_FAIL_STEP:-}" = "$CURRENT" ] && return "${DRY_FAIL_RC:-1}"; sleep "${DRY_SLEEP:-0}"; return 0
  fi
  if [ "$MODE" = local ]; then
    env PYTHONPATH="$SRC" NUMBA_CACHE_DIR="$CACHE_DIR" REFACTOR_ROOT="$ROOT" REFACTOR_SOURCE_ROOT="$ROOT" \
      REFACTOR_OVERLAY_ROOT="$OVERLAY_ROOT" REFACTOR_EXPORT="$EXPORT" REFACTOR_BUNDLE_SHA="${BUNDLE_SHA:-}" \
      REFACTOR_CHECKPOINT="$ROOT/$CKPT_REL" "$PY" "$@"
    return
  fi
  apptainer exec --bind "$REF:$ORIG:ro" --bind "$SRC:/lab:ro" --bind "$LAB:$LAB" --pwd /lab \
    --env PYTHONPATH=/lab,NUMBA_CACHE_DIR="$CACHE_DIR",REFACTOR_ROOT="$ORIG",REFACTOR_SOURCE_ROOT="$ORIG",REFACTOR_EXPORT="$EXPORT",REFACTOR_BUNDLE_SHA="${BUNDLE_SHA:-}",REFACTOR_CHECKPOINT="$ORIG/$CKPT_REL" \
    /share/apps/images/ubuntu-24.04.4.sif "$PY" "$@"
}
export -f py
export DRY_RUN REF ORIG SRC LAB CKPT_REL OUT MODE ROOT PY OVERLAY_ROOT EXPORT TESTS
capped() {   # local: budget_run.py (wall + tree RSS <= 12 GiB, own session); torch: coreutils timeout
  local s=$1; shift
  if [ "$MODE" = local ] && [ "$DRY_RUN" != 1 ]; then
    python3 "$SRC/refactor_lab/verification/budget_run.py" --seconds "$s" --max-rss-gib "$MAX_RSS_GIB" \
      --report "$OUT/$CURRENT.budget.json" -- bash -c '"$@"' _ "$@"
  elif command -v timeout >/dev/null; then timeout -k 30 "$s" bash -c '"$@"' _ "$@"
  else perl -e 'alarm shift; exec @ARGV' "$s" bash -c '"$@"' _ "$@"; fi
}
# step NAME [PER_STEP_CAP] -- CMD...: runs within the remaining phase budget; stops the phase on failure.
step() {
  local name=$1 limit=$2; shift 3
  CURRENT=$name; CACHE_DIR="${STEP_CACHE:-$LAB/numba_cache}"
  if [ -n "${STEP_CACHE:-}" ] && [ -n "$(ls -A "$STEP_CACHE" 2>/dev/null)" ] && [ "${WARM_CACHE_OK:-0}" != 1 ]; then
    echo "cold-comparison cache is not fresh: $STEP_CACHE (set WARM_CACHE_OK=1 only for a labelled warm run)" >&2; exit 2; fi
  mkdir -p "$CACHE_DIR"
  printf '%s\tbefore\t%s\n' "$name" "$(find "$CACHE_DIR" -type f 2>/dev/null | wc -l | tr -d ' ')" >> "$OUT/cache_state.tsv"
  export CURRENT CACHE_DIR BUNDLE_SHA DRY_FAIL_STEP DRY_SLEEP DRY_FAIL_RC; echo "$name" > "$OUT/.current"; progress running
  local remaining=$((CAP - ($(date +%s) - T_START)))
  [ "$limit" -gt 0 ] && [ "$limit" -lt "$remaining" ] && remaining=$limit
  if [ "$remaining" -le 0 ]; then echo "phase cap reached before $name" >&2; progress failed_timeout; exit 124; fi
  local t0 t1 rc; t0=$(now_s)
  if capped "$remaining" "$@" > "$OUT/$name.log" 2>&1; then rc=0; else rc=$?; fi
  t1=$(now_s)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$t0" "$t1" "$(perl -e "printf '%.3f', $t1-$t0")" "$rc" "$remaining" >> "$OUT/steps.tsv"
  printf '%s\tafter\t%s\n' "$name" "$(find "$CACHE_DIR" -type f 2>/dev/null | wc -l | tr -d ' ')" >> "$OUT/cache_state.tsv"
  COMPLETED="$COMPLETED${COMPLETED:+, }{\"step\": $(json_str "$name"), \"rc\": $rc}"
  printf '{"step": %s, "rc": %s, "cap_seconds": %s}\n' "$(json_str "$name")" "$rc" "$remaining" > "$OUT/latest_completed.json"
  printf '{"phase": %s, "steps": [%s], "elapsed_seconds": %s}\n' "$(json_str "$PHASE")" "$COMPLETED" "$(( $(date +%s) - T_START ))" > "$OUT/summary.json"
  if [ "$rc" -ne 0 ]; then
    if [ "$rc" -eq 124 ] || [ "$rc" -eq 142 ]; then progress failed_timeout
    elif [ "$rc" -eq 4 ]; then progress failed_budget
    elif [ "$rc" -eq 125 ]; then progress failed_memory; else progress failed; fi
    echo "step $name failed rc=$rc; phase stopped" >&2; exit "$rc"
  fi
  progress running
}
BUNDLE="$EXPORT/inputs"
INPUT_ARGS=(--reference-root "$ORIG" --bundle "$BUNDLE" --bundle-sha "${BUNDLE_SHA:-}")
case "$PHASE" in
  export)
    write_plan "export" "minutes; one checkpoint read"
    step export 0 -- py -m refactor_lab.export_inputs --reference-root "$ORIG" --checkpoint "$ORIG/$CKPT_REL" \
         --oracle-code "$ORIG/code/model" --out "$EXPORT"
    [ "$DRY_RUN" = 1 ] || sha_check "$BUNDLE/bundle.json" | tee "$OUT/bundle_sha.txt" ;;
  acceptance)
    need_var BUNDLE_SHA
    write_plan "pytest" "under 10 min; no solve"
    step pytest 0 -- py -m pytest "$TESTS" -q -p no:cacheprovider --junitxml="$OUT/acceptance.xml" ;;
  fixed-price)
    need_var BUNDLE_SHA
    write_plan "lab_fixed_price(2 solves), compare_rep1, compare_rep2, compare_reps, oracle(2 certificates)" \
               "about 2x107 s observed solves plus JIT/serialization; oracle time unmeasured (uncertain)"
    STEP_CACHE="$OUT/numba_cache_lab" step lab_fixed_price 0 -- py -m refactor_lab.run fixed-price --credit reference --repeat 2 "${INPUT_ARGS[@]}" --out "$OUT/lab"
    [ "$MODE" = local ] && need_var COMPARISON_REFERENCE COMPARISON_PIN   # local exactness only vs a pinned same-machine baseline
    if [ -n "${COMPARISON_REFERENCE:-}" ]; then   # same-machine gate first; Torch checkpoint = diagnostic
      need_var COMPARISON_PIN
      for k in 1 2; do
        step compare_rep${k}_vs_baseline 0 -- py -m refactor_lab.verification.compare --lab "$OUT/lab/rep$k/solution_arrays.npz" \
             --baseline-dir "$COMPARISON_REFERENCE" --baseline-pin "$COMPARISON_PIN" --root "$ORIG" --out "$OUT/comparison_rep${k}_vs_baseline.json"
        step compare_rep${k}_vs_torch_diag 0 -- py -m refactor_lab.verification.compare --record-only --lab "$OUT/lab/rep$k/solution_arrays.npz" \
             --reference "$EXPORT/verification/reference_solution.npz" --out "$OUT/diagnostic_rep${k}_vs_torch.json"
      done
    else
      for k in 1 2; do
        step compare_rep$k 0 -- py -m refactor_lab.verification.compare --lab "$OUT/lab/rep$k/solution_arrays.npz" \
             --reference "$EXPORT/verification/reference_solution.npz" --out "$OUT/comparison_rep$k.json"
      done
    fi
    step compare_reps 0 -- py -m refactor_lab.verification.compare --lab "$OUT/lab/rep2/solution_arrays.npz" \
         --reference "$OUT/lab/rep1/solution_arrays.npz" --out "$OUT/comparison_rep2_vs_rep1.json"
    CMP=(); [ -n "${COMPARISON_REFERENCE:-}" ] && { need_var COMPARISON_PIN; CMP=(--comparison-reference "$COMPARISON_REFERENCE" --comparison-pin "$COMPARISON_PIN"); }
    STEP_CACHE="$OUT/numba_cache_oracle" step oracle 0 -- py -m refactor_lab.verification.acceptance_oracle fixed-price --root "$ORIG" ${CMP[@]+"${CMP[@]}"} \
         --lab-solution "$OUT/lab/rep1/acceptance_solution.pkl.gz" "$OUT/lab/rep2/acceptance_solution.pkl.gz" --out "$OUT/oracle" ;;
  old-fixed-price)
    write_plan "old_fixed_price (ORIGINAL engine, 1 solve, gates, 14/31, 17 plots, census vs Torch checkpoint)" \
               "one lifecycle (~65-75 s on Torch; this machine unmeasured) plus oracle"
    STEP_CACHE="$OUT/numba_cache_old" step old_fixed_price 0 -- py -m refactor_lab.verification.acceptance_oracle old-fixed-price --root "$ORIG" --out "$OUT/old"
    [ "$DRY_RUN" = 1 ] || sha_check "$OUT/old/baseline_receipt.json" | tee "$OUT/baseline_receipt_sha.txt"
    echo "COMPARISON_REFERENCE=$OUT/old" | tee -a "$OUT/baseline_receipt_sha.txt" ;;
  oracle-fp)
    need_var FP_RESULT
    write_plan "oracle(2 certificates) on saved $FP_RESULT" "oracle only; no lifecycle solve"
    step oracle 0 -- py -m refactor_lab.verification.acceptance_oracle fixed-price --root "$ORIG" \
         --lab-solution "$FP_RESULT/lab/rep1/acceptance_solution.pkl.gz" "$FP_RESULT/lab/rep2/acceptance_solution.pkl.gz" --out "$OUT/oracle" ;;
  ge)
    need_var BUNDLE_SHA PRICE_FACTOR RENEWAL_TOL   # RENEWAL_TOL: lead-supplied (live transition value 1e-6)
    [ "$PRICE_FACTOR" = "1.05" ] || { echo "lead-selected start is 1.05" >&2; exit 2; }
    BUDGET=(--count-calls --max-lifecycle "$MAX_LIFECYCLE" --solve-deadline "$SOLVE_DEADLINE")
    write_plan "lab_ge, lab_ge_certify, old_ge(+certify, lab/old parameter identity), compare_lab_old" \
               "UNCERTAIN: ~107 s per full lifecycle observed (fixed price); historical 886.5 s GE (6 lifecycle calls) used a different credit contract; each engine hard-stopped at 18 at-price calls or 900 s of solve stage (profiled, same instrumentation both engines)"
    STEP_CACHE="$OUT/numba_cache_lab" step lab_ge "$ENGINE_CAP" -- py -m refactor_lab.run equilibrium --credit reference "${BUDGET[@]}" \
         --initial-price-factor "$PRICE_FACTOR" "${INPUT_ARGS[@]}" --out "$OUT/lab"
    STEP_CACHE="$OUT/numba_cache_oracle" step lab_ge_certify 0 -- py -m refactor_lab.verification.acceptance_oracle ge-certify --root "$ORIG" --renewal-tolerance "$RENEWAL_TOL" \
         --lab-solution "$OUT/lab/ge/acceptance_solution.pkl.gz" --out "$OUT/lab_certificate"
    STEP_CACHE="$OUT/numba_cache_old" step old_ge "$ENGINE_CAP" -- py -m refactor_lab.verification.acceptance_oracle ge-old --root "$ORIG" "${BUDGET[@]}" \
         --renewal-tolerance "$RENEWAL_TOL" --price-factor "$PRICE_FACTOR" \
         --lab-solution "$OUT/lab/ge/acceptance_solution.pkl.gz" --out "$OUT/old"
    step compare_lab_old 0 -- py -m refactor_lab.verification.compare --lab "$OUT/lab/ge/solution_arrays.npz" \
         --reference "$OUT/old/solution_arrays.npz" --out "$OUT/comparison_lab_vs_old.json" ;;
  smoke)
    # Exercises this driver's real sequencing with DRY_RUN stubs: pass, fail-fast, timeout.
    self="$0"; base="$OUT"; export MODE
    if [ "$MODE" = local ]; then   # local fixed-price must take the pinned-baseline branch (fake stub values)
      FPENV=(COMPARISON_REFERENCE=/nonexistent/fake_baseline COMPARISON_PIN=fakepin)
      FP_STEPS="lab_fixed_price compare_rep1_vs_baseline compare_rep1_vs_torch_diag compare_rep2_vs_baseline compare_rep2_vs_torch_diag compare_reps oracle "
      FP_FAIL=compare_rep1_vs_baseline
      if OUT_OVERRIDE="$base/fp_nobaseline" DRY_RUN=1 PHASE=fixed-price BUNDLE_SHA=x LAB="$LAB" bash "$self" 2>/dev/null; then
        echo "local fixed-price without baseline did not fail" >&2; exit 1; fi
    else
      FPENV=(SMOKE_TORCH=1)
      FP_STEPS="lab_fixed_price compare_rep1 compare_rep2 compare_reps oracle "
      FP_FAIL=compare_rep1
    fi
    env "${FPENV[@]}" OUT_OVERRIDE="$base/pass" DRY_RUN=1 PHASE=fixed-price BUNDLE_SHA=x LAB="$LAB" bash "$self"
    [ "$(cut -f1 "$base/pass/steps.tsv" | tr '\n' ' ')" = "$FP_STEPS" ]
    if env "${FPENV[@]}" OUT_OVERRIDE="$base/fail" DRY_RUN=1 DRY_FAIL_STEP=$FP_FAIL PHASE=fixed-price BUNDLE_SHA=x LAB="$LAB" bash "$self"; then
      echo "fail-fast smoke did not fail" >&2; exit 1; fi
    [ "$(cut -f1 "$base/fail/steps.tsv" | tr '\n' ' ')" = "lab_fixed_price $FP_FAIL " ]
    if OUT_OVERRIDE="$base/timeout" DRY_RUN=1 DRY_SLEEP=3 SMOKE_CAP=1 PHASE=acceptance BUNDLE_SHA=x LAB="$LAB" bash "$self"; then
      echo "timeout smoke did not fail" >&2; exit 1; fi
    grep -q failed_timeout "$base/timeout/progress.json"
    OUT_OVERRIDE="$base/ge_pass" DRY_RUN=1 PHASE=ge BUNDLE_SHA=x PRICE_FACTOR=1.05 RENEWAL_TOL=1e-6 LAB="$LAB" bash "$self"
    [ "$(cut -f1 "$base/ge_pass/steps.tsv" | tr '\n' ' ')" = "lab_ge lab_ge_certify old_ge compare_lab_old " ]
    if OUT_OVERRIDE="$base/ge_budget" DRY_RUN=1 DRY_FAIL_STEP=old_ge DRY_FAIL_RC=4 PHASE=ge BUNDLE_SHA=x PRICE_FACTOR=1.05 \
         RENEWAL_TOL=1e-6 LAB="$LAB" bash "$self"; then echo "budget smoke did not fail" >&2; exit 1; fi
    [ "$(cut -f1 "$base/ge_budget/steps.tsv" | tr '\n' ' ')" = "lab_ge lab_ge_certify old_ge " ]
    grep -q failed_budget "$base/ge_budget/progress.json"
    if OUT_OVERRIDE="$base/ge_factor" DRY_RUN=1 PHASE=ge BUNDLE_SHA=x PRICE_FACTOR=1.10 RENEWAL_TOL=1e-6 LAB="$LAB" bash "$self"; then
      echo "start-factor guard did not fail" >&2; exit 1; fi
    if OUT_OVERRIDE="$base/ge_notol" DRY_RUN=1 PHASE=ge BUNDLE_SHA=x PRICE_FACTOR=1.05 LAB="$LAB" bash "$self" 2>/dev/null; then
      echo "missing renewal tolerance did not fail" >&2; exit 1; fi
    echo "smoke passed" | tee "$OUT/smoke_result.txt" ;;
  credit-d0)
    echo "Not run here: corrected D=0 runtime is owned by the other chat; run.py keeps the census path." >&2; exit 2 ;;
  *) echo "unknown PHASE" >&2; exit 2 ;;
esac
progress passed
