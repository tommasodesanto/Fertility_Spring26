#!/usr/bin/env bash
set -euo pipefail

# Torch-side launcher. It stages each author .do file into a temporary copy so
# source/output roots are portable; each arm runs in a fresh Stata process
# because the author files intentionally begin with clear all.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE=""
OUTROOT=""
STAGE="first"
STATA_BIN="${STATA_BIN:-stata-mp}"
AUTHOR_SOURCE="/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
AUTHOR_OUTPUT="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
KNOWN_SHA="f1c5d48d5ef5357c40e743895dd25c73cd45789f779d1d575408f671fe637029"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

usage() {
  echo "Usage: $0 --source /remote/PSIDSHELF_MOBILITY.dta --out /remote/output --stage first|first_rooms|first_ownership|aligned_first_ownership|second|all [--stata /path/stata-mp]" >&2
  exit 2
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source) SOURCE="$2"; shift 2 ;;
    --out) OUTROOT="$2"; shift 2 ;;
    --stage) STAGE="$2"; shift 2 ;;
    --stata) STATA_BIN="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$SOURCE" && -n "$OUTROOT" ]] || usage
[[ -f "$SOURCE" ]] || { echo "Missing PSID shelf: $SOURCE" >&2; exit 601; }
command -v "$STATA_BIN" >/dev/null 2>&1 || { echo "StataMP unavailable: $STATA_BIN" >&2; exit 127; }
case "$STAGE" in first|first_rooms|first_ownership|aligned_first_ownership|second|all) ;; *) echo "Invalid stage: $STAGE" >&2; exit 2 ;; esac
if [[ -e "$OUTROOT" && -n "$(find "$OUTROOT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Refusing non-fresh output root: $OUTROOT" >&2
  exit 73
fi
mkdir -p "$OUTROOT"
REMOTE_SHA="$(sha256sum "$SOURCE" | awk '{print $1}')"
[[ "$REMOTE_SHA" == "$KNOWN_SHA" ]] || { echo "PSID source checksum mismatch: $REMOTE_SHA" >&2; exit 74; }
MANIFEST="$OUTROOT/run_manifest.csv"
printf '%s\n' 'arm,status,source,weighting,notes' > "$MANIFEST"
printf '%s\n' "source_preflight,passed,$SOURCE,known SHA-256,$(date -u +%FT%TZ)" >> "$MANIFEST"

stage_author_script() {
  local script="$1"
  local variant="${2:-}"
  local sentinel="$3"
  local temp_do
  temp_do="$(mktemp "${TMPDIR:-/tmp}/psid_author_XXXXXX.do")"
  python3 - "$script" "$temp_do" "$AUTHOR_SOURCE" "$SOURCE" "$AUTHOR_OUTPUT" "$OUTROOT" "$PROJECT_DIR" <<'PY'
from pathlib import Path
import sys
src, dst, old_source, new_source, old_output, new_output, project_dir = sys.argv[1:]
text = Path(src).read_text()
text = text.replace(old_source, new_source).replace(old_output, new_output)
text = text.replace("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26", project_dir)
name = Path(src).name
if name == "sa_replication_own_only.do":
    hook = '''

* Portable export hook: preserve the exact unweighted author regression while
* exporting the full covariance and requested +3 minus -1 contrast.
local hook_ip = colnumb(b, "L3event")
local hook_im = colnumb(b, "F1event")
assert !missing(`hook_ip') & !missing(`hook_im') & `hook_ip' > 0 & `hook_im' > 0
local hook_v = var[`hook_ip',`hook_ip'] + var[`hook_im',`hook_im'] - 2*var[`hook_ip',`hook_im']
assert `hook_v' >= 0
preserve
    clear
    svmat2 var, names(col) rnames(row_name)
    export delimited using "`outdir'/own_f_c_y_all_repl_covariance.csv", replace
restore
preserve
    clear
    set obs 1
    gen double estimate = b[1,`hook_ip'] - b[1,`hook_im']
    gen double standard_error = sqrt(`hook_v')
    gen double ci_lo = estimate - 1.96*standard_error
    gen double ci_hi = estimate + 1.96*standard_error
    export delimited using "`outdir'/own_f_c_y_all_repl_contrast.csv", replace
restore
'''
    text = text.replace("matrix var = e(V_iw)", "matrix var = e(V_iw)" + hook, 1)
elif name == "sa_rooms_first_birth_household_aligned_v1.do":
    hook = '''

preserve
    clear
    svmat2 V, names(col) rnames(row_name)
    export delimited using "`outdir'/event_study_covariance.csv", replace
restore
'''
    text = text.replace("matrix V = e(V_iw)", "matrix V = e(V_iw)" + hook, 1)
Path(dst).write_text(text)
PY
  if [[ -n "$variant" ]]; then
    "$STATA_BIN" -b do "$temp_do" "$variant"
  else
    "$STATA_BIN" -b do "$temp_do"
  fi
  rm -f "$temp_do"
  [[ -f "$3" ]] || { echo "Expected completion sentinel missing: $3" >&2; return 75; }
}

run_custom() {
  local arm="$1"
  local variant="${2:-all}"
  "$STATA_BIN" -b do "$SCRIPT_DIR/psid_housing_fullsample_driver.do" "$SOURCE" "$OUTROOT" "$arm" "$variant"
  [[ -f "$OUTROOT/$(if [[ "$arm" == "aligned_first_ownership" ]]; then echo first_birth_aligned_ownership/contrast.csv; else echo second_birth_ownership/${variant}/contrast.csv; fi)" ]] || { echo "Custom completion sentinel missing" >&2; return 76; }
}

if [[ "$STAGE" == first || "$STAGE" == first_rooms || "$STAGE" == all ]]; then
  stage_author_script "$SCRIPT_DIR/sa_rooms_first_birth_household_aligned_v1.do" "" "$OUTROOT/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv"
  [[ -f "$OUTROOT/sa_rooms_first_birth_household_aligned_v1/event_study_covariance.csv" ]] || { echo "First-birth rooms covariance export missing" >&2; exit 77; }
  printf '%s\n' "first_birth_rooms,passed,$SOURCE,IW pweight,exact corrected household-aligned eventstudyinteract" >> "$MANIFEST"
fi
if [[ "$STAGE" == first || "$STAGE" == first_ownership || "$STAGE" == all ]]; then
  stage_author_script "$SCRIPT_DIR/sa_replication_own_only.do" "" "$OUTROOT/sa_replication/own_f_c_y_all_repl_estimates.dta"
  [[ -f "$OUTROOT/sa_replication/own_f_c_y_all_repl_covariance.csv" && -f "$OUTROOT/sa_replication/own_f_c_y_all_repl_contrast.csv" ]] || { echo "First-birth ownership covariance/contrast export missing" >&2; exit 78; }
  printf '%s\n' "first_birth_ownership,passed,$SOURCE,unweighted,exact author command intentionally has no [pw=IW]" >> "$MANIFEST"
fi
if [[ "$STAGE" == first || "$STAGE" == aligned_first_ownership || "$STAGE" == all ]]; then
  run_custom aligned_first_ownership
  printf '%s\n' "first_birth_aligned_ownership,passed,$SOURCE,IW pweight,new corrected HH-year sensitivity; ownership missingness retained" >> "$MANIFEST"
fi
if [[ "$STAGE" == second || "$STAGE" == all ]]; then
  for variant in all no_third_by3 no_third_by3_gap5; do
    stage_author_script "$SCRIPT_DIR/sa_rooms_second_birth_with_onechild_controls_v1.do" "$variant" "$OUTROOT/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_${variant}_summary.csv"
    run_custom second_ownership "$variant"
  done
  printf '%s\n' "second_birth_rooms,passed,$SOURCE,unweighted,exact legacy direct-ACTUALROOMS_ arm; variants all/no_third_by3/no_third_by3_gap5" >> "$MANIFEST"
  printf '%s\n' "second_birth_ownership,passed,$SOURCE,unweighted,new same-clock ownership extension; not author-original" >> "$MANIFEST"
fi
printf '%s\n' 'source_sha256_known,recorded,f1c5d48d5ef5357c40e743895dd25c73cd45789f779d1d575408f671fe637029,,from corrected first-birth audit metadata' >> "$MANIFEST"
echo "PSID full-sample staging complete: $MANIFEST"
