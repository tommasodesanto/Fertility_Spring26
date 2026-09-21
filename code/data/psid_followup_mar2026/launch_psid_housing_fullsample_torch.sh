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
STATA_PLUS="${STATA_PLUS:-}"
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
STATA_PLUS="${STATA_PLUS:-$OUTROOT/ado}"
[[ -f "$SOURCE" ]] || { echo "Missing PSID shelf: $SOURCE" >&2; exit 601; }
command -v "$STATA_BIN" >/dev/null 2>&1 || { echo "StataMP unavailable: $STATA_BIN" >&2; exit 127; }
case "$STAGE" in first|first_rooms|first_ownership|aligned_first_ownership|second|all) ;; *) echo "Invalid stage: $STAGE" >&2; exit 2 ;; esac
if [[ -e "$OUTROOT" && -n "$(find "$OUTROOT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Refusing non-fresh output root: $OUTROOT" >&2
  exit 73
fi
mkdir -p "$OUTROOT"
mkdir -p "$STATA_PLUS"
REMOTE_SHA="$(sha256sum "$SOURCE" | awk '{print $1}')"
[[ "$REMOTE_SHA" == "$KNOWN_SHA" ]] || { echo "PSID source checksum mismatch: $REMOTE_SHA" >&2; exit 74; }
MANIFEST="$OUTROOT/run_manifest.csv"
printf '%s\n' 'arm,status,source,weighting,notes' > "$MANIFEST"
printf '%s\n' "source_preflight,passed,$SOURCE,known SHA-256,$(date -u +%FT%TZ)" >> "$MANIFEST"

stage_author_script() {
  local script="$1"
  local variant="${2:-}"
  local sentinel="$3"
  local stem
  stem="$(basename "$script" .do)"
  local suffix="${variant:-default}"
  local generated_do="$OUTROOT/generated_${stem}_${suffix}.do"
  local batch_log="$OUTROOT/batch_${stem}_${suffix}.log"
  python3 - "$script" "$generated_do" "$AUTHOR_SOURCE" "$SOURCE" "$AUTHOR_OUTPUT" "$OUTROOT" "$PROJECT_DIR" "$STATA_PLUS" <<'PY'
from pathlib import Path
import sys
src, dst, old_source, new_source, old_output, new_output, project_dir, stata_plus = sys.argv[1:]
text = Path(src).read_text()
name = Path(src).name
if name in {"sa_replication_own_only.do", "sa_rooms_second_birth_with_onechild_controls_v1.do"}:
    source_forms = ["local dta  \"`root'/PSID/PSIDSHELF_MOBILITY.dta\""]
else:
    source_forms = [f'local source  "{old_source}"']
source_hits = sum(text.count(form) for form in source_forms)
if source_hits != 1:
    raise SystemExit(f"expected one anchored source assignment in {name}, found {source_hits}")
for form in source_forms:
    text = text.replace(form, 'local dta  ' + f'"{new_source}"' if 'local dta' in form else 'local source  ' + f'"{new_source}"')
output_forms = [
    'local outroot "`project\'/code/data/psid_followup_mar2026/output"',
    f'local out_root "{old_output}"',
]
output_hits = sum(text.count(form) for form in output_forms)
if output_hits != 1:
    raise SystemExit(f"expected one anchored output assignment in {name}, found {output_hits}")
for form in output_forms:
    if form in text:
        text = text.replace(form, form.split('"')[0] + f'"{new_output}"')
text = text.replace("/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26", project_dir)
text = text.replace('local root "/Users/tommasodesanto/Desktop/Projects/Fertility"', 'local root "' + str(Path(project_dir).parent) + '"')
text = text.replace("/Users/tommasodesanto/Desktop/Projects/Fertility", str(Path(project_dir).parent))
startup_anchor = "clear all\n"
if text.count(startup_anchor) != 1:
    raise SystemExit(f"expected exactly one clear-all startup anchor in {name}")
text = text.replace(startup_anchor, startup_anchor + f'sysdir set PLUS "{stata_plus}"\nmata: mata mlib index\n', 1)
# The author template requests eight processors before opening its log.  Keep
# that template untouched, but cap only the staged copy to one processor so a
# small allocation or license cannot fail before any diagnostic is written.
processor_lines = [line for line in text.splitlines() if line.strip().startswith("set processors ")]
if len(processor_lines) > 1:
    raise SystemExit(f"multiple processor settings in {name}")
if processor_lines:
    text = text.replace(processor_lines[0], "set processors 1", 1)
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
    hook_anchor = "matrix var = e(V_iw)"
    if text.count(hook_anchor) != 1:
        raise SystemExit(f"expected exactly one covariance hook in {name}")
    text = text.replace(hook_anchor, hook_anchor + hook, 1)
elif name == "sa_rooms_first_birth_household_aligned_v1.do":
    hook = '''

preserve
    clear
    svmat2 V, names(col) rnames(row_name)
    export delimited using "`outdir'/event_study_covariance.csv", replace
restore
'''
    hook_anchor = "matrix V = e(V_iw)"
    if text.count(hook_anchor) != 1:
        raise SystemExit(f"expected exactly one covariance hook in {name}")
    text = text.replace(hook_anchor, hook_anchor + hook, 1)
elif name == "sa_rooms_second_birth_with_onechild_controls_v1.do":
    hook = '''

* Portable export hook: preserve the exact unweighted legacy rooms regression
* while exporting full covariance and the requested +3 minus -1 contrast.
local hook_ip = colnumb(var, "L3event")
local hook_im = colnumb(var, "F1event")
local hook_bp = colnumb(b, "L3event")
local hook_bm = colnumb(b, "F1event")
assert !missing(`hook_ip') & !missing(`hook_im') & `hook_ip' > 0 & `hook_im' > 0
assert `hook_bp' == `hook_ip' & `hook_bm' == `hook_im'
local hook_v = var[`hook_ip',`hook_ip'] + var[`hook_im',`hook_im'] - 2*var[`hook_ip',`hook_im']
assert `hook_v' >= 0 & `hook_v' < .
preserve
    clear
    svmat2 var, names(col) rnames(row_name)
    export delimited using "`outdir'/rooms_s_c_y_`variant'_covariance.csv", replace
restore
preserve
    clear
    set obs 1
    gen double estimate = b[1,`hook_bp'] - b[1,`hook_bm']
    gen double standard_error = sqrt(`hook_v')
    gen double ci_lo = estimate - 1.96*standard_error
    gen double ci_hi = estimate + 1.96*standard_error
    export delimited using "`outdir'/rooms_s_c_y_`variant'_contrast.csv", replace
restore
'''
    hook_anchor = "matrix var = e(V_iw)"
    if text.count(hook_anchor) != 1:
        raise SystemExit(f"expected exactly one covariance hook in {name}")
    text = text.replace(hook_anchor, hook_anchor + hook, 1)
completion_dir = {
    "sa_rooms_first_birth_household_aligned_v1.do": Path(new_output) / "sa_rooms_first_birth_household_aligned_v1",
    "sa_replication_own_only.do": Path(new_output) / "sa_replication",
    "sa_rooms_second_birth_with_onechild_controls_v1.do": Path(new_output) / "sa_rooms_second_birth_with_onechild_controls_v1",
}[name]
completion_path = completion_dir / "STATA_COMPLETE"
if name == "sa_rooms_second_birth_with_onechild_controls_v1.do":
    completion_path = completion_dir / "STATA_COMPLETE_`variant'"
text += f'''\n\ncapture file close _psid_complete
file open _psid_complete using "{completion_path}", write replace
file write _psid_complete "PASS: generated PSID arm completed" _n
file close _psid_complete
'''
Path(dst).write_text(text)
PY
  if [[ -n "$variant" ]]; then
    set +e
    "$STATA_BIN" -b do "$generated_do" "$variant" >"$batch_log" 2>&1
    stata_rc=$?
    set -e
  else
    set +e
    "$STATA_BIN" -b do "$generated_do" >"$batch_log" 2>&1
    stata_rc=$?
    set -e
  fi
  if [[ "$stata_rc" -ne 0 ]]; then
    cat "$batch_log" >&2
    return "$stata_rc"
  fi
  [[ -f "$3" ]] || { echo "Expected completion sentinel missing: $3" >&2; return 75; }
  case "$(basename "$script")" in
    sa_rooms_first_birth_household_aligned_v1.do) marker="$OUTROOT/sa_rooms_first_birth_household_aligned_v1/STATA_COMPLETE" ;;
    sa_replication_own_only.do) marker="$OUTROOT/sa_replication/STATA_COMPLETE" ;;
    sa_rooms_second_birth_with_onechild_controls_v1.do) marker="$OUTROOT/sa_rooms_second_birth_with_onechild_controls_v1/STATA_COMPLETE_${variant}" ;;
  esac
  [[ -f "$marker" ]] || { echo "Generated Stata completion marker missing: $marker" >&2; return 79; }
}

run_custom() {
  local arm="$1"
  local variant="${2:-all}"
  local temp_do
  temp_do="$(mktemp "${TMPDIR:-/tmp}/psid_driver_XXXXXX.do")"
  python3 - "$SCRIPT_DIR/psid_housing_fullsample_driver.do" "$temp_do" "$STATA_PLUS" <<'PY'
from pathlib import Path
import sys
src, dst, stata_plus = sys.argv[1:]
text = Path(src).read_text()
startup_anchor = "clear all\n"
if text.count(startup_anchor) != 1:
    raise SystemExit("expected exactly one clear-all startup anchor in custom driver")
Path(dst).write_text(text.replace(startup_anchor, startup_anchor + f'sysdir set PLUS "{stata_plus}"\nmata: mata mlib index\n', 1))
PY
  "$STATA_BIN" -b do "$temp_do" "$SOURCE" "$OUTROOT" "$arm" "$variant" "$STATA_PLUS"
  rm -f "$temp_do"
  [[ -f "$OUTROOT/$(if [[ "$arm" == "aligned_first_ownership" ]]; then echo first_birth_aligned_ownership/contrast.csv; else echo second_birth_ownership/${variant}/contrast.csv; fi)" ]] || { echo "Custom completion sentinel missing" >&2; return 76; }
}

if [[ "$STAGE" == first || "$STAGE" == first_rooms || "$STAGE" == all ]]; then
  stage_author_script "$SCRIPT_DIR/sa_rooms_first_birth_household_aligned_v1.do" "" "$OUTROOT/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv"
  [[ -f "$OUTROOT/sa_rooms_first_birth_household_aligned_v1/event_study_covariance.csv" ]] || { echo "First-birth rooms covariance export missing" >&2; exit 77; }
  printf '%s\n' "first_birth_rooms,passed,$SOURCE,IW pweight,exact corrected household-aligned eventstudyinteract" >> "$MANIFEST"
fi
if [[ "$STAGE" == first || "$STAGE" == first_ownership || "$STAGE" == all ]]; then
  stage_author_script "$SCRIPT_DIR/sa_replication_own_only.do" "" "$OUTROOT/sa_replication/own_f_c_y_all_repl_estimates.dta"
  [[ -f "$OUTROOT/sa_replication/own_f_c_y_all_repl_covariance.csv" && -f "$OUTROOT/sa_replication/own_f_c_y_all_repl_contrast.csv" && -f "$OUTROOT/sa_replication/own_f_c_y_all_repl.png" ]] || { echo "First-birth ownership covariance/contrast/graph export missing" >&2; exit 78; }
  printf '%s\n' "first_birth_ownership,passed,$SOURCE,unweighted,exact author command intentionally has no [pw=IW]" >> "$MANIFEST"
fi
if [[ "$STAGE" == first || "$STAGE" == aligned_first_ownership || "$STAGE" == all ]]; then
  run_custom aligned_first_ownership
  printf '%s\n' "first_birth_aligned_ownership,passed,$SOURCE,IW pweight,new corrected HH-year sensitivity; ownership missingness retained" >> "$MANIFEST"
fi
if [[ "$STAGE" == second || "$STAGE" == all ]]; then
  for variant in all no_third_by3 no_third_by3_gap5; do
    stage_author_script "$SCRIPT_DIR/sa_rooms_second_birth_with_onechild_controls_v1.do" "$variant" "$OUTROOT/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_${variant}_summary.csv"
    [[ -f "$OUTROOT/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_${variant}_covariance.csv" && -f "$OUTROOT/sa_rooms_second_birth_with_onechild_controls_v1/rooms_s_c_y_${variant}_contrast.csv" ]] || { echo "Second-birth rooms covariance/contrast export missing" >&2; exit 80; }
    run_custom second_ownership "$variant"
  done
  printf '%s\n' "second_birth_rooms,passed,$SOURCE,unweighted,exact legacy direct-ACTUALROOMS_ arm; variants all/no_third_by3/no_third_by3_gap5" >> "$MANIFEST"
  printf '%s\n' "second_birth_ownership,passed,$SOURCE,unweighted,corrected ownership extension; includes F6 unlike legacy rooms arm; not author-original" >> "$MANIFEST"
fi
printf '%s\n' 'source_sha256_known,recorded,f1c5d48d5ef5357c40e743895dd25c73cd45789f779d1d575408f671fe637029,,from corrected first-birth audit metadata' >> "$MANIFEST"
echo "PSID full-sample staging complete: $MANIFEST"
