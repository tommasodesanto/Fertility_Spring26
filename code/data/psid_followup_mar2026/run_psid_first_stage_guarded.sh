#!/usr/bin/env bash
set -euo pipefail

BASE="${1:?staged task root required}"
SOURCE="${2:?verified remote PSID source required}"
STATA_BIN="${3:-stata-mp}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ADO_ROOT="$BASE/ado"
SMOKE_OUT="$BASE/smoke_output"
FIRST_OUT="$BASE/first_stage_output"

[[ -f "$SOURCE" && ! -e "$SMOKE_OUT" && ! -e "$FIRST_OUT" ]] || {
  echo "guarded first-stage roots must be fresh and source must exist" >&2
  exit 73
}
mkdir -p "$ADO_ROOT"

"$SCRIPT_DIR/run_psid_housing_fullsample_stata_smoke.sh" \
  "$SCRIPT_DIR/psid_housing_fullsample_driver.do" \
  "$BASE/synthetic_fixture.dta" "$SMOKE_OUT" "$ADO_ROOT"
[[ -s "$SMOKE_OUT/SMOKE_PASS.txt" ]] || { echo "synthetic smoke PASS missing" >&2; exit 75; }

STATA_PLUS="$ADO_ROOT" "$SCRIPT_DIR/launch_psid_housing_fullsample_torch.sh" \
  --source "$SOURCE" --out "$FIRST_OUT" --stage first --stata "$STATA_BIN"

for required in \
  "$FIRST_OUT/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv" \
  "$FIRST_OUT/sa_rooms_first_birth_household_aligned_v1/event_study_covariance.csv" \
  "$FIRST_OUT/sa_rooms_first_birth_household_aligned_v1/STATA_COMPLETE" \
  "$FIRST_OUT/sa_replication/own_f_c_y_all_repl_estimates.dta" \
  "$FIRST_OUT/sa_replication/own_f_c_y_all_repl_covariance.csv" \
  "$FIRST_OUT/sa_replication/own_f_c_y_all_repl_contrast.csv" \
  "$FIRST_OUT/sa_replication/own_f_c_y_all_repl.png" \
  "$FIRST_OUT/sa_replication/STATA_COMPLETE" \
  "$FIRST_OUT/first_birth_aligned_ownership/contrast.csv" \
  "$FIRST_OUT/first_birth_aligned_ownership/STATA_COMPLETE"; do
  [[ -s "$required" ]] || { echo "required first-stage artifact missing or empty: $required" >&2; exit 76; }
done

printf '%s\n' "PASS: synthetic smoke and first three PSID arms completed with covariance, contrast, graph, and Stata markers" \
  > "$BASE/FIRST_STAGE_PASS.txt"
tail -n 1 "$BASE/FIRST_STAGE_PASS.txt"
