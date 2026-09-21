#!/usr/bin/env bash
set -euo pipefail

DRIVER="${1:?driver path required}"
FIXTURE="${2:?fixture path required}"
OUTROOT="${3:?fresh output root required}"
ADO_ROOT="${4:-$(dirname "$DRIVER")/ado}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATA_BIN="${STATA_BIN:-stata-mp}"

[[ -f "$DRIVER" && ! -e "$OUTROOT" ]] || {
  echo "driver must exist and output root must be fresh" >&2
  exit 73
}
mkdir -p "$OUTROOT"
mkdir -p "$ADO_ROOT"
cd "$OUTROOT"
required_ado=(
  require.ado
  eventstudyinteract.ado
  svmat2.ado
  a/avar.ado
  f/ftools.ado
  f/ftools.mata
  i/ivreg2.ado
  r/reghdfe.ado
  l/livreg2.mlib
  l/lmoremata.mlib
)
for rel in "${required_ado[@]}"; do
  [[ -s "$ADO_ROOT/$rel" ]] || { echo "required dependency missing or empty: $ADO_ROOT/$rel" >&2; exit 78; }
done
{
  printf '%s\n' 'PSID synthetic smoke dependency inventory'
  for rel in "${required_ado[@]}"; do
    version_line="$(grep -m1 '^*!' "$ADO_ROOT/$rel" || true)"
    printf '%s\t%s\n' "$rel" "$version_line"
    shasum -a 256 "$ADO_ROOT/$rel"
  done
} > "$OUTROOT/DEPENDENCY_INVENTORY.txt"
set +e
"$STATA_BIN" -b do "$SCRIPT_DIR/test_psid_housing_fullsample_stata_smoke.do" \
  "$DRIVER" "$FIXTURE" "$OUTROOT" "$ADO_ROOT" > "$OUTROOT/stata_stdout.log" 2>&1
stata_rc=$?
set -e
if [[ "$stata_rc" -ne 0 ]]; then
  cat "$OUTROOT/stata_stdout.log" >&2
  exit "$stata_rc"
fi
[[ -s "$ADO_ROOT/l/lftools.mlib" ]] || {
  echo "ftools compile did not produce lftools.mlib" >&2
  exit 79
}
printf '%s\n' 'l/lftools.mlib compiled during smoke' >> "$OUTROOT/DEPENDENCY_INVENTORY.txt"
shasum -a 256 "$ADO_ROOT/l/lftools.mlib" >> "$OUTROOT/DEPENDENCY_INVENTORY.txt"

[[ -f "$OUTROOT/first_birth_aligned_ownership/contrast.csv" ]] || {
  echo "contrast export missing" >&2
  exit 75
}
[[ -f "$OUTROOT/first_birth_aligned_ownership/event_study_covariance.csv" ]] || {
  echo "full covariance export missing" >&2
  exit 76
}
[[ -s "$OUTROOT/first_birth_aligned_ownership/STATA_COMPLETE" ]] || {
  echo "Stata completion marker missing" >&2
  exit 76
}
[[ -s "$OUTROOT/first_birth_aligned_ownership/aligned_first_ownership.log" ]] || {
  echo "Stata log missing" >&2
  exit 77
}
python3 - "$OUTROOT" <<'PY'
import csv
import math
import sys
from pathlib import Path

root = Path(sys.argv[1]) / "first_birth_aligned_ownership"
with (root / "contrast.csv").open(newline="") as fh:
    row = next(csv.DictReader(fh))
for key in ("contrast_l3_minus_f1", "contrast_se", "contrast_ci_lo", "contrast_ci_hi"):
    value = float(row[key])
    assert math.isfinite(value), (key, value)
assert int(row["estimation_observations"]) > 0
with (root / "event_study_covariance.csv").open(newline="") as fh:
    rows = list(csv.DictReader(fh))
assert rows
numeric = []
for row in rows:
    for key, value in row.items():
        if key == "row_name":
            continue
        assert value not in ("", "."), (key, value)
        numeric.append(float(value))
assert numeric and all(math.isfinite(value) for value in numeric)
PY
printf '%s\n' "PASS: synthetic dependency, cohort regression, full covariance, and export smoke" \
  > "$OUTROOT/SMOKE_PASS.txt"
tail -n 1 "$OUTROOT/SMOKE_PASS.txt"
