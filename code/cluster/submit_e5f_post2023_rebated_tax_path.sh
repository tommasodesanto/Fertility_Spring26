#!/bin/bash
# Short-run rebated property-tax paths under the maintained national supply elasticity.
# Smoke: 2023--2027. Production: 2023--2035, gated on the exact smoke receipt.
#SBATCH --job-name=e5ftaxpath
#SBATCH --output=logs/slurm_e5ftaxpath_%j.out
#SBATCH --error=logs/slurm_e5ftaxpath_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
MODEL_DIR="$PROJECT_ROOT/code/model"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

MODE="${E5F_TAX_PATH_MODE:?E5F_TAX_PATH_MODE must be smoke or production}"
RUN_TAG="${E5F_TAX_PATH_RUN_TAG:?E5F_TAX_PATH_RUN_TAG is required}"
DRIVER="$MODEL_DIR/tools/run_e5f_post2023_coven_property_tax_path.py"
EXPECTED_DRIVER_SHA256="${E5F_TAX_PATH_EXPECTED_DRIVER_SHA256:?required}"
EXPECTED_POLICY_BUNDLE_SHA256="${E5F_TAX_PATH_EXPECTED_POLICY_BUNDLE_SHA256:?required}"
SELECTED_REPORT="${E5F_TAX_PATH_SELECTED_REPORT:?required}"
SELECTED_TRANSITION="${E5F_TAX_PATH_SELECTED_TRANSITION:?required}"
SOURCE="${E5F_TAX_PATH_SOURCE:?required}"
ELASTICITY="1.75"

case "$MODE" in
    smoke) POST_PERIODS=1; SUFFIX=smoke ;;
    production) POST_PERIODS=3; SUFFIX=production ;;
    *) echo "E5F_TAX_PATH_MODE must be smoke or production" >&2; exit 2 ;;
esac

ACTUAL_DRIVER_SHA256="$(sha256sum "$DRIVER" | awk '{print $1}')"
if [ "$ACTUAL_DRIVER_SHA256" != "$EXPECTED_DRIVER_SHA256" ]; then
    echo "Driver hash mismatch: actual=$ACTUAL_DRIVER_SHA256 expected=$EXPECTED_DRIVER_SHA256" >&2
    exit 2
fi

POLICY_BUNDLE_FILES=(
    "$MODEL_DIR/tools/run_dynamic_population_transition.py"
    "$MODEL_DIR/tools/run_e5f_open_population_transition.py"
    "$MODEL_DIR/tools/run_e5f_post2023_no_policy_continuations.py"
    "$MODEL_DIR/tools/run_e5f_post2023_policy_mechanisms.py"
    "$MODEL_DIR/tools/run_e5f_post2023_rebated_property_tax_smoke.py"
    "$MODEL_DIR/tools/run_e5f_post2023_coven_property_tax_smoke.py"
    "$DRIVER"
)
ACTUAL_POLICY_BUNDLE_SHA256="$({
    for path in "${POLICY_BUNDLE_FILES[@]}"; do
        test -f "$path"
        sha256sum "$path" | awk '{print $1}'
    done
} | sha256sum | awk '{print $1}')"
if [ "$ACTUAL_POLICY_BUNDLE_SHA256" != "$EXPECTED_POLICY_BUNDLE_SHA256" ]; then
    echo "Policy bundle mismatch: actual=$ACTUAL_POLICY_BUNDLE_SHA256 expected=$EXPECTED_POLICY_BUNDLE_SHA256" >&2
    exit 2
fi

if [ "$MODE" = production ]; then
    SMOKE_RECEIPT="${E5F_TAX_PATH_SMOKE_RECEIPT:?production requires smoke receipt}"
    EXPECTED_SMOKE_RECEIPT_SHA256="${E5F_TAX_PATH_EXPECTED_SMOKE_RECEIPT_SHA256:?required}"
    ACTUAL_SMOKE_RECEIPT_SHA256="$(sha256sum "$SMOKE_RECEIPT" | awk '{print $1}')"
    if [ "$ACTUAL_SMOKE_RECEIPT_SHA256" != "$EXPECTED_SMOKE_RECEIPT_SHA256" ]; then
        echo "Smoke receipt hash mismatch" >&2
        exit 2
    fi
    "$PYTHON_BIN" - "$SMOKE_RECEIPT" "$EXPECTED_DRIVER_SHA256" "$EXPECTED_POLICY_BUNDLE_SHA256" <<'PY'
import csv, hashlib, json, math, sys
from pathlib import Path

receipt_path = Path(sys.argv[1])
d = json.loads(receipt_path.read_text())
if d.get("status") != "complete_smoke":
    raise SystemExit("Smoke receipt is incomplete")
if int(d.get("post_2023_periods", -1)) != 1:
    raise SystemExit("Smoke did not run the one-period loop")
if float(d.get("housing_supply_elasticity", -1.0)) != 1.75:
    raise SystemExit("Smoke used the wrong supply elasticity")
if d.get("driver_sha256") != sys.argv[2]:
    raise SystemExit("Smoke driver hash mismatch")
if d.get("policy_bundle_sha256") != sys.argv[3]:
    raise SystemExit("Smoke policy-bundle hash mismatch")
summary_path = receipt_path.parent / "summary.json"
if hashlib.sha256(summary_path.read_bytes()).hexdigest() != d.get("summary_sha256"):
    raise SystemExit("Smoke summary hash mismatch")
summary = json.loads(summary_path.read_text())
rows = list(csv.DictReader((receipt_path.parent / "policy_paths.csv").open()))
baseline = [
    row for row in rows
    if row["policy_case"] == "rebated-tax1-baseline"
    and int(float(row["calendar_year"])) == 2023
]
if len(baseline) != 1:
    raise SystemExit("Smoke has no unique rebated-1% 2023 row")
anchor = summary["supply_reanchor"]
if not math.isclose(float(baseline[0]["asset_price"]), float(anchor["inherited_2023_asset_price"]), rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit("Smoke baseline price does not nest the inherited 2023 price")
if not math.isclose(float(baseline[0]["equal_transfer_period_units"]), float(anchor["fixed_price_equal_transfer"]), rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit("Smoke baseline transfer does not nest the fixed-price fiscal root")
PY
fi

OUTDIR="$PROJECT_ROOT/output/model/e5f_post2023_rebated_tax_path_${RUN_TAG}_${SUFFIX}"
if [ -d "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
    echo "Refusing to overwrite nonempty output: $OUTDIR" >&2
    exit 2
fi

PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools" "$PYTHON_BIN" "$DRIVER" \
    --selected-report "$SELECTED_REPORT" \
    --selected-case-transition "$SELECTED_TRANSITION" \
    --source "$SOURCE" \
    --output-dir "$OUTDIR" \
    --expected-report-sha256 "${E5F_TAX_PATH_EXPECTED_REPORT_SHA256:?required}" \
    --expected-case-transition-sha256 "${E5F_TAX_PATH_EXPECTED_TRANSITION_SHA256:?required}" \
    --expected-source-sha256 "${E5F_TAX_PATH_EXPECTED_SOURCE_SHA256:?required}" \
    --expected-target-fingerprint "${E5F_TAX_PATH_EXPECTED_TARGET_FINGERPRINT:?required}" \
    --expected-code-bundle-sha256 "${E5F_TAX_PATH_EXPECTED_CODE_BUNDLE_SHA256:?required}" \
    --expected-renewal-contract-sha256 "${E5F_TAX_PATH_EXPECTED_RENEWAL_SHA256:?required}" \
    --expected-scientific-contract-sha256 "${E5F_TAX_PATH_EXPECTED_SCIENTIFIC_SHA256:?required}" \
    --expected-selection-sha256 "${E5F_TAX_PATH_EXPECTED_SELECTION_SHA256:?required}" \
    --housing-supply-elasticity "$ELASTICITY" \
    --post-2023-periods "$POST_PERIODS"

"$PYTHON_BIN" - "$OUTDIR/summary.json" "$OUTDIR/launch_receipt.json" "$POST_PERIODS" "$EXPECTED_DRIVER_SHA256" "$EXPECTED_POLICY_BUNDLE_SHA256" "$MODE" <<'PY'
import csv, hashlib, json, math, sys
from pathlib import Path

summary_path = Path(sys.argv[1])
receipt_path = Path(sys.argv[2])
periods = int(sys.argv[3])
d = json.loads(summary_path.read_text())
if d.get("status") != "complete_exact_rebated_property_tax_path":
    raise SystemExit("Path summary is incomplete")
if int(d.get("post_2023_periods", -1)) != periods:
    raise SystemExit("Post-period count mismatch")
if d.get("calendar_years") != [2023 + 4 * i for i in range(periods + 1)]:
    raise SystemExit("Calendar sequence mismatch")
if float(d.get("housing_supply_elasticity", -1.0)) != 1.75:
    raise SystemExit("Supply elasticity mismatch")
if len(d.get("effects_by_date", [])) != periods + 1:
    raise SystemExit("Effect-row count mismatch")
if abs(float(d.get("maximum_pre_2043_population_percent_gap", math.inf))) > 1e-10:
    raise SystemExit("Twenty-year population lag failed")
for rows in d.get("case_gates", {}).values():
    if len(rows) != periods + 1:
        raise SystemExit("Case gate-row count mismatch")
    for row in rows:
        if float(row["market_residual"]) > 2e-4:
            raise SystemExit("Market gate failed")
        if float(row["fiscal_residual"]) > 2.5e-5:
            raise SystemExit("Fiscal gate failed")
        if float(row["mass_residual"]) > 2e-8:
            raise SystemExit("Mass gate failed")
        if int(row["nonfinite_distribution_count"]) != 0:
            raise SystemExit("Nonfinite distribution")
for name in ("policy_paths.csv", "effects_by_date.csv", "rebated_property_tax_path.png", "rebated_property_tax_path.pdf"):
    artifact = summary_path.parent / name
    if not artifact.is_file() or artifact.stat().st_size == 0:
        raise SystemExit(f"Missing artifact: {artifact}")
rows = list(csv.DictReader((summary_path.parent / "policy_paths.csv").open()))
baseline = [
    row for row in rows
    if row["policy_case"] == "rebated-tax1-baseline"
    and int(float(row["calendar_year"])) == 2023
]
if len(baseline) != 1:
    raise SystemExit("No unique rebated-1% 2023 row")
anchor = d["supply_reanchor"]
if not math.isclose(float(baseline[0]["asset_price"]), float(anchor["inherited_2023_asset_price"]), rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit("Baseline price does not nest the inherited 2023 price")
if not math.isclose(float(baseline[0]["equal_transfer_period_units"]), float(anchor["fixed_price_equal_transfer"]), rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit("Baseline transfer does not nest the fixed-price fiscal root")
receipt = {
    "status": "complete_smoke" if sys.argv[6] == "smoke" else "complete_production",
    "post_2023_periods": periods,
    "housing_supply_elasticity": 1.75,
    "driver_sha256": sys.argv[4],
    "policy_bundle_sha256": sys.argv[5],
    "summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
}
receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print(f"REBATED_TAX_PATH_LAUNCHER_PASS mode={sys.argv[6]} summary={summary_path}")
PY
