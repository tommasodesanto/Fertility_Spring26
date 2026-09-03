#!/usr/bin/env bash
# Five-case post-2023 housing-policy mechanism array.  Smoke uses one later
# date; production uses ten later dates (through 2063) and requires the exact
# completed four-case smoke comparison receipt.
#SBATCH --job-name=e5fpol
#SBATCH --output=logs/slurm_e5fpol_%A_%a.out
#SBATCH --error=logs/slurm_e5fpol_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:40:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-5%5

set -euo pipefail

SCRIPT_NAME="submit_e5f_post2023_policy_mechanisms.sh"
LOCAL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REQUIRED=(
    E5F_POLICY_MODE
    E5F_POLICY_RUN_TAG
    E5F_POLICY_REPORT
    E5F_POLICY_EXPECTED_REPORT_SHA256
    E5F_POLICY_SELECTED_TRANSITION
    E5F_POLICY_EXPECTED_SELECTED_TRANSITION_SHA256
    E5F_POLICY_SOURCE
    E5F_POLICY_EXPECTED_SOURCE_SHA256
    E5F_POLICY_EXPECTED_TARGET_FINGERPRINT
    E5F_POLICY_EXPECTED_CODE_BUNDLE_SHA256
    E5F_POLICY_EXPECTED_RENEWAL_CONTRACT_SHA256
    E5F_POLICY_EXPECTED_SCIENTIFIC_CONTRACT_SHA256
    E5F_POLICY_EXPECTED_SELECTION_SHA256
    E5F_POLICY_EXPECTED_DRIVER_SHA256
    E5F_POLICY_EXPECTED_LAUNCHER_SHA256
)

require_environment() {
    local name
    for name in "${REQUIRED[@]}"; do
        if [ -z "${!name:-}" ]; then
            echo "$name is required" >&2
            exit 2
        fi
    done
    case "$E5F_POLICY_MODE" in
        smoke|production) ;;
        *) echo "E5F_POLICY_MODE must be smoke or production" >&2; exit 2 ;;
    esac
    if [ -n "${E5F_POLICY_SELECTED_TASK_SUMMARY:-}${E5F_POLICY_SELECTED_CASE_DIR:-}${E5F_POLICY_EXPECTED_TASK_SUMMARY_SHA256:-}" ]; then
        : "${E5F_POLICY_SELECTED_TASK_SUMMARY:?selected task summary is required for a collector-style report}"
        : "${E5F_POLICY_SELECTED_CASE_DIR:?selected case directory is required for a collector-style report}"
        : "${E5F_POLICY_EXPECTED_TASK_SUMMARY_SHA256:?selected task-summary hash is required for a collector-style report}"
    fi
    if [ "$E5F_POLICY_MODE" = "production" ]; then
        : "${E5F_POLICY_SMOKE_ROOT:?E5F_POLICY_SMOKE_ROOT is required in production}"
        : "${E5F_POLICY_SMOKE_COMPARISON_SUMMARY:?smoke comparison summary is required}"
        : "${E5F_POLICY_EXPECTED_SMOKE_COMPARISON_SHA256:?smoke comparison hash is required}"
        : "${E5F_POLICY_SMOKE_COMPARISON_MANIFEST:?smoke comparison manifest is required}"
        : "${E5F_POLICY_EXPECTED_SMOKE_MANIFEST_SHA256:?smoke manifest hash is required}"
    fi
}

if [ "${1:-}" = "--submit" ]; then
    shift
    require_environment
    if [ "$E5F_POLICY_MODE" = "smoke" ]; then
        TIME_LIMIT="${E5F_POLICY_TIME_LIMIT:-00:40:00}"
    else
        TIME_LIMIT="${E5F_POLICY_TIME_LIMIT:-00:55:00}"
    fi
    cd "$LOCAL_SCRIPT_DIR"
    exec sbatch --time="$TIME_LIMIT" --mem="${E5F_POLICY_MEMORY:-8G}" \
        --array=1-5%5 --export=ALL "$LOCAL_SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Run with --submit or sbatch}"
: "${SLURM_ARRAY_TASK_ID:?Five-case array task id is required}"
require_environment

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$LOCAL_SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
MODEL_DIR="$PROJECT_ROOT/code/model"
DRIVER="$MODEL_DIR/tools/run_e5f_post2023_policy_mechanisms.py"
SCRIPT_PATH="$CLUSTER_DIR/$SCRIPT_NAME"
mkdir -p "$CLUSTER_DIR/logs"

case "$SLURM_ARRAY_TASK_ID" in
    1) POLICY_CASE="baseline" ;;
    2) POLICY_CASE="supply-plus-20" ;;
    3) POLICY_CASE="dependent-child-ltv95" ;;
    4) POLICY_CASE="combined" ;;
    5) POLICY_CASE="property-tax-2pct-no-rebate" ;;
    *) echo "Array task must lie in 1..5" >&2; exit 2 ;;
esac

resolve_path() {
    if [[ "$1" = /* ]]; then printf '%s\n' "$1"; else printf '%s\n' "$PROJECT_ROOT/$1"; fi
}
hash_file() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi
}
check_hash() {
    if [ "$2" != "$3" ]; then echo "$1 hash mismatch: actual=$2 expected=$3" >&2; exit 2; fi
}

REPORT="$(resolve_path "$E5F_POLICY_REPORT")"
TRANSITION="$(resolve_path "$E5F_POLICY_SELECTED_TRANSITION")"
SOURCE="$(resolve_path "$E5F_POLICY_SOURCE")"
CALIBRATION_ARGS=()
if [ -n "${E5F_POLICY_SELECTED_TASK_SUMMARY:-}" ]; then
    TASK_SUMMARY="$(resolve_path "$E5F_POLICY_SELECTED_TASK_SUMMARY")"
    CASE_DIR="$(resolve_path "$E5F_POLICY_SELECTED_CASE_DIR")"
    [ -s "$TASK_SUMMARY" ] || { echo "Missing selected task summary: $TASK_SUMMARY" >&2; exit 2; }
    [ -d "$CASE_DIR" ] || { echo "Missing selected case directory: $CASE_DIR" >&2; exit 2; }
    check_hash selected-task-summary "$(hash_file "$TASK_SUMMARY")" "$E5F_POLICY_EXPECTED_TASK_SUMMARY_SHA256"
    CALIBRATION_ARGS+=(
        --selected-task-summary "$TASK_SUMMARY"
        --selected-case-dir "$CASE_DIR"
        --expected-task-summary-sha256 "$E5F_POLICY_EXPECTED_TASK_SUMMARY_SHA256"
    )
fi
RUN_ROOT="${E5F_POLICY_OUT_ROOT:-$PROJECT_ROOT/output/model/e5f_post2023_policy_mechanisms_${E5F_POLICY_RUN_TAG}}"
RUN_ROOT="$(resolve_path "$RUN_ROOT")"
OUTDIR="$RUN_ROOT/$POLICY_CASE"
case "$OUTDIR" in "$PROJECT_ROOT"/output/model/*) ;; *) echo "Unsafe output path: $OUTDIR" >&2; exit 2 ;; esac
if [ -e "$OUTDIR" ]; then echo "Refusing existing output: $OUTDIR" >&2; exit 2; fi
for path in "$REPORT" "$TRANSITION" "$SOURCE" "$DRIVER" "$SCRIPT_PATH"; do
    [ -s "$path" ] || { echo "Missing input: $path" >&2; exit 2; }
done

check_hash report "$(hash_file "$REPORT")" "$E5F_POLICY_EXPECTED_REPORT_SHA256"
check_hash selected-transition "$(hash_file "$TRANSITION")" "$E5F_POLICY_EXPECTED_SELECTED_TRANSITION_SHA256"
check_hash source "$(hash_file "$SOURCE")" "$E5F_POLICY_EXPECTED_SOURCE_SHA256"
check_hash driver "$(hash_file "$DRIVER")" "$E5F_POLICY_EXPECTED_DRIVER_SHA256"
check_hash launcher "$(hash_file "$SCRIPT_PATH")" "$E5F_POLICY_EXPECTED_LAUNCHER_SHA256"

if [ "$E5F_POLICY_MODE" = "production" ]; then
    SMOKE_ROOT="$(resolve_path "$E5F_POLICY_SMOKE_ROOT")"
    SMOKE_SUMMARY="$(resolve_path "$E5F_POLICY_SMOKE_COMPARISON_SUMMARY")"
    SMOKE_MANIFEST="$(resolve_path "$E5F_POLICY_SMOKE_COMPARISON_MANIFEST")"
    check_hash smoke-comparison "$(hash_file "$SMOKE_SUMMARY")" "$E5F_POLICY_EXPECTED_SMOKE_COMPARISON_SHA256"
    check_hash smoke-manifest "$(hash_file "$SMOKE_MANIFEST")" "$E5F_POLICY_EXPECTED_SMOKE_MANIFEST_SHA256"
    python3 - \
        "$SMOKE_ROOT" "$SMOKE_SUMMARY" "$E5F_POLICY_EXPECTED_DRIVER_SHA256" \
        "$E5F_POLICY_EXPECTED_REPORT_SHA256" \
        "$E5F_POLICY_EXPECTED_SELECTED_TRANSITION_SHA256" \
        "$E5F_POLICY_EXPECTED_SOURCE_SHA256" <<'PY'
import json, sys
from pathlib import Path
root = Path(sys.argv[1])
comparison = json.loads(Path(sys.argv[2]).read_text())
driver_sha, report_sha, transition_sha, source_sha = sys.argv[3:]
cases = ("baseline", "supply-plus-20", "dependent-child-ltv95", "combined", "property-tax-2pct-no-rebate")
if comparison.get("status") != "complete_post2023_policy_mechanism_comparison":
    raise SystemExit("smoke comparison is incomplete")
if comparison.get("cases") != list(cases) or int(comparison.get("post_2023_periods", -1)) != 1:
    raise SystemExit("smoke comparison has the wrong case set or horizon")
hashes = comparison.get("input_hashes") or {}
for field, expected in (
    ("selected_report_sha256", report_sha),
    ("selected_case_transition_sha256", transition_sha),
    ("source_sha256", source_sha),
):
    if str(hashes.get(field)) != expected:
        raise SystemExit(f"smoke input mismatch: {field}")
for case in cases:
    status = json.loads((root / case / "launcher_status.json").read_text())
    if status.get("status") != "complete" or status.get("mode") != "smoke":
        raise SystemExit(f"incomplete smoke launcher status: {case}")
    if status.get("policy_case") != case or int(status.get("post_2023_periods", -1)) != 1:
        raise SystemExit(f"wrong smoke case contract: {case}")
    if status.get("driver_sha256") != driver_sha:
        raise SystemExit(f"smoke driver differs for {case}")
PY
fi

if [ "$E5F_POLICY_MODE" = "smoke" ]; then POST_PERIODS=1; else POST_PERIODS=10; fi
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"

"$PYTHON_BIN" "$DRIVER" \
    --selected-report "$REPORT" \
    "${CALIBRATION_ARGS[@]}" \
    --selected-case-transition "$TRANSITION" \
    --source "$SOURCE" \
    --output-dir "$OUTDIR" \
    --policy-case "$POLICY_CASE" \
    --expected-report-sha256 "$E5F_POLICY_EXPECTED_REPORT_SHA256" \
    --expected-case-transition-sha256 "$E5F_POLICY_EXPECTED_SELECTED_TRANSITION_SHA256" \
    --expected-source-sha256 "$E5F_POLICY_EXPECTED_SOURCE_SHA256" \
    --expected-target-fingerprint "$E5F_POLICY_EXPECTED_TARGET_FINGERPRINT" \
    --expected-code-bundle-sha256 "$E5F_POLICY_EXPECTED_CODE_BUNDLE_SHA256" \
    --expected-renewal-contract-sha256 "$E5F_POLICY_EXPECTED_RENEWAL_CONTRACT_SHA256" \
    --expected-scientific-contract-sha256 "$E5F_POLICY_EXPECTED_SCIENTIFIC_CONTRACT_SHA256" \
    --expected-selection-sha256 "$E5F_POLICY_EXPECTED_SELECTION_SHA256" \
    --post-2023-periods "$POST_PERIODS"

"$PYTHON_BIN" - "$OUTDIR" "$POLICY_CASE" "$E5F_POLICY_MODE" "$POST_PERIODS" "$DRIVER" "$SCRIPT_PATH" <<'PY'
import hashlib, json, math, sys
from pathlib import Path
outdir, case, mode, periods, driver, launcher = sys.argv[1:]
outdir = Path(outdir)
def read(path): return json.loads(Path(path).read_text())
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
summary = read(outdir / "summary.json")
manifest = read(outdir / "manifest.json")
if summary.get("status") != "complete_post2023_policy_mechanism_case": raise SystemExit("incomplete summary")
if (summary.get("policy") or {}).get("case") != case: raise SystemExit("wrong policy case")
if int(summary.get("post_2023_periods", -1)) != int(periods): raise SystemExit("wrong horizon")
for name, expected in (manifest.get("artifacts") or {}).items():
    if sha(outdir / name) != expected: raise SystemExit(f"artifact hash mismatch: {name}")
g = summary.get("numerical_gates") or {}
if float(g.get("maximum_market_residual", math.inf)) > 2e-4: raise SystemExit("market gate failed")
if abs(float(g.get("maximum_mass_residual", math.inf))) > 2e-8: raise SystemExit("mass gate failed")
status = {
    "status": "complete",
    "mode": mode,
    "policy_case": case,
    "post_2023_periods": int(periods),
    "summary_sha256": sha(outdir / "summary.json"),
    "manifest_sha256": sha(outdir / "manifest.json"),
    "driver_sha256": sha(driver),
    "launcher_sha256": sha(launcher),
}
(outdir / "launcher_status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
print(f"POLICY_LAUNCHER_COMPLETE case={case} mode={mode}")
PY
