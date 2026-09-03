#!/usr/bin/env bash
# Exact current-calibration 2023 rebated-tax impact and Shapley decomposition.
# The one-date tax solve must pass every gate before the eight-cell decomposition.
# All scientific inputs and both Python drivers are hash-pinned; outputs are immutable.
#SBATCH --job-name=e5fshap
#SBATCH --output=logs/slurm_e5fshap_%j.out
#SBATCH --error=logs/slurm_e5fshap_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_rebated_tax_shapley_current.sh"
LOCAL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

required=(
    E5F_SHAPLEY_RUN_TAG
    E5F_SHAPLEY_REPORT
    E5F_SHAPLEY_TASK_SUMMARY
    E5F_SHAPLEY_CASE_DIR
    E5F_SHAPLEY_TRANSITION
    E5F_SHAPLEY_SOURCE
    E5F_SHAPLEY_EXPECTED_REPORT_SHA256
    E5F_SHAPLEY_EXPECTED_TASK_SUMMARY_SHA256
    E5F_SHAPLEY_EXPECTED_TRANSITION_SHA256
    E5F_SHAPLEY_EXPECTED_SOURCE_SHA256
    E5F_SHAPLEY_EXPECTED_TARGET_FINGERPRINT
    E5F_SHAPLEY_EXPECTED_CODE_BUNDLE_SHA256
    E5F_SHAPLEY_EXPECTED_RENEWAL_CONTRACT_SHA256
    E5F_SHAPLEY_EXPECTED_IMPACT_DRIVER_SHA256
    E5F_SHAPLEY_EXPECTED_SHAPLEY_DRIVER_SHA256
    E5F_SHAPLEY_EXPECTED_LAUNCHER_SHA256
)

require_environment() {
    local name
    for name in "${required[@]}"; do
        if [ -z "${!name:-}" ]; then
            echo "$name is required" >&2
            exit 2
        fi
    done
    if [[ ! "$E5F_SHAPLEY_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
        echo "invalid run tag: $E5F_SHAPLEY_RUN_TAG" >&2
        exit 2
    fi
}

if [ "${1:-}" = "--submit" ]; then
    shift
    require_environment
    cd "$LOCAL_SCRIPT_DIR"
    exec sbatch --export=ALL "$LOCAL_SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Run with --submit or sbatch}"
require_environment

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$LOCAL_SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
MODEL_DIR="$PROJECT_ROOT/code/model"
IMPACT_DRIVER="$MODEL_DIR/tools/run_e5f_post2023_rebated_property_tax_smoke.py"
SHAPLEY_DRIVER="$MODEL_DIR/tools/run_e5f_rebated_tax_shapley_diagnosis.py"
SCRIPT_PATH="$CLUSTER_DIR/$SCRIPT_NAME"
mkdir -p "$CLUSTER_DIR/logs"

resolve_path() {
    if [[ "$1" = /* ]]; then printf '%s\n' "$1"; else printf '%s\n' "$PROJECT_ROOT/$1"; fi
}
hash_file() {
    sha256sum "$1" | awk '{print $1}'
}
check_hash() {
    if [ "$2" != "$3" ]; then
        echo "$1 hash mismatch: actual=$2 expected=$3" >&2
        exit 2
    fi
}

REPORT="$(resolve_path "$E5F_SHAPLEY_REPORT")"
TASK_SUMMARY="$(resolve_path "$E5F_SHAPLEY_TASK_SUMMARY")"
CASE_DIR="$(resolve_path "$E5F_SHAPLEY_CASE_DIR")"
TRANSITION="$(resolve_path "$E5F_SHAPLEY_TRANSITION")"
SOURCE="$(resolve_path "$E5F_SHAPLEY_SOURCE")"
for path in "$REPORT" "$TASK_SUMMARY" "$TRANSITION" "$SOURCE" "$IMPACT_DRIVER" "$SHAPLEY_DRIVER" "$SCRIPT_PATH"; do
    [ -s "$path" ] || { echo "missing input: $path" >&2; exit 2; }
done
[ -d "$CASE_DIR" ] || { echo "missing case directory: $CASE_DIR" >&2; exit 2; }
check_hash report "$(hash_file "$REPORT")" "$E5F_SHAPLEY_EXPECTED_REPORT_SHA256"
check_hash task-summary "$(hash_file "$TASK_SUMMARY")" "$E5F_SHAPLEY_EXPECTED_TASK_SUMMARY_SHA256"
check_hash transition "$(hash_file "$TRANSITION")" "$E5F_SHAPLEY_EXPECTED_TRANSITION_SHA256"
check_hash source "$(hash_file "$SOURCE")" "$E5F_SHAPLEY_EXPECTED_SOURCE_SHA256"
check_hash impact-driver "$(hash_file "$IMPACT_DRIVER")" "$E5F_SHAPLEY_EXPECTED_IMPACT_DRIVER_SHA256"
check_hash shapley-driver "$(hash_file "$SHAPLEY_DRIVER")" "$E5F_SHAPLEY_EXPECTED_SHAPLEY_DRIVER_SHA256"
check_hash launcher "$(hash_file "$SCRIPT_PATH")" "$E5F_SHAPLEY_EXPECTED_LAUNCHER_SHA256"

IMPACT_OUT="$PROJECT_ROOT/output/model/e5f_post2023_rebated_tax_impact_${E5F_SHAPLEY_RUN_TAG}"
SHAPLEY_OUT="$PROJECT_ROOT/output/model/e5f_rebated_tax_shapley_diagnosis_${E5F_SHAPLEY_RUN_TAG}"
for path in "$IMPACT_OUT" "$SHAPLEY_OUT"; do
    case "$path" in "$PROJECT_ROOT"/output/model/*) ;; *) echo "unsafe output: $path" >&2; exit 2 ;; esac
    if [ -e "$path" ]; then echo "refusing existing output: $path" >&2; exit 2; fi
done

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

common_args=(
    --selected-report "$REPORT"
    --selected-task-summary "$TASK_SUMMARY"
    --selected-case-dir "$CASE_DIR"
    --selected-case-transition "$TRANSITION"
    --source "$SOURCE"
    --expected-report-sha256 "$E5F_SHAPLEY_EXPECTED_REPORT_SHA256"
    --expected-task-summary-sha256 "$E5F_SHAPLEY_EXPECTED_TASK_SUMMARY_SHA256"
    --expected-case-transition-sha256 "$E5F_SHAPLEY_EXPECTED_TRANSITION_SHA256"
    --expected-source-sha256 "$E5F_SHAPLEY_EXPECTED_SOURCE_SHA256"
    --expected-target-fingerprint "$E5F_SHAPLEY_EXPECTED_TARGET_FINGERPRINT"
    --expected-code-bundle-sha256 "$E5F_SHAPLEY_EXPECTED_CODE_BUNDLE_SHA256"
    --expected-renewal-contract-sha256 "$E5F_SHAPLEY_EXPECTED_RENEWAL_CONTRACT_SHA256"
    --expected-scientific-contract-sha256 collector-style-not-applicable
    --expected-selection-sha256 collector-style-not-applicable
)

"$PYTHON_BIN" "$IMPACT_DRIVER" "${common_args[@]}" --output-dir "$IMPACT_OUT"
"$PYTHON_BIN" - "$IMPACT_OUT/summary.json" <<'PY'
import json, math, sys
p=json.load(open(sys.argv[1], encoding="utf-8"))
if p.get("status") != "complete_exact_rebated_property_tax_smoke":
    raise SystemExit("rebated-tax impact is incomplete")
if list(p.get("cases", [])) != ["status-quo-tax1-unrebated", "tax1-equal-rebate", "tax2-equal-rebate"]:
    raise SystemExit("rebated-tax impact case menu mismatch")
for name, gates in (p.get("case_gates") or {}).items():
    if float(gates.get("market_residual", math.inf)) > 2e-4:
        raise SystemExit(f"market gate failed for {name}")
    if float(gates.get("mass_residual", math.inf)) > 2e-8:
        raise SystemExit(f"mass gate failed for {name}")
    fiscal=gates.get("fiscal_residual")
    if fiscal is not None and float(fiscal) > 2.5e-5:
        raise SystemExit(f"fiscal gate failed for {name}")
PY

"$PYTHON_BIN" "$SHAPLEY_DRIVER" "${common_args[@]}" \
    --rebated-path-csv "$IMPACT_OUT/impact_results.csv" \
    --output-dir "$SHAPLEY_OUT"

"$PYTHON_BIN" - "$IMPACT_OUT" "$SHAPLEY_OUT" "$IMPACT_DRIVER" "$SHAPLEY_DRIVER" "$SCRIPT_PATH" <<'PY'
import hashlib, json, sys
from pathlib import Path
impact, shapley, impact_driver, shapley_driver, launcher = map(Path, sys.argv[1:])
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
s=json.loads((shapley / "summary.json").read_text(encoding="utf-8"))
if s.get("status") != "complete" or int(s.get("model_solve_count", -1)) < 8:
    raise SystemExit("Shapley decomposition is incomplete")
required=(
    impact / "summary.json",
    impact / "impact_results.csv",
    shapley / "summary.json",
    shapley / "component_cells.csv",
    shapley / "shapley_decomposition.csv",
    shapley / "rebated_tax_shapley.pdf",
    shapley / "rebated_tax_shapley.png",
)
for path in required:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"missing output: {path}")
status={
    "status": "complete",
    "impact_summary_sha256": sha(impact / "summary.json"),
    "impact_results_sha256": sha(impact / "impact_results.csv"),
    "shapley_summary_sha256": sha(shapley / "summary.json"),
    "shapley_csv_sha256": sha(shapley / "shapley_decomposition.csv"),
    "impact_driver_sha256": sha(impact_driver),
    "shapley_driver_sha256": sha(shapley_driver),
    "launcher_sha256": sha(launcher),
}
(shapley / "launcher_status.json").write_text(json.dumps(status, indent=2, sort_keys=True)+"\n", encoding="utf-8")
print("CURRENT_REBATED_TAX_SHAPLEY_COMPLETE")
PY
