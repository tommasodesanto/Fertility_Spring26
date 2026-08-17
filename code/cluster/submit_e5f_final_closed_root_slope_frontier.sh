#!/bin/bash
# Fixed-parameter fertility-slope frontier around the final certified E5F ridge
# estimate.  This launcher changes no scientific code and performs no
# optimization.  Array tasks 1--5 map exactly to multipliers
# 1, 0.5, 0.25, 0.1, and 0.05 on both fertility-logit scales.
#
# E5F_FINAL_SLOPE_MODE=factor (default): submit --array=1-5%5.
# E5F_FINAL_SLOPE_MODE=smoke: submit --array=2-2%1 for the exact factor-0.5 loop.
# E5F_FINAL_SLOPE_MODE=collect: submit --array=1-1%1 after all five factors finish.
#
# Every mode requires the explicit paths and hashes listed below.  No value is
# inferred from an old status note or diagnostic packet.
#SBATCH --job-name=e5frootslope
#SBATCH --output=logs/slurm_e5frootslope_%A_%a.out
#SBATCH --error=logs/slurm_e5frootslope_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-5%5

set -euo pipefail

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

MODE="${E5F_FINAL_SLOPE_MODE:-factor}"
RUN_TAG="${E5F_FINAL_SLOPE_RUN_TAG:?E5F_FINAL_SLOPE_RUN_TAG is required}"
RIDGE_REPORT="${E5F_FINAL_SLOPE_RIDGE_REPORT:?E5F_FINAL_SLOPE_RIDGE_REPORT is required}"
TASK_SUMMARY="${E5F_FINAL_SLOPE_TASK_SUMMARY:?E5F_FINAL_SLOPE_TASK_SUMMARY is required}"
SOURCE="${E5F_FINAL_SLOPE_SOURCE:?E5F_FINAL_SLOPE_SOURCE is required}"
EXPECTED_REPORT_SHA="${E5F_FINAL_SLOPE_EXPECTED_REPORT_SHA256:?expected report hash is required}"
EXPECTED_TASK_SHA="${E5F_FINAL_SLOPE_EXPECTED_TASK_SHA256:?expected task hash is required}"
EXPECTED_SOURCE_SHA="${E5F_FINAL_SLOPE_EXPECTED_SOURCE_SHA256:?expected source hash is required}"
EXPECTED_TARGET_SET="${E5F_FINAL_SLOPE_EXPECTED_TARGET_SET:?expected target set is required}"
EXPECTED_TARGET_FP="${E5F_FINAL_SLOPE_EXPECTED_TARGET_FINGERPRINT:?expected target fingerprint is required}"
EXPECTED_CODE_SHA="${E5F_FINAL_SLOPE_EXPECTED_CODE_BUNDLE_SHA256:?expected code bundle is required}"
EXPECTED_SCIENTIFIC_SHA="${E5F_FINAL_SLOPE_EXPECTED_SCIENTIFIC_CONTRACT_SHA256:?expected scientific-contract hash is required}"
EXPECTED_SELECTION_SHA="${E5F_FINAL_SLOPE_EXPECTED_SELECTION_SHA256:?expected selection hash is required}"
EXPECTED_CANDIDATE_SHA="${E5F_FINAL_SLOPE_EXPECTED_SELECTED_CANDIDATE_SHA256:?expected candidate hash is required}"
EXPECTED_CANDIDATE_ID="${E5F_FINAL_SLOPE_EXPECTED_SELECTED_CANDIDATE_ID:?expected candidate id is required}"
EXPECTED_RENEWAL_SHA="${E5F_FINAL_SLOPE_EXPECTED_RENEWAL_CONTRACT_SHA256:?expected renewal-contract hash is required}"
EXPECTED_DATED_SHA="${E5F_FINAL_SLOPE_EXPECTED_DATED_CONTRACT_SHA256:?expected dated-contract hash is required}"
EXPECTED_FRONTIER_DRIVER_SHA="${E5F_FINAL_SLOPE_EXPECTED_FRONTIER_DRIVER_SHA256:?expected frontier-driver hash is required}"
EXPECTED_CONTINUATION_DRIVER_SHA="${E5F_FINAL_SLOPE_EXPECTED_CONTINUATION_DRIVER_SHA256:?expected continuation-driver hash is required}"
EXPECTED_LAUNCHER_SHA="${E5F_FINAL_SLOPE_EXPECTED_LAUNCHER_SHA256:?expected launcher hash is required}"
EXPECTED_DIAGNOSTIC_BUNDLE_SHA="${E5F_FINAL_SLOPE_EXPECTED_DIAGNOSTIC_BUNDLE_SHA256:?expected diagnostic-bundle hash is required}"
SMOKE_RECEIPT="${E5F_FINAL_SLOPE_SMOKE_RECEIPT:-}"
EXPECTED_SMOKE_RECEIPT_SHA="${E5F_FINAL_SLOPE_EXPECTED_SMOKE_RECEIPT_SHA256:-}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if [ "$MODE" != "factor" ] && [ "$MODE" != "smoke" ] && [ "$MODE" != "collect" ]; then
    echo "E5F_FINAL_SLOPE_MODE must be factor, smoke, or collect" >&2
    exit 2
fi
if [ "$MODE" = "smoke" ] && [ "$TASK_ID" != "2" ]; then
    echo "exact smoke must be submitted with --array=2-2%1" >&2
    exit 2
fi
if [ "$MODE" = "collect" ] && [ "$TASK_ID" != "1" ]; then
    echo "collection must be submitted with --array=1-1%1" >&2
    exit 2
fi
if [ "$MODE" = "factor" ] && { [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 5 ]; }; then
    echo "factor mode requires an array task in 1..5" >&2
    exit 2
fi
if [[ ! "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]] || [ "$RUN_TAG" = "." ] || [ "$RUN_TAG" = ".." ]; then
    echo "E5F_FINAL_SLOPE_RUN_TAG may contain only letters, digits, '.', '_', and '-'" >&2
    exit 2
fi
if [ "$MODE" = "smoke" ]; then
    if [ -n "$SMOKE_RECEIPT" ] || [ -n "$EXPECTED_SMOKE_RECEIPT_SHA" ]; then
        echo "smoke mode forbids a recursive smoke receipt" >&2
        exit 2
    fi
else
    if [ -z "$SMOKE_RECEIPT" ] || [ -z "$EXPECTED_SMOKE_RECEIPT_SHA" ]; then
        echo "$MODE mode requires a completed exact-loop smoke receipt and hash" >&2
        exit 2
    fi
fi

FRONTIER_DRIVER="$MODEL_DIR/tools/run_e5f_final_closed_root_slope_frontier.py"
CONTINUATION_DRIVER="$MODEL_DIR/tools/run_e5f_post2023_no_policy_continuations.py"
LAUNCHER_PATH="$CLUSTER_DIR/submit_e5f_final_closed_root_slope_frontier.sh"
for artifact in "$RIDGE_REPORT" "$TASK_SUMMARY" "$SOURCE" "$FRONTIER_DRIVER" "$CONTINUATION_DRIVER" "$LAUNCHER_PATH"; do
    if [ ! -f "$artifact" ]; then
        echo "required artifact is unavailable: $artifact" >&2
        exit 2
    fi
done
ACTUAL_REPORT_SHA="$(sha256sum "$RIDGE_REPORT" | awk '{print $1}')"
ACTUAL_TASK_SHA="$(sha256sum "$TASK_SUMMARY" | awk '{print $1}')"
ACTUAL_SOURCE_SHA="$(sha256sum "$SOURCE" | awk '{print $1}')"
ACTUAL_FRONTIER_DRIVER_SHA="$(sha256sum "$FRONTIER_DRIVER" | awk '{print $1}')"
ACTUAL_CONTINUATION_DRIVER_SHA="$(sha256sum "$CONTINUATION_DRIVER" | awk '{print $1}')"
ACTUAL_LAUNCHER_SHA="$(sha256sum "$LAUNCHER_PATH" | awk '{print $1}')"
if [ "$ACTUAL_REPORT_SHA" != "$EXPECTED_REPORT_SHA" ]; then
    echo "ridge-report hash mismatch: $ACTUAL_REPORT_SHA" >&2
    exit 2
fi
if [ "$ACTUAL_TASK_SHA" != "$EXPECTED_TASK_SHA" ]; then
    echo "selected-task hash mismatch: $ACTUAL_TASK_SHA" >&2
    exit 2
fi
if [ "$ACTUAL_SOURCE_SHA" != "$EXPECTED_SOURCE_SHA" ]; then
    echo "source hash mismatch: $ACTUAL_SOURCE_SHA" >&2
    exit 2
fi
if [ "$ACTUAL_FRONTIER_DRIVER_SHA" != "$EXPECTED_FRONTIER_DRIVER_SHA" ]; then
    echo "frontier-driver hash mismatch: $ACTUAL_FRONTIER_DRIVER_SHA" >&2
    exit 2
fi
if [ "$ACTUAL_CONTINUATION_DRIVER_SHA" != "$EXPECTED_CONTINUATION_DRIVER_SHA" ]; then
    echo "continuation-driver hash mismatch: $ACTUAL_CONTINUATION_DRIVER_SHA" >&2
    exit 2
fi
if [ "$ACTUAL_LAUNCHER_SHA" != "$EXPECTED_LAUNCHER_SHA" ]; then
    echo "launcher hash mismatch: $ACTUAL_LAUNCHER_SHA" >&2
    exit 2
fi
ACTUAL_DIAGNOSTIC_BUNDLE_SHA="$("$PYTHON_BIN" - "$ACTUAL_FRONTIER_DRIVER_SHA" "$ACTUAL_CONTINUATION_DRIVER_SHA" "$ACTUAL_LAUNCHER_SHA" <<'PY'
import hashlib
import json
import sys

payload = {
    "frontier_driver_sha256": sys.argv[1],
    "continuation_driver_sha256": sys.argv[2],
    "launcher_sha256": sys.argv[3],
}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
print(hashlib.sha256(encoded).hexdigest())
PY
)"
if [ "$ACTUAL_DIAGNOSTIC_BUNDLE_SHA" != "$EXPECTED_DIAGNOSTIC_BUNDLE_SHA" ]; then
    echo "diagnostic-bundle hash mismatch: $ACTUAL_DIAGNOSTIC_BUNDLE_SHA" >&2
    exit 2
fi

MODEL_OUTPUT_ROOT="$PROJECT_ROOT/output/model"
OUTPUT_ROOT="$MODEL_OUTPUT_ROOT/e5f_final_closed_root_slope_frontier_${RUN_TAG}"
if [ "$(dirname "$OUTPUT_ROOT")" != "$MODEL_OUTPUT_ROOT" ]; then
    echo "resolved output root is not a direct output/model child: $OUTPUT_ROOT" >&2
    exit 2
fi
if [ "$MODE" != "smoke" ]; then
    if [ ! -f "$SMOKE_RECEIPT" ]; then
        echo "smoke receipt is unavailable: $SMOKE_RECEIPT" >&2
        exit 2
    fi
    ACTUAL_SMOKE_RECEIPT_SHA="$(sha256sum "$SMOKE_RECEIPT" | awk '{print $1}')"
    if [ "$ACTUAL_SMOKE_RECEIPT_SHA" != "$EXPECTED_SMOKE_RECEIPT_SHA" ]; then
        echo "smoke-receipt hash mismatch: $ACTUAL_SMOKE_RECEIPT_SHA" >&2
        exit 2
    fi
    SMOKE_RESULTS_ROOT="$(cd "$(dirname "$SMOKE_RECEIPT")/.." && pwd)"
    if [ "$SMOKE_RESULTS_ROOT" = "$OUTPUT_ROOT" ]; then
        echo "production and smoke must use distinct output roots" >&2
        exit 2
    fi
fi
COMMON_ARGS=(
    --ridge-report "$RIDGE_REPORT"
    --selected-task-summary "$TASK_SUMMARY"
    --source "$SOURCE"
    --expected-ridge-report-sha256 "$EXPECTED_REPORT_SHA"
    --expected-task-summary-sha256 "$EXPECTED_TASK_SHA"
    --expected-source-sha256 "$EXPECTED_SOURCE_SHA"
    --expected-target-set "$EXPECTED_TARGET_SET"
    --expected-target-fingerprint "$EXPECTED_TARGET_FP"
    --expected-code-bundle-sha256 "$EXPECTED_CODE_SHA"
    --expected-scientific-contract-sha256 "$EXPECTED_SCIENTIFIC_SHA"
    --expected-selection-sha256 "$EXPECTED_SELECTION_SHA"
    --expected-selected-candidate-sha256 "$EXPECTED_CANDIDATE_SHA"
    --expected-selected-candidate-id "$EXPECTED_CANDIDATE_ID"
    --expected-renewal-contract-sha256 "$EXPECTED_RENEWAL_SHA"
    --expected-dated-contract-sha256 "$EXPECTED_DATED_SHA"
    --expected-frontier-driver-sha256 "$EXPECTED_FRONTIER_DRIVER_SHA"
    --expected-continuation-driver-sha256 "$EXPECTED_CONTINUATION_DRIVER_SHA"
    --expected-launcher-sha256 "$EXPECTED_LAUNCHER_SHA"
    --expected-diagnostic-bundle-sha256 "$EXPECTED_DIAGNOSTIC_BUNDLE_SHA"
    --run-stage "$MODE"
)
if [ "$MODE" != "smoke" ]; then
    COMMON_ARGS+=(
        --smoke-receipt "$SMOKE_RECEIPT"
        --expected-smoke-receipt-sha256 "$EXPECTED_SMOKE_RECEIPT_SHA"
    )
fi

export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [ "$MODE" = "collect" ]; then
    OUTDIR="$OUTPUT_ROOT/report"
    if [ -e "$OUTDIR" ]; then
        echo "refusing to overwrite collection output: $OUTDIR" >&2
        exit 2
    fi
    exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_e5f_final_closed_root_slope_frontier.py" \
        "${COMMON_ARGS[@]}" \
        --collect \
        --results-root "$OUTPUT_ROOT" \
        --output-dir "$OUTDIR"
fi

OUTDIR="$OUTPUT_ROOT/factor_$(printf '%03d' "$TASK_ID")"
if [ -e "$OUTDIR" ]; then
    echo "refusing to overwrite factor output: $OUTDIR" >&2
    exit 2
fi
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_e5f_final_closed_root_slope_frontier.py" \
    "${COMMON_ARGS[@]}" \
    --factor-index "$TASK_ID" \
    --output-dir "$OUTDIR"
