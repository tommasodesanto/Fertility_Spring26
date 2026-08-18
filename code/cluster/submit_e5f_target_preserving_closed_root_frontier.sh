#!/bin/bash
# Target-preserving fertility-slope profile and closed-root diagnostic.
#
# Smoke: E5F_PROFILE_MODE=smoke, submit --array=3-3%1.
# Production: E5F_PROFILE_MODE=factor, submit --array=1-7%7 and require the
# completed smoke receipt.  Collection is one task with E5F_PROFILE_MODE=collect.
# Every factor re-solves the 2007 preference intercept and the 2007--2023
# preference change before auditing the terminal 25-price renewal schedule.
#SBATCH --job-name=e5frootprofile
#SBATCH --output=logs/slurm_e5frootprofile_%A_%a.out
#SBATCH --error=logs/slurm_e5frootprofile_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-7%7

set -euo pipefail

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

MODE="${E5F_PROFILE_MODE:-factor}"
RUN_TAG="${E5F_PROFILE_RUN_TAG:?E5F_PROFILE_RUN_TAG is required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
RIDGE_REPORT="${E5F_PROFILE_RIDGE_REPORT:?ridge report is required}"
TASK_SUMMARY="${E5F_PROFILE_TASK_SUMMARY:?selected task summary is required}"
SOURCE="${E5F_PROFILE_SOURCE:?stationary source is required}"
SMOKE_RECEIPT="${E5F_PROFILE_SMOKE_RECEIPT:-}"
EXPECTED_SMOKE_SHA="${E5F_PROFILE_EXPECTED_SMOKE_RECEIPT_SHA256:-}"

required=(
  E5F_PROFILE_EXPECTED_REPORT_SHA256
  E5F_PROFILE_EXPECTED_TASK_SHA256
  E5F_PROFILE_EXPECTED_SOURCE_SHA256
  E5F_PROFILE_EXPECTED_TARGET_SET
  E5F_PROFILE_EXPECTED_TARGET_FINGERPRINT
  E5F_PROFILE_EXPECTED_CODE_BUNDLE_SHA256
  E5F_PROFILE_EXPECTED_SCIENTIFIC_CONTRACT_SHA256
  E5F_PROFILE_EXPECTED_SELECTION_SHA256
  E5F_PROFILE_EXPECTED_SELECTED_CANDIDATE_SHA256
  E5F_PROFILE_EXPECTED_SELECTED_CANDIDATE_ID
  E5F_PROFILE_EXPECTED_RENEWAL_CONTRACT_SHA256
  E5F_PROFILE_EXPECTED_DATED_CONTRACT_SHA256
  E5F_PROFILE_EXPECTED_FIXED_FRONTIER_DRIVER_SHA256
  E5F_PROFILE_EXPECTED_CONTINUATION_DRIVER_SHA256
  E5F_PROFILE_EXPECTED_FIXED_FRONTIER_LAUNCHER_SHA256
  E5F_PROFILE_EXPECTED_FIXED_DIAGNOSTIC_BUNDLE_SHA256
  E5F_PROFILE_EXPECTED_DRIVER_SHA256
  E5F_PROFILE_EXPECTED_LAUNCHER_SHA256
  E5F_PROFILE_EXPECTED_BUNDLE_SHA256
)
for name in "${required[@]}"; do
  if [ -z "${!name:-}" ]; then
    echo "missing required contract variable: $name" >&2
    exit 2
  fi
done

if [ "$MODE" = "smoke" ]; then
  if [ "$TASK_ID" != "3" ] || [ -n "$SMOKE_RECEIPT" ] || [ -n "$EXPECTED_SMOKE_SHA" ]; then
    echo "smoke requires --array=3-3%1 and forbids a smoke receipt" >&2
    exit 2
  fi
elif [ "$MODE" = "factor" ]; then
  if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 7 ] || [ -z "$SMOKE_RECEIPT" ] || [ -z "$EXPECTED_SMOKE_SHA" ]; then
    echo "factor mode requires task 1..7 and a completed smoke receipt" >&2
    exit 2
  fi
elif [ "$MODE" = "collect" ]; then
  if [ "$TASK_ID" != "1" ] || [ -z "$SMOKE_RECEIPT" ] || [ -z "$EXPECTED_SMOKE_SHA" ]; then
    echo "collect mode requires --array=1-1%1 and a completed smoke receipt" >&2
    exit 2
  fi
else
  echo "E5F_PROFILE_MODE must be smoke, factor, or collect" >&2
  exit 2
fi

DRIVER="$MODEL_DIR/tools/run_e5f_target_preserving_closed_root_frontier.py"
COMMON=(
  --ridge-report "$RIDGE_REPORT"
  --selected-task-summary "$TASK_SUMMARY"
  --source "$SOURCE"
  --expected-ridge-report-sha256 "$E5F_PROFILE_EXPECTED_REPORT_SHA256"
  --expected-task-summary-sha256 "$E5F_PROFILE_EXPECTED_TASK_SHA256"
  --expected-source-sha256 "$E5F_PROFILE_EXPECTED_SOURCE_SHA256"
  --expected-target-set "$E5F_PROFILE_EXPECTED_TARGET_SET"
  --expected-target-fingerprint "$E5F_PROFILE_EXPECTED_TARGET_FINGERPRINT"
  --expected-code-bundle-sha256 "$E5F_PROFILE_EXPECTED_CODE_BUNDLE_SHA256"
  --expected-scientific-contract-sha256 "$E5F_PROFILE_EXPECTED_SCIENTIFIC_CONTRACT_SHA256"
  --expected-selection-sha256 "$E5F_PROFILE_EXPECTED_SELECTION_SHA256"
  --expected-selected-candidate-sha256 "$E5F_PROFILE_EXPECTED_SELECTED_CANDIDATE_SHA256"
  --expected-selected-candidate-id "$E5F_PROFILE_EXPECTED_SELECTED_CANDIDATE_ID"
  --expected-renewal-contract-sha256 "$E5F_PROFILE_EXPECTED_RENEWAL_CONTRACT_SHA256"
  --expected-dated-contract-sha256 "$E5F_PROFILE_EXPECTED_DATED_CONTRACT_SHA256"
  --expected-frontier-driver-sha256 "$E5F_PROFILE_EXPECTED_FIXED_FRONTIER_DRIVER_SHA256"
  --expected-continuation-driver-sha256 "$E5F_PROFILE_EXPECTED_CONTINUATION_DRIVER_SHA256"
  --expected-launcher-sha256 "$E5F_PROFILE_EXPECTED_FIXED_FRONTIER_LAUNCHER_SHA256"
  --expected-diagnostic-bundle-sha256 "$E5F_PROFILE_EXPECTED_FIXED_DIAGNOSTIC_BUNDLE_SHA256"
  --expected-profile-driver-sha256 "$E5F_PROFILE_EXPECTED_DRIVER_SHA256"
  --expected-profile-launcher-sha256 "$E5F_PROFILE_EXPECTED_LAUNCHER_SHA256"
  --expected-profile-bundle-sha256 "$E5F_PROFILE_EXPECTED_BUNDLE_SHA256"
  --run-stage "$MODE"
)
if [ "$MODE" != "smoke" ]; then
  COMMON+=(--smoke-receipt "$SMOKE_RECEIPT" --expected-smoke-receipt-sha256 "$EXPECTED_SMOKE_SHA")
fi

OUTPUT_ROOT="$PROJECT_ROOT/output/model/e5f_target_preserving_closed_root_frontier_${RUN_TAG}"
export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

if [ "$MODE" = "collect" ]; then
  exec "$PYTHON_BIN" "$DRIVER" "${COMMON[@]}" --collect --results-root "$OUTPUT_ROOT" --output-dir "$OUTPUT_ROOT/report"
fi

OUTDIR="$OUTPUT_ROOT/factor_$(printf '%03d' "$TASK_ID")"
if [ -e "$OUTDIR" ]; then
  echo "refusing to overwrite factor output: $OUTDIR" >&2
  exit 2
fi
exec "$PYTHON_BIN" "$DRIVER" "${COMMON[@]}" --factor-index "$TASK_ID" --output-dir "$OUTDIR"
