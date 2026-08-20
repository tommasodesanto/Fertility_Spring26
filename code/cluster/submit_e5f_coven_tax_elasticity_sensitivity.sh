#!/bin/bash
# Exact two-point housing-supply-elasticity sensitivity for the rebated-tax impact benchmark.
# Task 1 uses Baum-Snow--Han's average U.S. urban floor-space elasticity (0.5).
# Task 2 uses Saiz's population-weighted U.S. metropolitan elasticity (1.75).
# Both are diagnostics around the same fitted 2023 state; neither is silently
# substituted for the California benchmark used by Coven et al. (2025).
#SBATCH --job-name=e5fcovelas
#SBATCH --output=logs/slurm_e5fcovelas_%A_%a.out
#SBATCH --error=logs/slurm_e5fcovelas_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-2%2

set -euo pipefail

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
MODEL_DIR="$PROJECT_ROOT/code/model"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

RUN_TAG="${E5F_COVEN_ELASTICITY_RUN_TAG:?E5F_COVEN_ELASTICITY_RUN_TAG is required}"
DRIVER="$MODEL_DIR/tools/run_e5f_post2023_coven_property_tax_smoke.py"
EXPECTED_DRIVER_SHA256="${E5F_COVEN_EXPECTED_DRIVER_SHA256:?E5F_COVEN_EXPECTED_DRIVER_SHA256 is required}"
EXPECTED_POLICY_BUNDLE_SHA256="${E5F_COVEN_EXPECTED_POLICY_BUNDLE_SHA256:?E5F_COVEN_EXPECTED_POLICY_BUNDLE_SHA256 is required}"
SELECTED_REPORT="${E5F_COVEN_SELECTED_REPORT:?E5F_COVEN_SELECTED_REPORT is required}"
SELECTED_TRANSITION="${E5F_COVEN_SELECTED_TRANSITION:?E5F_COVEN_SELECTED_TRANSITION is required}"
SOURCE="${E5F_COVEN_SOURCE:?E5F_COVEN_SOURCE is required}"

case "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}" in
    1) ELASTICITY="0.5"; LABEL="baum_snow_han_floor_space" ;;
    2) ELASTICITY="1.75"; LABEL="saiz_population_weighted_metro" ;;
    *) echo "Array task must be 1 or 2" >&2; exit 2 ;;
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
    "$DRIVER"
)
for path in "${POLICY_BUNDLE_FILES[@]}"; do
    if [ ! -f "$path" ]; then
        echo "Policy bundle file is missing: $path" >&2
        exit 2
    fi
done
ACTUAL_POLICY_BUNDLE_SHA256="$({
    for path in "${POLICY_BUNDLE_FILES[@]}"; do
        sha256sum "$path" | awk '{print $1}'
    done
} | sha256sum | awk '{print $1}')"
if [ "$ACTUAL_POLICY_BUNDLE_SHA256" != "$EXPECTED_POLICY_BUNDLE_SHA256" ]; then
    echo "Policy bundle hash mismatch: actual=$ACTUAL_POLICY_BUNDLE_SHA256 expected=$EXPECTED_POLICY_BUNDLE_SHA256" >&2
    exit 2
fi

OUTDIR="$PROJECT_ROOT/output/model/e5f_post2023_coven_tax_elasticity_${RUN_TAG}/${LABEL}"
if [ -d "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
    echo "Refusing to overwrite nonempty output: $OUTDIR" >&2
    exit 2
fi

PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools" "$PYTHON_BIN" "$DRIVER" \
    --selected-report "$SELECTED_REPORT" \
    --selected-case-transition "$SELECTED_TRANSITION" \
    --source "$SOURCE" \
    --output-dir "$OUTDIR" \
    --expected-report-sha256 "${E5F_COVEN_EXPECTED_REPORT_SHA256:?required}" \
    --expected-case-transition-sha256 "${E5F_COVEN_EXPECTED_TRANSITION_SHA256:?required}" \
    --expected-source-sha256 "${E5F_COVEN_EXPECTED_SOURCE_SHA256:?required}" \
    --expected-target-fingerprint "${E5F_COVEN_EXPECTED_TARGET_FINGERPRINT:?required}" \
    --expected-code-bundle-sha256 "${E5F_COVEN_EXPECTED_CODE_BUNDLE_SHA256:?required}" \
    --expected-renewal-contract-sha256 "${E5F_COVEN_EXPECTED_RENEWAL_SHA256:?required}" \
    --expected-scientific-contract-sha256 "${E5F_COVEN_EXPECTED_SCIENTIFIC_SHA256:?required}" \
    --expected-selection-sha256 "${E5F_COVEN_EXPECTED_SELECTION_SHA256:?required}" \
    --housing-supply-elasticity "$ELASTICITY"

"$PYTHON_BIN" - "$OUTDIR/summary.json" "$ELASTICITY" <<'PY'
import json
import math
import sys
from pathlib import Path

path = Path(sys.argv[1])
expected = float(sys.argv[2])
payload = json.loads(path.read_text(encoding="utf-8"))
if payload.get("status") != "complete_exact_coven_property_tax_impact_smoke":
    raise SystemExit("Coven sensitivity summary is incomplete")
actual = float(payload.get("housing_supply_elasticity", math.nan))
if not math.isfinite(actual) or actual != expected:
    raise SystemExit(f"elasticity mismatch: {actual} versus {expected}")
if len(payload.get("impact_results", [])) != 2:
    raise SystemExit("Coven sensitivity must contain exactly two impact rows")
print(f"COVEN_ELASTICITY_SENSITIVITY_COMPLETE elasticity={actual} summary={path}")
PY
