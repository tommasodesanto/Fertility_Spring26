#!/bin/bash
# Dated-transition calibration panel for the sequential child-room-floor model.
# Each array task evaluates one reproducible candidate on the exact 2007--2023
# calendar loop and writes a complete twelve-row target table.  The baseline
# has nine transition parameters; the measured-income repair has ten.
#
# Required production launch variables:
#   E5F_TRANSITION_RUN_TAG, E5F_TRANSITION_PANEL_SIZE,
#   E5F_TRANSITION_PANEL_SEED, E5F_TRANSITION_MODEL_PROFILE,
#   E5F_TRANSITION_EXPECTED_CODE_BUNDLE_SHA256.
#   Refinement rounds may additionally set
#   E5F_TRANSITION_PANEL_CENTER_JSON and E5F_TRANSITION_LOCAL_RADIUS.
#SBATCH --job-name=e5ftrcal
#SBATCH --output=logs/slurm_e5ftrcal_%A_%a.out
#SBATCH --error=logs/slurm_e5ftrcal_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-2%2

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
mkdir -p "$SCRIPT_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

RUN_TAG="${E5F_TRANSITION_RUN_TAG:?E5F_TRANSITION_RUN_TAG is required}"
PANEL_SIZE="${E5F_TRANSITION_PANEL_SIZE:?E5F_TRANSITION_PANEL_SIZE is required}"
PANEL_SEED="${E5F_TRANSITION_PANEL_SEED:?E5F_TRANSITION_PANEL_SEED is required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
LOCAL_RADIUS="${E5F_TRANSITION_LOCAL_RADIUS:-0.18}"
PANEL_DESIGN="${E5F_TRANSITION_PANEL_DESIGN:-mixed}"
POST_2023_PERIODS="${E5F_TRANSITION_POST_2023_PERIODS:-0}"
POLICY_CASE="${E5F_TRANSITION_POLICY_CASE:-none}"
MODEL_PROFILE="${E5F_TRANSITION_MODEL_PROFILE:?E5F_TRANSITION_MODEL_PROFILE is required}"
EXPECTED_CODE_BUNDLE_SHA256="${E5F_TRANSITION_EXPECTED_CODE_BUNDLE_SHA256:?E5F_TRANSITION_EXPECTED_CODE_BUNDLE_SHA256 is required}"
if [ "$POST_2023_PERIODS" != "0" ]; then
    echo "calibration panels must stop at 2023; run continuations separately" >&2
    exit 2
fi
if [ "$POLICY_CASE" != "none" ]; then
    echo "calibration panels cannot contain a policy experiment" >&2
    exit 2
fi
if [ "$MODEL_PROFILE" != "e5f-income-entry" ]; then
    echo "production transition calibration requires e5f-income-entry" >&2
    exit 2
fi
SOURCE="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
EXPECTED_SOURCE_SHA256="0afcb82d4735bd15aaa143ea04e3105a5d43df152122d02b983372102f20eef6"
EXPECTED_TARGET_SET="e5_fullhistory_roomsfix_h1_20260817"
EXPECTED_TARGET_FINGERPRINT="3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497"
ACTUAL_SOURCE_SHA256="$(sha256sum "$SOURCE" | awk '{print $1}')"
if [ "$ACTUAL_SOURCE_SHA256" != "$EXPECTED_SOURCE_SHA256" ]; then
    echo "source hash mismatch: $ACTUAL_SOURCE_SHA256" >&2
    exit 2
fi

OUTDIR="$PROJECT_ROOT/output/model/e5f_transition_calibration_${RUN_TAG}/task_$(printf '%03d' "$TASK_ID")"
if [ -d "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
    echo "refusing to overwrite nonempty task output: $OUTDIR" >&2
    exit 2
fi
ARGS=(
    --source "$SOURCE"
    --expected-source-sha256 "$EXPECTED_SOURCE_SHA256"
    --expected-target-set "$EXPECTED_TARGET_SET"
    --expected-target-fingerprint "$EXPECTED_TARGET_FINGERPRINT"
    --expected-code-bundle-sha256 "$EXPECTED_CODE_BUNDLE_SHA256"
    --outdir "$OUTDIR"
    --panel-task-id "$TASK_ID"
    --panel-size "$PANEL_SIZE"
    --panel-seed "$PANEL_SEED"
    --panel-local-radius "$LOCAL_RADIUS"
    --panel-design "$PANEL_DESIGN"
    --post-2023-periods "$POST_2023_PERIODS"
    --policy-case "$POLICY_CASE"
    --model-profile "$MODEL_PROFILE"
    --replacement-fertility 2.1
    --old-completed-fertility-target 2.1
    --outside-origin-entry-share 0.169
)
if [ -n "${E5F_TRANSITION_PANEL_CENTER_JSON:-}" ]; then
    ARGS+=(--panel-center-json "$E5F_TRANSITION_PANEL_CENTER_JSON")
fi

export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec "$PYTHON_BIN" "$MODEL_DIR/tools/run_e5f_transition_calibration.py" "${ARGS[@]}"
