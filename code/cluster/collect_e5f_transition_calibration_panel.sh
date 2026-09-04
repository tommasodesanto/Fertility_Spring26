#!/bin/bash
# Collect one hashed E5F transition-calibration panel on Torch.
# The scientific task summaries remain the source of truth; this job only
# verifies their common contract and selects the lowest valid loss.
#SBATCH --job-name=e5ftrcollect
#SBATCH --output=logs/slurm_e5ftrcollect_%j.out
#SBATCH --error=logs/slurm_e5ftrcollect_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SUBMIT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
if [ -d "$SUBMIT_DIR/../model" ] && [ -d "$SUBMIT_DIR/../../output" ]; then
    CLUSTER_DIR="$SUBMIT_DIR"
elif [ -d "$SUBMIT_DIR/code/model" ] && [ -d "$SUBMIT_DIR/code/cluster" ]; then
    CLUSTER_DIR="$SUBMIT_DIR/code/cluster"
else
    echo "cannot resolve the project checkout from SLURM_SUBMIT_DIR=$SUBMIT_DIR" >&2
    exit 2
fi
MODEL_DIR="$(cd "${CLUSTER_DIR}/../model" && pwd)"
PROJECT_ROOT="$(cd "${CLUSTER_DIR}/../.." && pwd)"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

RUN_TAG="${E5F_TRANSITION_RUN_TAG:?E5F_TRANSITION_RUN_TAG is required}"
EXPECTED_TASKS="${E5F_TRANSITION_PANEL_SIZE:?E5F_TRANSITION_PANEL_SIZE is required}"
EXPECTED_SEED="${E5F_TRANSITION_PANEL_SEED:?E5F_TRANSITION_PANEL_SEED is required}"
EXPECTED_CODE="${E5F_TRANSITION_EXPECTED_CODE_BUNDLE_SHA256:?E5F_TRANSITION_EXPECTED_CODE_BUNDLE_SHA256 is required}"
EXPECTED_PROFILE="${E5F_TRANSITION_MODEL_PROFILE:?E5F_TRANSITION_MODEL_PROFILE is required}"
EXPECTED_CENTER_SHA="${E5F_TRANSITION_EXPECTED_CENTER_SHA256:?E5F_TRANSITION_EXPECTED_CENTER_SHA256 is required}"
EXPECTED_DESIGN="${E5F_TRANSITION_PANEL_DESIGN:-mixed}"
EXPECTED_RADIUS="${E5F_TRANSITION_LOCAL_RADIUS:-0.18}"
EXPECTED_ETA="${E5F_TRANSITION_HOUSING_SUPPLY_ELASTICITY:-}"
EXPECTED_KAPPA="${E5F_TRANSITION_FIXED_TENURE_CHOICE_KAPPA:-}"
TARGET_PROFILE="${E5F_TRANSITION_TARGET_PROFILE:-baseline}"
case "$TARGET_PROFILE" in
    baseline)
        EXPECTED_TARGET_SET="e5_fullhistory_roomsfix_h1_20260817"
        EXPECTED_TARGET_FINGERPRINT="3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497"
        ;;
    young-ownership-overid-v1)
        EXPECTED_TARGET_SET="e5_fullhistory_roomsfix_h1_20260817_plus_young_ownership_v1"
        EXPECTED_TARGET_FINGERPRINT="186d692167d2d0905b2621f5dc31a9f9edca2fe3e7f9c9a3c20d35201442f0ac"
        ;;
    *)
        echo "unsupported target profile: $TARGET_PROFILE" >&2
        exit 2
        ;;
esac

ARGS=(
    --results-dir "$PROJECT_ROOT/output/model/e5f_transition_calibration_${RUN_TAG}"
    --expected-tasks "$EXPECTED_TASKS"
    --expected-source-sha256 0afcb82d4735bd15aaa143ea04e3105a5d43df152122d02b983372102f20eef6
    --expected-target-set "$EXPECTED_TARGET_SET"
    --expected-target-fingerprint "$EXPECTED_TARGET_FINGERPRINT"
    --expected-code-bundle-sha256 "$EXPECTED_CODE"
    --expected-model-profile "$EXPECTED_PROFILE"
    --expected-panel-seed "$EXPECTED_SEED"
    --expected-center-sha256 "$EXPECTED_CENTER_SHA"
    --expected-panel-design "$EXPECTED_DESIGN"
    --expected-local-radius "$EXPECTED_RADIUS"
    --require-complete
)
if [ -n "$EXPECTED_ETA" ]; then
    ARGS+=(--expected-housing-supply-elasticity "$EXPECTED_ETA")
fi
if [ -n "$EXPECTED_KAPPA" ]; then
    ARGS+=(--expected-tenure-choice-kappa "$EXPECTED_KAPPA")
fi

export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec "$PYTHON_BIN" "$MODEL_DIR/tools/collect_e5f_transition_calibration.py" "${ARGS[@]}"
