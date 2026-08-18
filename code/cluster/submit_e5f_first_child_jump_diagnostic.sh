#!/usr/bin/env bash
# Bounded diagnostic for a one-time first-child housing-service floor.
#
# The six predeclared values are external diagnostic coordinates, not
# calibrated parameters.  Every task re-solves the old steady state and the
# full 2007--2023 transition from the same byte-pinned center.  Smoke runs only
# task 3 (jump 0.20); production runs tasks 1--6 and requires that smoke receipt.
#
# Required: E5F_JUMP_MODE, E5F_JUMP_RUN_TAG, E5F_JUMP_CENTER,
# E5F_JUMP_EXPECTED_CENTER_SHA256, E5F_JUMP_SOURCE,
# E5F_JUMP_EXPECTED_SOURCE_SHA256, E5F_JUMP_EXPECTED_TARGET_SET,
# E5F_JUMP_EXPECTED_TARGET_FINGERPRINT, E5F_JUMP_EXPECTED_CODE_BUNDLE_SHA256,
# E5F_JUMP_EXPECTED_LAUNCHER_SHA256.
# Production additionally requires E5F_JUMP_SMOKE_RECEIPT and its SHA-256.
#SBATCH --job-name=e5fjump
#SBATCH --output=logs/slurm_e5fjump_%A_%a.out
#SBATCH --error=logs/slurm_e5fjump_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-6%6

set -euo pipefail

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

MODE="${E5F_JUMP_MODE:?E5F_JUMP_MODE is required}"
RUN_TAG="${E5F_JUMP_RUN_TAG:?E5F_JUMP_RUN_TAG is required}"
CENTER="${E5F_JUMP_CENTER:?E5F_JUMP_CENTER is required}"
SOURCE="${E5F_JUMP_SOURCE:?E5F_JUMP_SOURCE is required}"
EXPECTED_CENTER_SHA="${E5F_JUMP_EXPECTED_CENTER_SHA256:?expected center SHA-256 is required}"
EXPECTED_SOURCE_SHA="${E5F_JUMP_EXPECTED_SOURCE_SHA256:?expected source SHA-256 is required}"
EXPECTED_TARGET_SET="${E5F_JUMP_EXPECTED_TARGET_SET:?expected target set is required}"
EXPECTED_TARGET_FP="${E5F_JUMP_EXPECTED_TARGET_FINGERPRINT:?expected target fingerprint is required}"
EXPECTED_CODE_SHA="${E5F_JUMP_EXPECTED_CODE_BUNDLE_SHA256:?expected code bundle is required}"
EXPECTED_LAUNCHER_SHA="${E5F_JUMP_EXPECTED_LAUNCHER_SHA256:?expected launcher SHA-256 is required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
SMOKE_RECEIPT="${E5F_JUMP_SMOKE_RECEIPT:-}"
EXPECTED_SMOKE_SHA="${E5F_JUMP_EXPECTED_SMOKE_RECEIPT_SHA256:-}"

if [ "$MODE" = "smoke" ]; then
    if [ "$TASK_ID" != "3" ] || [ -n "$SMOKE_RECEIPT" ] || [ -n "$EXPECTED_SMOKE_SHA" ]; then
        echo "smoke requires --array=3-3%1 and forbids a prior smoke receipt" >&2
        exit 2
    fi
elif [ "$MODE" = "factor" ]; then
    if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 6 ] || [ -z "$SMOKE_RECEIPT" ] || [ -z "$EXPECTED_SMOKE_SHA" ]; then
        echo "factor mode requires task 1..6 and a completed smoke receipt" >&2
        exit 2
    fi
else
    echo "E5F_JUMP_MODE must be smoke or factor" >&2
    exit 2
fi

if [[ ! "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $RUN_TAG" >&2
    exit 2
fi
if [ ! -s "$CENTER" ] || [ ! -s "$SOURCE" ]; then
    echo "center or source is missing" >&2
    exit 2
fi

hash_file() {
    sha256sum "$1" | awk '{print $1}'
}
if [ "$(hash_file "$CENTER")" != "$EXPECTED_CENTER_SHA" ]; then
    echo "center hash mismatch" >&2
    exit 2
fi
if [ "$(hash_file "$SOURCE")" != "$EXPECTED_SOURCE_SHA" ]; then
    echo "source hash mismatch" >&2
    exit 2
fi
SCRIPT_PATH="$CLUSTER_DIR/submit_e5f_first_child_jump_diagnostic.sh"
if [ "$(hash_file "$SCRIPT_PATH")" != "$EXPECTED_LAUNCHER_SHA" ]; then
    echo "launcher hash mismatch" >&2
    exit 2
fi
if [ "$MODE" = "factor" ]; then
    if [ ! -s "$SMOKE_RECEIPT" ] || [ "$(hash_file "$SMOKE_RECEIPT")" != "$EXPECTED_SMOKE_SHA" ]; then
        echo "smoke receipt missing or hash mismatch" >&2
        exit 2
    fi
    "$PYTHON_BIN" - "$SMOKE_RECEIPT" "$EXPECTED_CENTER_SHA" "$EXPECTED_SOURCE_SHA" \
        "$EXPECTED_TARGET_SET" "$EXPECTED_TARGET_FP" "$EXPECTED_CODE_SHA" \
        "$EXPECTED_LAUNCHER_SHA" <<'PY'
import json
import sys
payload = json.load(open(sys.argv[1], encoding="utf-8"))
if payload.get("schema") != "e5f_first_child_jump_task_contract_v1":
    raise SystemExit("wrong smoke-receipt schema")
if payload.get("status") != "complete" or payload.get("mode") != "smoke":
    raise SystemExit("smoke receipt is not complete smoke mode")
if abs(float(payload.get("hbar_first_child_jump")) - 0.2) > 1e-14:
    raise SystemExit("smoke receipt did not evaluate jump 0.2")
for key, expected in zip(
    (
        "center_sha256",
        "source_sha256",
        "target_set",
        "target_fingerprint",
        "code_bundle_sha256",
        "launcher_sha256",
    ),
    sys.argv[2:],
    strict=True,
):
    if str(payload.get(key)) != str(expected):
        raise SystemExit(f"smoke receipt {key} mismatch")
PY
fi

JUMPS=(0.00 0.10 0.20 0.25 0.35 0.50)
JUMP="${JUMPS[$((TASK_ID - 1))]}"
OUTROOT="$PROJECT_ROOT/output/model/e5f_first_child_jump_diagnostic_${RUN_TAG}"
OUTDIR="$OUTROOT/factor_$(printf '%03d' "$TASK_ID")"
if [ -e "$OUTDIR" ]; then
    echo "refusing to overwrite output: $OUTDIR" >&2
    exit 2
fi

export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$PYTHON_BIN" "$MODEL_DIR/tools/run_e5f_transition_calibration.py" \
    --source "$SOURCE" \
    --expected-source-sha256 "$EXPECTED_SOURCE_SHA" \
    --expected-target-set "$EXPECTED_TARGET_SET" \
    --expected-target-fingerprint "$EXPECTED_TARGET_FP" \
    --expected-code-bundle-sha256 "$EXPECTED_CODE_SHA" \
    --outdir "$OUTDIR" \
    --model-profile e5f-income-entry \
    --fixed-first-child-room-jump "$JUMP" \
    --panel-task-id 1 \
    --panel-size 1 \
    --panel-seed 2026081801 \
    --panel-design mixed \
    --panel-center-json "$CENTER" \
    --post-2023-periods 0 \
    --policy-case none \
    --replacement-fertility 2.1 \
    --old-completed-fertility-target 2.1 \
    --outside-origin-entry-share 0.169

"$PYTHON_BIN" - \
    "$OUTDIR/summary.json" "$OUTDIR/jump_task_contract.json" "$MODE" "$JUMP" \
    "$EXPECTED_CENTER_SHA" "$EXPECTED_SOURCE_SHA" "$EXPECTED_TARGET_SET" \
    "$EXPECTED_TARGET_FP" "$EXPECTED_CODE_SHA" "$EXPECTED_LAUNCHER_SHA" <<'PY'
import hashlib
import json
import math
import sys
from pathlib import Path

(
    summary_raw,
    contract_raw,
    mode,
    jump_raw,
    center_sha,
    source_sha,
    target_set,
    target_fp,
    code_sha,
    launcher_sha,
) = sys.argv[1:]
summary_path = Path(summary_raw)
summary = json.loads(summary_path.read_text(encoding="utf-8"))
if summary.get("status") != "complete_transition_calibration_panel_task":
    raise SystemExit("scientific driver did not finish")
checks = (
    (summary.get("source_sha256"), source_sha, "source"),
    (summary.get("target_set"), target_set, "target set"),
    (summary.get("target_fingerprint"), target_fp, "target fingerprint"),
    (summary["code_fingerprints"]["bundle_sha256"], code_sha, "code bundle"),
    (summary["panel_design"]["center_sha256"], center_sha, "center"),
)
for actual, expected, label in checks:
    if str(actual) != str(expected):
        raise SystemExit(f"{label} mismatch")
jump = float(jump_raw)
profile = summary.get("model_profile", {})
if profile.get("first_child_room_jump_status") != "externally fixed diagnostic; not included in the free-parameter count":
    raise SystemExit("missing external-jump status")
if not math.isclose(float(profile.get("first_child_room_jump")), jump, rel_tol=0.0, abs_tol=1e-14):
    raise SystemExit("jump mismatch")
contract = {
    "schema": "e5f_first_child_jump_task_contract_v1",
    "status": "complete",
    "mode": mode,
    "hbar_first_child_jump": jump,
    "summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
    "center_sha256": center_sha,
    "source_sha256": source_sha,
    "target_set": target_set,
    "target_fingerprint": target_fp,
    "code_bundle_sha256": code_sha,
    "launcher_sha256": launcher_sha,
}
Path(contract_raw).write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
