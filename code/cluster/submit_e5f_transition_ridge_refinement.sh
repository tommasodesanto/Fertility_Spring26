#!/bin/bash
# Validate a hashed E5F ridge-refinement plan with the existing scientific driver.
#
# This file performs no optimization.  In validation mode (the default), array
# tasks 1--7 each run one of the seven predeclared centers.  After the
# validation-only collector writes report/selection.json, submit a separate
# repeat array with E5F_RIDGE_MODE=repeat and --array=1-2%2.  Each scientific
# invocation is deliberately the driver's one-candidate anchor mode:
# panel-task-id=1, panel-size=1, panel-design=mixed.
#
# Required in both modes:
#   E5F_RIDGE_RUN_TAG
#   E5F_RIDGE_PLAN
#   E5F_RIDGE_EXPECTED_PLAN_SHA256
#   E5F_RIDGE_EXPECTED_SOURCE_SHA256
#   E5F_RIDGE_EXPECTED_TARGET_SET
#   E5F_RIDGE_EXPECTED_TARGET_FINGERPRINT
#   E5F_RIDGE_EXPECTED_CODE_BUNDLE_SHA256
#   E5F_RIDGE_EXPECTED_MODEL_PROFILE
#   E5F_RIDGE_EXPECTED_RENEWAL_CONTRACT_SHA256
#   E5F_RIDGE_EXPECTED_DATED_CONTRACT_SHA256
# Repeat mode additionally requires:
#   E5F_RIDGE_SELECTION
#   E5F_RIDGE_EXPECTED_SELECTION_SHA256
#
# The final timing-measurement repair changes both the target fingerprint and
# scientific bundle hash.  Supply the active values; this launcher contains no
# provisional target, source, candidate, or code hash.
#SBATCH --job-name=e5ftrridge
#SBATCH --output=logs/slurm_e5ftrridge_%A_%a.out
#SBATCH --error=logs/slurm_e5ftrridge_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-7%7

set -euo pipefail

# Slurm executes a private copy of this file from its spool directory, so
# BASH_SOURCE does not identify the project checkout inside a batch job.
CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

RUN_TAG="${E5F_RIDGE_RUN_TAG:?E5F_RIDGE_RUN_TAG is required}"
PLAN="${E5F_RIDGE_PLAN:?E5F_RIDGE_PLAN is required}"
EXPECTED_PLAN_SHA256="${E5F_RIDGE_EXPECTED_PLAN_SHA256:?E5F_RIDGE_EXPECTED_PLAN_SHA256 is required}"
EXPECTED_SOURCE_SHA256="${E5F_RIDGE_EXPECTED_SOURCE_SHA256:?E5F_RIDGE_EXPECTED_SOURCE_SHA256 is required}"
EXPECTED_TARGET_SET="${E5F_RIDGE_EXPECTED_TARGET_SET:?E5F_RIDGE_EXPECTED_TARGET_SET is required}"
EXPECTED_TARGET_FINGERPRINT="${E5F_RIDGE_EXPECTED_TARGET_FINGERPRINT:?E5F_RIDGE_EXPECTED_TARGET_FINGERPRINT is required}"
EXPECTED_CODE_BUNDLE_SHA256="${E5F_RIDGE_EXPECTED_CODE_BUNDLE_SHA256:?E5F_RIDGE_EXPECTED_CODE_BUNDLE_SHA256 is required}"
EXPECTED_MODEL_PROFILE="${E5F_RIDGE_EXPECTED_MODEL_PROFILE:?E5F_RIDGE_EXPECTED_MODEL_PROFILE is required}"
EXPECTED_RENEWAL_SHA256="${E5F_RIDGE_EXPECTED_RENEWAL_CONTRACT_SHA256:?E5F_RIDGE_EXPECTED_RENEWAL_CONTRACT_SHA256 is required}"
EXPECTED_DATED_SHA256="${E5F_RIDGE_EXPECTED_DATED_CONTRACT_SHA256:?E5F_RIDGE_EXPECTED_DATED_CONTRACT_SHA256 is required}"
MODE="${E5F_RIDGE_MODE:-validation}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
SELECTION="${E5F_RIDGE_SELECTION:-}"
EXPECTED_SELECTION_SHA256="${E5F_RIDGE_EXPECTED_SELECTION_SHA256:-}"

if [ "$MODE" != "validation" ] && [ "$MODE" != "repeat" ]; then
    echo "E5F_RIDGE_MODE must be validation or repeat" >&2
    exit 2
fi
if [ "$MODE" = "validation" ]; then
    if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 7 ]; then
        echo "validation mode requires array task 1..7" >&2
        exit 2
    fi
    if [ -n "$SELECTION" ] || [ -n "$EXPECTED_SELECTION_SHA256" ]; then
        echo "validation mode refuses selection inputs" >&2
        exit 2
    fi
else
    if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt 2 ]; then
        echo "repeat mode requires an explicit --array=1-2%2 submission override" >&2
        exit 2
    fi
    : "${SELECTION:?E5F_RIDGE_SELECTION is required in repeat mode}"
    : "${EXPECTED_SELECTION_SHA256:?E5F_RIDGE_EXPECTED_SELECTION_SHA256 is required in repeat mode}"
fi

mapfile -t PLAN_FIELDS < <(
    "$PYTHON_BIN" - \
        "$PLAN" \
        "$EXPECTED_PLAN_SHA256" \
        "$MODE" \
        "$TASK_ID" \
        "$SELECTION" \
        "$EXPECTED_SELECTION_SHA256" \
        "$EXPECTED_SOURCE_SHA256" \
        "$EXPECTED_TARGET_SET" \
        "$EXPECTED_TARGET_FINGERPRINT" \
        "$EXPECTED_CODE_BUNDLE_SHA256" \
        "$EXPECTED_MODEL_PROFILE" \
        "$EXPECTED_RENEWAL_SHA256" \
        "$EXPECTED_DATED_SHA256" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

(
    plan_raw,
    expected_plan_sha,
    mode,
    task_raw,
    selection_raw,
    expected_selection_sha,
    expected_source_sha,
    expected_target_set,
    expected_target_fingerprint,
    expected_code_sha,
    expected_profile,
    expected_renewal_sha,
    expected_dated_sha,
) = sys.argv[1:]

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def emit(values):
    for value in values:
        text = str(value)
        if "\n" in text or "\r" in text:
            raise SystemExit("plan field contains a newline")
        print(text)

plan_path = Path(plan_raw).resolve()
if digest(plan_path) != expected_plan_sha:
    raise SystemExit("refinement plan hash mismatch")
plan = json.loads(plan_path.read_text(encoding="utf-8"))
if plan.get("schema") != "e5f_transition_ridge_refinement_plan_v1":
    raise SystemExit("wrong refinement-plan schema")
contract = plan["scientific_contract"]
checks = (
    (contract["source_sha256"], expected_source_sha, "source hash"),
    (contract["target_set"], expected_target_set, "target set"),
    (contract["target_fingerprint"], expected_target_fingerprint, "target fingerprint"),
    (contract["code_fingerprints"]["bundle_sha256"], expected_code_sha, "code bundle"),
    (contract["model_profile"]["name"], expected_profile, "model profile"),
    (plan["renewal_contract_sha256"], expected_renewal_sha, "renewal contract"),
    (plan["dated_contract_sha256"], expected_dated_sha, "dated contract"),
)
for actual, expected, label in checks:
    if str(actual) != str(expected):
        raise SystemExit(f"{label} mismatch: {actual} versus {expected}")
if len(plan.get("candidates", [])) != 7:
    raise SystemExit("plan does not contain seven candidates")

outer_task = int(task_raw)
selection_sha = ""
if mode == "validation":
    candidate_id = outer_task
else:
    selection_path = Path(selection_raw).resolve()
    if digest(selection_path) != expected_selection_sha:
        raise SystemExit("selection receipt hash mismatch")
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if selection.get("schema") != "e5f_transition_ridge_selection_v1":
        raise SystemExit("wrong selection-receipt schema")
    if selection.get("plan_sha256") != expected_plan_sha:
        raise SystemExit("selection receipt belongs to a different plan")
    candidate_id = int(selection["selected_candidate_id"])
    selection_sha = expected_selection_sha

candidate = next(
    (row for row in plan["candidates"] if int(row["candidate_id"]) == candidate_id),
    None,
)
if candidate is None:
    raise SystemExit("selected candidate id is absent from plan")
candidate_path = (plan_path.parent / candidate["candidate_file"]).resolve()
if digest(candidate_path) != candidate["candidate_sha256"]:
    raise SystemExit("candidate file hash mismatch")
if mode == "repeat":
    if selection["selected_candidate_sha256"] != candidate["candidate_sha256"]:
        raise SystemExit("selection candidate hash differs from plan")

emit(
    (
        candidate_path,
        candidate["candidate_sha256"],
        candidate_id,
        contract["source"],
        contract["source_sha256"],
        contract["target_set"],
        contract["target_fingerprint"],
        contract["code_fingerprints"]["bundle_sha256"],
        contract["model_profile"]["name"],
        plan["coordinate_panel_contract"]["panel_seed"],
        contract["outside_origin_entry_share"],
        contract["old_completed_fertility_reference"],
        plan["scientific_contract_sha256"],
        plan["renewal_contract_sha256"],
        plan["dated_contract_sha256"],
        selection_sha,
    )
)
PY
)

if [ "${#PLAN_FIELDS[@]}" -ne 16 ]; then
    echo "failed to extract the exact refinement contract" >&2
    exit 2
fi
CANDIDATE_PATH="${PLAN_FIELDS[0]}"
CANDIDATE_SHA256="${PLAN_FIELDS[1]}"
CANDIDATE_ID="${PLAN_FIELDS[2]}"
SOURCE="${PLAN_FIELDS[3]}"
SOURCE_SHA256="${PLAN_FIELDS[4]}"
TARGET_SET="${PLAN_FIELDS[5]}"
TARGET_FINGERPRINT="${PLAN_FIELDS[6]}"
CODE_BUNDLE_SHA256="${PLAN_FIELDS[7]}"
MODEL_PROFILE="${PLAN_FIELDS[8]}"
PANEL_SEED="${PLAN_FIELDS[9]}"
OUTSIDE_ORIGIN_ENTRY_SHARE="${PLAN_FIELDS[10]}"
OLD_FERTILITY_REFERENCE="${PLAN_FIELDS[11]}"
SCIENTIFIC_CONTRACT_SHA256="${PLAN_FIELDS[12]}"
RENEWAL_SHA256="${PLAN_FIELDS[13]}"
DATED_SHA256="${PLAN_FIELDS[14]}"
SELECTION_SHA256="${PLAN_FIELDS[15]}"

if [ ! -f "$SOURCE" ]; then
    echo "plan source is unavailable on Torch: $SOURCE" >&2
    exit 2
fi
ACTUAL_SOURCE_SHA256="$(sha256sum "$SOURCE" | awk '{print $1}')"
if [ "$ACTUAL_SOURCE_SHA256" != "$SOURCE_SHA256" ]; then
    echo "source hash mismatch: actual=$ACTUAL_SOURCE_SHA256 expected=$SOURCE_SHA256" >&2
    exit 2
fi

OUTPUT_ROOT="$PROJECT_ROOT/output/model/e5f_transition_ridge_refinement_${RUN_TAG}"
if [ "$MODE" = "validation" ]; then
    OUTDIR="$OUTPUT_ROOT/task_$(printf '%03d' "$TASK_ID")"
    REPLICATE_ID=""
else
    OUTDIR="$OUTPUT_ROOT/repeat_$(printf '%03d' "$TASK_ID")"
    REPLICATE_ID="$TASK_ID"
fi
if [ -d "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
    echo "refusing to overwrite nonempty task output: $OUTDIR" >&2
    exit 2
fi

export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
"$PYTHON_BIN" "$MODEL_DIR/tools/run_e5f_transition_calibration.py" \
    --source "$SOURCE" \
    --expected-source-sha256 "$SOURCE_SHA256" \
    --expected-target-set "$TARGET_SET" \
    --expected-target-fingerprint "$TARGET_FINGERPRINT" \
    --expected-code-bundle-sha256 "$CODE_BUNDLE_SHA256" \
    --outdir "$OUTDIR" \
    --panel-task-id 1 \
    --panel-size 1 \
    --panel-seed "$PANEL_SEED" \
    --panel-local-radius 0.02 \
    --panel-design mixed \
    --panel-center-json "$CANDIDATE_PATH" \
    --post-2023-periods 0 \
    --policy-case none \
    --model-profile "$MODEL_PROFILE" \
    --replacement-fertility "$OLD_FERTILITY_REFERENCE" \
    --old-completed-fertility-target "$OLD_FERTILITY_REFERENCE" \
    --outside-origin-entry-share "$OUTSIDE_ORIGIN_ENTRY_SHARE"

if [ ! -s "$OUTDIR/summary.json" ] || [ ! -s "$OUTDIR/target_fit_long.csv" ]; then
    echo "scientific driver returned without complete core artifacts" >&2
    exit 2
fi

EXECUTION_IDENTITY="${SLURM_JOB_ID:-local}:${SLURM_ARRAY_JOB_ID:-local}:${TASK_ID}:${HOSTNAME:-unknown}:$MODE"
"$PYTHON_BIN" - \
    "$OUTDIR/refinement_task_contract.json" \
    "$MODE" \
    "$TASK_ID" \
    "$REPLICATE_ID" \
    "$CANDIDATE_ID" \
    "$CANDIDATE_SHA256" \
    "$EXPECTED_PLAN_SHA256" \
    "$SCIENTIFIC_CONTRACT_SHA256" \
    "$SOURCE_SHA256" \
    "$TARGET_SET" \
    "$TARGET_FINGERPRINT" \
    "$CODE_BUNDLE_SHA256" \
    "$RENEWAL_SHA256" \
    "$DATED_SHA256" \
    "$SELECTION_SHA256" \
    "$EXECUTION_IDENTITY" <<'PY'
import json
import sys
from pathlib import Path

(
    output_raw,
    mode,
    outer_task,
    replicate,
    candidate_id,
    candidate_sha,
    plan_sha,
    scientific_sha,
    source_sha,
    target_set,
    target_fingerprint,
    code_sha,
    renewal_sha,
    dated_sha,
    selection_sha,
    execution_identity,
) = sys.argv[1:]
payload = {
    "schema": "e5f_transition_ridge_task_contract_v1",
    "mode": mode,
    "outer_task_id": int(outer_task),
    "replicate_id": int(replicate) if replicate else None,
    "candidate_id": int(candidate_id),
    "candidate_sha256": candidate_sha,
    "plan_sha256": plan_sha,
    "scientific_contract_sha256": scientific_sha,
    "source_sha256": source_sha,
    "target_set": target_set,
    "target_fingerprint": target_fingerprint,
    "code_bundle_sha256": code_sha,
    "renewal_contract_sha256": renewal_sha,
    "dated_contract_sha256": dated_sha,
    "selection_sha256": selection_sha or None,
    "execution_identity": execution_identity,
}
path = Path(output_raw)
temporary = path.with_suffix(path.suffix + ".tmp")
temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
temporary.replace(path)
PY

echo "E5F_RIDGE_TASK_COMPLETE mode=$MODE outer_task=$TASK_ID candidate=$CANDIDATE_ID output=$OUTDIR"
