#!/usr/bin/env bash
# Run the final paired closed/open no-policy continuation from a certified E5F
# ridge report.  The same launcher has two modes:
#
#   smoke       one post-2023 date, but the complete two-branch, closed-root,
#               open-endpoint, figure, manifest, and final-audit loop;
#   production  forty post-2023 dates (through 2183).  Production refuses to
#               start without a hashed complete receipt from smoke mode under
#               the identical input and launcher contract.
#
# No report, repeat, source, target, model-bundle, selection, or driver hash is
# embedded here.  Every final artifact is supplied and byte-pinned at launch.
# See the REQUIRED array below for the common contract.
#SBATCH --job-name=e5fcont
#SBATCH --output=logs/slurm_e5fcont_%j.out
#SBATCH --error=logs/slurm_e5fcont_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

# Slurm executes a private spool copy named ``slurm_script``.  Keep the
# repository filename explicit so the batch process hashes and validates the
# submitted source rather than looking for the transient spool name under
# SLURM_SUBMIT_DIR.
SCRIPT_NAME="submit_e5f_post2023_no_policy_continuations.sh"
LOCAL_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REQUIRED=(
    E5F_CONT_MODE
    E5F_CONT_RUN_TAG
    E5F_CONT_REPORT
    E5F_CONT_EXPECTED_REPORT_SHA256
    E5F_CONT_SELECTION
    E5F_CONT_EXPECTED_SELECTION_SHA256
    E5F_CONT_REPEAT1_SUMMARY
    E5F_CONT_EXPECTED_REPEAT1_SUMMARY_SHA256
    E5F_CONT_REPEAT1_CONTRACT
    E5F_CONT_EXPECTED_REPEAT1_CONTRACT_SHA256
    E5F_CONT_REPEAT2_SUMMARY
    E5F_CONT_EXPECTED_REPEAT2_SUMMARY_SHA256
    E5F_CONT_REPEAT2_CONTRACT
    E5F_CONT_EXPECTED_REPEAT2_CONTRACT_SHA256
    E5F_CONT_SELECTED_TRANSITION
    E5F_CONT_EXPECTED_SELECTED_TRANSITION_SHA256
    E5F_CONT_SOURCE
    E5F_CONT_EXPECTED_SOURCE_SHA256
    E5F_CONT_EXPECTED_TARGET_SET
    E5F_CONT_EXPECTED_TARGET_FINGERPRINT
    E5F_CONT_EXPECTED_CODE_BUNDLE_SHA256
    E5F_CONT_EXPECTED_RENEWAL_CONTRACT_SHA256
    E5F_CONT_EXPECTED_SCIENTIFIC_CONTRACT_SHA256
    E5F_CONT_EXPECTED_DRIVER_SHA256
    E5F_CONT_EXPECTED_LAUNCHER_SHA256
)

require_environment() {
    local name
    for name in "${REQUIRED[@]}"; do
        if [ -z "${!name:-}" ]; then
            echo "$name is required" >&2
            exit 2
        fi
    done
    if [ "$E5F_CONT_MODE" != "smoke" ] && [ "$E5F_CONT_MODE" != "production" ]; then
        echo "E5F_CONT_MODE must be smoke or production" >&2
        exit 2
    fi
    if [ "$E5F_CONT_MODE" = "production" ]; then
        : "${E5F_CONT_SMOKE_RECEIPT:?E5F_CONT_SMOKE_RECEIPT is required in production}"
        : "${E5F_CONT_EXPECTED_SMOKE_RECEIPT_SHA256:?expected smoke-receipt SHA-256 is required in production}"
    elif [ -n "${E5F_CONT_SMOKE_RECEIPT:-}" ] || [ -n "${E5F_CONT_EXPECTED_SMOKE_RECEIPT_SHA256:-}" ]; then
        echo "smoke mode refuses production smoke-receipt inputs" >&2
        exit 2
    fi
}

if [ "${1:-}" = "--submit" ]; then
    shift
    require_environment
    if [ "$E5F_CONT_MODE" = "smoke" ]; then
        TIME_LIMIT="${E5F_CONT_TIME_LIMIT:-00:40:00}"
        MEMORY="${E5F_CONT_MEMORY:-8G}"
    else
        TIME_LIMIT="${E5F_CONT_TIME_LIMIT:-02:00:00}"
        MEMORY="${E5F_CONT_MEMORY:-8G}"
    fi
    cd "$LOCAL_SCRIPT_DIR"
    exec sbatch --time="$TIME_LIMIT" --mem="$MEMORY" --export=ALL \
        "$LOCAL_SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Run with --submit or submit this file through sbatch}"
require_environment

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$LOCAL_SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
MODEL_DIR="$PROJECT_ROOT/code/model"
SCRIPT_PATH="$CLUSTER_DIR/$SCRIPT_NAME"
DRIVER="$MODEL_DIR/tools/run_e5f_post2023_no_policy_continuations.py"
mkdir -p "$CLUSTER_DIR/logs"

RUN_TAG="$E5F_CONT_RUN_TAG"
if [[ ! "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $RUN_TAG" >&2
    exit 2
fi

resolve_path() {
    local value="$1"
    if [[ "$value" = /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$PROJECT_ROOT/$value"
    fi
}

REPORT="$(resolve_path "$E5F_CONT_REPORT")"
SELECTION="$(resolve_path "$E5F_CONT_SELECTION")"
REPEAT1_SUMMARY="$(resolve_path "$E5F_CONT_REPEAT1_SUMMARY")"
REPEAT1_CONTRACT="$(resolve_path "$E5F_CONT_REPEAT1_CONTRACT")"
REPEAT2_SUMMARY="$(resolve_path "$E5F_CONT_REPEAT2_SUMMARY")"
REPEAT2_CONTRACT="$(resolve_path "$E5F_CONT_REPEAT2_CONTRACT")"
SELECTED_TRANSITION="$(resolve_path "$E5F_CONT_SELECTED_TRANSITION")"
SOURCE="$(resolve_path "$E5F_CONT_SOURCE")"
if [ -n "${E5F_CONT_OUTDIR:-}" ]; then
    OUTDIR="$(resolve_path "$E5F_CONT_OUTDIR")"
else
    OUTDIR="$PROJECT_ROOT/output/model/e5f_post2023_no_policy_continuations_${RUN_TAG}"
fi
if [ "$E5F_CONT_MODE" = "production" ]; then
    SMOKE_RECEIPT="$(resolve_path "$E5F_CONT_SMOKE_RECEIPT")"
    POST_PERIODS=40
else
    SMOKE_RECEIPT=""
    POST_PERIODS=1
fi

HEARTBEAT_SECONDS="${E5F_CONT_HEARTBEAT_SECONDS:-300}"
if [[ ! "$HEARTBEAT_SECONDS" =~ ^[0-9]+$ ]] || [ "$HEARTBEAT_SECONDS" -lt 30 ]; then
    echo "E5F_CONT_HEARTBEAT_SECONDS must be an integer of at least 30" >&2
    exit 2
fi

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
OUTDIR="$($PYTHON_BIN -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve())' "$OUTDIR")"
case "$OUTDIR" in
    "$PROJECT_ROOT"/output/model/*) ;;
    *)
        echo "output must be a named child of $PROJECT_ROOT/output/model" >&2
        exit 2
        ;;
esac
hash_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

INPUT_PATHS=(
    "$REPORT" "$SELECTION" "$REPEAT1_SUMMARY" "$REPEAT1_CONTRACT"
    "$REPEAT2_SUMMARY" "$REPEAT2_CONTRACT" "$SELECTED_TRANSITION" "$SOURCE"
    "$DRIVER" "$SCRIPT_PATH"
)
for path in "${INPUT_PATHS[@]}"; do
    if [ ! -s "$path" ]; then
        echo "required input is missing or empty: $path" >&2
        exit 2
    fi
done

REPORT_SHA="$(hash_file "$REPORT")"
SELECTION_SHA="$(hash_file "$SELECTION")"
REPEAT1_SUMMARY_SHA="$(hash_file "$REPEAT1_SUMMARY")"
REPEAT1_CONTRACT_SHA="$(hash_file "$REPEAT1_CONTRACT")"
REPEAT2_SUMMARY_SHA="$(hash_file "$REPEAT2_SUMMARY")"
REPEAT2_CONTRACT_SHA="$(hash_file "$REPEAT2_CONTRACT")"
SELECTED_TRANSITION_SHA="$(hash_file "$SELECTED_TRANSITION")"
SOURCE_SHA="$(hash_file "$SOURCE")"
DRIVER_SHA="$(hash_file "$DRIVER")"
LAUNCHER_SHA="$(hash_file "$SCRIPT_PATH")"

check_hash() {
    local label="$1" actual="$2" expected="$3"
    if [ "$actual" != "$expected" ]; then
        echo "$label hash mismatch: actual=$actual expected=$expected" >&2
        exit 2
    fi
}
check_hash report "$REPORT_SHA" "$E5F_CONT_EXPECTED_REPORT_SHA256"
check_hash selection "$SELECTION_SHA" "$E5F_CONT_EXPECTED_SELECTION_SHA256"
check_hash repeat1-summary "$REPEAT1_SUMMARY_SHA" "$E5F_CONT_EXPECTED_REPEAT1_SUMMARY_SHA256"
check_hash repeat1-contract "$REPEAT1_CONTRACT_SHA" "$E5F_CONT_EXPECTED_REPEAT1_CONTRACT_SHA256"
check_hash repeat2-summary "$REPEAT2_SUMMARY_SHA" "$E5F_CONT_EXPECTED_REPEAT2_SUMMARY_SHA256"
check_hash repeat2-contract "$REPEAT2_CONTRACT_SHA" "$E5F_CONT_EXPECTED_REPEAT2_CONTRACT_SHA256"
check_hash selected-transition "$SELECTED_TRANSITION_SHA" "$E5F_CONT_EXPECTED_SELECTED_TRANSITION_SHA256"
check_hash source "$SOURCE_SHA" "$E5F_CONT_EXPECTED_SOURCE_SHA256"
check_hash continuation-driver "$DRIVER_SHA" "$E5F_CONT_EXPECTED_DRIVER_SHA256"
check_hash continuation-launcher "$LAUNCHER_SHA" "$E5F_CONT_EXPECTED_LAUNCHER_SHA256"

# Re-run the driver's complete read-only input contract and independently tie
# the two byte-pinned repeat summaries/wrappers to the final ridge report.
export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"
"$PYTHON_BIN" - \
    "$PROJECT_ROOT" "$REPORT" "$SELECTION" \
    "$REPEAT1_SUMMARY" "$REPEAT1_CONTRACT" \
    "$REPEAT2_SUMMARY" "$REPEAT2_CONTRACT" \
    "$SELECTED_TRANSITION" "$SOURCE" \
    "$REPORT_SHA" "$SELECTION_SHA" "$SELECTED_TRANSITION_SHA" "$SOURCE_SHA" \
    "$E5F_CONT_EXPECTED_TARGET_SET" "$E5F_CONT_EXPECTED_TARGET_FINGERPRINT" \
    "$E5F_CONT_EXPECTED_CODE_BUNDLE_SHA256" \
    "$E5F_CONT_EXPECTED_RENEWAL_CONTRACT_SHA256" \
    "$E5F_CONT_EXPECTED_SCIENTIFIC_CONTRACT_SHA256" <<'PY'
import json
import math
import sys
from pathlib import Path

(root_raw, report_raw, selection_raw, r1_summary_raw, r1_contract_raw,
 r2_summary_raw, r2_contract_raw, transition_raw, source_raw, report_sha,
 selection_sha, transition_sha, source_sha, target_set, target_fingerprint,
 code_sha, renewal_sha, scientific_sha) = sys.argv[1:]
root = Path(root_raw).resolve()
sys.path.insert(0, str(root / "code/model"))
sys.path.insert(0, str(root / "code/model/tools"))
import run_e5f_open_population_transition as transition
import run_e5f_post2023_no_policy_continuations as driver

def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

report = read(report_raw)
selection = read(selection_raw)
repeat_summaries = [read(r1_summary_raw), read(r2_summary_raw)]
repeat_contracts = [read(r1_contract_raw), read(r2_contract_raw)]
if report.get("schema") != "e5f_transition_ridge_refinement_report_v1":
    raise SystemExit("final report has the wrong schema")
if report.get("status") != "complete_refinement_with_two_independent_identity_repeats":
    raise SystemExit("final report is not complete with two repeats")
if not report.get("promotion_eligible"):
    raise SystemExit("final report is not promotion eligible")
if report.get("selection_sha256") != selection_sha:
    raise SystemExit("selection does not belong to the final report")
if selection.get("schema") != "e5f_transition_ridge_selection_v1":
    raise SystemExit("selection has the wrong schema")
scientific = report.get("scientific_contract") or {}
checks = (
    (scientific.get("target_set"), target_set, "target set"),
    (scientific.get("target_fingerprint"), target_fingerprint, "target fingerprint"),
    ((scientific.get("code_fingerprints") or {}).get("bundle_sha256"), code_sha, "code bundle"),
    (scientific.get("source_sha256"), source_sha, "source"),
    (report.get("renewal_contract_sha256"), renewal_sha, "renewal contract"),
    (report.get("scientific_contract_sha256"), scientific_sha, "scientific contract"),
    (report.get("selected_candidate_sha256"), selection.get("selected_candidate_sha256"), "selected candidate"),
)
for actual, expected, label in checks:
    if str(actual) != str(expected):
        raise SystemExit(f"{label} mismatch: {actual} versus {expected}")

best = report.get("best_candidate") or {}
identities = list((report.get("repeat_gate") or {}).get("execution_identities") or [])
if len(identities) != 2 or len(set(identities)) != 2:
    raise SystemExit("final report does not name two distinct repeat executions")
for index, (summary, wrapper) in enumerate(zip(repeat_summaries, repeat_contracts), start=1):
    if summary.get("status") != "complete_transition_calibration_panel_task":
        raise SystemExit(f"repeat {index} summary is incomplete")
    if wrapper.get("schema") != "e5f_transition_ridge_task_contract_v1":
        raise SystemExit(f"repeat {index} wrapper has the wrong schema")
    if wrapper.get("mode") != "repeat" or int(wrapper.get("replicate_id", -1)) != index:
        raise SystemExit(f"repeat {index} wrapper has the wrong repeat identity")
    if wrapper.get("execution_identity") != identities[index - 1]:
        raise SystemExit(f"repeat {index} execution identity differs from the final report")
    if wrapper.get("selection_sha256") != selection_sha:
        raise SystemExit(f"repeat {index} selection hash mismatch")
    if wrapper.get("candidate_sha256") != report.get("selected_candidate_sha256"):
        raise SystemExit(f"repeat {index} candidate hash mismatch")
    for field, expected in (
        ("source_sha256", source_sha),
        ("target_fingerprint", target_fingerprint),
        ("code_bundle_sha256", code_sha),
        ("renewal_contract_sha256", renewal_sha),
        ("scientific_contract_sha256", scientific_sha),
    ):
        if str(wrapper.get(field)) != str(expected):
            raise SystemExit(f"repeat {index} {field} mismatch")
    candidate = summary.get("best_candidate") or {}
    if (summary.get("panel_design") or {}).get("center_sha256") != report.get("selected_candidate_sha256"):
        raise SystemExit(f"repeat {index} did not evaluate the selected center")
    if set(candidate.get("theta") or {}) != set(best.get("theta") or {}):
        raise SystemExit(f"repeat {index} theta keys differ")
    for name, value in (best.get("theta") or {}).items():
        if not math.isclose(float(candidate["theta"][name]), float(value), rel_tol=0.0, abs_tol=5e-12):
            raise SystemExit(f"repeat {index} theta differs at {name}")
    if not math.isclose(float(candidate.get("transition_loss")), float(best.get("transition_loss")), rel_tol=0.0, abs_tol=1e-9):
        raise SystemExit(f"repeat {index} loss differs")

chain, model = transition.configure_sequential_model()
driver.validate_input_contracts(
    report_path=Path(report_raw),
    task_summary_path=None,
    case_dir=None,
    case_transition_path=Path(transition_raw),
    source_path=Path(source_raw),
    expected_report_sha256=report_sha,
    expected_task_summary_sha256=None,
    expected_case_transition_sha256=transition_sha,
    expected_source_sha256=source_sha,
    expected_target_fingerprint=target_fingerprint,
    expected_code_bundle_sha256=code_sha,
    expected_renewal_contract_sha256=renewal_sha,
    expected_scientific_contract_sha256=scientific_sha,
    expected_selection_sha256=selection_sha,
    chain=chain,
    model=model,
)
PY

INPUT_CONTRACT_JSON="$(mktemp "$PROJECT_ROOT/output/model/.e5f_cont_input_contract.XXXXXX")"
cleanup_input_contract() {
    if [ -n "${INPUT_CONTRACT_JSON:-}" ] && [ -f "$INPUT_CONTRACT_JSON" ]; then
        rm -f "$INPUT_CONTRACT_JSON"
    fi
}
trap cleanup_input_contract EXIT
"$PYTHON_BIN" - \
    "$INPUT_CONTRACT_JSON" "$REPORT_SHA" "$SELECTION_SHA" \
    "$REPEAT1_SUMMARY_SHA" "$REPEAT1_CONTRACT_SHA" \
    "$REPEAT2_SUMMARY_SHA" "$REPEAT2_CONTRACT_SHA" \
    "$SELECTED_TRANSITION_SHA" "$SOURCE_SHA" \
    "$E5F_CONT_EXPECTED_TARGET_SET" "$E5F_CONT_EXPECTED_TARGET_FINGERPRINT" \
    "$E5F_CONT_EXPECTED_CODE_BUNDLE_SHA256" \
    "$E5F_CONT_EXPECTED_RENEWAL_CONTRACT_SHA256" \
    "$E5F_CONT_EXPECTED_SCIENTIFIC_CONTRACT_SHA256" \
    "$DRIVER_SHA" "$LAUNCHER_SHA" <<'PY'
import json, sys
path, *values = sys.argv[1:]
names = (
    "report_sha256", "selection_sha256", "repeat1_summary_sha256",
    "repeat1_contract_sha256", "repeat2_summary_sha256",
    "repeat2_contract_sha256", "selected_transition_sha256", "source_sha256",
    "target_set", "target_fingerprint", "code_bundle_sha256",
    "renewal_contract_sha256", "scientific_contract_sha256", "driver_sha256",
    "launcher_sha256",
)
with open(path, "w", encoding="utf-8") as handle:
    json.dump(dict(zip(names, values)), handle, sort_keys=True)
PY

if [ "$E5F_CONT_MODE" = "production" ]; then
    if [ ! -s "$SMOKE_RECEIPT" ]; then
        echo "smoke receipt is missing or empty: $SMOKE_RECEIPT" >&2
        exit 2
    fi
    SMOKE_RECEIPT_SHA="$(hash_file "$SMOKE_RECEIPT")"
    check_hash smoke-receipt "$SMOKE_RECEIPT_SHA" "$E5F_CONT_EXPECTED_SMOKE_RECEIPT_SHA256"
    "$PYTHON_BIN" - "$SMOKE_RECEIPT" "$INPUT_CONTRACT_JSON" "$SLURM_JOB_ID" <<'PY'
import hashlib, json, sys
from pathlib import Path
receipt_path, contract_path, job_id = sys.argv[1:]
receipt = json.loads(Path(receipt_path).read_text(encoding="utf-8"))
contract = json.loads(Path(contract_path).read_text(encoding="utf-8"))
if receipt.get("status") != "complete" or receipt.get("mode") != "smoke":
    raise SystemExit("smoke receipt is not a completed smoke run")
if receipt.get("input_contract") != contract:
    raise SystemExit("smoke receipt belongs to a different scientific or launcher contract")
if int(receipt.get("post_2023_periods", -1)) != 1 or not receipt.get("artifact_gate_passed"):
    raise SystemExit("smoke receipt did not complete the exact one-date loop and artifact gate")
if str(receipt.get("slurm_job_id")) == str(job_id):
    raise SystemExit("production and smoke must have distinct Slurm executions")
manifest = Path(str(receipt.get("driver_manifest", "")))
if not manifest.is_file():
    raise SystemExit("smoke receipt's driver manifest is unavailable")
actual = hashlib.sha256(manifest.read_bytes()).hexdigest()
if actual != receipt.get("driver_manifest_sha256"):
    raise SystemExit("smoke driver manifest changed after certification")
PY
else
    SMOKE_RECEIPT_SHA=""
fi

if [ -e "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    echo "refusing to overwrite nonempty output directory: $OUTDIR" >&2
    exit 2
fi
mkdir -p "$OUTDIR/checkpoints"
mv "$INPUT_CONTRACT_JSON" "$OUTDIR/input_contract.json"
INPUT_CONTRACT_JSON=""

REQUESTED_TIME="${E5F_CONT_TIME_LIMIT:-$([ "$E5F_CONT_MODE" = smoke ] && echo 00:40:00 || echo 02:00:00)}"
REQUESTED_MEMORY="${E5F_CONT_MEMORY:-8G}"
"$PYTHON_BIN" - \
    "$OUTDIR/launch_contract.json" "$OUTDIR/input_contract.json" \
    "$E5F_CONT_MODE" "$RUN_TAG" "$POST_PERIODS" "$SLURM_JOB_ID" \
    "$REQUESTED_TIME" "$REQUESTED_MEMORY" "$SMOKE_RECEIPT" "$SMOKE_RECEIPT_SHA" <<'PY'
import json, os, sys
from datetime import datetime, timezone
path, contract_path, mode, tag, periods, job, walltime, memory, smoke, smoke_sha = sys.argv[1:]
with open(contract_path, encoding="utf-8") as handle:
    contract = json.load(handle)
p = int(periods)
payload = {
    "status": "launched",
    "launched_utc": datetime.now(timezone.utc).isoformat(),
    "mode": mode,
    "run_tag": tag,
    "slurm_job_id": job,
    "input_contract": contract,
    "post_2023_periods": p,
    "terminal_calendar_year": 2023 + 4 * p,
    "dynamic_market_evaluations": 5 + 2 * p,
    "paired_output_rows": 2 * (5 + p),
    "closed_schedule_declared_grid_points": 25,
    "run_size": f"one task; 5 shared dates through 2023; two no-policy branches x {p} later dates; full closed and open endpoint searches",
    "rough_walltime": "production: roughly 35--50 minutes from prior comparable 40-date work; smoke determines whether that estimate remains credible",
    "requested_time": walltime,
    "requested_memory": memory,
    "slurm_mem_per_node": os.environ.get("SLURM_MEM_PER_NODE"),
    "heartbeat_seconds": int(os.environ.get("E5F_CONT_HEARTBEAT_SECONDS", "300")),
    "smoke_receipt": smoke or None,
    "smoke_receipt_sha256": smoke_sha or None,
    "stop_criteria": "complete both no-policy paths and both endpoint audits; reproduce the selected five-date history; pass input, identity, market, mass, finite, artifact, and manifest hashes",
}
tmp = path + ".tmp"
with open(tmp, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(tmp, path)
PY

export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
DRIVER_LOG="$OUTDIR/driver.log"
ARGS=(
    --selected-report "$REPORT"
    --selected-case-transition "$SELECTED_TRANSITION"
    --source "$SOURCE"
    --output-dir "$OUTDIR"
    --expected-report-sha256 "$REPORT_SHA"
    --expected-case-transition-sha256 "$SELECTED_TRANSITION_SHA"
    --expected-source-sha256 "$SOURCE_SHA"
    --expected-target-fingerprint "$E5F_CONT_EXPECTED_TARGET_FINGERPRINT"
    --expected-code-bundle-sha256 "$E5F_CONT_EXPECTED_CODE_BUNDLE_SHA256"
    --expected-renewal-contract-sha256 "$E5F_CONT_EXPECTED_RENEWAL_CONTRACT_SHA256"
    --expected-scientific-contract-sha256 "$E5F_CONT_EXPECTED_SCIENTIFIC_CONTRACT_SHA256"
    --expected-selection-sha256 "$SELECTION_SHA"
    --post-2023-periods "$POST_PERIODS"
    --market-tol 2e-4
    --market-max-iter 30
    --closed-price-min-ratio 1e-4
    --closed-price-max-ratio 3.0
    --closed-grid-points 25
)

snapshot_progress() {
    "$PYTHON_BIN" - "$OUTDIR" "$SLURM_JOB_ID" <<'PY'
import hashlib, json, os, shutil, sys
from datetime import datetime, timezone
outdir, job = sys.argv[1:]
checkpoints = os.path.join(outdir, "checkpoints")
os.makedirs(checkpoints, exist_ok=True)

def read_json(path):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except Exception as error:
        return {"snapshot_error": repr(error), "path": path}

latest_period_path = os.path.join(outdir, "latest_completed_period.json")
latest_endpoint_path = os.path.join(outdir, "latest_endpoint_search.json")
period = read_json(latest_period_path) if os.path.isfile(latest_period_path) else None
endpoint = read_json(latest_endpoint_path) if os.path.isfile(latest_endpoint_path) else None
if period and "snapshot_error" not in period:
    token = str(period.get("calendar_year") or period.get("terminal_calendar_year") or "state")
    destination = os.path.join(checkpoints, f"continuation_{token}.json")
    if not os.path.exists(destination):
        shutil.copy2(latest_period_path, destination)
if endpoint and "snapshot_error" not in endpoint:
    encoded = json.dumps(endpoint, sort_keys=True).encode()
    token = hashlib.sha256(encoded).hexdigest()[:12]
    destination = os.path.join(checkpoints, f"endpoint_{token}.json")
    if not os.path.exists(destination):
        shutil.copy2(latest_endpoint_path, destination)
payload = {
    "status": "running",
    "utc": datetime.now(timezone.utc).isoformat(),
    "slurm_job_id": job,
    "latest_completed_period": period,
    "latest_endpoint_search": endpoint,
    "driver_log_bytes": os.path.getsize(os.path.join(outdir, "driver.log")) if os.path.isfile(os.path.join(outdir, "driver.log")) else 0,
}
path = os.path.join(outdir, "heartbeat.json")
tmp = path + ".tmp"
with open(tmp, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(tmp, path)
PY
}

CHILD_PID=""
terminate_child() {
    echo "termination signal received at $(date -u +%Y-%m-%dT%H:%M:%SZ)" >&2
    if [ -n "$CHILD_PID" ] && kill -0 "$CHILD_PID" 2>/dev/null; then
        kill -TERM "$CHILD_PID" 2>/dev/null || true
    fi
    snapshot_progress || true
}
trap terminate_child TERM INT

"$PYTHON_BIN" "$DRIVER" "${ARGS[@]}" >"$DRIVER_LOG" 2>&1 &
CHILD_PID=$!
echo "no-policy continuation started: mode=$E5F_CONT_MODE job=$SLURM_JOB_ID pid=$CHILD_PID outdir=$OUTDIR"
while kill -0 "$CHILD_PID" 2>/dev/null; do
    snapshot_progress
    echo "heartbeat $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    tail -n 8 "$DRIVER_LOG" 2>/dev/null || true
    elapsed=0
    while [ "$elapsed" -lt "$HEARTBEAT_SECONDS" ] && kill -0 "$CHILD_PID" 2>/dev/null; do
        sleep 1
        elapsed=$((elapsed + 1))
    done
done
set +e
wait "$CHILD_PID"
DRIVER_STATUS=$?
set -e
snapshot_progress || true
if [ "$DRIVER_STATUS" -ne 0 ]; then
    echo "no-policy continuation failed with status $DRIVER_STATUS" >&2
    tail -n 100 "$DRIVER_LOG" >&2 || true
    exit "$DRIVER_STATUS"
fi

# Independent packet audit.  A missing closed root is an admissible scientific
# result; a root that fails its declared numerical gates is not.
"$PYTHON_BIN" - \
    "$PROJECT_ROOT" "$OUTDIR" "$POST_PERIODS" "$DRIVER_SHA" \
    "$REPORT_SHA" "$SELECTION_SHA" "$REPEAT1_SUMMARY_SHA" "$REPEAT1_CONTRACT_SHA" \
    "$REPEAT2_SUMMARY_SHA" "$REPEAT2_CONTRACT_SHA" "$SELECTED_TRANSITION_SHA" \
    "$SOURCE_SHA" "$E5F_CONT_EXPECTED_TARGET_SET" "$E5F_CONT_EXPECTED_TARGET_FINGERPRINT" \
    "$E5F_CONT_EXPECTED_CODE_BUNDLE_SHA256" "$E5F_CONT_EXPECTED_RENEWAL_CONTRACT_SHA256" \
    "$E5F_CONT_EXPECTED_SCIENTIFIC_CONTRACT_SHA256" <<'PY'
import csv, hashlib, json, math, os, sys
from pathlib import Path

(root_raw, outdir_raw, periods_raw, driver_sha, report_sha, selection_sha,
 r1s, r1c, r2s, r2c, transition_sha, source_sha, target_set,
 target_fingerprint, code_sha, renewal_sha, scientific_sha) = sys.argv[1:]
root, outdir, periods = Path(root_raw), Path(outdir_raw), int(periods_raw)
sys.path.insert(0, str(root / "code/model/tools"))
import run_e5f_post2023_no_policy_continuations as driver

def read_json(name):
    return json.loads((outdir / name).read_text(encoding="utf-8"))
def read_csv(name):
    with (outdir / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def finite(value, label):
    result = float(value)
    if not math.isfinite(result):
        raise SystemExit(f"nonfinite {label}")
    return result

summary = read_json("summary.json")
manifest = read_json("manifest.json")
latest = read_json("latest_completed_period.json")
input_contract = read_json("input_contract.json")
expected_input_contract = {
    "report_sha256": report_sha,
    "selection_sha256": selection_sha,
    "repeat1_summary_sha256": r1s,
    "repeat1_contract_sha256": r1c,
    "repeat2_summary_sha256": r2s,
    "repeat2_contract_sha256": r2c,
    "selected_transition_sha256": transition_sha,
    "source_sha256": source_sha,
    "target_set": target_set,
    "target_fingerprint": target_fingerprint,
    "code_bundle_sha256": code_sha,
    "renewal_contract_sha256": renewal_sha,
    "scientific_contract_sha256": scientific_sha,
    "driver_sha256": driver_sha,
}
for key, expected in expected_input_contract.items():
    if str(input_contract.get(key)) != str(expected):
        raise SystemExit(f"saved launcher input contract mismatch at {key}")
if summary.get("status") != "complete_no_policy_post2023_continuation_pair":
    raise SystemExit("continuation summary is incomplete")
if summary.get("paper_scope") != "new equilibrium and transition validation; no policy experiment":
    raise SystemExit("packet has the wrong paper scope")
if summary.get("policy_case") != "none" or summary.get("fiscal_change") != "none":
    raise SystemExit("packet contains a policy or fiscal experiment")
if int(summary.get("post_2023_periods", -1)) != periods:
    raise SystemExit("packet has the wrong continuation horizon")
if int(summary.get("terminal_calendar_year", -1)) != 2023 + 4 * periods:
    raise SystemExit("packet has the wrong terminal year")
if summary.get("target_set") != target_set or summary.get("target_fingerprint") != target_fingerprint:
    raise SystemExit("packet target contract mismatch")
expected_provenance = {
    "selected_report_sha256": report_sha,
    "selected_case_transition_sha256": transition_sha,
    "source_sha256": source_sha,
    "target_fingerprint": target_fingerprint,
    "code_bundle_sha256": code_sha,
    "renewal_contract_sha256": renewal_sha,
    "scientific_contract_sha256": scientific_sha,
    "selection_sha256": selection_sha,
}
provenance = summary.get("provenance") or {}
for key, expected in expected_provenance.items():
    if str(provenance.get(key)) != str(expected):
        raise SystemExit(f"packet provenance mismatch at {key}")
if latest.get("status") != "complete":
    raise SystemExit("continuation loop did not write a complete checkpoint")

history = summary.get("history_reproduction_audit") or {}
if history.get("status") != "passed" or finite(history.get("maximum_absolute_gap"), "history gap") > 5e-10:
    raise SystemExit("selected 2007--2023 history did not reproduce")
initial = summary.get("paired_initial_state_audit") or {}
if initial.get("status") != "passed" or finite(initial.get("shared_2007_2023_history_maximum_absolute_gap"), "paired state gap") != 0.0:
    raise SystemExit("paired paths do not share the exact matched state")
outside = summary.get("outside_share_invariance_audit") or {}
if outside.get("status") != "passed" or finite(outside.get("maximum_absolute_identity_residual"), "outside-share identity") > 2e-10:
    raise SystemExit("outside-share accounting gate failed")

rows = read_csv("paired_continuation_path.csv")
scenarios = ("closed_national_benchmark", "open_cbsa_sensitivity")
expected_years = list(range(2007, 2023 + 4 * periods + 1, 4))
if len(rows) != 2 * len(expected_years):
    raise SystemExit("paired path has the wrong row count")
by_scenario = {label: [] for label in scenarios}
for row in rows:
    label = row.get("scenario")
    if label not in by_scenario:
        raise SystemExit(f"unexpected scenario {label}")
    if row.get("policy_case") != "none":
        raise SystemExit("paired path contains a policy row")
    by_scenario[label].append(row)
    if abs(finite(row["relative_market_residual"], "market residual")) > 2e-4:
        raise SystemExit("row-level market gate failed")
    if abs(finite(row["mass_accounting_residual"], "mass residual")) > 2e-8:
        raise SystemExit("row-level mass gate failed")
    if finite(row["feasibility_frontier_projection_mass"], "projection mass") > 1e-6:
        raise SystemExit("row-level feasibility gate failed")
    if finite(row["min_distribution_mass"], "minimum mass") < -1e-14:
        raise SystemExit("row-level minimum-mass gate failed")
    if int(row["nonfinite_distribution_count"]) != 0:
        raise SystemExit("row-level finite-distribution gate failed")
for label in scenarios:
    scenario_rows = sorted(by_scenario[label], key=lambda row: int(row["calendar_year"]))
    if [int(row["calendar_year"]) for row in scenario_rows] != expected_years:
        raise SystemExit(f"{label} has the wrong calendar sequence")
    gates = ((summary.get("paths") or {}).get(label) or {}).get("path_gates") or {}
    if finite(gates.get("maximum_market_residual"), "summary market gate") > 2e-4:
        raise SystemExit(f"{label} summary market gate failed")
    if finite(gates.get("maximum_mass_residual"), "summary mass gate") > 2e-8:
        raise SystemExit(f"{label} summary mass gate failed")
    if finite(gates.get("maximum_feasibility_projection_mass"), "summary projection gate") > 1e-6:
        raise SystemExit(f"{label} summary projection gate failed")
    if finite(gates.get("minimum_distribution_mass"), "summary minimum mass") < -1e-14:
        raise SystemExit(f"{label} summary minimum-mass gate failed")
    if int(gates.get("maximum_nonfinite_distribution_count", -1)) != 0:
        raise SystemExit(f"{label} summary finite gate failed")
closed_2023 = next(row for row in by_scenario[scenarios[0]] if int(row["calendar_year"]) == 2023)
open_2023 = next(row for row in by_scenario[scenarios[1]] if int(row["calendar_year"]) == 2023)
for field in driver.HISTORY_STATE_COLUMNS:
    if float(closed_2023[field]) != float(open_2023[field]):
        raise SystemExit(f"2023 branch state differs at {field}")

closed = summary.get("closed_stationary_endpoint") or {}
closed_status = closed.get("status")
allowed = {
    "no_positive_renewal_root_on_audited_grid",
    "multiple_positive_renewal_roots_on_audited_grid",
    "complete_usable_positive_root",
}
if closed_status not in allowed:
    raise SystemExit(f"closed endpoint has a failed or unknown status: {closed_status}")
usable_closed = bool(closed.get("usable_closed_root"))
if usable_closed != bool(closed.get("between_steady_states_label_allowed")):
    raise SystemExit("closed-root and between-steady-states labels disagree")
if usable_closed:
    if closed_status != "complete_usable_positive_root":
        raise SystemExit("usable closed endpoint has the wrong status")
    if finite(closed.get("renewal_root_absolute_residual"), "closed renewal root") > 2.5e-5:
        raise SystemExit("closed renewal root gate failed")
    if finite(closed.get("housing_clearing_absolute_residual"), "closed housing gap") > 1e-10:
        raise SystemExit("closed endpoint housing gate failed")
elif "not described as a transition between steady states" not in str(summary.get("between_steady_states_statement")):
    raise SystemExit("packet overstates a closed finite-horizon benchmark")
schedule = read_csv("closed_stationary_schedule.csv")
if len(schedule) < 25:
    raise SystemExit("closed stationary schedule has fewer than 25 declared grid points")
for row in schedule:
    if row.get("policy_case") != "none" or row.get("closure") != "closed_national_benchmark":
        raise SystemExit("closed schedule has the wrong scope")
    if finite(row.get("asset_price"), "closed price") <= 0.0:
        raise SystemExit("closed schedule has a nonpositive price")
    finite(row.get("queue_B_over_E"), "closed renewal ratio")

open_endpoint = summary.get("open_stationary_endpoint") or {}
if open_endpoint.get("status") != "complete" or not open_endpoint.get("usable_open_endpoint"):
    raise SystemExit("open endpoint is incomplete or unusable")
if finite(open_endpoint.get("stationary_population_scale"), "open scale") <= 0.0:
    raise SystemExit("open endpoint has a nonpositive scale")
if abs(finite(open_endpoint.get("renewal_identity_residual"), "open renewal identity")) > 2e-10:
    raise SystemExit("open endpoint renewal gate failed")
if abs(finite(open_endpoint.get("fixed_stock_relative_market_gap"), "open market gap")) > 2.5e-5:
    raise SystemExit("open endpoint market gate failed")

if manifest.get("status") != "complete_no_policy_post2023_continuation_manifest":
    raise SystemExit("driver manifest is incomplete")
if manifest.get("driver_sha256") != driver_sha:
    raise SystemExit("driver manifest hash differs from launched driver")
manifest_provenance = manifest.get("input_provenance") or {}
for key, expected in expected_provenance.items():
    if str(manifest_provenance.get(key)) != str(expected):
        raise SystemExit(f"manifest provenance mismatch at {key}")
for name, expected in (manifest.get("artifacts") or {}).items():
    path = outdir / name
    if not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"manifest artifact is missing or empty: {name}")
    if digest(path) != expected:
        raise SystemExit(f"manifest artifact hash mismatch: {name}")
if "NO_POLICY_CONTINUATIONS_COMPLETE" not in (outdir / "driver.log").read_text(encoding="utf-8"):
    raise SystemExit("driver log omits its completion marker")

audit = {
    "status": "passed",
    "artifact_gate_passed": True,
    "post_2023_periods": periods,
    "terminal_calendar_year": 2023 + 4 * periods,
    "paired_path_rows": len(rows),
    "closed_endpoint_status": closed_status,
    "usable_closed_root": usable_closed,
    "usable_open_endpoint": True,
    "driver_manifest_sha256": digest(outdir / "manifest.json"),
}
path = outdir / "launcher_audit.json"
temporary = path.with_suffix(path.suffix + ".tmp")
temporary.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
temporary.replace(path)
PY

"$PYTHON_BIN" - \
    "$OUTDIR/launcher_status.json" "$OUTDIR/heartbeat.json" \
    "$OUTDIR/input_contract.json" "$OUTDIR/launcher_audit.json" \
    "$OUTDIR/manifest.json" "$E5F_CONT_MODE" "$POST_PERIODS" "$SLURM_JOB_ID" <<'PY'
import hashlib, json, os, sys
from datetime import datetime, timezone
status_path, heartbeat_path, contract_path, audit_path, manifest_path, mode, periods, job = sys.argv[1:]
with open(contract_path, encoding="utf-8") as handle:
    contract = json.load(handle)
with open(audit_path, encoding="utf-8") as handle:
    audit = json.load(handle)
manifest_sha = hashlib.sha256(open(manifest_path, "rb").read()).hexdigest()
payload = {
    "status": "complete",
    "completed_utc": datetime.now(timezone.utc).isoformat(),
    "mode": mode,
    "slurm_job_id": job,
    "post_2023_periods": int(periods),
    "input_contract": contract,
    "artifact_gate_passed": bool(audit.get("artifact_gate_passed")),
    "driver_manifest": os.path.realpath(manifest_path),
    "driver_manifest_sha256": manifest_sha,
    "launcher_audit": os.path.realpath(audit_path),
}
for path in (status_path, heartbeat_path):
    temporary = path + ".tmp"
    with open(temporary, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)
PY
STATUS_SHA="$(hash_file "$OUTDIR/launcher_status.json")"
printf '%s  launcher_status.json\n' "$STATUS_SHA" >"$OUTDIR/launcher_status.sha256"
echo "E5F_NO_POLICY_CONTINUATION_COMPLETE mode=$E5F_CONT_MODE periods=$POST_PERIODS receipt_sha256=$STATUS_SHA output=$OUTDIR"
