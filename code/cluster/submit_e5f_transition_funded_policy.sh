#!/usr/bin/env bash
# Exact four-case policy transition from the selected dated calibration.
#
# Preferred submission (from Torch's code/cluster directory):
#   E5F_FUNDED_RUN_TAG=... \
#   E5F_FUNDED_REPORT=output/model/.../summary.json \
#   E5F_FUNDED_EXPECTED_REPORT_SHA256=... \
#   E5F_FUNDED_EXPECTED_SOURCE_SHA256=... \
#   E5F_FUNDED_EXPECTED_TARGET_SET=... \
#   E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT=... \
#   E5F_FUNDED_EXPECTED_MODEL_PROFILE=e5f-income-entry \
#   E5F_FUNDED_OUTSIDE_ORIGIN_SHARE=... \
#   E5F_FUNDED_EXPECTED_CODE_BUNDLE_SHA256=... \
#   E5F_FUNDED_TIME_LIMIT=03:55:00 E5F_FUNDED_MEMORY=32G \
#   bash submit_e5f_transition_funded_policy.sh --submit
#
# The --submit front end passes explicit walltime and memory to sbatch. Direct
# sbatch submission remains possible and uses the conservative defaults below.
#SBATCH --job-name=e5ftrfund
#SBATCH --output=logs/slurm_e5ftrfund_%j.out
#SBATCH --error=logs/slurm_e5ftrfund_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}" && pwd)"
SCRIPT_PATH="$SCRIPT_DIR/$SCRIPT_NAME"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MODEL_DIR="$PROJECT_ROOT/code/model"
mkdir -p "$SCRIPT_DIR/logs"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_FUNDED_RUN_TAG:?E5F_FUNDED_RUN_TAG is required}"
    : "${E5F_FUNDED_REPORT:?E5F_FUNDED_REPORT is required}"
    : "${E5F_FUNDED_EXPECTED_REPORT_SHA256:?expected report SHA-256 is required}"
    : "${E5F_FUNDED_EXPECTED_SOURCE_SHA256:?expected source SHA-256 is required}"
    : "${E5F_FUNDED_EXPECTED_TARGET_SET:?expected target set is required}"
    : "${E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT:?expected target fingerprint is required}"
    : "${E5F_FUNDED_EXPECTED_MODEL_PROFILE:?expected model profile is required}"
    : "${E5F_FUNDED_OUTSIDE_ORIGIN_SHARE:?outside-origin entry share is required}"
    : "${E5F_FUNDED_EXPECTED_CODE_BUNDLE_SHA256:?expected code-bundle SHA-256 is required}"
    TIME_LIMIT="${E5F_FUNDED_TIME_LIMIT:-03:55:00}"
    MEMORY="${E5F_FUNDED_MEMORY:-32G}"
    cd "$SCRIPT_DIR"
    exec sbatch \
        --time="$TIME_LIMIT" \
        --mem="$MEMORY" \
        --export=ALL \
        "$SCRIPT_PATH" "$@"
fi

: "${SLURM_JOB_ID:?Run with --submit or submit this script through sbatch}"
: "${E5F_FUNDED_RUN_TAG:?E5F_FUNDED_RUN_TAG is required}"
: "${E5F_FUNDED_REPORT:?E5F_FUNDED_REPORT is required}"
: "${E5F_FUNDED_EXPECTED_REPORT_SHA256:?expected report SHA-256 is required}"
: "${E5F_FUNDED_EXPECTED_SOURCE_SHA256:?expected source SHA-256 is required}"
: "${E5F_FUNDED_EXPECTED_TARGET_SET:?expected target set is required}"
: "${E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT:?expected target fingerprint is required}"
: "${E5F_FUNDED_EXPECTED_MODEL_PROFILE:?expected model profile is required}"
: "${E5F_FUNDED_OUTSIDE_ORIGIN_SHARE:?outside-origin entry share is required}"
: "${E5F_FUNDED_EXPECTED_CODE_BUNDLE_SHA256:?expected code-bundle SHA-256 is required}"

RUN_TAG="$E5F_FUNDED_RUN_TAG"
if [[ ! "$RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $RUN_TAG" >&2
    exit 2
fi
REPORT="$E5F_FUNDED_REPORT"
if [[ "$REPORT" != /* ]]; then
    REPORT="$PROJECT_ROOT/$REPORT"
fi
SOURCE="${E5F_FUNDED_SOURCE:-$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json}"
if [[ "$SOURCE" != /* ]]; then
    SOURCE="$PROJECT_ROOT/$SOURCE"
fi
OUTDIR="${E5F_FUNDED_OUTDIR:-$PROJECT_ROOT/output/model/e5f_transition_funded_policy_${RUN_TAG}}"
if [[ "$OUTDIR" != /* ]]; then
    OUTDIR="$PROJECT_ROOT/$OUTDIR"
fi
case "$OUTDIR" in
    "$PROJECT_ROOT"/output/model/*) ;;
    *)
        echo "output directory must be a named child of $PROJECT_ROOT/output/model" >&2
        exit 2
        ;;
esac
HEARTBEAT_SECONDS="${E5F_FUNDED_HEARTBEAT_SECONDS:-300}"
if [[ ! "$HEARTBEAT_SECONDS" =~ ^[0-9]+$ ]] || [ "$HEARTBEAT_SECONDS" -lt 30 ]; then
    echo "E5F_FUNDED_HEARTBEAT_SECONDS must be an integer of at least 30" >&2
    exit 2
fi
if [ ! -f "$REPORT" ] || [ ! -f "$SOURCE" ]; then
    echo "missing report or source: report=$REPORT source=$SOURCE" >&2
    exit 2
fi
if [ -e "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    echo "refusing to overwrite nonempty output directory: $OUTDIR" >&2
    exit 2
fi
mkdir -p "$OUTDIR/checkpoints"

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
hash_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}
ACTUAL_REPORT_SHA256="$(hash_file "$REPORT")"
ACTUAL_SOURCE_SHA256="$(hash_file "$SOURCE")"
if [ "$ACTUAL_REPORT_SHA256" != "$E5F_FUNDED_EXPECTED_REPORT_SHA256" ]; then
    echo "report hash mismatch: actual=$ACTUAL_REPORT_SHA256 expected=$E5F_FUNDED_EXPECTED_REPORT_SHA256" >&2
    exit 2
fi
if [ "$ACTUAL_SOURCE_SHA256" != "$E5F_FUNDED_EXPECTED_SOURCE_SHA256" ]; then
    echo "source hash mismatch: actual=$ACTUAL_SOURCE_SHA256 expected=$E5F_FUNDED_EXPECTED_SOURCE_SHA256" >&2
    exit 2
fi

# The outside-origin share is an explicit externally fixed closure object.  A
# production job may not inherit the driver's provisional diagnostic default.
"$PYTHON_BIN" - "$E5F_FUNDED_OUTSIDE_ORIGIN_SHARE" <<'PY'
import math
import sys

try:
    value = float(sys.argv[1])
except ValueError as error:
    raise SystemExit("outside-origin entry share must be numeric") from error
if not math.isfinite(value) or not 0.0 < value < 1.0:
    raise SystemExit("outside-origin entry share must lie strictly between zero and one")
PY

# Match the driver's deterministic bundle hash exactly: ordered relative path,
# NUL, file bytes, NUL for each model-critical source file.
ACTUAL_CODE_BUNDLE_SHA256="$("$PYTHON_BIN" - "$PROJECT_ROOT" <<'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root / "code/model"))
sys.path.insert(0, str(root / "code/model/tools"))
import run_e5f_transition_funded_policy as funded

if funded.ROOT.resolve() != root:
    raise SystemExit(f"funded driver resolved the wrong project root: {funded.ROOT}")
print(funded.code_bundle_sha256(root))
PY
)"
if [ "$ACTUAL_CODE_BUNDLE_SHA256" != "$E5F_FUNDED_EXPECTED_CODE_BUNDLE_SHA256" ]; then
    echo "code-bundle hash mismatch: actual=$ACTUAL_CODE_BUNDLE_SHA256 expected=$E5F_FUNDED_EXPECTED_CODE_BUNDLE_SHA256" >&2
    exit 2
fi

# Gate the report's own contract before spending a single model solve.
"$PYTHON_BIN" - \
    "$REPORT" \
    "$SOURCE" \
    "$ACTUAL_REPORT_SHA256" \
    "$ACTUAL_SOURCE_SHA256" \
    "$E5F_FUNDED_EXPECTED_TARGET_SET" \
    "$E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT" \
    "$E5F_FUNDED_EXPECTED_MODEL_PROFILE" \
    "$E5F_FUNDED_OUTSIDE_ORIGIN_SHARE" \
    "$PROJECT_ROOT" <<'PY'
import sys
from pathlib import Path

(report, source, report_hash, source_hash, target_set, target_fingerprint,
 expected_profile, outside_share, project_root) = sys.argv[1:]
root = Path(project_root)
sys.path.insert(0, str(root / "code/model"))
sys.path.insert(0, str(root / "code/model/tools"))
import run_e5f_open_population_transition as transition
import run_e5f_transition_calibration as calibration
import run_e5f_transition_funded_policy as funded

_, model = transition.configure_sequential_model()
live_code = calibration.code_fingerprint_contract(model)
funded.load_transition_contract(
    Path(report),
    Path(source),
    expected_report_sha256=report_hash,
    expected_source_sha256=source_hash,
    expected_target_set=target_set,
    expected_target_fingerprint=target_fingerprint,
    expected_model_profile=expected_profile,
    outside_origin_entry_share=float(outside_share),
    live_calibration_code_fingerprints=live_code,
)
PY

REQUESTED_TIME="${E5F_FUNDED_TIME_LIMIT:-03:55:00}"
REQUESTED_MEMORY="${E5F_FUNDED_MEMORY:-32G}"
"$PYTHON_BIN" - \
    "$OUTDIR/launch_contract.json" "$RUN_TAG" "$REPORT" "$ACTUAL_REPORT_SHA256" \
    "$SOURCE" "$ACTUAL_SOURCE_SHA256" "$E5F_FUNDED_EXPECTED_TARGET_SET" \
    "$E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT" "$E5F_FUNDED_EXPECTED_MODEL_PROFILE" \
    "$E5F_FUNDED_OUTSIDE_ORIGIN_SHARE" "$ACTUAL_CODE_BUNDLE_SHA256" \
    "$REQUESTED_TIME" "$REQUESTED_MEMORY" "$SLURM_JOB_ID" "$PROJECT_ROOT" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone

(path, run_tag, report, report_hash, source, source_hash, target_set,
 target_fingerprint, profile, outside_share, code_bundle_hash, requested_time,
 requested_memory, job_id, project_root) = sys.argv[1:]
root = os.path.realpath(project_root)
sys.path.insert(0, os.path.join(root, "code/model"))
sys.path.insert(0, os.path.join(root, "code/model/tools"))
import run_e5f_transition_funded_policy as funded

payload = {
    "status": "launched",
    "launched_utc": datetime.now(timezone.utc).isoformat(),
    "run_tag": run_tag,
    "slurm_job_id": job_id,
    "report": report,
    "report_sha256": report_hash,
    "source": source,
    "source_sha256": source_hash,
    "target_set": target_set,
    "target_fingerprint": target_fingerprint,
    "model_profile": profile,
    "outside_origin_entry_share": float(outside_share),
    "code_bundle_sha256": code_bundle_hash,
    "code_bundle_files": list(funded.CODE_BUNDLE_FILES),
    "policy_cases": [
        "unrebated_tax1_status_quo",
        "rebated_tax1_baseline",
        "rebated_tax2",
        "rebated_tax2_grant0p4_Hge6",
    ],
    "post_2023_periods": 40,
    "total_calendar_points": 180,
    "run_size": "4 cases x 45 dates; 20 shared dates through 2023; 40 active dates for each of 3 funded cases",
    "rough_walltime": "floor smoke extrapolation about 40 minutes; bounded by the 3:55 short-partition limit for the 15-income-state profile",
    "requested_time": requested_time,
    "requested_memory": requested_memory,
    "slurm_mem_per_node": os.environ.get("SLURM_MEM_PER_NODE"),
    "stop_criteria": "all 4 cases complete 45 dates and pass timing, market, fiscal, mass, finite, profile, outside-share, and fingerprint gates",
}
tmp = path + ".tmp"
with open(tmp, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(tmp, path)
PY

export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
DRIVER_LOG="$OUTDIR/driver.log"
ARGS=(
    --source "$SOURCE"
    --transition-report "$REPORT"
    --expected-report-sha256 "$ACTUAL_REPORT_SHA256"
    --expected-source-sha256 "$ACTUAL_SOURCE_SHA256"
    --expected-target-set "$E5F_FUNDED_EXPECTED_TARGET_SET"
    --expected-target-fingerprint "$E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT"
    --expected-model-profile "$E5F_FUNDED_EXPECTED_MODEL_PROFILE"
    --expected-code-bundle-sha256 "$ACTUAL_CODE_BUNDLE_SHA256"
    --outside-origin-entry-share "$E5F_FUNDED_OUTSIDE_ORIGIN_SHARE"
    --output-dir "$OUTDIR"
    --post-2023-periods 40
)

snapshot_progress() {
    "$PYTHON_BIN" - "$OUTDIR" "$SLURM_JOB_ID" <<'PY'
import glob
import json
import os
import shutil
import sys
from datetime import datetime, timezone

outdir, job_id = sys.argv[1:]
checkpoints = os.path.join(outdir, "checkpoints")
os.makedirs(checkpoints, exist_ok=True)
latest = []
for path in sorted(glob.glob(os.path.join(outdir, "cases", "*", "latest_completed_period.json"))):
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
        case = str(payload.get("scenario") or os.path.basename(os.path.dirname(path)))
        period = int(payload["completed_period"])
        destination = os.path.join(checkpoints, f"{case}_period_{period:03d}.json")
        if not os.path.exists(destination):
            shutil.copy2(path, destination)
        latest.append({"case": case, "period": period, "path": path})
    except Exception as error:
        latest.append({"path": path, "snapshot_error": repr(error)})
case_path = os.path.join(outdir, "latest_completed_case.json")
case_payload = None
if os.path.isfile(case_path):
    try:
        with open(case_path, encoding="utf-8") as handle:
            case_payload = json.load(handle)
        shutil.copy2(case_path, os.path.join(checkpoints, "latest_completed_case.json"))
    except Exception as error:
        case_payload = {"snapshot_error": repr(error)}
heartbeat = {
    "status": "running",
    "utc": datetime.now(timezone.utc).isoformat(),
    "slurm_job_id": job_id,
    "latest_periods": latest,
    "latest_case": case_payload,
}
path = os.path.join(outdir, "heartbeat.json")
tmp = path + ".tmp"
with open(tmp, "w", encoding="utf-8") as handle:
    json.dump(heartbeat, handle, indent=2, sort_keys=True)
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

"$PYTHON_BIN" "$MODEL_DIR/tools/run_e5f_transition_funded_policy.py" "${ARGS[@]}" \
    >"$DRIVER_LOG" 2>&1 &
CHILD_PID=$!
echo "funded transition started: job=$SLURM_JOB_ID pid=$CHILD_PID outdir=$OUTDIR"

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
    echo "funded transition failed with status $DRIVER_STATUS" >&2
    tail -n 80 "$DRIVER_LOG" >&2 || true
    exit "$DRIVER_STATUS"
fi

# Final production gates: exactly four complete 45-date paths and pinned metadata.
"$PYTHON_BIN" - \
    "$OUTDIR" "$ACTUAL_REPORT_SHA256" "$ACTUAL_SOURCE_SHA256" \
    "$E5F_FUNDED_EXPECTED_TARGET_SET" "$E5F_FUNDED_EXPECTED_TARGET_FINGERPRINT" \
    "$E5F_FUNDED_EXPECTED_MODEL_PROFILE" "$E5F_FUNDED_OUTSIDE_ORIGIN_SHARE" \
    "$ACTUAL_CODE_BUNDLE_SHA256" "$PROJECT_ROOT" <<'PY'
import csv
import json
import math
import os
import sys

(outdir, report_hash, source_hash, target_set, target_fingerprint, profile,
 outside_share, code_bundle_hash, project_root) = sys.argv[1:]
sys.path.insert(0, os.path.join(project_root, "code/model"))
sys.path.insert(0, os.path.join(project_root, "code/model/tools"))
import run_e5f_transition_funded_policy as funded

with open(os.path.join(outdir, "metadata.json"), encoding="utf-8") as handle:
    metadata = json.load(handle)
if metadata.get("status") != "complete" or int(metadata.get("periods", -1)) != 45:
    raise SystemExit("driver metadata is not a complete 45-date packet")
contract = metadata.get("contract") or {}
if contract.get("report_sha256") != report_hash or contract.get("source_sha256") != source_hash:
    raise SystemExit("completed packet fingerprint mismatch")
if contract.get("target_set") != target_set or contract.get("target_fingerprint") != target_fingerprint:
    raise SystemExit("completed packet target-contract mismatch")
if (contract.get("model_profile") or {}).get("name") != profile:
    raise SystemExit("completed packet model-profile mismatch")
if (metadata.get("active_model_profile") or {}).get("name") != profile:
    raise SystemExit("solved model-profile mismatch")
solved_outside_share = float(metadata.get("outside_origin_entry_share", float("nan")))
if not math.isfinite(solved_outside_share) or abs(solved_outside_share - float(outside_share)) > 1e-15:
    raise SystemExit("completed packet outside-origin share mismatch")
if metadata.get("code_bundle_sha256") != code_bundle_hash:
    raise SystemExit("completed packet code-bundle fingerprint mismatch")
if metadata.get("code_bundle_files") != list(funded.CODE_BUNDLE_FILES):
    raise SystemExit("completed packet code-bundle file contract mismatch")
provenance = contract.get("provenance") or {}
reported_outside = provenance.get("outside_origin_entry_share_reported")
if profile == "e5f-income-entry" and reported_outside is None:
    raise SystemExit("repaired packet used a legacy outside-share exception")
if reported_outside is not None and abs(float(reported_outside) - solved_outside_share) > 1e-15:
    raise SystemExit("report and solved outside-origin shares disagree")
reported_calibration_code = provenance.get("calibration_code_fingerprints_reported") or {}
live_calibration_code = provenance.get("calibration_code_fingerprints_live") or {}
if profile == "e5f-income-entry" and not reported_calibration_code:
    raise SystemExit("repaired packet used a legacy calibration-code exception")
if reported_calibration_code and reported_calibration_code.get("bundle_sha256") != live_calibration_code.get("bundle_sha256"):
    raise SystemExit("report and live calibration code bundles disagree")
if metadata.get("relative_effect_baseline") != "unrebated_tax1_status_quo":
    raise SystemExit("completed packet uses the wrong relative-effect baseline")
if float(metadata.get("prepolicy_path_identity_max_abs_gap", 1.0)) > 1e-12:
    raise SystemExit("completed policy paths are not identical through 2023")
timing = metadata.get("matched_2023_and_policy_timing_gates") or {}
if max(
    float(timing.get("matched_2023_population_gap", 1.0)),
    float(timing.get("matched_2023_housing_cost_gap", 1.0)),
) > 2e-8:
    raise SystemExit("completed packet does not preserve the matched 2023 economy")
endpoint = metadata.get("stationary_open_endpoint") or {}
if endpoint.get("status") != "complete":
    raise SystemExit("completed packet omits the exact status-quo stationary endpoint")
endpoint_values = {
    name: float(endpoint.get(name, float("nan")))
    for name in (
        "asset_price",
        "stationary_population_scale",
        "renewal_denominator",
        "fixed_stock_relative_market_gap",
        "psi_child",
    )
}
if not all(math.isfinite(value) for value in endpoint_values.values()):
    raise SystemExit("stationary endpoint contains a nonfinite or missing value")
if endpoint_values["asset_price"] <= 0.0 or endpoint_values["stationary_population_scale"] <= 0.0:
    raise SystemExit("stationary endpoint has a nonpositive price or population")
if endpoint_values["renewal_denominator"] <= 0.0:
    raise SystemExit("stationary endpoint has a nonpositive renewal denominator")
if abs(endpoint_values["fixed_stock_relative_market_gap"]) > 2.5e-5:
    raise SystemExit("stationary endpoint housing market does not clear")
if abs(endpoint_values["psi_child"] - float(contract["new_psi_child"])) > 1e-12:
    raise SystemExit("stationary endpoint does not use the selected terminal preference")
method = endpoint.get("market_clearing_method")
if method not in {
    "pure-price root",
    "persistent entry-type convexification at a discrete household-policy jump",
}:
    raise SystemExit("stationary endpoint omits its market-clearing semantics")
if method.startswith("persistent entry-type"):
    weight = float(
        endpoint.get("persistent_type_population_share_positive_branch", float("nan"))
    )
    width = float(endpoint.get("indifference_price_bracket_width", float("nan")))
    if not math.isfinite(weight) or not 0.0 < weight < 1.0:
        raise SystemExit("persistent endpoint has an invalid population-type share")
    if not math.isfinite(width) or not 0.0 < width <= 1e-6 * max(1.0, endpoint_values["asset_price"]):
        raise SystemExit("persistent endpoint has an invalid indifference-price bracket")
    d_pos = float(endpoint.get("positive_branch_renewal_denominator", float("nan")))
    d_neg = float(endpoint.get("negative_branch_renewal_denominator", float("nan")))
    outside_pos = float(endpoint.get("positive_branch_outside_entry_flow", float("nan")))
    outside_neg = float(endpoint.get("negative_branch_outside_entry_flow", float("nan")))
    share_pos = float(
        endpoint.get("implied_outside_entrant_share_positive_branch", float("nan"))
    )
    share_neg = float(
        endpoint.get("implied_outside_entrant_share_negative_branch", float("nan"))
    )
    identity = float(endpoint.get("outside_entry_flow_identity_residual", float("nan")))
    positive_identity = float(
        endpoint.get("positive_branch_renewal_identity_residual", float("nan"))
    )
    negative_identity = float(
        endpoint.get("negative_branch_renewal_identity_residual", float("nan"))
    )
    persistent_values = (
        d_pos, d_neg, outside_pos, outside_neg, share_pos, share_neg,
        identity, positive_identity, negative_identity,
    )
    if not all(math.isfinite(value) for value in persistent_values):
        raise SystemExit("persistent endpoint accounting contains a nonfinite value")
    if min(d_pos, d_neg, outside_pos, outside_neg, share_pos, share_neg) <= 0.0:
        raise SystemExit("persistent endpoint has a nonpositive branch flow or denominator")
    if max(share_pos, share_neg) >= 1.0 or abs(share_pos + share_neg - 1.0) > 1e-12:
        raise SystemExit("persistent endpoint outside-entrant shares are invalid")
    if abs(outside_pos + outside_neg - float((metadata.get("renewal") or {})["outside_flow"])) > 1e-12:
        raise SystemExit("persistent endpoint branch flows do not reproduce outside inflow")
    if max(abs(identity), abs(positive_identity), abs(negative_identity)) > 1e-12:
        raise SystemExit("persistent endpoint renewal identity failed")
pure_endpoint_path = os.path.join(
    outdir,
    "stationary_endpoint_status_quo",
    "stationary_pure_branch_endpoint.csv",
)
search_path = os.path.join(
    outdir,
    "stationary_endpoint_status_quo",
    "stationary_endpoint_search.csv",
)
if not os.path.isfile(pure_endpoint_path) or os.path.getsize(pure_endpoint_path) <= 0:
    raise SystemExit("stationary pure-branch endpoint diagnostic is missing")
if not os.path.isfile(search_path) or os.path.getsize(search_path) <= 0:
    raise SystemExit("stationary endpoint search diagnostic is missing")
diagnostics = metadata.get("diagnostic_figure") or {}
for suffix in ("png", "pdf"):
    expected_figure = os.path.join(outdir, f"funded_policy_transition_paths.{suffix}")
    if os.path.realpath(str(diagnostics.get(suffix, ""))) != os.path.realpath(expected_figure):
        raise SystemExit(f"diagnostic {suffix} path does not match the stable packet contract")
    if not os.path.isfile(expected_figure) or os.path.getsize(expected_figure) <= 0:
        raise SystemExit(f"diagnostic {suffix} is missing or empty")
if profile == "e5f-income-entry":
    if int(metadata.get("income_state_count", -1)) != 15:
        raise SystemExit("income-entry packet did not solve the measured 15-state process")
    estimated_cost = float((contract.get("theta") or {}).get("first_birth_fixed_cost", float("nan")))
    solved_cost = float(metadata.get("first_birth_fixed_cost", float("nan")))
    if not (math.isfinite(estimated_cost) and math.isfinite(solved_cost)) or abs(solved_cost - estimated_cost) > 1e-12:
        raise SystemExit("income-entry packet did not preserve the estimated first-birth cost")
    declared_profile = (contract.get("model_profile") or {}).get("report_metadata") or {}
    live_profile = metadata.get("active_model_profile") or {}
    if int(declared_profile.get("income_state_count", -1)) != int(live_profile.get("income_state_count", -2)):
        raise SystemExit("report/live income-state-count semantics mismatch")
    declared_variance = float(
        declared_profile.get("permanent_income_log_variance", float("nan"))
    )
    live_variance = float(
        live_profile.get("permanent_income_log_variance", float("nan"))
    )
    if not (math.isfinite(declared_variance) and math.isfinite(live_variance)):
        raise SystemExit("report/live permanent-income variance is missing or nonfinite")
    if abs(declared_variance - live_variance) > 1e-15:
        raise SystemExit("report/live permanent-income variance mismatch")
    for semantic_field in (
        "first_birth_fixed_cost",
        "first_birth_fixed_cost_semantics",
    ):
        if declared_profile.get(semantic_field) != live_profile.get(semantic_field):
            raise SystemExit(f"report/live {semantic_field} semantics mismatch")
expected = [
    "unrebated_tax1_status_quo",
    "rebated_tax1_baseline",
    "rebated_tax2",
    "rebated_tax2_grant0p4_Hge6",
]
with open(os.path.join(outdir, "summary.csv"), newline="", encoding="utf-8") as handle:
    summary = list(csv.DictReader(handle))
if [row.get("case") for row in summary] != expected:
    raise SystemExit("summary does not contain the exact ordered four-case contract")
for case in expected:
    path = os.path.join(outdir, "cases", case, "transition_path.csv")
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 45 or [int(row["period"]) for row in rows] != list(range(45)):
        raise SystemExit(f"incomplete transition path for {case}")
    if int(rows[4]["year"]) != 2023 or rows[4]["policy_active"].lower() != "false":
        raise SystemExit(f"{case} alters or mislabels the matched 2023 date")
    expected_active = case != "unrebated_tax1_status_quo"
    if int(rows[5]["year"]) != 2027 or (rows[5]["policy_active"].lower() == "true") != expected_active:
        raise SystemExit(f"{case} has the wrong first-reform-date timing")
gates = metadata.get("numerical_gates") or {}
gate_values = {
    "maximum_market_residual": float(gates.get("maximum_market_residual", float("nan"))),
    "market_tolerance": float(gates.get("market_tolerance", float("nan"))),
    "maximum_fiscal_residual": float(gates.get("maximum_fiscal_residual", float("nan"))),
    "fiscal_absolute_tolerance": float(gates.get("fiscal_absolute_tolerance", float("nan"))),
    "maximum_mass_residual": float(gates.get("maximum_mass_residual", float("nan"))),
    "mass_tolerance": float(gates.get("mass_tolerance", float("nan"))),
}
if not all(math.isfinite(value) for value in gate_values.values()):
    raise SystemExit("final numerical gates contain a nonfinite or missing value")
if gate_values["maximum_market_residual"] > gate_values["market_tolerance"]:
    raise SystemExit("final market gate failed")
if gate_values["maximum_fiscal_residual"] > gate_values["fiscal_absolute_tolerance"]:
    raise SystemExit("final fiscal gate failed")
if gate_values["maximum_mass_residual"] > gate_values["mass_tolerance"]:
    raise SystemExit("final mass gate failed")
if int(gates.get("maximum_nonfinite_reported_objects", -1)) != 0:
    raise SystemExit("final reported-object nonfinite gate failed")
if int(gates.get("maximum_nonfinite_distribution_count", -1)) != 0:
    raise SystemExit("final distribution nonfinite gate failed")
minimum_mass = float(gates.get("minimum_distribution_mass", float("nan")))
maximum_projection = float(
    gates.get("maximum_feasibility_frontier_projection_mass", float("nan"))
)
if not math.isfinite(minimum_mass) or minimum_mass < -1e-14:
    raise SystemExit("final minimum-distribution-mass gate failed")
if not math.isfinite(maximum_projection) or maximum_projection > 1e-6:
    raise SystemExit("final feasibility-projection gate failed")
PY

"$PYTHON_BIN" - "$OUTDIR/launcher_status.json" "$OUTDIR/heartbeat.json" "$SLURM_JOB_ID" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone

path, heartbeat_path, job_id = sys.argv[1:]
payload = {
    "status": "complete",
    "completed_utc": datetime.now(timezone.utc).isoformat(),
    "slurm_job_id": job_id,
}
tmp = path + ".tmp"
with open(tmp, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(tmp, path)
heartbeat = dict(payload)
heartbeat["final_packet"] = os.path.dirname(path)
tmp = heartbeat_path + ".tmp"
with open(tmp, "w", encoding="utf-8") as handle:
    json.dump(heartbeat, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(tmp, heartbeat_path)
PY
echo "funded transition complete: $OUTDIR"
