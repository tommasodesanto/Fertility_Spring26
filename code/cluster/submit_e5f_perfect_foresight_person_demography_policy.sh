#!/usr/bin/env bash
#SBATCH --job-name=e5f_pf_person_policy
#SBATCH --output=logs/slurm_e5f_pf_person_policy_%j.out
#SBATCH --error=logs/slurm_e5f_pf_person_policy_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_perfect_foresight_person_demography_policy.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_PF_PERSON_POLICY_RUN_TAG:?run tag is required}"
    : "${E5F_PF_PERSON_POLICY_MODE:?mode is required}"
    case "$E5F_PF_PERSON_POLICY_MODE" in
        smoke) REQUESTED_TIME="${E5F_PF_PERSON_POLICY_TIME_LIMIT:-00:45:00}" ;;
        convergence) REQUESTED_TIME="${E5F_PF_PERSON_POLICY_TIME_LIMIT:-04:00:00}" ;;
        *) echo "mode must be smoke or convergence" >&2; exit 2 ;;
    esac
    SBATCH_OVERRIDES=()
    if [ -n "${E5F_PF_PERSON_POLICY_PARTITION:-}" ]; then
        [[ "$E5F_PF_PERSON_POLICY_PARTITION" =~ ^[A-Za-z0-9._-]+$ ]] || {
            echo "invalid partition: $E5F_PF_PERSON_POLICY_PARTITION" >&2
            exit 2
        }
        SBATCH_OVERRIDES+=(--partition="$E5F_PF_PERSON_POLICY_PARTITION")
    fi
    if [ -n "${E5F_PF_PERSON_POLICY_QOS:-}" ]; then
        [[ "$E5F_PF_PERSON_POLICY_QOS" =~ ^[A-Za-z0-9._-]+$ ]] || {
            echo "invalid QOS: $E5F_PF_PERSON_POLICY_QOS" >&2
            exit 2
        }
        SBATCH_OVERRIDES+=(--qos="$E5F_PF_PERSON_POLICY_QOS")
    fi
    cd "$SCRIPT_DIR"
    exec sbatch --time="$REQUESTED_TIME" \
        --mem="${E5F_PF_PERSON_POLICY_MEMORY:-12G}" --cpus-per-task=1 \
        "${SBATCH_OVERRIDES[@]}" --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_PF_PERSON_POLICY_RUN_TAG:?run tag is required}"
: "${E5F_PF_PERSON_POLICY_MODE:?mode is required}"
: "${E5F_PF_PERSON_POLICY_CASE:?case is required}"
: "${E5F_PF_PERSON_POLICY_TERMINAL_SUMMARY:?terminal summary is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_TERMINAL_SUMMARY_SHA256:?terminal-summary hash is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_DRIVER_SHA256:?path-driver hash is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_PERSON_DRIVER_SHA256:?person-demography hash is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_FUNDED_DRIVER_SHA256:?funded-policy hash is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_PF_DRIVER_SHA256:?PF-driver hash is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_SOLVER_SHA256:?solver hash is required}"
: "${E5F_PF_PERSON_POLICY_EXPECTED_CONTINUATION_SHA256:?continuation hash is required}"

case "$E5F_PF_PERSON_POLICY_CASE" in
    rebated-tax1-baseline|rebated-tax2-reform) ;;
    *) echo "invalid case: $E5F_PF_PERSON_POLICY_CASE" >&2; exit 2 ;;
esac
case "$E5F_PF_PERSON_POLICY_MODE" in
    smoke|convergence) ;;
    *) echo "invalid mode: $E5F_PF_PERSON_POLICY_MODE" >&2; exit 2 ;;
esac
if [[ ! "$E5F_PF_PERSON_POLICY_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $E5F_PF_PERSON_POLICY_RUN_TAG" >&2
    exit 2
fi
SAME_HORIZON_POLISH="${E5F_PF_PERSON_POLICY_SAME_HORIZON_POLISH:-0}"
case "$SAME_HORIZON_POLISH" in
    0|1) ;;
    *) echo "same-horizon polish flag must be 0 or 1" >&2; exit 2 ;;
esac

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography_policy.py"
PERSON_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography.py"
FUNDED_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_rebated_property_tax.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
SOURCE_DIR="$PROJECT_ROOT/output/model/e5f_persons_demographic_satellite_20260826b/source_data"
HEADSHIP_DIR="$PROJECT_ROOT/output/model/e5f_headship_demographic_bridge_20260826a"
TERMINAL_SUMMARY="$E5F_PF_PERSON_POLICY_TERMINAL_SUMMARY"
[[ "$TERMINAL_SUMMARY" = /* ]] || TERMINAL_SUMMARY="$PROJECT_ROOT/$TERMINAL_SUMMARY"
OUTDIR="$PROJECT_ROOT/output/model/e5f_pf_person_policy_${E5F_PF_PERSON_POLICY_CASE}_${E5F_PF_PERSON_POLICY_RUN_TAG}_${E5F_PF_PERSON_POLICY_MODE}"
mkdir -p "$CLUSTER_DIR/logs"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    if [ "$(hash_file "$path")" != "$expected" ]; then
        echo "$label hash mismatch" >&2
        exit 2
    fi
}
check_hash "$DRIVER" "$E5F_PF_PERSON_POLICY_EXPECTED_DRIVER_SHA256" "path driver"
check_hash "$PERSON_DRIVER" "$E5F_PF_PERSON_POLICY_EXPECTED_PERSON_DRIVER_SHA256" "person-demography driver"
check_hash "$FUNDED_DRIVER" "$E5F_PF_PERSON_POLICY_EXPECTED_FUNDED_DRIVER_SHA256" "funded-policy driver"
check_hash "$PF_DRIVER" "$E5F_PF_PERSON_POLICY_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$SOLVER" "$E5F_PF_PERSON_POLICY_EXPECTED_SOLVER_SHA256" "solver"
check_hash "$CONTINUATION" "$E5F_PF_PERSON_POLICY_EXPECTED_CONTINUATION_SHA256" "continuation"
check_hash "$TERMINAL_SUMMARY" "$E5F_PF_PERSON_POLICY_EXPECTED_TERMINAL_SUMMARY_SHA256" "terminal summary"
for required in \
    "$SOURCE_DIR/population_mid.csv" \
    "$SOURCE_DIR/births_mid.csv" \
    "$SOURCE_DIR/survival.csv" \
    "$SOURCE_DIR/vintage_2025_age_sex.csv" \
    "$HEADSHIP_DIR/acs_headship_profiles.csv"; do
    [ -s "$required" ] || { echo "missing demographic input: $required" >&2; exit 2; }
done
if [ -e "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    echo "refusing to overwrite nonempty output: $OUTDIR" >&2
    exit 2
fi

if [ "$E5F_PF_PERSON_POLICY_MODE" = "smoke" ]; then
    HORIZON="${E5F_PF_PERSON_POLICY_HORIZON:-2}"
    MAX_ITERATIONS="${E5F_PF_PERSON_POLICY_MAX_ITERATIONS:-1}"
else
    : "${E5F_PF_PERSON_POLICY_HORIZON:?convergence horizon is required}"
    : "${E5F_PF_PERSON_POLICY_SMOKE_SUMMARY:?convergence requires smoke summary}"
    HORIZON="$E5F_PF_PERSON_POLICY_HORIZON"
    MAX_ITERATIONS="${E5F_PF_PERSON_POLICY_MAX_ITERATIONS:-35}"
    SMOKE_SUMMARY="$E5F_PF_PERSON_POLICY_SMOKE_SUMMARY"
    [[ "$SMOKE_SUMMARY" = /* ]] || SMOKE_SUMMARY="$PROJECT_ROOT/$SMOKE_SUMMARY"
    python3 - "$SMOKE_SUMMARY" "$E5F_PF_PERSON_POLICY_EXPECTED_DRIVER_SHA256" \
        "$E5F_PF_PERSON_POLICY_CASE" <<'PY'
import json, sys
from pathlib import Path
summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
    raise SystemExit("person-demography policy smoke is incomplete")
if summary.get("case") != sys.argv[3]:
    raise SystemExit("smoke used a different policy case")
if (summary.get("source_hashes") or {}).get("driver") != sys.argv[2]:
    raise SystemExit("smoke used a different path-driver hash")
if not all((summary.get("accounting_gates") or {}).values()):
    raise SystemExit("smoke accounting gates failed")
if summary.get("path_status") not in {"converged", "maximum_iterations_reached"}:
    raise SystemExit("smoke did not execute the path loop")
PY
    if [ -n "${E5F_PF_PERSON_POLICY_INITIAL_PATH:-}" ]; then
        : "${E5F_PF_PERSON_POLICY_INITIAL_SUMMARY:?seeded convergence requires initial summary}"
        : "${E5F_PF_PERSON_POLICY_EXPECTED_INITIAL_SUMMARY_SHA256:?seeded convergence requires initial-summary hash}"
        : "${E5F_PF_PERSON_POLICY_EXPECTED_INITIAL_PATH_SHA256:?seeded convergence requires initial-path hash}"
        INITIAL_PATH="$E5F_PF_PERSON_POLICY_INITIAL_PATH"
        [[ "$INITIAL_PATH" = /* ]] || INITIAL_PATH="$PROJECT_ROOT/$INITIAL_PATH"
        INITIAL_SUMMARY="$E5F_PF_PERSON_POLICY_INITIAL_SUMMARY"
        [[ "$INITIAL_SUMMARY" = /* ]] || INITIAL_SUMMARY="$PROJECT_ROOT/$INITIAL_SUMMARY"
        [ -s "$INITIAL_PATH" ] || { echo "missing initial path: $INITIAL_PATH" >&2; exit 2; }
        [ -s "$INITIAL_SUMMARY" ] || { echo "missing initial summary: $INITIAL_SUMMARY" >&2; exit 2; }
        check_hash "$INITIAL_PATH" "$E5F_PF_PERSON_POLICY_EXPECTED_INITIAL_PATH_SHA256" "initial path"
        check_hash "$INITIAL_SUMMARY" "$E5F_PF_PERSON_POLICY_EXPECTED_INITIAL_SUMMARY_SHA256" "initial summary"
        python3 - "$INITIAL_SUMMARY" "$INITIAL_PATH" "$TERMINAL_SUMMARY" \
            "$E5F_PF_PERSON_POLICY_CASE" "$HORIZON" "$SAME_HORIZON_POLISH" <<'PY'
import csv, json, math, os, sys
from pathlib import Path

summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
rows = list(csv.DictReader(Path(sys.argv[2]).open(newline="", encoding="utf-8")))
terminal = json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
case = sys.argv[4]
requested_horizon = int(sys.argv[5])
same_horizon_polish = sys.argv[6] == "1"
if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
    raise SystemExit("initial summary is incomplete")
if summary.get("case") != case:
    raise SystemExit("initial summary uses a different policy case")
if not all((summary.get("accounting_gates") or {}).values()):
    raise SystemExit("initial path accounting gates failed")
seed_horizon = int(summary.get("horizon", 0))
if same_horizon_polish:
    if seed_horizon != requested_horizon:
        raise SystemExit("same-horizon polish requires an equal-horizon seed")
    if summary.get("path_converged") or summary.get("path_status") != "maximum_iterations_reached":
        raise SystemExit("same-horizon polish requires a bounded nonconverged seed")
else:
    if not summary.get("path_converged") or summary.get("path_status") != "converged":
        raise SystemExit("initial path has not passed its equilibrium gates")
    if seed_horizon < 2 or seed_horizon >= requested_horizon:
        raise SystemExit("initial-path horizon must be shorter than requested horizon")
if len(rows) != seed_horizon:
    raise SystemExit("initial-path row count differs from its summary horizon")
if [int(row["period"]) for row in rows] != list(range(seed_horizon)):
    raise SystemExit("initial-path periods are not contiguous from zero")
try:
    market = max(abs(float(row["relative_market_residual"])) for row in rows)
    fiscal = max(abs(float(row["government_budget_residual"])) for row in rows)
except (KeyError, TypeError, ValueError) as exc:
    raise SystemExit("initial path lacks finite equilibrium residuals") from exc
if not math.isfinite(market) or not math.isfinite(fiscal):
    raise SystemExit("initial path has nonfinite equilibrium residuals")
if not (
    math.isclose(market, float(summary["maximum_market_residual"]), rel_tol=1e-12, abs_tol=1e-14)
    and math.isclose(fiscal, float(summary["maximum_fiscal_residual"]), rel_tol=1e-12, abs_tol=1e-14)
):
    raise SystemExit("initial summary residuals differ from the saved path")
recomputed_converged = market <= 2.0e-4 and fiscal <= 2.5e-5
if bool(summary.get("path_converged")) != recomputed_converged:
    raise SystemExit("initial summary convergence flag differs from the saved path")
if summary.get("terminal_root") != terminal:
    raise SystemExit("initial path used a different terminal root")
expected = {
    "driver": os.environ["E5F_PF_PERSON_POLICY_EXPECTED_DRIVER_SHA256"],
    "person_demography": os.environ[
        "E5F_PF_PERSON_POLICY_EXPECTED_PERSON_DRIVER_SHA256"
    ],
    "funded_policy": os.environ[
        "E5F_PF_PERSON_POLICY_EXPECTED_FUNDED_DRIVER_SHA256"
    ],
    "perfect_foresight": os.environ[
        "E5F_PF_PERSON_POLICY_EXPECTED_PF_DRIVER_SHA256"
    ],
    "solver": os.environ["E5F_PF_PERSON_POLICY_EXPECTED_SOLVER_SHA256"],
}
for name, digest in expected.items():
    if (summary.get("source_hashes") or {}).get(name) != digest:
        raise SystemExit(f"initial path used a different {name} hash")
PY
    fi
fi

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
mkdir -p "$OUTDIR"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

"$PYTHON_BIN" - "$OUTDIR/launch_contract.json" <<'PY'
import json, os, sys
path = sys.argv[1]
payload = {
    "status": "launched",
    "run_tag": os.environ["E5F_PF_PERSON_POLICY_RUN_TAG"],
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "mode": os.environ["E5F_PF_PERSON_POLICY_MODE"],
    "case": os.environ["E5F_PF_PERSON_POLICY_CASE"],
    "horizon": int(os.environ.get("E5F_PF_PERSON_POLICY_HORIZON", "2")),
    "maximum_path_iterations": int(os.environ.get("E5F_PF_PERSON_POLICY_MAX_ITERATIONS", "1")),
    "scheduler": {
        "partition": os.environ.get("SLURM_JOB_PARTITION"),
        "qos": os.environ.get("SLURM_JOB_QOS"),
        "requested_time_limit": os.environ.get("E5F_PF_PERSON_POLICY_TIME_LIMIT"),
        "requested_memory": os.environ.get("E5F_PF_PERSON_POLICY_MEMORY"),
    },
    "terminal_summary": os.environ["E5F_PF_PERSON_POLICY_TERMINAL_SUMMARY"],
    "initial_path": os.environ.get("E5F_PF_PERSON_POLICY_INITIAL_PATH"),
    "initial_summary": os.environ.get("E5F_PF_PERSON_POLICY_INITIAL_SUMMARY"),
    "same_horizon_polish": os.environ.get(
        "E5F_PF_PERSON_POLICY_SAME_HORIZON_POLISH", "0"
    ) == "1",
    "hashes": {
        key: value for key, value in os.environ.items()
        if key.startswith("E5F_PF_PERSON_POLICY_EXPECTED_")
    },
}
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY

ARGS=(
    --case "$E5F_PF_PERSON_POLICY_CASE"
    --terminal-summary "$TERMINAL_SUMMARY"
    --horizon "$HORIZON"
    --source-dir "$SOURCE_DIR"
    --headship-dir "$HEADSHIP_DIR"
    --output-dir "$OUTDIR"
    --maximum-path-iterations "$MAX_ITERATIONS"
    --price-damping "${E5F_PF_PERSON_POLICY_PRICE_DAMPING:-0.25}"
    --transfer-damping "${E5F_PF_PERSON_POLICY_TRANSFER_DAMPING:-0.50}"
    --maximum-log-price-step "${E5F_PF_PERSON_POLICY_MAXIMUM_LOG_PRICE_STEP:-0.12}"
    --maximum-transfer-step "${E5F_PF_PERSON_POLICY_MAXIMUM_TRANSFER_STEP:-0.08}"
)
if [ -n "${E5F_PF_PERSON_POLICY_INITIAL_PATH:-}" ]; then
    INITIAL_PATH="${INITIAL_PATH:-$E5F_PF_PERSON_POLICY_INITIAL_PATH}"
    [[ "$INITIAL_PATH" = /* ]] || INITIAL_PATH="$PROJECT_ROOT/$INITIAL_PATH"
    [ -s "$INITIAL_PATH" ] || { echo "missing initial path: $INITIAL_PATH" >&2; exit 2; }
    ARGS+=(--initial-path "$INITIAL_PATH")
fi

"$PYTHON_BIN" "$DRIVER" "${ARGS[@]}" >"$OUTDIR/driver.log" 2>&1 &
CHILD_PID=$!
while kill -0 "$CHILD_PID" 2>/dev/null; do
    "$PYTHON_BIN" - "$OUTDIR" <<'PY'
import json, os, sys
from datetime import datetime, timezone
outdir = sys.argv[1]
latest = None
path = os.path.join(outdir, "latest_iteration.json")
if os.path.isfile(path):
    try:
        with open(path, encoding="utf-8") as handle:
            latest = json.load(handle)
    except Exception:
        pass
payload = {
    "status": "running",
    "utc": datetime.now(timezone.utc).isoformat(),
    "latest_iteration": latest,
    "driver_log_bytes": os.path.getsize(os.path.join(outdir, "driver.log")),
}
target = os.path.join(outdir, "heartbeat.json")
with open(target + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(target + ".tmp", target)
PY
    sleep 300
done
wait "$CHILD_PID"

"$PYTHON_BIN" - "$OUTDIR/summary.json" <<'PY'
import json, sys
from pathlib import Path
summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if summary.get("status") != "complete_unpromoted_person_demography_policy_path":
    raise SystemExit("person-demography policy packet is incomplete")
if not all((summary.get("accounting_gates") or {}).values()):
    raise SystemExit("person-demography policy accounting gates failed")
PY
