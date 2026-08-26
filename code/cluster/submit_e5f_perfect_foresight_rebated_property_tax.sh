#!/usr/bin/env bash
#SBATCH --job-name=e5f_pf_tax
#SBATCH --output=logs/slurm_e5f_pf_tax_%j.out
#SBATCH --error=logs/slurm_e5f_pf_tax_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_perfect_foresight_rebated_property_tax.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_PF_TAX_MODE:?E5F_PF_TAX_MODE is required}"
    : "${E5F_PF_TAX_RUN_TAG:?E5F_PF_TAX_RUN_TAG is required}"
    if [ "$E5F_PF_TAX_MODE" = "smoke" ]; then
        REQUESTED_TIME="${E5F_PF_TAX_TIME_LIMIT:-01:00:00}"
    elif [ "$E5F_PF_TAX_MODE" = "convergence" ]; then
        : "${E5F_PF_TAX_SMOKE_SUMMARY:?convergence requires a smoke summary}"
        : "${E5F_PF_TAX_HORIZONS:?convergence requires horizons}"
        REQUESTED_TIME="${E5F_PF_TAX_TIME_LIMIT:-04:00:00}"
    else
        echo "E5F_PF_TAX_MODE must be smoke or convergence" >&2
        exit 2
    fi
    REQUESTED_MEMORY="${E5F_PF_TAX_MEMORY:-8G}"
    REQUESTED_CPUS="${E5F_PF_TAX_CPUS_PER_TASK:-1}"
    cd "$SCRIPT_DIR"
    exec sbatch --time="$REQUESTED_TIME" --mem="$REQUESTED_MEMORY" \
        --cpus-per-task="$REQUESTED_CPUS" --export=ALL \
        "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_PF_TAX_MODE:?E5F_PF_TAX_MODE is required}"
: "${E5F_PF_TAX_RUN_TAG:?E5F_PF_TAX_RUN_TAG is required}"
: "${E5F_PF_TAX_CASES:?E5F_PF_TAX_CASES is required}"
: "${E5F_PF_TAX_EXPECTED_DRIVER_SHA256:?policy driver hash is required}"
: "${E5F_PF_TAX_EXPECTED_PF_DRIVER_SHA256:?PF driver hash is required}"
: "${E5F_PF_TAX_EXPECTED_SOLVER_SHA256:?Bellman solver hash is required}"
: "${E5F_PF_TAX_EXPECTED_CONTINUATION_SHA256:?continuation hash is required}"
: "${E5F_PF_TAX_EXPECTED_TAX_HELPER_SHA256:?tax helper hash is required}"
: "${E5F_PF_TAX_EXPECTED_IMPACT_HELPER_SHA256:?impact helper hash is required}"

if [[ ! "$E5F_PF_TAX_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $E5F_PF_TAX_RUN_TAG" >&2
    exit 2
fi

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_rebated_property_tax.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
TAX_HELPER="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_rebated_property_tax_smoke.py"
IMPACT_HELPER="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_coven_property_tax_smoke.py"
REPORT="$PROJECT_ROOT/output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/summary.json"
SELECTED_TRANSITION="$PROJECT_ROOT/output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/selected_transition_path.csv"
SOURCE="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
OUTDIR="$PROJECT_ROOT/output/model/e5f_perfect_foresight_rebated_tax_${E5F_PF_TAX_RUN_TAG}_${E5F_PF_TAX_MODE}"
mkdir -p "$CLUSTER_DIR/logs"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    if [ "$(hash_file "$path")" != "$expected" ]; then
        echo "$label hash mismatch" >&2
        exit 2
    fi
}
check_hash "$DRIVER" "$E5F_PF_TAX_EXPECTED_DRIVER_SHA256" "policy driver"
check_hash "$PF_DRIVER" "$E5F_PF_TAX_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$SOLVER" "$E5F_PF_TAX_EXPECTED_SOLVER_SHA256" "Bellman solver"
check_hash "$CONTINUATION" "$E5F_PF_TAX_EXPECTED_CONTINUATION_SHA256" "continuation"
check_hash "$TAX_HELPER" "$E5F_PF_TAX_EXPECTED_TAX_HELPER_SHA256" "tax helper"
check_hash "$IMPACT_HELPER" "$E5F_PF_TAX_EXPECTED_IMPACT_HELPER_SHA256" "impact helper"
for required in "$REPORT" "$SELECTED_TRANSITION" "$SOURCE"; do
    [ -s "$required" ] || { echo "missing input: $required" >&2; exit 2; }
done
if [ -e "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    echo "refusing to overwrite nonempty output: $OUTDIR" >&2
    exit 2
fi

IFS=',' read -r -a CASES <<< "$E5F_PF_TAX_CASES"
for case_name in "${CASES[@]}"; do
    case "$case_name" in
        rebated-tax1-baseline|rebated-tax2-reform) ;;
        *) echo "invalid policy case: $case_name" >&2; exit 2 ;;
    esac
done

if [ "$E5F_PF_TAX_MODE" = "smoke" ]; then
    HORIZONS=(2)
    MAX_ITERATIONS="${E5F_PF_TAX_MAX_ITERATIONS:-1}"
    ENDPOINT_GRID="${E5F_PF_TAX_ENDPOINT_GRID_POINTS:-5}"
else
    IFS=',' read -r -a HORIZONS <<< "$E5F_PF_TAX_HORIZONS"
    MAX_ITERATIONS="${E5F_PF_TAX_MAX_ITERATIONS:-8}"
    ENDPOINT_GRID="${E5F_PF_TAX_ENDPOINT_GRID_POINTS:-25}"
    SMOKE_SUMMARY="$E5F_PF_TAX_SMOKE_SUMMARY"
    [[ "$SMOKE_SUMMARY" = /* ]] || SMOKE_SUMMARY="$PROJECT_ROOT/$SMOKE_SUMMARY"
    python3 - "$SMOKE_SUMMARY" "$E5F_PF_TAX_EXPECTED_DRIVER_SHA256" \
        "$E5F_PF_TAX_EXPECTED_PF_DRIVER_SHA256" \
        "${E5F_PF_TAX_TERMINAL_PSI_CHILD:-}" "${CASES[@]}" <<'PY'
import json, math, sys
from pathlib import Path
summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
driver_sha, pf_driver_sha = sys.argv[2:4]
requested_terminal_psi = sys.argv[4]
if summary.get("status") != "complete_unpromoted_perfect_foresight_funded_policy_diagnostic":
    raise SystemExit("policy smoke summary is incomplete")
provenance = summary.get("provenance") or {}
if provenance.get("driver_sha256") != driver_sha:
    raise SystemExit("smoke used a different funded-policy driver")
if provenance.get("perfect_foresight_driver_sha256") != pf_driver_sha:
    raise SystemExit("smoke used a different perfect-foresight driver")
if requested_terminal_psi and not math.isclose(
    float(summary.get("terminal_psi_child")),
    float(requested_terminal_psi),
    rel_tol=0.0,
    abs_tol=1e-13,
):
    raise SystemExit("smoke used a different terminal preference value")
for case in sys.argv[5:]:
    terminal = (summary.get("terminal_steady_states") or {}).get(case) or {}
    if (terminal.get("fixed_point_gates") or {}).get("status") != "passed":
        raise SystemExit(f"smoke terminal gate failed for {case}")
    horizon = ((summary.get("solutions") or {}).get(case) or {}).get("2") or {}
    if horizon.get("status") not in {"converged", "maximum_iterations_reached"}:
        raise SystemExit(f"smoke did not execute the two-date loop for {case}")
PY
fi

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
mkdir -p "$OUTDIR"
COMPUTE_THREADS="${E5F_PF_TAX_NUMBA_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export NUMBA_NUM_THREADS="$COMPUTE_THREADS"
export OMP_NUM_THREADS="$COMPUTE_THREADS"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

"$PYTHON_BIN" - "$OUTDIR/launch_contract.json" "${HORIZONS[@]}" <<'PY'
import json, os, sys
path, *horizons = sys.argv[1:]
payload = {
    "status": "launched",
    "mode": os.environ["E5F_PF_TAX_MODE"],
    "run_tag": os.environ["E5F_PF_TAX_RUN_TAG"],
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "cases": os.environ["E5F_PF_TAX_CASES"].split(","),
    "horizons": [int(value) for value in horizons],
    "maximum_path_iterations": int(os.environ.get("E5F_PF_TAX_MAX_ITERATIONS", "1")),
    "initial_paths": os.environ.get("E5F_PF_TAX_INITIAL_PATHS"),
    "terminal_psi_child": os.environ.get("E5F_PF_TAX_TERMINAL_PSI_CHILD"),
    "price_damping": os.environ.get("E5F_PF_TAX_PRICE_DAMPING", "0.25"),
    "transfer_damping": os.environ.get("E5F_PF_TAX_TRANSFER_DAMPING", "0.50"),
    "maximum_log_price_step": os.environ.get("E5F_PF_TAX_MAXIMUM_LOG_PRICE_STEP", "0.12"),
    "maximum_transfer_step": os.environ.get("E5F_PF_TAX_MAXIMUM_TRANSFER_STEP", "0.08"),
    "hashes": {key: value for key, value in os.environ.items() if key.startswith("E5F_PF_TAX_EXPECTED_")},
}
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY

ARGS=(
    --selected-report "$REPORT"
    --selected-transition "$SELECTED_TRANSITION"
    --source "$SOURCE"
    --output-dir "$OUTDIR"
    --horizons "${HORIZONS[@]}"
    --cases "${CASES[@]}"
    --psi-persistence 0.60
    --housing-supply-elasticity 1.75
    --endpoint-price-grid-points "$ENDPOINT_GRID"
    --market-tol 2e-4
    --fiscal-tol 2.5e-5
    --max-path-iterations "$MAX_ITERATIONS"
    --price-damping "${E5F_PF_TAX_PRICE_DAMPING:-0.25}"
    --transfer-damping "${E5F_PF_TAX_TRANSFER_DAMPING:-0.50}"
    --maximum-log-price-step "${E5F_PF_TAX_MAXIMUM_LOG_PRICE_STEP:-0.12}"
    --maximum-transfer-step "${E5F_PF_TAX_MAXIMUM_TRANSFER_STEP:-0.08}"
)
if [ -n "${E5F_PF_TAX_TERMINAL_PSI_CHILD:-}" ]; then
    ARGS+=(--terminal-psi-child "$E5F_PF_TAX_TERMINAL_PSI_CHILD")
fi
if [ -n "${E5F_PF_TAX_INITIAL_PATHS:-}" ]; then
    IFS=';' read -r -a INITIAL_PATHS <<< "$E5F_PF_TAX_INITIAL_PATHS"
    for specification in "${INITIAL_PATHS[@]}"; do
        case_name="${specification%%=*}"
        source_path="${specification#*=}"
        [[ "$source_path" = /* ]] || source_path="$PROJECT_ROOT/$source_path"
        [ -s "$source_path" ] || { echo "missing initial path: $source_path" >&2; exit 2; }
        ARGS+=(--initial-path "$case_name=$source_path")
    done
fi

DRIVER_LOG="$OUTDIR/driver.log"
"$PYTHON_BIN" "$DRIVER" "${ARGS[@]}" >"$DRIVER_LOG" 2>&1 &
CHILD_PID=$!
while kill -0 "$CHILD_PID" 2>/dev/null; do
    "$PYTHON_BIN" - "$OUTDIR" <<'PY'
import json, os, sys
from datetime import datetime, timezone
outdir = sys.argv[1]
latest = []
for root, _, files in os.walk(outdir):
    if "latest_iteration.json" in files:
        path = os.path.join(root, "latest_iteration.json")
        try:
            with open(path, encoding="utf-8") as handle:
                item = json.load(handle)
            item["path"] = os.path.relpath(path, outdir)
            latest.append(item)
        except Exception:
            pass
payload = {
    "status": "running",
    "utc": datetime.now(timezone.utc).isoformat(),
    "latest_iterations": latest,
    "driver_log_bytes": os.path.getsize(os.path.join(outdir, "driver.log")),
}
path = os.path.join(outdir, "heartbeat.json")
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY
    sleep 300
done
wait "$CHILD_PID"

"$PYTHON_BIN" - "$OUTDIR/summary.json" "${CASES[@]}" <<'PY'
import json, sys
from pathlib import Path
summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if summary.get("status") != "complete_unpromoted_perfect_foresight_funded_policy_diagnostic":
    raise SystemExit("policy packet is incomplete")
for case in sys.argv[2:]:
    terminal = (summary.get("terminal_steady_states") or {}).get(case) or {}
    if (terminal.get("fixed_point_gates") or {}).get("status") != "passed":
        raise SystemExit(f"terminal fixed-point gate failed for {case}")
PY
