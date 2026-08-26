#!/usr/bin/env bash
#SBATCH --job-name=e5f_pf
#SBATCH --output=logs/slurm_e5f_pf_%j.out
#SBATCH --error=logs/slurm_e5f_pf_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_perfect_foresight_transition.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_PF_MODE:?E5F_PF_MODE is required}"
    : "${E5F_PF_RUN_TAG:?E5F_PF_RUN_TAG is required}"
    if [ "$E5F_PF_MODE" = "smoke" ]; then
        REQUESTED_TIME="${E5F_PF_TIME_LIMIT:-00:45:00}"
        REQUESTED_MEMORY="${E5F_PF_MEMORY:-8G}"
    elif [ "$E5F_PF_MODE" = "production" ]; then
        : "${E5F_PF_SMOKE_SUMMARY:?production requires E5F_PF_SMOKE_SUMMARY}"
        REQUESTED_TIME="${E5F_PF_TIME_LIMIT:-08:00:00}"
        REQUESTED_MEMORY="${E5F_PF_MEMORY:-8G}"
    elif [ "$E5F_PF_MODE" = "convergence" ]; then
        : "${E5F_PF_SMOKE_SUMMARY:?convergence requires E5F_PF_SMOKE_SUMMARY}"
        : "${E5F_PF_HORIZONS:?convergence requires E5F_PF_HORIZONS}"
        REQUESTED_TIME="${E5F_PF_TIME_LIMIT:-04:00:00}"
        REQUESTED_MEMORY="${E5F_PF_MEMORY:-8G}"
    else
        echo "E5F_PF_MODE must be smoke, production, or convergence" >&2
        exit 2
    fi
    REQUESTED_CPUS="${E5F_PF_CPUS_PER_TASK:-1}"
    if [[ ! "$REQUESTED_CPUS" =~ ^[0-9]+$ ]] || [ "$REQUESTED_CPUS" -lt 1 ]; then
        echo "E5F_PF_CPUS_PER_TASK must be a positive integer" >&2
        exit 2
    fi
    cd "$SCRIPT_DIR"
    exec sbatch --time="$REQUESTED_TIME" --mem="$REQUESTED_MEMORY" \
        --cpus-per-task="$REQUESTED_CPUS" \
        --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_PF_MODE:?E5F_PF_MODE is required}"
: "${E5F_PF_RUN_TAG:?E5F_PF_RUN_TAG is required}"
: "${E5F_PF_EXPECTED_DRIVER_SHA256:?driver hash is required}"
: "${E5F_PF_EXPECTED_SOLVER_SHA256:?solver hash is required}"
: "${E5F_PF_EXPECTED_CONTINUATION_SHA256:?closed-endpoint solver hash is required}"

if [[ ! "$E5F_PF_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $E5F_PF_RUN_TAG" >&2
    exit 2
fi

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
REPORT="$PROJECT_ROOT/output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/summary.json"
SELECTED_TRANSITION="$PROJECT_ROOT/output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/selected_transition_path.csv"
SOURCE="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
OUTDIR="$PROJECT_ROOT/output/model/e5f_perfect_foresight_${E5F_PF_RUN_TAG}_${E5F_PF_MODE}"
mkdir -p "$CLUSTER_DIR/logs"

INITIAL_PRICE_PATH=""
if [ -n "${E5F_PF_INITIAL_PRICE_PATH:-}" ]; then
    INITIAL_PRICE_PATH="$E5F_PF_INITIAL_PRICE_PATH"
    if [[ "$INITIAL_PRICE_PATH" != /* ]]; then
        INITIAL_PRICE_PATH="$PROJECT_ROOT/$INITIAL_PRICE_PATH"
    fi
    if [ ! -s "$INITIAL_PRICE_PATH" ]; then
        echo "missing initial price seed: $INITIAL_PRICE_PATH" >&2
        exit 2
    fi
    export E5F_PF_INITIAL_PRICE_PATH="$INITIAL_PRICE_PATH"
fi

hash_file() {
    sha256sum "$1" | awk '{print $1}'
}

if [ "$(hash_file "$DRIVER")" != "$E5F_PF_EXPECTED_DRIVER_SHA256" ]; then
    echo "perfect-foresight driver hash mismatch" >&2
    exit 2
fi
if [ "$(hash_file "$SOLVER")" != "$E5F_PF_EXPECTED_SOLVER_SHA256" ]; then
    echo "calendar Bellman solver hash mismatch" >&2
    exit 2
fi
if [ "$(hash_file "$CONTINUATION")" != "$E5F_PF_EXPECTED_CONTINUATION_SHA256" ]; then
    echo "closed-endpoint solver hash mismatch" >&2
    exit 2
fi
for required in "$REPORT" "$SELECTED_TRANSITION" "$SOURCE"; do
    if [ ! -s "$required" ]; then
        echo "missing required input: $required" >&2
        exit 2
    fi
done
if [ -e "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    echo "refusing to overwrite nonempty output directory: $OUTDIR" >&2
    exit 2
fi

if [ "$E5F_PF_MODE" = "production" ] || [ "$E5F_PF_MODE" = "convergence" ]; then
    SMOKE_SUMMARY="$E5F_PF_SMOKE_SUMMARY"
    if [[ "$SMOKE_SUMMARY" != /* ]]; then
        SMOKE_SUMMARY="$PROJECT_ROOT/$SMOKE_SUMMARY"
    fi
    python3 - "$SMOKE_SUMMARY" "$E5F_PF_EXPECTED_DRIVER_SHA256" \
        "$E5F_PF_EXPECTED_SOLVER_SHA256" \
        "$E5F_PF_EXPECTED_CONTINUATION_SHA256" <<'PY'
import json, math, os, sys
from pathlib import Path
summary_path, driver_sha, solver_sha, continuation_sha = sys.argv[1:]
summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
if summary.get("status") != "complete_isolated_perfect_foresight_diagnostic":
    raise SystemExit("smoke summary is incomplete")
if (summary.get("zero_shock_test") or {}).get("status") != "passed":
    raise SystemExit("smoke zero-shock gate did not pass")
terminal = summary.get("terminal_steady_state") or {}
if (terminal.get("fixed_point_gates") or {}).get("status") != "passed":
    raise SystemExit("smoke terminal fixed-point gate did not pass")
declared_terminal = os.environ.get("E5F_PF_TERMINAL_PSI_CHILD")
if declared_terminal and not math.isclose(
    float(terminal.get("psi_child")), float(declared_terminal),
    rel_tol=0.0, abs_tol=1e-14,
):
    raise SystemExit("smoke used a different terminal preference regime")
horizons = summary.get("horizons") or {}
if "2" not in horizons:
    raise SystemExit("smoke did not execute the exact two-date full-state loop")
provenance = summary.get("provenance") or {}
if provenance.get("perfect_foresight_solver_sha256") != driver_sha:
    raise SystemExit("smoke used a different perfect-foresight driver")
if provenance.get("calendar_bellman_sha256") != solver_sha:
    raise SystemExit("smoke used a different calendar Bellman solver")
if provenance.get("closed_endpoint_solver_sha256") != continuation_sha:
    raise SystemExit("smoke used a different closed-endpoint solver")
PY
fi

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
mkdir -p "$OUTDIR"

COMPUTE_THREADS="${E5F_PF_NUMBA_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
if [[ ! "$COMPUTE_THREADS" =~ ^[0-9]+$ ]] || [ "$COMPUTE_THREADS" -lt 1 ]; then
    echo "E5F_PF_NUMBA_THREADS must be a positive integer" >&2
    exit 2
fi
if [ "$COMPUTE_THREADS" -gt "${SLURM_CPUS_PER_TASK:-1}" ]; then
    echo "E5F_PF_NUMBA_THREADS cannot exceed SLURM_CPUS_PER_TASK" >&2
    exit 2
fi

if [ "$E5F_PF_MODE" = "smoke" ]; then
    HORIZONS=(2)
    MAX_ITERATIONS=1
elif [ "$E5F_PF_MODE" = "convergence" ]; then
    IFS=',' read -r -a HORIZONS <<< "$E5F_PF_HORIZONS"
    if [ "${#HORIZONS[@]}" -lt 1 ]; then
        echo "convergence mode requires at least one horizon" >&2
        exit 2
    fi
    for horizon in "${HORIZONS[@]}"; do
        if [[ ! "$horizon" =~ ^[0-9]+$ ]] || [ "$horizon" -lt 2 ]; then
            echo "invalid convergence horizon: $horizon" >&2
            exit 2
        fi
    done
    MAX_ITERATIONS="${E5F_PF_MAX_PRICE_ITERATIONS:-35}"
else
    HORIZONS=(8 12)
    MAX_ITERATIONS="${E5F_PF_MAX_PRICE_ITERATIONS:-30}"
fi

"$PYTHON_BIN" - "$OUTDIR/launch_contract.json" "$E5F_PF_MODE" \
    "$E5F_PF_RUN_TAG" "$SLURM_JOB_ID" "$MAX_ITERATIONS" \
    "$COMPUTE_THREADS" "${SLURM_CPUS_PER_TASK:-1}" \
    "${HORIZONS[@]}" <<'PY'
import json, os, sys
path, mode, tag, job, iterations, threads, allocated_cpus, *horizons = sys.argv[1:]
payload = {
    "status": "launched",
    "mode": mode,
    "run_tag": tag,
    "slurm_job_id": job,
    "horizons": [int(value) for value in horizons],
    "maximum_price_iterations": int(iterations),
    "numba_threads": int(threads),
    "allocated_cpus_per_task": int(allocated_cpus),
    "driver_sha256": os.environ["E5F_PF_EXPECTED_DRIVER_SHA256"],
    "solver_sha256": os.environ["E5F_PF_EXPECTED_SOLVER_SHA256"],
    "closed_endpoint_solver_sha256": os.environ[
        "E5F_PF_EXPECTED_CONTINUATION_SHA256"
    ],
    "terminal_psi_child": (
        float(os.environ["E5F_PF_TERMINAL_PSI_CHILD"])
        if os.environ.get("E5F_PF_TERMINAL_PSI_CHILD")
        else None
    ),
    "initial_price_path": os.environ.get("E5F_PF_INITIAL_PRICE_PATH"),
    "requested_memory": os.environ.get("SLURM_MEM_PER_NODE"),
}
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY

export NUMBA_NUM_THREADS="$COMPUTE_THREADS"
export OMP_NUM_THREADS="$COMPUTE_THREADS"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
DRIVER_LOG="$OUTDIR/driver.log"
ARGS=(
    --selected-report "$REPORT"
    --selected-transition "$SELECTED_TRANSITION"
    --source "$SOURCE"
    --output-dir "$OUTDIR"
    --horizons "${HORIZONS[@]}"
    --psi-persistence 0.60
    --market-tol 2e-4
    --max-price-iterations "$MAX_ITERATIONS"
    --price-damping 0.25
    --maximum-log-price-step 0.12
)
if [ -n "${E5F_PF_TERMINAL_PSI_CHILD:-}" ]; then
    ARGS+=(--terminal-psi-child "$E5F_PF_TERMINAL_PSI_CHILD")
fi
if [ -n "$INITIAL_PRICE_PATH" ]; then
    ARGS+=(--initial-price-path "$INITIAL_PRICE_PATH")
fi

"$PYTHON_BIN" "$DRIVER" "${ARGS[@]}" >"$DRIVER_LOG" 2>&1 &
CHILD_PID=$!
while kill -0 "$CHILD_PID" 2>/dev/null; do
    "$PYTHON_BIN" - "$OUTDIR" "$SLURM_JOB_ID" <<'PY'
import json, os, sys
from datetime import datetime, timezone
outdir, job = sys.argv[1:]
latest = None
for root, _, files in os.walk(outdir):
    if "latest_iteration.json" in files:
        path = os.path.join(root, "latest_iteration.json")
        try:
            with open(path, encoding="utf-8") as handle:
                candidate = json.load(handle)
            if latest is None or candidate.get("elapsed_seconds", 0) > latest.get("elapsed_seconds", 0):
                latest = candidate
        except Exception:
            pass
payload = {
    "status": "running",
    "utc": datetime.now(timezone.utc).isoformat(),
    "slurm_job_id": job,
    "latest_iteration": latest,
    "driver_log_bytes": os.path.getsize(os.path.join(outdir, "driver.log")) if os.path.isfile(os.path.join(outdir, "driver.log")) else 0,
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

"$PYTHON_BIN" - "$OUTDIR/summary.json" "$OUTDIR/launch_contract.json" <<'PY'
import json, sys
from pathlib import Path
path, launch_path = sys.argv[1:]
summary = json.loads(Path(path).read_text(encoding="utf-8"))
launch = json.loads(Path(launch_path).read_text(encoding="utf-8"))
if summary.get("status") != "complete_isolated_perfect_foresight_diagnostic":
    raise SystemExit("perfect-foresight packet is incomplete")
if (summary.get("zero_shock_test") or {}).get("status") != "passed":
    raise SystemExit("zero-shock gate failed")
if (
    ((summary.get("terminal_steady_state") or {}).get("fixed_point_gates") or {})
    .get("status") != "passed"
):
    raise SystemExit("terminal fixed-point gate failed")
expected = {str(value) for value in launch.get("horizons", [])}
if set(summary.get("horizons") or {}) != expected:
    raise SystemExit("perfect-foresight packet has the wrong horizons")
for horizon, diagnostics in (summary.get("horizons") or {}).items():
    terminal = diagnostics.get("terminal_convergence") or {}
    if terminal.get("status") not in {"passed", "not_converged"}:
        raise SystemExit(f"horizon {horizon} lacks a terminal convergence verdict")
PY
