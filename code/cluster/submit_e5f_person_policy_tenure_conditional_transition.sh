#!/usr/bin/env bash
#SBATCH --job-name=e5f_pf_cond_tenure
#SBATCH --output=logs/slurm_e5f_pf_cond_tenure_%j.out
#SBATCH --error=logs/slurm_e5f_pf_cond_tenure_%j.err
#SBATCH --partition=cs
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_person_policy_tenure_conditional_transition.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_TENURE_CONDITIONAL_RUN_TAG:?run tag is required}"
    cd "$SCRIPT_DIR"
    exec sbatch --time="${E5F_TENURE_CONDITIONAL_TIME_LIMIT:-16:00:00}" \
        --mem="${E5F_TENURE_CONDITIONAL_MEMORY:-16G}" \
        --partition="${E5F_TENURE_CONDITIONAL_PARTITION:-cs}" \
        --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_TENURE_CONDITIONAL_RUN_TAG:?run tag is required}"
: "${E5F_TENURE_CONDITIONAL_TERMINAL_SUMMARY:?terminal summary is required}"
: "${E5F_TENURE_CONDITIONAL_INITIAL_PATH:?initial path is required}"
: "${E5F_TENURE_CONDITIONAL_OWNER_SHARE:?owner share is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_DRIVER_SHA256:?conditional driver hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_ORACLE_SHA256:?oracle hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_KINK_SHA256:?kink helper hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_PRODUCTION_DRIVER_SHA256:?production driver hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_PERSON_DRIVER_SHA256:?person driver hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_PF_DRIVER_SHA256:?PF driver hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_FUNDED_DRIVER_SHA256:?funded-policy hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_CONTINUATION_SHA256:?continuation hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_SOLVER_SHA256:?solver hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_TERMINAL_SHA256:?terminal hash is required}"
: "${E5F_TENURE_CONDITIONAL_EXPECTED_INITIAL_PATH_SHA256:?initial path hash is required}"

[[ "$E5F_TENURE_CONDITIONAL_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "invalid run tag" >&2
    exit 2
}

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_person_policy_tenure_conditional_transition.py"
ORACLE="$PROJECT_ROOT/code/model/tools/run_e5f_person_policy_tenure_complementarity_oracle.py"
KINK="$PROJECT_ROOT/code/model/tools/run_e5f_person_policy_tenure_kink_diagnostic.py"
PRODUCTION_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography_policy.py"
PERSON_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
FUNDED_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_rebated_property_tax.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
TERMINAL="$E5F_TENURE_CONDITIONAL_TERMINAL_SUMMARY"
INITIAL="$E5F_TENURE_CONDITIONAL_INITIAL_PATH"
[[ "$TERMINAL" = /* ]] || TERMINAL="$PROJECT_ROOT/$TERMINAL"
[[ "$INITIAL" = /* ]] || INITIAL="$PROJECT_ROOT/$INITIAL"
OUTDIR="$PROJECT_ROOT/output/model/e5f_pf_person_policy_tenure_conditional_${E5F_TENURE_CONDITIONAL_RUN_TAG}"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    [ -s "$path" ] || { echo "missing $label: $path" >&2; exit 2; }
    [ "$(hash_file "$path")" = "$expected" ] || {
        echo "$label hash mismatch" >&2
        exit 2
    }
}
check_hash "$DRIVER" "$E5F_TENURE_CONDITIONAL_EXPECTED_DRIVER_SHA256" "conditional driver"
check_hash "$ORACLE" "$E5F_TENURE_CONDITIONAL_EXPECTED_ORACLE_SHA256" "oracle helper"
check_hash "$KINK" "$E5F_TENURE_CONDITIONAL_EXPECTED_KINK_SHA256" "kink helper"
check_hash "$PRODUCTION_DRIVER" "$E5F_TENURE_CONDITIONAL_EXPECTED_PRODUCTION_DRIVER_SHA256" "production driver"
check_hash "$PERSON_DRIVER" "$E5F_TENURE_CONDITIONAL_EXPECTED_PERSON_DRIVER_SHA256" "person driver"
check_hash "$PF_DRIVER" "$E5F_TENURE_CONDITIONAL_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$FUNDED_DRIVER" "$E5F_TENURE_CONDITIONAL_EXPECTED_FUNDED_DRIVER_SHA256" "funded-policy driver"
check_hash "$CONTINUATION" "$E5F_TENURE_CONDITIONAL_EXPECTED_CONTINUATION_SHA256" "continuation driver"
check_hash "$SOLVER" "$E5F_TENURE_CONDITIONAL_EXPECTED_SOLVER_SHA256" "solver"
check_hash "$TERMINAL" "$E5F_TENURE_CONDITIONAL_EXPECTED_TERMINAL_SHA256" "terminal summary"
check_hash "$INITIAL" "$E5F_TENURE_CONDITIONAL_EXPECTED_INITIAL_PATH_SHA256" "initial path"
if [ -e "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]; then
    echo "refusing to overwrite nonempty output: $OUTDIR" >&2
    exit 2
fi

mkdir -p "$CLUSTER_DIR/logs" "$OUTDIR"
python3 - "$OUTDIR/launch_contract.json" <<'PY'
import json, os, sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
    "status": "launched",
    "job_id": os.environ["SLURM_JOB_ID"],
    "run_tag": os.environ["E5F_TENURE_CONDITIONAL_RUN_TAG"],
    "horizon": int(os.environ.get("E5F_TENURE_CONDITIONAL_HORIZON", "128")),
    "mix_calendar_year": int(os.environ.get("E5F_TENURE_CONDITIONAL_CALENDAR_YEAR", "2351")),
    "owner_share": float(os.environ["E5F_TENURE_CONDITIONAL_OWNER_SHARE"]),
    "maximum_path_iterations": int(os.environ.get("E5F_TENURE_CONDITIONAL_MAX_ITERATIONS", "35")),
    "price_damping": float(os.environ.get("E5F_TENURE_CONDITIONAL_PRICE_DAMPING", "0.25")),
    "transfer_damping": float(os.environ.get("E5F_TENURE_CONDITIONAL_TRANSFER_DAMPING", "0.50")),
    "maximum_log_price_step": float(os.environ.get("E5F_TENURE_CONDITIONAL_MAX_LOG_PRICE_STEP", "0.12")),
    "maximum_transfer_step": float(os.environ.get("E5F_TENURE_CONDITIONAL_MAX_TRANSFER_STEP", "0.08")),
    "terminal_summary_sha256": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_TERMINAL_SHA256"],
    "initial_path_sha256": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_INITIAL_PATH_SHA256"],
    "source_hashes": {
        "conditional_driver": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_DRIVER_SHA256"],
        "oracle_helper": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_ORACLE_SHA256"],
        "kink_helper": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_KINK_SHA256"],
        "production_driver": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_PRODUCTION_DRIVER_SHA256"],
        "person_demography": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_PERSON_DRIVER_SHA256"],
        "perfect_foresight": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_PF_DRIVER_SHA256"],
        "funded_policy": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_FUNDED_DRIVER_SHA256"],
        "continuation": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_CONTINUATION_SHA256"],
        "solver": os.environ["E5F_TENURE_CONDITIONAL_EXPECTED_SOLVER_SHA256"],
    },
    "scheduler": {
        "partition": os.environ.get("SLURM_JOB_PARTITION"),
        "qos": os.environ.get("SLURM_JOB_QOS"),
        "requested_time_limit": os.environ.get("E5F_TENURE_CONDITIONAL_TIME_LIMIT", "16:00:00"),
        "requested_memory": os.environ.get("E5F_TENURE_CONDITIONAL_MEMORY", "16G"),
    },
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

"$PYTHON_BIN" "$DRIVER" \
    --terminal-summary "$TERMINAL" \
    --initial-path "$INITIAL" \
    --output-dir "$OUTDIR" \
    --horizon "${E5F_TENURE_CONDITIONAL_HORIZON:-128}" \
    --mix-calendar-year "${E5F_TENURE_CONDITIONAL_CALENDAR_YEAR:-2351}" \
    --owner-share "$E5F_TENURE_CONDITIONAL_OWNER_SHARE" \
    --maximum-path-iterations "${E5F_TENURE_CONDITIONAL_MAX_ITERATIONS:-35}" \
    --price-damping "${E5F_TENURE_CONDITIONAL_PRICE_DAMPING:-0.25}" \
    --transfer-damping "${E5F_TENURE_CONDITIONAL_TRANSFER_DAMPING:-0.50}" \
    --maximum-log-price-step "${E5F_TENURE_CONDITIONAL_MAX_LOG_PRICE_STEP:-0.12}" \
    --maximum-transfer-step "${E5F_TENURE_CONDITIONAL_MAX_TRANSFER_STEP:-0.08}"
