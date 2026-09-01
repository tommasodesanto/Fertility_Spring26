#!/usr/bin/env bash
#SBATCH --job-name=e5f_pf_fert_mech
#SBATCH --output=logs/slurm_e5f_pf_fert_mech_%j.out
#SBATCH --error=logs/slurm_e5f_pf_fert_mech_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_person_policy_fertility_mechanism_diagnostic.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_FERT_MECH_RUN_TAG:?run tag is required}"
    cd "$SCRIPT_DIR"
    exec sbatch --time="${E5F_FERT_MECH_TIME_LIMIT:-01:00:00}" \
        --mem="${E5F_FERT_MECH_MEMORY:-16G}" \
        --partition="${E5F_FERT_MECH_PARTITION:-cpu_short}" \
        --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_FERT_MECH_RUN_TAG:?run tag is required}"
: "${E5F_FERT_MECH_CASE:?case is required}"
: "${E5F_FERT_MECH_TERMINAL_SUMMARY:?terminal summary is required}"
: "${E5F_FERT_MECH_INITIAL_PATH:?initial path is required}"
: "${E5F_FERT_MECH_TENURE_KAPPA:?positive tenure kappa is required}"
: "${E5F_FERT_MECH_EXPECTED_WRAPPER_SHA256:?wrapper hash is required}"
: "${E5F_FERT_MECH_EXPECTED_PRODUCTION_DRIVER_SHA256:?production driver hash is required}"
: "${E5F_FERT_MECH_EXPECTED_PERSON_DRIVER_SHA256:?person driver hash is required}"
: "${E5F_FERT_MECH_EXPECTED_PF_DRIVER_SHA256:?PF driver hash is required}"
: "${E5F_FERT_MECH_EXPECTED_SOLVER_SHA256:?solver hash is required}"
: "${E5F_FERT_MECH_EXPECTED_TENURE_SENSITIVITY_SHA256:?tenure sensitivity hash is required}"
: "${E5F_FERT_MECH_EXPECTED_TERMINAL_SHA256:?terminal hash is required}"
: "${E5F_FERT_MECH_EXPECTED_INITIAL_PATH_SHA256:?initial-path hash is required}"

[[ "$E5F_FERT_MECH_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "invalid run tag" >&2
    exit 2
}
case "$E5F_FERT_MECH_CASE" in
    rebated-tax1-baseline|rebated-tax2-reform) ;;
    *) echo "invalid policy case" >&2; exit 2 ;;
esac

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
WRAPPER="$PROJECT_ROOT/code/model/tools/run_e5f_person_policy_fertility_mechanism_diagnostic.py"
PRODUCTION_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography_policy.py"
PERSON_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
TENURE_SENSITIVITY="$PROJECT_ROOT/code/model/tools/e5f_post2023_tenure_sensitivity.py"
TERMINAL="$E5F_FERT_MECH_TERMINAL_SUMMARY"
INITIAL="$E5F_FERT_MECH_INITIAL_PATH"
[[ "$TERMINAL" = /* ]] || TERMINAL="$PROJECT_ROOT/$TERMINAL"
[[ "$INITIAL" = /* ]] || INITIAL="$PROJECT_ROOT/$INITIAL"
OUTDIR="$PROJECT_ROOT/output/model/e5f_pf_person_policy_fertility_mechanism_${E5F_FERT_MECH_RUN_TAG}"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    [ -s "$path" ] || { echo "missing $label: $path" >&2; exit 2; }
    [ "$(hash_file "$path")" = "$expected" ] || {
        echo "$label hash mismatch" >&2
        exit 2
    }
}
check_hash "$WRAPPER" "$E5F_FERT_MECH_EXPECTED_WRAPPER_SHA256" "wrapper"
check_hash "$PRODUCTION_DRIVER" "$E5F_FERT_MECH_EXPECTED_PRODUCTION_DRIVER_SHA256" "production driver"
check_hash "$PERSON_DRIVER" "$E5F_FERT_MECH_EXPECTED_PERSON_DRIVER_SHA256" "person driver"
check_hash "$PF_DRIVER" "$E5F_FERT_MECH_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$SOLVER" "$E5F_FERT_MECH_EXPECTED_SOLVER_SHA256" "solver"
check_hash "$TENURE_SENSITIVITY" "$E5F_FERT_MECH_EXPECTED_TENURE_SENSITIVITY_SHA256" "tenure sensitivity"
check_hash "$TERMINAL" "$E5F_FERT_MECH_EXPECTED_TERMINAL_SHA256" "terminal summary"
check_hash "$INITIAL" "$E5F_FERT_MECH_EXPECTED_INITIAL_PATH_SHA256" "initial path"
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
    "run_tag": os.environ["E5F_FERT_MECH_RUN_TAG"],
    "case": os.environ["E5F_FERT_MECH_CASE"],
    "horizon": int(os.environ.get("E5F_FERT_MECH_HORIZON", "128")),
    "diagnostic_years": os.environ.get("E5F_FERT_MECH_DIAGNOSTIC_YEARS", "2023,2035,2051,2079,2103,2355,2531"),
    "post2023_tenure_choice_kappa": float(os.environ["E5F_FERT_MECH_TENURE_KAPPA"]),
    "initial_path_sha256": os.environ["E5F_FERT_MECH_EXPECTED_INITIAL_PATH_SHA256"],
    "terminal_summary_sha256": os.environ["E5F_FERT_MECH_EXPECTED_TERMINAL_SHA256"],
    "wrapper_sha256": os.environ["E5F_FERT_MECH_EXPECTED_WRAPPER_SHA256"],
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

"$PYTHON_BIN" "$WRAPPER" \
    --case "$E5F_FERT_MECH_CASE" \
    --terminal-summary "$TERMINAL" \
    --initial-path "$INITIAL" \
    --output-dir "$OUTDIR" \
    --horizon "${E5F_FERT_MECH_HORIZON:-128}" \
    --diagnostic-years "${E5F_FERT_MECH_DIAGNOSTIC_YEARS:-2023,2035,2051,2079,2103,2355,2531}" \
    --post2023-tenure-choice-kappa "$E5F_FERT_MECH_TENURE_KAPPA"
