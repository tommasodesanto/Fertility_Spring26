#!/usr/bin/env bash
#SBATCH --job-name=e5f_person_root
#SBATCH --output=logs/slurm_e5f_person_root_%j.out
#SBATCH --error=logs/slurm_e5f_person_root_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_person_demography_terminal_root.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_PERSON_ROOT_RUN_TAG:?run tag is required}"
    cd "$SCRIPT_DIR"
    exec sbatch --time="${E5F_PERSON_ROOT_TIME_LIMIT:-04:00:00}" \
        --mem="${E5F_PERSON_ROOT_MEMORY:-8G}" --cpus-per-task=1 \
        --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_PERSON_ROOT_RUN_TAG:?run tag is required}"
: "${E5F_PERSON_ROOT_CASE:?policy case is required}"
: "${E5F_PERSON_ROOT_ASSET_PRICE:?initial asset price is required}"
: "${E5F_PERSON_ROOT_EQUAL_TRANSFER:?initial transfer is required}"
: "${E5F_PERSON_ROOT_PSI_CHILD:?psi_child is required}"
: "${E5F_PERSON_ROOT_EXPECTED_DRIVER_SHA256:?terminal-root driver hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_PERSON_DRIVER_SHA256:?person-demography hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_POINT_DRIVER_SHA256:?terminal-point driver hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_FUNDED_DRIVER_SHA256:?funded-policy hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_PF_DRIVER_SHA256:?PF driver hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_SOLVER_SHA256:?solver hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_CONTINUATION_SHA256:?continuation hash is required}"
: "${E5F_PERSON_ROOT_EXPECTED_TENURE_SENSITIVITY_SHA256:?tenure-sensitivity hash is required}"

case "$E5F_PERSON_ROOT_CASE" in
    rebated-tax1-baseline|rebated-tax2-reform) ;;
    *) echo "invalid case: $E5F_PERSON_ROOT_CASE" >&2; exit 2 ;;
esac
if [[ ! "$E5F_PERSON_ROOT_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $E5F_PERSON_ROOT_RUN_TAG" >&2
    exit 2
fi

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/solve_e5f_person_demography_terminal_root.py"
PERSON_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography.py"
POINT_DRIVER="$PROJECT_ROOT/code/model/tools/build_e5f_person_demography_terminal_point.py"
FUNDED_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_rebated_property_tax.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
TENURE_SENSITIVITY="$PROJECT_ROOT/code/model/tools/e5f_post2023_tenure_sensitivity.py"
SOURCE_DIR="$PROJECT_ROOT/output/model/e5f_persons_demographic_satellite_20260826b/source_data"
HEADSHIP_DIR="$PROJECT_ROOT/output/model/e5f_headship_demographic_bridge_20260826a"
OUTDIR="$PROJECT_ROOT/output/model/e5f_person_demography_terminal_root_${E5F_PERSON_ROOT_CASE}_${E5F_PERSON_ROOT_RUN_TAG}"
mkdir -p "$CLUSTER_DIR/logs"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    if [ "$(hash_file "$path")" != "$expected" ]; then
        echo "$label hash mismatch" >&2
        exit 2
    fi
}
check_hash "$DRIVER" "$E5F_PERSON_ROOT_EXPECTED_DRIVER_SHA256" "terminal-root driver"
check_hash "$PERSON_DRIVER" "$E5F_PERSON_ROOT_EXPECTED_PERSON_DRIVER_SHA256" "person-demography driver"
check_hash "$POINT_DRIVER" "$E5F_PERSON_ROOT_EXPECTED_POINT_DRIVER_SHA256" "terminal-point driver"
check_hash "$FUNDED_DRIVER" "$E5F_PERSON_ROOT_EXPECTED_FUNDED_DRIVER_SHA256" "funded-policy driver"
check_hash "$PF_DRIVER" "$E5F_PERSON_ROOT_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$SOLVER" "$E5F_PERSON_ROOT_EXPECTED_SOLVER_SHA256" "solver"
check_hash "$CONTINUATION" "$E5F_PERSON_ROOT_EXPECTED_CONTINUATION_SHA256" "continuation"
check_hash "$TENURE_SENSITIVITY" "$E5F_PERSON_ROOT_EXPECTED_TENURE_SENSITIVITY_SHA256" "tenure sensitivity"
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
    "run_tag": os.environ["E5F_PERSON_ROOT_RUN_TAG"],
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "case": os.environ["E5F_PERSON_ROOT_CASE"],
    "initial_asset_price": float(os.environ["E5F_PERSON_ROOT_ASSET_PRICE"]),
    "initial_equal_transfer": float(os.environ["E5F_PERSON_ROOT_EQUAL_TRANSFER"]),
    "psi_child": float(os.environ["E5F_PERSON_ROOT_PSI_CHILD"]),
    "housing_supply_elasticity": float(
        os.environ.get("E5F_PERSON_ROOT_HOUSING_SUPPLY_ELASTICITY", "1.75")
    ),
    "post2023_tenure_choice_kappa": (
        float(os.environ["E5F_PERSON_ROOT_POST2023_TENURE_CHOICE_KAPPA"])
        if os.environ.get("E5F_PERSON_ROOT_POST2023_TENURE_CHOICE_KAPPA")
        else None
    ),
    "maximum_root_evaluations": int(os.environ.get("E5F_PERSON_ROOT_MAXIMUM_EVALUATIONS", "24")),
    "maximum_root_iterations": int(os.environ.get("E5F_PERSON_ROOT_MAXIMUM_ITERATIONS", "16")),
    "log_price_difference_step": float(os.environ.get("E5F_PERSON_ROOT_LOG_PRICE_DIFFERENCE_STEP", "0.01")),
    "transfer_difference_step": float(os.environ.get("E5F_PERSON_ROOT_TRANSFER_DIFFERENCE_STEP", "0.01")),
    "maximum_inner_iterations": int(os.environ.get("E5F_PERSON_ROOT_MAXIMUM_INNER_ITERATIONS", "250")),
    "inner_damping": float(os.environ.get("E5F_PERSON_ROOT_INNER_DAMPING", "0.50")),
    "interpretation": "joint housing-market, rebate-budget, household, and person-demography terminal root",
    "hashes": {
        key: value
        for key, value in os.environ.items()
        if key.startswith("E5F_PERSON_ROOT_EXPECTED_")
    },
}
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY

ARGS=(
    --case "$E5F_PERSON_ROOT_CASE"
    --psi-child "$E5F_PERSON_ROOT_PSI_CHILD"
    --initial-asset-price "$E5F_PERSON_ROOT_ASSET_PRICE"
    --initial-equal-transfer "$E5F_PERSON_ROOT_EQUAL_TRANSFER"
    --source-dir "$SOURCE_DIR"
    --headship-dir "$HEADSHIP_DIR"
    --output-dir "$OUTDIR/result"
    --maximum-root-evaluations "${E5F_PERSON_ROOT_MAXIMUM_EVALUATIONS:-24}"
    --maximum-root-iterations "${E5F_PERSON_ROOT_MAXIMUM_ITERATIONS:-16}"
    --log-price-difference-step "${E5F_PERSON_ROOT_LOG_PRICE_DIFFERENCE_STEP:-0.01}"
    --transfer-difference-step "${E5F_PERSON_ROOT_TRANSFER_DIFFERENCE_STEP:-0.01}"
    --maximum-inner-iterations "${E5F_PERSON_ROOT_MAXIMUM_INNER_ITERATIONS:-250}"
    --inner-damping "${E5F_PERSON_ROOT_INNER_DAMPING:-0.50}"
    --housing-supply-elasticity "${E5F_PERSON_ROOT_HOUSING_SUPPLY_ELASTICITY:-1.75}"
)
if [ -n "${E5F_PERSON_ROOT_POST2023_TENURE_CHOICE_KAPPA:-}" ]; then
    ARGS+=(--post2023-tenure-choice-kappa "$E5F_PERSON_ROOT_POST2023_TENURE_CHOICE_KAPPA")
fi

"$PYTHON_BIN" "$DRIVER" "${ARGS[@]}" \
    >"$OUTDIR/driver.log" 2>&1

cp "$OUTDIR/result/summary.json" "$OUTDIR/summary.json"
