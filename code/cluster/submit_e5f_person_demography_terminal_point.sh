#!/usr/bin/env bash
#SBATCH --job-name=e5f_person_ss
#SBATCH --output=logs/slurm_e5f_person_ss_%j.out
#SBATCH --error=logs/slurm_e5f_person_ss_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_person_demography_terminal_point.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_PERSON_SS_RUN_TAG:?run tag is required}"
    cd "$SCRIPT_DIR"
    exec sbatch --time="${E5F_PERSON_SS_TIME_LIMIT:-00:45:00}" \
        --mem="${E5F_PERSON_SS_MEMORY:-8G}" --cpus-per-task=1 \
        --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_PERSON_SS_RUN_TAG:?run tag is required}"
: "${E5F_PERSON_SS_CASE:?policy case is required}"
: "${E5F_PERSON_SS_ASSET_PRICE:?asset price is required}"
: "${E5F_PERSON_SS_EQUAL_TRANSFER:?equal transfer is required}"
: "${E5F_PERSON_SS_PSI_CHILD:?psi_child is required}"
: "${E5F_PERSON_SS_EXPECTED_DRIVER_SHA256:?terminal-point driver hash is required}"
: "${E5F_PERSON_SS_EXPECTED_PERSON_DRIVER_SHA256:?person-demography hash is required}"
: "${E5F_PERSON_SS_EXPECTED_FUNDED_DRIVER_SHA256:?funded-policy hash is required}"
: "${E5F_PERSON_SS_EXPECTED_PF_DRIVER_SHA256:?PF driver hash is required}"
: "${E5F_PERSON_SS_EXPECTED_SOLVER_SHA256:?solver hash is required}"
: "${E5F_PERSON_SS_EXPECTED_CONTINUATION_SHA256:?continuation hash is required}"

case "$E5F_PERSON_SS_CASE" in
    rebated-tax1-baseline|rebated-tax2-reform) ;;
    *) echo "invalid case: $E5F_PERSON_SS_CASE" >&2; exit 2 ;;
esac
if [[ ! "$E5F_PERSON_SS_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $E5F_PERSON_SS_RUN_TAG" >&2
    exit 2
fi

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/build_e5f_person_demography_terminal_point.py"
PERSON_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography.py"
FUNDED_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_rebated_property_tax.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
SOURCE_DIR="$PROJECT_ROOT/output/model/e5f_persons_demographic_satellite_20260826b/source_data"
HEADSHIP_DIR="$PROJECT_ROOT/output/model/e5f_headship_demographic_bridge_20260826a"
OUTDIR="$PROJECT_ROOT/output/model/e5f_person_demography_terminal_${E5F_PERSON_SS_CASE}_${E5F_PERSON_SS_RUN_TAG}"
mkdir -p "$CLUSTER_DIR/logs"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    if [ "$(hash_file "$path")" != "$expected" ]; then
        echo "$label hash mismatch" >&2
        exit 2
    fi
}
check_hash "$DRIVER" "$E5F_PERSON_SS_EXPECTED_DRIVER_SHA256" "terminal-point driver"
check_hash "$PERSON_DRIVER" "$E5F_PERSON_SS_EXPECTED_PERSON_DRIVER_SHA256" "person-demography driver"
check_hash "$FUNDED_DRIVER" "$E5F_PERSON_SS_EXPECTED_FUNDED_DRIVER_SHA256" "funded-policy driver"
check_hash "$PF_DRIVER" "$E5F_PERSON_SS_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$SOLVER" "$E5F_PERSON_SS_EXPECTED_SOLVER_SHA256" "solver"
check_hash "$CONTINUATION" "$E5F_PERSON_SS_EXPECTED_CONTINUATION_SHA256" "continuation"
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
    "run_tag": os.environ["E5F_PERSON_SS_RUN_TAG"],
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "case": os.environ["E5F_PERSON_SS_CASE"],
    "asset_price": float(os.environ["E5F_PERSON_SS_ASSET_PRICE"]),
    "equal_transfer": float(os.environ["E5F_PERSON_SS_EQUAL_TRANSFER"]),
    "psi_child": float(os.environ["E5F_PERSON_SS_PSI_CHILD"]),
    "interpretation": "fixed policy point; inner demographic-household root only",
    "hashes": {
        key: value
        for key, value in os.environ.items()
        if key.startswith("E5F_PERSON_SS_EXPECTED_")
    },
}
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY

"$PYTHON_BIN" "$DRIVER" \
    --case "$E5F_PERSON_SS_CASE" \
    --asset-price "$E5F_PERSON_SS_ASSET_PRICE" \
    --equal-transfer "$E5F_PERSON_SS_EQUAL_TRANSFER" \
    --psi-child "$E5F_PERSON_SS_PSI_CHILD" \
    --source-dir "$SOURCE_DIR" \
    --headship-dir "$HEADSHIP_DIR" \
    --output-dir "$OUTDIR/result" \
    --maximum-iterations "${E5F_PERSON_SS_MAXIMUM_ITERATIONS:-250}" \
    --damping "${E5F_PERSON_SS_DAMPING:-0.50}" \
    >"$OUTDIR/driver.log" 2>&1

cp "$OUTDIR/result/summary.json" "$OUTDIR/summary.json"
