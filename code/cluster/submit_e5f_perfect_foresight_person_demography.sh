#!/usr/bin/env bash
#SBATCH --job-name=e5f_pf_person
#SBATCH --output=logs/slurm_e5f_pf_person_%j.out
#SBATCH --error=logs/slurm_e5f_pf_person_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_NAME="submit_e5f_perfect_foresight_person_demography.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--submit" ]; then
    shift
    : "${E5F_PF_PERSON_RUN_TAG:?run tag is required}"
    cd "$SCRIPT_DIR"
    exec sbatch --time="${E5F_PF_PERSON_TIME_LIMIT:-00:45:00}" \
        --mem="${E5F_PF_PERSON_MEMORY:-8G}" --cpus-per-task=1 \
        --export=ALL "$SCRIPT_DIR/$SCRIPT_NAME" "$@"
fi

: "${SLURM_JOB_ID:?Submit with --submit or sbatch}"
: "${E5F_PF_PERSON_RUN_TAG:?run tag is required}"
: "${E5F_PF_PERSON_EXPECTED_DRIVER_SHA256:?person-demography driver hash is required}"
: "${E5F_PF_PERSON_EXPECTED_PF_DRIVER_SHA256:?PF driver hash is required}"
: "${E5F_PF_PERSON_EXPECTED_SOLVER_SHA256:?Bellman solver hash is required}"
: "${E5F_PF_PERSON_EXPECTED_CONTINUATION_SHA256:?continuation hash is required}"
: "${E5F_PF_PERSON_EXPECTED_PERSON_LAW_SHA256:?person-law hash is required}"
: "${E5F_PF_PERSON_EXPECTED_COUPLING_SHA256:?coupling hash is required}"
: "${E5F_PF_PERSON_EXPECTED_DEMOGRAPHIC_BUILDER_SHA256:?demographic builder hash is required}"
: "${E5F_PF_PERSON_EXPECTED_SATELLITE_HELPER_SHA256:?demographic satellite helper hash is required}"
: "${E5F_PF_PERSON_EXPECTED_TEST_SHA256:?terminal-mapping test hash is required}"

if [[ ! "$E5F_PF_PERSON_RUN_TAG" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "invalid run tag: $E5F_PF_PERSON_RUN_TAG" >&2
    exit 2
fi

CLUSTER_DIR="$(cd "${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_person_demography.py"
PF_DRIVER="$PROJECT_ROOT/code/model/tools/run_e5f_perfect_foresight_transition.py"
SOLVER="$PROJECT_ROOT/code/model/intergen_eqscale_seq_optimized/solver.py"
CONTINUATION="$PROJECT_ROOT/code/model/tools/run_e5f_post2023_no_policy_continuations.py"
PERSON_LAW="$PROJECT_ROOT/code/model/demographic_transition/person_cohort_law.py"
COUPLING="$PROJECT_ROOT/code/model/demographic_transition/household_person_coupling.py"
DEMOGRAPHIC_BUILDER="$PROJECT_ROOT/code/model/tools/build_e5f_coherent_person_cohort_path.py"
SATELLITE_HELPER="$PROJECT_ROOT/code/model/tools/build_e5f_persons_demographic_satellite.py"
TERMINAL_TEST="$PROJECT_ROOT/code/model/tools/test_run_e5f_perfect_foresight_person_demography.py"
SOURCE_DIR="$PROJECT_ROOT/output/model/e5f_persons_demographic_satellite_20260826b/source_data"
HEADSHIP_DIR="$PROJECT_ROOT/output/model/e5f_headship_demographic_bridge_20260826a"
OUTDIR="$PROJECT_ROOT/output/model/e5f_pf_person_demography_${E5F_PF_PERSON_RUN_TAG}_smoke"
mkdir -p "$CLUSTER_DIR/logs"

hash_file() { sha256sum "$1" | awk '{print $1}'; }
check_hash() {
    local path="$1" expected="$2" label="$3"
    if [ "$(hash_file "$path")" != "$expected" ]; then
        echo "$label hash mismatch" >&2
        exit 2
    fi
}
check_hash "$DRIVER" "$E5F_PF_PERSON_EXPECTED_DRIVER_SHA256" "person-demography driver"
check_hash "$PF_DRIVER" "$E5F_PF_PERSON_EXPECTED_PF_DRIVER_SHA256" "PF driver"
check_hash "$SOLVER" "$E5F_PF_PERSON_EXPECTED_SOLVER_SHA256" "Bellman solver"
check_hash "$CONTINUATION" "$E5F_PF_PERSON_EXPECTED_CONTINUATION_SHA256" "continuation"
check_hash "$PERSON_LAW" "$E5F_PF_PERSON_EXPECTED_PERSON_LAW_SHA256" "person law"
check_hash "$COUPLING" "$E5F_PF_PERSON_EXPECTED_COUPLING_SHA256" "household-person coupling"
check_hash "$DEMOGRAPHIC_BUILDER" "$E5F_PF_PERSON_EXPECTED_DEMOGRAPHIC_BUILDER_SHA256" "demographic builder"
check_hash "$SATELLITE_HELPER" "$E5F_PF_PERSON_EXPECTED_SATELLITE_HELPER_SHA256" "demographic satellite helper"
check_hash "$TERMINAL_TEST" "$E5F_PF_PERSON_EXPECTED_TEST_SHA256" "terminal-mapping test"

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
    "run_tag": os.environ["E5F_PF_PERSON_RUN_TAG"],
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "purpose": "one-period coherent-person-law accounting smoke, not equilibrium",
    "hashes": {
        key: value
        for key, value in os.environ.items()
        if key.startswith("E5F_PF_PERSON_EXPECTED_")
    },
}
with open(path + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(path + ".tmp", path)
PY

PYTHONPATH="$PROJECT_ROOT/code/model" "$PYTHON_BIN" -m unittest discover \
    -s "$PROJECT_ROOT/code/model/demographic_transition/tests" \
    -p 'test_*.py' -v >"$OUTDIR/unit_tests.log" 2>&1
PYTHONPATH="$PROJECT_ROOT/code/model:$PROJECT_ROOT/code/model/tools" \
    "$PYTHON_BIN" "$TERMINAL_TEST" -v \
    >"$OUTDIR/terminal_mapping_test.log" 2>&1

"$PYTHON_BIN" "$DRIVER" \
    --source-dir "$SOURCE_DIR" \
    --headship-dir "$HEADSHIP_DIR" \
    --output-dir "$OUTDIR/accounting" \
    --periods 1 >"$OUTDIR/driver.log" 2>&1

"$PYTHON_BIN" - "$OUTDIR/summary.json" "$OUTDIR/accounting/summary.json" <<'PY'
import json, os, sys
output, accounting = sys.argv[1:]
with open(accounting, encoding="utf-8") as handle:
    result = json.load(handle)
if result.get("status") != "passed_accounting_smoke_not_equilibrium":
    raise SystemExit("accounting smoke did not pass")
payload = {
    "status": "complete_accounting_smoke_not_equilibrium",
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "accounting_summary": result,
}
with open(output + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(output + ".tmp", output)
PY
