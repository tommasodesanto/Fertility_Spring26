#!/bin/bash
# Four-group no-date-specific-reweighting coordinate panel.
# One task independently normalizes the old steady state and solves the exact
# five-date household/housing transition for one predeclared candidate.
#SBATCH --job-name=e5fh4c
#SBATCH --output=logs/slurm_e5fh4c_%A_%a.out
#SBATCH --error=logs/slurm_e5fh4c_%A_%a.err
#SBATCH --partition=cpu_short
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general
#SBATCH --array=1-23%12

set -euo pipefail

SUBMIT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
if [ -d "$SUBMIT_DIR/../model" ] && [ -d "$SUBMIT_DIR/../../output" ]; then
    CLUSTER_DIR="$SUBMIT_DIR"
elif [ -d "$SUBMIT_DIR/code/model" ] && [ -d "$SUBMIT_DIR/code/cluster" ]; then
    CLUSTER_DIR="$SUBMIT_DIR/code/cluster"
else
    echo "cannot resolve project checkout from SLURM_SUBMIT_DIR=$SUBMIT_DIR" >&2
    exit 2
fi
MODEL_DIR="$(cd "$CLUSTER_DIR/../model" && pwd)"
PROJECT_ROOT="$(cd "$CLUSTER_DIR/../.." && pwd)"
mkdir -p "$CLUSTER_DIR/logs"
module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || true
PYTHON_BIN="$(command -v python3 || command -v python)"

RUN_TAG="${E5F_HIST4_RUN_TAG:?E5F_HIST4_RUN_TAG is required}"
PLAN_SHA="${E5F_HIST4_PLAN_SHA256:?E5F_HIST4_PLAN_SHA256 is required}"
SMOKE_SHA="${E5F_HIST4_SMOKE_SUMMARY_SHA256:?E5F_HIST4_SMOKE_SUMMARY_SHA256 is required}"
DRIVER_SHA="${E5F_HIST4_DRIVER_SHA256:?E5F_HIST4_DRIVER_SHA256 is required}"
SMOKE_HELPER_SHA="${E5F_HIST4_SMOKE_HELPER_SHA256:?E5F_HIST4_SMOKE_HELPER_SHA256 is required}"
DESIGN_HELPER_SHA="${E5F_HIST4_DESIGN_HELPER_SHA256:?E5F_HIST4_DESIGN_HELPER_SHA256 is required}"
AGE_PATH_SHA="${E5F_HIST4_AGE_PATH_SHA256:?E5F_HIST4_AGE_PATH_SHA256 is required}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
PLAN_DIR="$PROJECT_ROOT/output/model/e5f_historical_four_group_coordinate_plan"
PLAN="$PLAN_DIR/coordinate_plan.json"
SMOKE="$PROJECT_ROOT/output/model/e5f_historical_four_group_coordinate_torch_smoke/summary.json"
DRIVER="$MODEL_DIR/tools/run_e5f_historical_demographic_closure_candidate.py"
SMOKE_HELPER="$MODEL_DIR/tools/run_e5f_historical_demographic_closure_smoke.py"
DESIGN_HELPER="$MODEL_DIR/tools/build_e5f_transition_design_feasibility.py"
AGE_PATH="$PROJECT_ROOT/code/data/Spatial_aggregate_withmicrodata/output/national_householder_age_path/national_householder_age_path.csv"
SOURCE="$PROJECT_ROOT/output/model/intergen_e5f_child_room_floor_psinneg_extended_20260806/report/results.json"
REPORT="$PROJECT_ROOT/output/model/e5f_transition_ridge_refinement_jump11_polish_r2_20260818/report/summary.json"
CANDIDATE="$PLAN_DIR/candidate_$(printf '%03d' "$TASK_ID").json"
OUTDIR="$PROJECT_ROOT/output/model/e5f_historical_four_group_coordinate_${RUN_TAG}/task_$(printf '%03d' "$TASK_ID")"

for specification in \
    "$PLAN:$PLAN_SHA" \
    "$SMOKE:$SMOKE_SHA" \
    "$DRIVER:$DRIVER_SHA" \
    "$SMOKE_HELPER:$SMOKE_HELPER_SHA" \
    "$DESIGN_HELPER:$DESIGN_HELPER_SHA" \
    "$AGE_PATH:$AGE_PATH_SHA"; do
    path="${specification%%:*}"
    expected="${specification##*:}"
    actual="$(sha256sum "$path" | awk '{print $1}')"
    if [ "$actual" != "$expected" ]; then
        echo "hash gate failed for $path: actual=$actual expected=$expected" >&2
        exit 2
    fi
done
if [ -d "$OUTDIR" ] && [ -n "$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -print -quit)" ]; then
    echo "refusing to overwrite nonempty task output: $OUTDIR" >&2
    exit 2
fi

"$PYTHON_BIN" - "$PLAN" "$CANDIDATE" "$TASK_ID" "$SMOKE" <<'PY'
import hashlib,json,math,sys
plan_path,candidate_path,task_id,smoke_path=sys.argv[1:]
plan=json.load(open(plan_path)); smoke=json.load(open(smoke_path)); task=int(task_id)
assert plan['status']=='complete' and plan['task_count']==23
row=plan['tasks'][task-1]
assert row['task_id']==task
actual=hashlib.sha256(open(candidate_path,'rb').read()).hexdigest()
assert actual==row['candidate_sha256']
assert smoke['status']=='complete_historical_four_group_candidate'
smoke_id=smoke['candidate']['candidate_id']
smoke_task=next(row for row in plan['tasks'] if row['candidate_id']==smoke_id)
assert smoke['candidate_sha256']==smoke_task['candidate_sha256']
assert math.isfinite(smoke['transition_loss_at_selected_parameters'])
PY

export PYTHONPATH="$MODEL_DIR:$MODEL_DIR/tools:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
exec "$PYTHON_BIN" "$DRIVER" \
    --candidate-json "$CANDIDATE" \
    --selected-report "$REPORT" \
    --source "$SOURCE" \
    --output-dir "$OUTDIR"
