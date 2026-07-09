#!/bin/bash
# True derivative-free local-polish panel for code/model/intergen_housing_fertility.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26/code/cluster
#
# Example:
#   INTERGEN_RUN_TAG=intergen_age18_local_polish_20260628 \
#   INTERGEN_SEED_THETA_JSON=/path/to/best.json \
#   sbatch --array=1-12%12 submit_intergen_housing_fertility_local_polish.sh

#SBATCH --job-name=ihfpol
#SBATCH --output=logs/slurm_ihf_polish_%A_%a.out
#SBATCH --error=logs/slurm_ihf_polish_%A_%a.err
#SBATCH --array=1-8%8
#SBATCH --partition=cpu_short
#SBATCH --time=02:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"

mkdir -p "${SCRIPT_DIR}/logs"

export INTERGEN_RUN_TAG="${INTERGEN_RUN_TAG:-intergen_local_polish_$(date +%Y%m%d_%H%M%S)}"
export INTERGEN_PROFILE="${INTERGEN_PROFILE:-intergen_july9_repair_v1}"
export INTERGEN_TARGET_SET="${INTERGEN_TARGET_SET:-candidate_replacement_post_audit_v1}"
export INTERGEN_LOCAL_MAX_EVALS="${INTERGEN_LOCAL_MAX_EVALS:-400}"
export INTERGEN_MINUTES="${INTERGEN_MINUTES:-115}"
export INTERGEN_SEED_BASE="${INTERGEN_SEED_BASE:-2026062800}"
export INTERGEN_SEED_THETA_JSON="${INTERGEN_SEED_THETA_JSON:-}"
export INTERGEN_J="${INTERGEN_J:-17}"
export INTERGEN_NB="${INTERGEN_NB:-120}"
export INTERGEN_INCOME_STATES="${INTERGEN_INCOME_STATES:-5}"
export INTERGEN_N_HOUSE="${INTERGEN_N_HOUSE:-5}"
export INTERGEN_MAX_ITER_EQ="${INTERGEN_MAX_ITER_EQ:-25}"
export INTERGEN_LOCAL_MIN_STEP="${INTERGEN_LOCAL_MIN_STEP:-0.003}"
export INTERGEN_LOCAL_SHRINK="${INTERGEN_LOCAL_SHRINK:-0.5}"
export INTERGEN_FIXED_BETA_ANNUAL="${INTERGEN_FIXED_BETA_ANNUAL:-}"
export INTERGEN_FIXED_THETA_N="${INTERGEN_FIXED_THETA_N:-}"
export INTERGEN_FIXED_CHI="${INTERGEN_FIXED_CHI:-}"
export INTERGEN_REPAIR_DESIGN="${INTERGEN_REPAIR_DESIGN:-0}"
export INTERGEN_RESULTS_DIR="${INTERGEN_RESULTS_DIR:-${SCRIPT_DIR}/results_intergen_housing_fertility_${INTERGEN_RUN_TAG}}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
SEED=$((INTERGEN_SEED_BASE + 1000 * TASK_ID))

ARM_LABEL="free_chi"
if [ "${INTERGEN_REPAIR_DESIGN}" = "1" ]; then
    INTERGEN_FIXED_BETA_ANNUAL="${INTERGEN_FIXED_BETA_ANNUAL:-0.94}"
    INTERGEN_FIXED_THETA_N="${INTERGEN_FIXED_THETA_N:-0.9811}"
    if [ "${TASK_ID}" -ge 3 ] && [ -z "${INTERGEN_FIXED_CHI}" ]; then
        CHI_RUNGS=(0.40 0.55 0.70 0.85 1.00 1.15)
        CHI_INDEX=$((TASK_ID - 3))
        if [ "${CHI_INDEX}" -ge "${#CHI_RUNGS[@]}" ]; then
            echo "Repair design supports array tasks 1-8 only." >&2
            exit 2
        fi
        INTERGEN_FIXED_CHI="${CHI_RUNGS[$CHI_INDEX]}"
    fi
    if [ -n "${INTERGEN_FIXED_CHI}" ]; then
        ARM_LABEL="fixed_chi_${INTERGEN_FIXED_CHI}"
    fi
fi
export INTERGEN_FIXED_BETA_ANNUAL INTERGEN_FIXED_THETA_N INTERGEN_FIXED_CHI
TASK_OUTDIR="${INTERGEN_RESULTS_DIR}/task_${TASK_ID}_${ARM_LABEL}"

if [ -z "${INTERGEN_SEED_THETA_JSON}" ]; then
    echo "INTERGEN_SEED_THETA_JSON is required for local polish." >&2
    exit 2
fi

STEP_OPTIONS=(0.025 0.040 0.060 0.085 0.115 0.150)
STEP_INDEX=$(( (TASK_ID - 1) % ${#STEP_OPTIONS[@]} ))
DEFAULT_STEP="${STEP_OPTIONS[$STEP_INDEX]}"
if [ -n "${INTERGEN_LOCAL_INITIAL_STEP:-}" ]; then
    LOCAL_STEP="${INTERGEN_LOCAL_INITIAL_STEP}"
else
    LOCAL_STEP="${DEFAULT_STEP}"
fi

if [ -n "${INTERGEN_LOCAL_METHOD:-}" ]; then
    LOCAL_METHOD="${INTERGEN_LOCAL_METHOD}"
else
    if [ $(( (TASK_ID - 1) / ${#STEP_OPTIONS[@]} % 2 )) -eq 0 ]; then
        LOCAL_METHOD="nelder-mead"
    else
        LOCAL_METHOD="pattern"
    fi
fi

if command -v module >/dev/null 2>&1; then
    module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || module load python/3.12 2>/dev/null || module load python/3.11 2>/dev/null || true
fi

if [ -n "${INTERGEN_PYTHON:-}" ]; then
    PYTHON_BIN="${INTERGEN_PYTHON}"
elif [ "$(uname)" = "Darwin" ] && [ -x "${MODEL_DIR}/.venv/bin/python" ]; then
    PYTHON_BIN="${MODEL_DIR}/.venv/bin/python"
else
    PYTHON_BIN="$(command -v python3 || command -v python)"
fi

export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}/intergen_housing_fertility}"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mkdir -p "${NUMBA_CACHE_DIR}" "${TASK_OUTDIR}"

echo "============================================"
echo "Intergenerational housing-fertility local polish"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array Task ID: ${TASK_ID}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Working dir: $(pwd)"
echo "Model code: ${MODEL_DIR}"
echo "Python: ${PYTHON_BIN}"
echo "Run tag: ${INTERGEN_RUN_TAG}"
echo "Task outdir: ${TASK_OUTDIR}"
echo "Target set: ${INTERGEN_TARGET_SET}"
echo "Production profile: ${INTERGEN_PROFILE}"
echo "Repair design: ${INTERGEN_REPAIR_DESIGN} arm=${ARM_LABEL}"
echo "method=${LOCAL_METHOD} initial_step=${LOCAL_STEP} min_step=${INTERGEN_LOCAL_MIN_STEP} shrink=${INTERGEN_LOCAL_SHRINK}"
echo "max_evals=${INTERGEN_LOCAL_MAX_EVALS} minutes=${INTERGEN_MINUTES}"
echo "seed=${SEED} J=${INTERGEN_J} Nb=${INTERGEN_NB} income_states=${INTERGEN_INCOME_STATES} n_house=${INTERGEN_N_HOUSE} max_iter_eq=${INTERGEN_MAX_ITER_EQ}"
echo "seed_theta_json=${INTERGEN_SEED_THETA_JSON}"
echo "fixed_beta_annual=${INTERGEN_FIXED_BETA_ANNUAL:-free} fixed_theta_n=${INTERGEN_FIXED_THETA_N:-free} fixed_chi=${INTERGEN_FIXED_CHI:-free}"
echo "Started: $(date)"
echo "============================================"

cd "${MODEL_DIR}"
FIXED_ARGS=()
if [ -n "${INTERGEN_FIXED_BETA_ANNUAL}" ]; then
    FIXED_ARGS+=(--fixed-beta-annual "${INTERGEN_FIXED_BETA_ANNUAL}")
fi
if [ -n "${INTERGEN_FIXED_THETA_N}" ]; then
    FIXED_ARGS+=(--fixed-theta-n "${INTERGEN_FIXED_THETA_N}")
fi
if [ -n "${INTERGEN_FIXED_CHI}" ]; then
    FIXED_ARGS+=(--fixed-chi "${INTERGEN_FIXED_CHI}")
fi
"${PYTHON_BIN}" -m intergen_housing_fertility.cli local-polish \
    --max-evals "${INTERGEN_LOCAL_MAX_EVALS}" \
    --seed "${SEED}" \
    --J "${INTERGEN_J}" \
    --Nb "${INTERGEN_NB}" \
    --income-states "${INTERGEN_INCOME_STATES}" \
    --n-house "${INTERGEN_N_HOUSE}" \
    --max-iter-eq "${INTERGEN_MAX_ITER_EQ}" \
    --minutes "${INTERGEN_MINUTES}" \
    --target-set "${INTERGEN_TARGET_SET}" \
    --profile "${INTERGEN_PROFILE}" \
    --seed-theta-json "${INTERGEN_SEED_THETA_JSON}" \
    --method "${LOCAL_METHOD}" \
    --initial-step "${LOCAL_STEP}" \
    --min-step "${INTERGEN_LOCAL_MIN_STEP}" \
    --shrink "${INTERGEN_LOCAL_SHRINK}" \
    "${FIXED_ARGS[@]}" \
    --outdir "${TASK_OUTDIR}"

echo "Finished: $(date)"
echo "Done: ${TASK_OUTDIR}"
