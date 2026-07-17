#!/bin/bash
# Two-hour global differential-evolution panel for code/model/intergen_housing_fertility.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26/code/cluster
#
# Example:
#   INTERGEN_RUN_TAG=intergen_candidate_no_timing_v0_globalde_3g_20260609 \
#   sbatch --array=1-64%64 submit_intergen_housing_fertility_global_de.sh

#SBATCH --job-name=ihfde
#SBATCH --output=logs/slurm_ihf_de_%A_%a.out
#SBATCH --error=logs/slurm_ihf_de_%A_%a.err
#SBATCH --array=1-8%8
#SBATCH --partition=cpu_short
#SBATCH --time=02:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"

mkdir -p "${SCRIPT_DIR}/logs"

export INTERGEN_RUN_TAG="${INTERGEN_RUN_TAG:-intergen_july9_repair_globalde_$(date +%Y%m%d_%H%M%S)}"
export INTERGEN_PROFILE="${INTERGEN_PROFILE:-intergen_july9_repair_v1}"
export INTERGEN_TARGET_SET="${INTERGEN_TARGET_SET:-candidate_replacement_post_audit_v1}"
export INTERGEN_GLOBAL_EVALS_PER_TASK="${INTERGEN_GLOBAL_EVALS_PER_TASK:-300}"
export INTERGEN_MINUTES="${INTERGEN_MINUTES:-115}"
export INTERGEN_GLOBAL_POP_SIZE="${INTERGEN_GLOBAL_POP_SIZE:-22}"
export INTERGEN_GLOBAL_MUTATION="${INTERGEN_GLOBAL_MUTATION:-0.85}"
export INTERGEN_GLOBAL_CROSSOVER="${INTERGEN_GLOBAL_CROSSOVER:-0.70}"
export INTERGEN_SEED_BASE="${INTERGEN_SEED_BASE:-2026060900}"
export INTERGEN_SEED_THETA_JSON="${INTERGEN_SEED_THETA_JSON:-}"
export INTERGEN_PREVIOUS_RESULTS_DIR="${INTERGEN_PREVIOUS_RESULTS_DIR:-}"
export INTERGEN_J="${INTERGEN_J:-17}"
export INTERGEN_NB="${INTERGEN_NB:-120}"
export INTERGEN_INCOME_STATES="${INTERGEN_INCOME_STATES:-5}"
export INTERGEN_INCOME_PROCESS="${INTERGEN_INCOME_PROCESS:-current}"
export INTERGEN_INCOME_ANNUAL_RHO="${INTERGEN_INCOME_ANNUAL_RHO:-0.90}"
export INTERGEN_INCOME_INNOVATION_SD="${INTERGEN_INCOME_INNOVATION_SD:-0.20}"
export INTERGEN_N_HOUSE="${INTERGEN_N_HOUSE:-5}"
export INTERGEN_MAX_ITER_EQ="${INTERGEN_MAX_ITER_EQ:-10}"
export INTERGEN_RESULTS_DIR="${INTERGEN_RESULTS_DIR:-${SCRIPT_DIR}/results_intergen_housing_fertility_${INTERGEN_RUN_TAG}}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
SEED=$((INTERGEN_SEED_BASE + 1000 * TASK_ID))
TASK_OUTDIR="${INTERGEN_RESULTS_DIR}/task_${TASK_ID}"

if [ -n "${INTERGEN_PREVIOUS_RESULTS_DIR}" ]; then
    PREVIOUS_BEST="${INTERGEN_PREVIOUS_RESULTS_DIR}/task_${TASK_ID}/best.json"
    if [ ! -s "${PREVIOUS_BEST}" ]; then
        echo "Missing previous-wave strict best: ${PREVIOUS_BEST}" >&2
        exit 2
    fi
    INTERGEN_SEED_THETA_JSON="${PREVIOUS_BEST}"
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
echo "Intergenerational housing-fertility global DE panel"
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
echo "max_evals=${INTERGEN_GLOBAL_EVALS_PER_TASK} minutes=${INTERGEN_MINUTES} pop_size=${INTERGEN_GLOBAL_POP_SIZE}"
echo "mutation=${INTERGEN_GLOBAL_MUTATION} crossover=${INTERGEN_GLOBAL_CROSSOVER}"
echo "seed=${SEED} J=${INTERGEN_J} Nb=${INTERGEN_NB} income_states=${INTERGEN_INCOME_STATES} income_process=${INTERGEN_INCOME_PROCESS} income_annual_rho=${INTERGEN_INCOME_ANNUAL_RHO} income_innovation_sd=${INTERGEN_INCOME_INNOVATION_SD} n_house=${INTERGEN_N_HOUSE} max_iter_eq=${INTERGEN_MAX_ITER_EQ}"
echo "seed_theta_json=${INTERGEN_SEED_THETA_JSON:-none}"
echo "Started: $(date)"
echo "============================================"

cd "${MODEL_DIR}"
EXTRA_ARGS=()
if [ -n "${INTERGEN_SEED_THETA_JSON}" ]; then
    EXTRA_ARGS+=(--seed-theta-json "${INTERGEN_SEED_THETA_JSON}")
fi
"${PYTHON_BIN}" -m intergen_housing_fertility.cli global-de-panel \
    --max-evals "${INTERGEN_GLOBAL_EVALS_PER_TASK}" \
    --seed "${SEED}" \
    --J "${INTERGEN_J}" \
    --Nb "${INTERGEN_NB}" \
    --income-states "${INTERGEN_INCOME_STATES}" \
    --income-process "${INTERGEN_INCOME_PROCESS}" \
    --income-annual-rho "${INTERGEN_INCOME_ANNUAL_RHO}" \
    --income-innovation-sd "${INTERGEN_INCOME_INNOVATION_SD}" \
    --n-house "${INTERGEN_N_HOUSE}" \
    --max-iter-eq "${INTERGEN_MAX_ITER_EQ}" \
    --minutes "${INTERGEN_MINUTES}" \
    --pop-size "${INTERGEN_GLOBAL_POP_SIZE}" \
    --mutation "${INTERGEN_GLOBAL_MUTATION}" \
    --crossover "${INTERGEN_GLOBAL_CROSSOVER}" \
    --target-set "${INTERGEN_TARGET_SET}" \
    --profile "${INTERGEN_PROFILE}" \
    --outdir "${TASK_OUTDIR}" \
    "${EXTRA_ARGS[@]}"

echo "Finished: $(date)"
echo "Done: ${TASK_OUTDIR}"
