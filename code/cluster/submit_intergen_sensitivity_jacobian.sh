#!/bin/bash
# Finite-difference identification audit for code/model/intergen_housing_fertility.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26_20260617_fast/code/cluster
#
# Example:
#   sbatch submit_intergen_sensitivity_jacobian.sh

#SBATCH --job-name=ihfsens
#SBATCH --output=logs/slurm_ihf_sens_%j.out
#SBATCH --error=logs/slurm_ihf_sens_%j.err
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MODEL_DIR="${PROJECT_DIR}/code/model"

mkdir -p "${SCRIPT_DIR}/logs"

if command -v module >/dev/null 2>&1; then
    module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || module load python/3.12 2>/dev/null || module load python/3.11 2>/dev/null || true
fi

if [ -n "${INTERGEN_PYTHON:-}" ]; then
    PYTHON_BIN="${INTERGEN_PYTHON}"
elif [ -x "/share/apps/anaconda3/2025.06/bin/python3" ]; then
    PYTHON_BIN="/share/apps/anaconda3/2025.06/bin/python3"
else
    PYTHON_BIN="$(command -v python3 || command -v python)"
fi

export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}/intergen_housing_fertility}"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

RUN_TAG="${INTERGEN_SENS_RUN_TAG:-intergen_sensitivity_jacobian_20260618}"
OUTDIR="${INTERGEN_SENS_OUTDIR:-${PROJECT_DIR}/output/model/${RUN_TAG}}"
FRONTIER_ROOT="${INTERGEN_FRONTIER_ROOT:-${PROJECT_DIR}/output/model/cluster_pulls/intergen_overnight_frontier_20260617}"
TARGET_SET="${INTERGEN_SENS_TARGET_SET:-candidate_no_timing_v0}"

POINT_ARGS=()
if [ -n "${INTERGEN_SENS_POINT_JSONS:-}" ]; then
    IFS=',' read -r -a POINT_JSON_ARRAY <<< "${INTERGEN_SENS_POINT_JSONS}"
    for point_json in "${POINT_JSON_ARRAY[@]}"; do
        POINT_ARGS+=(--point-json "${point_json}")
    done
elif [ -n "${INTERGEN_SENS_POINT_LABELS:-}" ]; then
    IFS=',' read -r -a POINT_LABEL_ARRAY <<< "${INTERGEN_SENS_POINT_LABELS}"
    for point_label in "${POINT_LABEL_ARRAY[@]}"; do
        POINT_ARGS+=(--point-label "${point_label}")
    done
else
    POINT_ARGS+=(--point-label core_feasibility_v1 --point-label roomcost_test_v1)
fi

mkdir -p "${NUMBA_CACHE_DIR}" "${OUTDIR}"

echo "============================================"
echo "Intergen one-market sensitivity Jacobian audit"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Working dir: $(pwd)"
echo "Project dir: ${PROJECT_DIR}"
echo "Model dir: ${MODEL_DIR}"
echo "Python: ${PYTHON_BIN}"
echo "Outdir: ${OUTDIR}"
echo "Frontier root: ${FRONTIER_ROOT}"
echo "Target set: ${TARGET_SET}"
echo "Point args: ${POINT_ARGS[*]}"
echo "Started: $(date)"
echo "============================================"

cd "${PROJECT_DIR}"
"${PYTHON_BIN}" code/model/tools/audit_intergen_sensitivity_jacobian.py \
    --frontier-root "${FRONTIER_ROOT}" \
    --outdir "${OUTDIR}" \
    "${POINT_ARGS[@]}" \
    --target-set "${TARGET_SET}" \
    --rel-step "${INTERGEN_SENS_REL_STEP:-0.01}" \
    --J "${INTERGEN_J:-17}" \
    --Nb "${INTERGEN_NB:-60}" \
    --income-states "${INTERGEN_INCOME_STATES:-5}" \
    --n-house "${INTERGEN_N_HOUSE:-5}" \
    --max-iter-eq "${INTERGEN_MAX_ITER_EQ:-3}"

echo "Finished: $(date)"
echo "Done: ${OUTDIR}"
