#!/bin/bash
# Focused repair-grid audit for the intergen mechanism near-misses.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26_20260617_fast/code/cluster
#
# Example:
#   INTERGEN_MECH_RUN_TAG=intergen_mechanism_repair_grid_20260622 \
#   sbatch --array=0-5%6 submit_intergen_mechanism_repair_grid.sh

#SBATCH --job-name=ihfrepair
#SBATCH --output=logs/slurm_ihf_repair_%A_%a.out
#SBATCH --error=logs/slurm_ihf_repair_%A_%a.err
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

RUN_TAG="${INTERGEN_MECH_RUN_TAG:-intergen_mechanism_repair_grid_20260622}"
ROOT_OUTDIR="${INTERGEN_MECH_OUTDIR:-${PROJECT_DIR}/output/model/${RUN_TAG}}"
POINT_DIR="${INTERGEN_MECH_POINT_DIR:-${SCRIPT_DIR}/audit_points_20260621_frontier}"
CASE_SLICES="${INTERGEN_MECH_CASE_SLICES:-2}"

POINTS=(
    "young_gap_no_old_young_old_own_w3"
    "young_old_no_gap_young_old_own_w3"
    "soft_joint_old_retention_w3"
)

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
TOTAL_TASKS="$((${#POINTS[@]} * CASE_SLICES))"
if [ "${TASK_ID}" -lt 0 ] || [ "${TASK_ID}" -ge "${TOTAL_TASKS}" ]; then
    echo "Invalid SLURM_ARRAY_TASK_ID=${TASK_ID}; expected 0..$((TOTAL_TASKS - 1))" >&2
    exit 2
fi
POINT_INDEX="$((TASK_ID / CASE_SLICES))"
SLICE_INDEX="$((TASK_ID % CASE_SLICES))"
POINT_LABEL="${POINTS[$POINT_INDEX]}"
POINT_JSON="${POINT_DIR}/${POINT_LABEL}.json"
OUTDIR="${ROOT_OUTDIR}/task_${TASK_ID}_${POINT_LABEL}_slice${SLICE_INDEX}"

mkdir -p "${NUMBA_CACHE_DIR}" "${OUTDIR}"

echo "============================================"
echo "Intergen focused mechanism repair grid"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array task: ${TASK_ID}"
echo "Point index: ${POINT_INDEX}"
echo "Point: ${POINT_LABEL}"
echo "Case slice: ${SLICE_INDEX}/${CASE_SLICES}"
echo "Working dir: $(pwd)"
echo "Project dir: ${PROJECT_DIR}"
echo "Model dir: ${MODEL_DIR}"
echo "Python: ${PYTHON_BIN}"
echo "Outdir: ${OUTDIR}"
echo "Point JSON: ${POINT_JSON}"
echo "Started: $(date)"
echo "============================================"

cd "${PROJECT_DIR}"
"${PYTHON_BIN}" code/model/tools/audit_intergen_mechanism_grid.py \
    --candidate-mode repair \
    --point-json "${POINT_LABEL}=${POINT_JSON}" \
    --outdir "${OUTDIR}" \
    --J "${INTERGEN_J:-17}" \
    --Nb "${INTERGEN_NB:-60}" \
    --income-states "${INTERGEN_INCOME_STATES:-5}" \
    --n-house "${INTERGEN_N_HOUSE:-5}" \
    --max-iter-eq "${INTERGEN_MAX_ITER_EQ:-3}" \
    --case-slices "${CASE_SLICES}" \
    --case-slice-index "${SLICE_INDEX}" \
    --max-cases "${INTERGEN_MECH_MAX_CASES:-0}"

echo "Finished: $(date)"
echo "Done: ${OUTDIR}/${POINT_LABEL}"
