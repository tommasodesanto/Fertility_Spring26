#!/bin/bash
# Diagnostic calibration for code/model/intergen_housing_fertility on Torch.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26/code/cluster
#
# Smoke:
#   INTERGEN_RUN_TAG=intergen_smoke_YYYYMMDD INTERGEN_CASES_PER_TASK=2 \
#   INTERGEN_J=12 INTERGEN_NB=30 INTERGEN_MAX_ITER_EQ=12 \
#   sbatch --array=1-2%2 submit_intergen_housing_fertility_calibration.sh
#
# Moderate pass:
#   INTERGEN_RUN_TAG=intergen_calib_YYYYMMDD INTERGEN_CASES_PER_TASK=24 \
#   INTERGEN_J=16 INTERGEN_NB=60 INTERGEN_N_HOUSE=6 INTERGEN_MAX_ITER_EQ=45 \
#   sbatch --array=1-24%24 submit_intergen_housing_fertility_calibration.sh
#
# Old-target, no-location overnight diagnostic pass:
#   INTERGEN_RUN_TAG=intergen_old_nonlocation_YYYYMMDD INTERGEN_TARGET_SET=old_nonlocation \
#   INTERGEN_CASES_PER_TASK=48 INTERGEN_J=16 INTERGEN_NB=70 \
#   INTERGEN_N_HOUSE=6 INTERGEN_MAX_ITER_EQ=60 \
#   sbatch --array=1-48%48 submit_intergen_housing_fertility_calibration.sh

#SBATCH --job-name=ihfcal
#SBATCH --output=logs/slurm_ihf_cal_%A_%a.out
#SBATCH --error=logs/slurm_ihf_cal_%A_%a.err
#SBATCH --array=1-24%24
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"

mkdir -p "${SCRIPT_DIR}/logs"

export INTERGEN_RUN_TAG="${INTERGEN_RUN_TAG:-intergen_calib_$(date +%Y%m%d_%H%M%S)}"
export INTERGEN_CASES_PER_TASK="${INTERGEN_CASES_PER_TASK:-24}"
export INTERGEN_SEED_BASE="${INTERGEN_SEED_BASE:-260604}"
export INTERGEN_J="${INTERGEN_J:-16}"
export INTERGEN_NB="${INTERGEN_NB:-60}"
export INTERGEN_N_HOUSE="${INTERGEN_N_HOUSE:-6}"
export INTERGEN_MAX_ITER_EQ="${INTERGEN_MAX_ITER_EQ:-45}"
export INTERGEN_TARGET_SET="${INTERGEN_TARGET_SET:-old_nonlocation}"
export INTERGEN_RESULTS_DIR="${INTERGEN_RESULTS_DIR:-${SCRIPT_DIR}/results_intergen_housing_fertility_${INTERGEN_RUN_TAG}}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
SEED=$((INTERGEN_SEED_BASE + TASK_ID))
TASK_OUTDIR="${INTERGEN_RESULTS_DIR}/task_${TASK_ID}"

if command -v module >/dev/null 2>&1; then
    module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || module load python/3.12 2>/dev/null || module load python/3.11 2>/dev/null || true
fi

PYTHON_BIN="${INTERGEN_PYTHON:-$(command -v python3 || command -v python)}"

export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}/intergen_housing_fertility}"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mkdir -p "${NUMBA_CACHE_DIR}" "${TASK_OUTDIR}"

echo "============================================"
echo "Intergenerational housing-fertility diagnostic calibration"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array Task ID: ${TASK_ID}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Working dir: $(pwd)"
echo "Model code: ${MODEL_DIR}"
echo "Python: ${PYTHON_BIN}"
echo "Run tag: ${INTERGEN_RUN_TAG}"
echo "Task outdir: ${TASK_OUTDIR}"
echo "cases=${INTERGEN_CASES_PER_TASK} seed=${SEED} J=${INTERGEN_J} Nb=${INTERGEN_NB} n_house=${INTERGEN_N_HOUSE} max_iter_eq=${INTERGEN_MAX_ITER_EQ} target_set=${INTERGEN_TARGET_SET}"
echo "============================================"

cd "${MODEL_DIR}"
"${PYTHON_BIN}" -m intergen_housing_fertility.cli calibrate-small \
    --cases "${INTERGEN_CASES_PER_TASK}" \
    --seed "${SEED}" \
    --J "${INTERGEN_J}" \
    --Nb "${INTERGEN_NB}" \
    --n-house "${INTERGEN_N_HOUSE}" \
    --max-iter-eq "${INTERGEN_MAX_ITER_EQ}" \
    --target-set "${INTERGEN_TARGET_SET}" \
    --outdir "${TASK_OUTDIR}"

echo "Done: ${TASK_OUTDIR}"
