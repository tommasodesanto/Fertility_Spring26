#!/bin/bash
# Direct no-inversion renewal-valve Python calibration on Torch.
# Submit from:
#   $SCRATCH/projects/Fertility_Spring26/code/cluster
#
# Pilot:
#   DT_DIRECT_RUN_TAG=py_direct_open_pilot DT_DIRECT_SETUP=fast DT_DIRECT_MAX_ITER_EQ=4 \
#   DT_DIRECT_MAX_EVALS=2 DT_DIRECT_BUDGET_SEC=1800 sbatch --array=1-2%2 submit_python_direct_geometry_overnight.sh
#
# Short-partition production:
#   DT_DIRECT_RUN_TAG=py_direct_open_wide_YYYYMMDD DT_DIRECT_BOUNDS=wide \
#   sbatch --array=1-40%40 submit_python_direct_geometry_overnight.sh

#SBATCH --job-name=dtpydir
#SBATCH --output=logs/slurm_py_direct_%A_%a.out
#SBATCH --error=logs/slurm_py_direct_%A_%a.err
#SBATCH --array=1-40%40
#SBATCH --partition=cpu_short
#SBATCH --time=03:55:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --account=torch_pr_570_general

set -euo pipefail

SCRIPT_DIR="$(cd "${SLURM_SUBMIT_DIR:-$(pwd)}" && pwd)"
MODEL_DIR="$(cd "${SCRIPT_DIR}/../model" && pwd)"

mkdir -p "${SCRIPT_DIR}/logs"

export DT_DIRECT_RUN_TAG="${DT_DIRECT_RUN_TAG:-py_direct_open_city_$(date +%Y%m%d_%H%M%S)}"
export DT_DIRECT_SETUP="${DT_DIRECT_SETUP:-benchmark}"
export DT_DIRECT_BOUNDS="${DT_DIRECT_BOUNDS:-scout}"
export DT_DIRECT_RESULTS_DIR="${DT_DIRECT_RESULTS_DIR:-${SCRIPT_DIR}/results_python_direct_geometry_${DT_DIRECT_RUN_TAG}}"
export DT_DIRECT_BUDGET_SEC="${DT_DIRECT_BUDGET_SEC:-13500}"
export DT_DIRECT_MAX_EVALS="${DT_DIRECT_MAX_EVALS:-1000000}"
export DT_DIRECT_GEO_WEIGHT="${DT_DIRECT_GEO_WEIGHT:-100}"
export DT_DIRECT_POPULATION_CLOSURE="${DT_DIRECT_POPULATION_CLOSURE:-renewal_valve_calibrated}"
export DT_DIRECT_SCALE_TARGET="${DT_DIRECT_SCALE_TARGET:-1.0}"
export DT_DIRECT_SCALE_WEIGHT="${DT_DIRECT_SCALE_WEIGHT:-100}"
export DT_DIRECT_OUTSIDE_VALUE_X0="${DT_DIRECT_OUTSIDE_VALUE_X0:--41.95}"
export DT_DIRECT_OUTSIDE_FLOW_X0="${DT_DIRECT_OUTSIDE_FLOW_X0:-0.003}"
export DT_DIRECT_RENEWAL_RETENTION="${DT_DIRECT_RENEWAL_RETENTION:-1.0}"
export DT_DIRECT_EQ_PENALTY_WEIGHT="${DT_DIRECT_EQ_PENALTY_WEIGHT:-0}"
export DT_DIRECT_MAX_ITER_EQ="${DT_DIRECT_MAX_ITER_EQ:-0}"
export DT_DIRECT_LOG_EVERY="${DT_DIRECT_LOG_EVERY:-5}"

mkdir -p "${DT_DIRECT_RESULTS_DIR}"

if command -v module >/dev/null 2>&1; then
    module load anaconda3/2025.06 2>/dev/null || module load anaconda3 2>/dev/null || module load python/3.11 2>/dev/null || module load python/3.10 2>/dev/null || true
fi

if [ -n "${DT_DIRECT_PYTHON:-}" ]; then
    PYTHON_BIN="${DT_DIRECT_PYTHON}"
elif [ "$(uname)" = "Darwin" ] && [ -x "${MODEL_DIR}/.venv/bin/python" ]; then
    PYTHON_BIN="${MODEL_DIR}/.venv/bin/python"
else
    PYTHON_BIN="$(command -v python3 || command -v python)"
fi

export PYTHONPATH="${MODEL_DIR}:${PYTHONPATH:-}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${SCRATCH:-/tmp}/fertility_numba_cache/${USER:-user}}"
export NUMBA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mkdir -p "${NUMBA_CACHE_DIR}"

echo "============================================"
echo "Direct no-inversion renewal-valve Python calibration"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-1}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Working dir: $(pwd)"
echo "Model code: ${MODEL_DIR}"
echo "Python: ${PYTHON_BIN}"
echo "Run tag: ${DT_DIRECT_RUN_TAG}"
echo "Results: ${DT_DIRECT_RESULTS_DIR}"
echo "Setup: ${DT_DIRECT_SETUP} | bounds: ${DT_DIRECT_BOUNDS} | closure: ${DT_DIRECT_POPULATION_CLOSURE}"
if [ "${DT_DIRECT_POPULATION_CLOSURE}" = "renewal_valve_calibrated" ]; then
    echo "Scale target: ${DT_DIRECT_SCALE_TARGET} imposed mechanically | scale weight inactive"
else
    echo "Scale target: ${DT_DIRECT_SCALE_TARGET} | scale weight: ${DT_DIRECT_SCALE_WEIGHT}"
fi
echo "Outside value x0: ${DT_DIRECT_OUTSIDE_VALUE_X0} | outside flow x0: ${DT_DIRECT_OUTSIDE_FLOW_X0} | renewal rho: ${DT_DIRECT_RENEWAL_RETENTION}"
echo "Budget: ${DT_DIRECT_BUDGET_SEC}s | max evals: ${DT_DIRECT_MAX_EVALS}"
echo "Started: $(date)"
echo "============================================"

"${PYTHON_BIN}" - <<'PY'
import importlib.util
import sys

missing = [pkg for pkg in ("numpy", "numba") if importlib.util.find_spec(pkg) is None]
if missing:
    raise SystemExit(
        "Missing required Python packages: "
        + ", ".join(missing)
        + ". Load a Python module with NumPy/Numba or set DT_DIRECT_PYTHON to a prepared environment."
    )
print("python", sys.version)
for pkg in ("numpy", "numba"):
    mod = __import__(pkg)
    print(pkg, mod.__version__)
PY

ARGS=(
    "${MODEL_DIR}/tools/calibrate_direct_geometry.py"
    --job-id "${SLURM_ARRAY_TASK_ID:-1}"
    --run-tag "${DT_DIRECT_RUN_TAG}"
    --setup "${DT_DIRECT_SETUP}"
    --bound-profile "${DT_DIRECT_BOUNDS}"
    --results-dir "${DT_DIRECT_RESULTS_DIR}"
    --max-evals "${DT_DIRECT_MAX_EVALS}"
    --budget-sec "${DT_DIRECT_BUDGET_SEC}"
    --geo-weight "${DT_DIRECT_GEO_WEIGHT}"
    --population-closure "${DT_DIRECT_POPULATION_CLOSURE}"
    --scale-target "${DT_DIRECT_SCALE_TARGET}"
    --scale-weight "${DT_DIRECT_SCALE_WEIGHT}"
    --outside-value-x0 "${DT_DIRECT_OUTSIDE_VALUE_X0}"
    --outside-flow-x0 "${DT_DIRECT_OUTSIDE_FLOW_X0}"
    --renewal-retention "${DT_DIRECT_RENEWAL_RETENTION}"
    --eq-penalty-weight "${DT_DIRECT_EQ_PENALTY_WEIGHT}"
    --log-every "${DT_DIRECT_LOG_EVERY}"
    --resume
)

if [ "${DT_DIRECT_MAX_ITER_EQ}" != "0" ]; then
    ARGS+=(--max-iter-eq "${DT_DIRECT_MAX_ITER_EQ}")
fi

"${PYTHON_BIN}" "${ARGS[@]}"

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
