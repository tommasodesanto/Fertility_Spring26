#!/bin/bash
# Direct no-inversion open-city Python calibration on Torch.
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
export DT_DIRECT_POPULATION_CLOSURE="${DT_DIRECT_POPULATION_CLOSURE:-outside_option_benchmark_normalized}"
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
echo "Direct no-inversion Python calibration"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-1}"
echo "Node: ${SLURM_NODELIST:-local}"
echo "Working dir: $(pwd)"
echo "Model code: ${MODEL_DIR}"
echo "Python: ${PYTHON_BIN}"
echo "Run tag: ${DT_DIRECT_RUN_TAG}"
echo "Results: ${DT_DIRECT_RESULTS_DIR}"
echo "Setup: ${DT_DIRECT_SETUP} | bounds: ${DT_DIRECT_BOUNDS} | closure: ${DT_DIRECT_POPULATION_CLOSURE}"
if [ "${DT_DIRECT_POPULATION_CLOSURE}" = "outside_option_benchmark_normalized" ]; then
    echo "Scale target: ${DT_DIRECT_SCALE_TARGET} imposed mechanically by benchmark outside-option accounting"
elif [ "${DT_DIRECT_POPULATION_CLOSURE}" = "renewal_valve_calibrated" ]; then
    echo "Scale target: ${DT_DIRECT_SCALE_TARGET} imposed mechanically | scale weight inactive"
else
    echo "Scale target: ${DT_DIRECT_SCALE_TARGET} | scale weight: ${DT_DIRECT_SCALE_WEIGHT}"
fi
echo "Outside value x0: ${DT_DIRECT_OUTSIDE_VALUE_X0} | outside flow x0: ${DT_DIRECT_OUTSIDE_FLOW_X0} | renewal rho: ${DT_DIRECT_RENEWAL_RETENTION}"
if [ -n "${DT_DIRECT_HR_MAX:-}" ]; then
    echo "Renter cap override hR_max: ${DT_DIRECT_HR_MAX}"
fi
if [ -n "${DT_DIRECT_OWNER_H_BAR_SCALE:-}" ]; then
    echo "Owner h_bar scale override: ${DT_DIRECT_OWNER_H_BAR_SCALE}"
fi
if [ -n "${DT_DIRECT_OWNER_SIZE_COST:-}" ]; then
    echo "Owner size cost: ${DT_DIRECT_OWNER_SIZE_COST} | ref: ${DT_DIRECT_OWNER_SIZE_COST_REF:-default} | power: ${DT_DIRECT_OWNER_SIZE_COST_POWER:-default}"
fi
if [ -n "${DT_DIRECT_TENURE_CHOICE_KAPPA:-}" ]; then
    echo "Tenure/product choice kappa: ${DT_DIRECT_TENURE_CHOICE_KAPPA}"
fi
if [ -n "${DT_DIRECT_WEIGHT_OVERRIDES:-}" ]; then
    echo "Weight overrides: ${DT_DIRECT_WEIGHT_OVERRIDES}"
fi
if [ -n "${DT_DIRECT_EXTRA_TARGETS:-}" ]; then
    echo "Extra targets: ${DT_DIRECT_EXTRA_TARGETS}"
fi
if [ -n "${DT_DIRECT_WARM_START_JSON:-}" ]; then
    echo "Warm start JSON: ${DT_DIRECT_WARM_START_JSON}"
fi
if [ -n "${DT_DIRECT_PARENT_DP_WAIVER:-}" ]; then
    echo "Parent DP waiver: ${DT_DIRECT_PARENT_DP_WAIVER} | phi: ${DT_DIRECT_PARENT_DP_WAIVER_PHI:-default}"
fi
if [ -n "${DT_DIRECT_H_OWN:-}" ]; then
    echo "Owner ladder override H_own: ${DT_DIRECT_H_OWN}"
fi
if [ -n "${DT_DIRECT_BOUND_OVERRIDES:-}" ]; then
    echo "Bound overrides: ${DT_DIRECT_BOUND_OVERRIDES}"
fi
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

if [ -n "${DT_DIRECT_HR_MAX:-}" ]; then
    ARGS+=(--hR-max "${DT_DIRECT_HR_MAX}")
fi

if [ -n "${DT_DIRECT_OWNER_H_BAR_SCALE:-}" ]; then
    ARGS+=(--owner-h-bar-scale "${DT_DIRECT_OWNER_H_BAR_SCALE}")
fi

if [ -n "${DT_DIRECT_OWNER_SIZE_COST:-}" ]; then
    ARGS+=(--owner-size-cost "${DT_DIRECT_OWNER_SIZE_COST}")
fi

if [ -n "${DT_DIRECT_OWNER_SIZE_COST_REF:-}" ]; then
    ARGS+=(--owner-size-cost-ref "${DT_DIRECT_OWNER_SIZE_COST_REF}")
fi

if [ -n "${DT_DIRECT_OWNER_SIZE_COST_POWER:-}" ]; then
    ARGS+=(--owner-size-cost-power "${DT_DIRECT_OWNER_SIZE_COST_POWER}")
fi

if [ -n "${DT_DIRECT_TENURE_CHOICE_KAPPA:-}" ]; then
    ARGS+=(--tenure-choice-kappa "${DT_DIRECT_TENURE_CHOICE_KAPPA}")
fi

if [ -n "${DT_DIRECT_WEIGHT_OVERRIDES:-}" ]; then
    ARGS+=(--weight-overrides "${DT_DIRECT_WEIGHT_OVERRIDES}")
fi

if [ -n "${DT_DIRECT_EXTRA_TARGETS:-}" ]; then
    ARGS+=(--extra-targets "${DT_DIRECT_EXTRA_TARGETS}")
fi

if [ -n "${DT_DIRECT_WARM_START_JSON:-}" ]; then
    ARGS+=(--warm-start-json "${DT_DIRECT_WARM_START_JSON}")
fi

if [ -n "${DT_DIRECT_PARENT_DP_WAIVER:-}" ]; then
    if [ "${DT_DIRECT_PARENT_DP_WAIVER}" = "1" ] || [ "${DT_DIRECT_PARENT_DP_WAIVER}" = "true" ] || [ "${DT_DIRECT_PARENT_DP_WAIVER}" = "yes" ]; then
        ARGS+=(--parent-dp-waiver)
    else
        ARGS+=(--no-parent-dp-waiver)
    fi
fi

if [ -n "${DT_DIRECT_PARENT_DP_WAIVER_PHI:-}" ]; then
    ARGS+=(--parent-dp-waiver-phi "${DT_DIRECT_PARENT_DP_WAIVER_PHI}")
fi

if [ -n "${DT_DIRECT_H_OWN:-}" ]; then
    ARGS+=(--H-own "${DT_DIRECT_H_OWN}")
fi

if [ -n "${DT_DIRECT_BOUND_OVERRIDES:-}" ]; then
    ARGS+=(--bound-overrides "${DT_DIRECT_BOUND_OVERRIDES}")
fi

if [ "${DT_DIRECT_MAX_ITER_EQ}" != "0" ]; then
    ARGS+=(--max-iter-eq "${DT_DIRECT_MAX_ITER_EQ}")
fi

"${PYTHON_BIN}" "${ARGS[@]}"

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
