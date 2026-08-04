#!/bin/bash
# Local two-chain E6a+E6b diagnostic under transparent block-equal relative weights.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MODEL_DIR="$(cd "$SCRIPT_DIR/../model" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PYTHON_BIN="$MODEL_DIR/.venv/bin/python"
RUN_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_plainvanilla_local_20260804/production"
LOG_ROOT="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_plainvanilla_local_20260804/logs"
CHAIN_IDS="${E6_LOCAL_CHAIN_IDS:-1 2}"
COLLECT_AFTER="${E6_LOCAL_COLLECT_AFTER:-1}"

mkdir -p "$RUN_ROOT" "$LOG_ROOT"

export E3_L4=1 E5=1 E6A=1 E6B=1
unset E6C
export E3_TFR_TOP_BIN_WEIGHT=3.602359422009
export E2_SEED_RECORD="$PROJECT_ROOT/output/model/eqscale_seq_e6ab_rescue_recalibration_20260728/report/results.json"
export PYTHONPATH="$MODEL_DIR:${PYTHONPATH:-}"
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

run_chain() {
    local task_id="$1"
    local start_mix="$2"
    "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/run_e1_chain.py" \
        --outdir "$RUN_ROOT/chain_$task_id" \
        --seed "$((2026080400 + task_id))" \
        --start-mix "$start_mix" \
        --minutes 465 \
        --max-evals 1200 \
        --weight-scheme target_relative_block_equal \
        >"$LOG_ROOT/chain_$task_id.log" 2>&1
}

pids=()
for task_id in $CHAIN_IDS; do
    case "$task_id" in
        1) start_mix=0.00 ;;
        2) start_mix=0.10 ;;
        3) start_mix=0.03 ;;
        4) start_mix=0.06 ;;
        5) start_mix=0.12 ;;
        6) start_mix=0.16 ;;
        7) start_mix=0.20 ;;
        8) start_mix=0.25 ;;
        *) echo "unsupported local chain id: $task_id" >&2; exit 2 ;;
    esac
    run_chain "$task_id" "$start_mix" &
    pids+=("$!")
done

exit_status=0
for pid in "${pids[@]}"; do
    wait "$pid" || exit_status=1
done

if [ "$COLLECT_AFTER" = "1" ]; then
    "$PYTHON_BIN" "$MODEL_DIR/intergen_eqscale_seq_optimized/collect_e1.py" \
        --results-root "$RUN_ROOT" \
        --outdir "$PROJECT_ROOT/output/model/eqscale_seq_e6ab_plainvanilla_local_20260804/report"
fi

exit "$exit_status"
