#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PYTHON="$ROOT/code/model/.venv/bin/python"
HISTORICAL_CONFIG="$ROOT/output/model/fixed_reference_economics_20260928/soft_timing_review_v1/explorer_cases.json"
LATEST_CONFIG="$ROOT/output/model/local_solution/latest/explorer_cases.json"
if [[ $# -gt 1 ]]; then
  echo "Usage: $0 [explorer_cases.json]" >&2
  exit 2
elif [[ $# -eq 1 ]]; then
  CONFIG="$1"
  if [[ "$CONFIG" != /* ]]; then
    CONFIG="$(pwd)/$CONFIG"
  fi
elif [[ -f "$LATEST_CONFIG" ]]; then
  CONFIG="$LATEST_CONFIG"
else
  CONFIG="$HISTORICAL_CONFIG"
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "Explorer case config not found: $CONFIG" >&2
  exit 1
fi
URL="http://127.0.0.1:8765"

if "$PYTHON" -c 'import urllib.request; urllib.request.urlopen("http://127.0.0.1:8765/api/meta", timeout=1)' >/dev/null 2>&1; then
  echo "The model explorer is already responding at $URL. No process was stopped."
  exit 0
fi

echo "Model explorer: $URL"
echo "Case config: $CONFIG"
echo "Saved solutions only; this server performs no model solves. Press Ctrl-C to stop."
cd "$ROOT"
exec env NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$PYTHON" "$ROOT/code/model/tools/economics_explorer.py" --config "$CONFIG"
