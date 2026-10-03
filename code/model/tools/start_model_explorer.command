#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PYTHON="$ROOT/code/model/.venv/bin/python"
CONFIG="$ROOT/output/model/fixed_reference_economics_20260928/soft_timing_review_v1/explorer_cases.json"
URL="http://127.0.0.1:8765"

if "$PYTHON" -c 'import urllib.request; urllib.request.urlopen("http://127.0.0.1:8765/api/meta", timeout=1)' >/dev/null 2>&1; then
  echo "The model explorer is already responding at $URL. No process was stopped."
  exit 0
fi

echo "Model explorer: $URL"
echo "Saved solutions only; this server performs no model solves. Press Ctrl-C to stop."
cd "$ROOT"
exec env NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$PYTHON" "$ROOT/code/model/tools/economics_explorer.py" --config "$CONFIG"
