#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PYTHON="$ROOT/code/model/.venv/bin/python"
usage() { echo "Usage: $0 [--params parameter_file.py | explorer_cases.json]"; }
if [[ $# -eq 1 && ( "$1" == --help || "$1" == -h ) ]]; then
  usage
  exit 0
elif [[ $# -eq 2 && "$1" == --params ]]; then
  PARAMS="$2"
  CONFIG="$("$PYTHON" -c 'import sys; from pathlib import Path; sys.path.insert(0, sys.argv[1]); from production.parameter_files import output_root_for, describe_saved_case; root = output_root_for(sys.argv[2]); print(describe_saved_case(root / "latest", sys.argv[2]), file=sys.stderr); print(root / "latest/explorer_cases.json")' "$ROOT/code/model" "$PARAMS")"
elif [[ $# -eq 0 ]]; then
  CONFIG="$("$PYTHON" -c 'import sys; sys.path.insert(0, sys.argv[1]); from run_model import PARAMETER_FILE; from production.parameter_files import output_root_for, describe_saved_case; root = output_root_for(PARAMETER_FILE); print(describe_saved_case(root / "latest", PARAMETER_FILE), file=sys.stderr); print(root / "latest/explorer_cases.json")' "$ROOT/code/model")"
elif [[ $# -eq 1 && "$1" != --* ]]; then
  CONFIG="$1"
  if [[ "$CONFIG" != /* ]]; then
    CONFIG="$(pwd)/$CONFIG"
  fi
else
  usage >&2
  exit 2
fi
if [[ ! -f "$CONFIG" ]]; then
  echo "No completed saved explorer config found: $CONFIG. Run run_model.py with the matching parameter file first." >&2
  exit 1
fi
# Resolve latest once, then compare the exact saved configuration with live servers.
# This probe uses only the standard library; it never loads or solves the model.
SELECTION="$("$PYTHON" - "$CONFIG" <<'PY'
import hashlib
import json
import socket
import sys
import urllib.request
from pathlib import Path

config = Path(sys.argv[1]).resolve(strict=True)
fingerprint = hashlib.sha256(config.read_bytes()).hexdigest()
first_free = None
for port in range(8765, 8786):
    with socket.socket() as probe:
        probe.settimeout(0.2)
        if probe.connect_ex(('127.0.0.1', port)) != 0:
            if first_free is None:
                first_free = port
            continue
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/meta', timeout=0.5) as response:
            meta = json.load(response)
        if (isinstance(meta, dict) and meta.get('config_path') == str(config)
                and meta.get('config_sha256') == fingerprint):
            print(config)
            print(port)
            print('reuse')
            break
    except (OSError, ValueError):
        pass
else:
    if first_free is None:
        sys.exit('No available explorer port between 8765 and 8785; no process was stopped.')
    print(config)
    print(first_free)
    print('start')
PY
)"
CONFIG="$(printf '%s\n' "$SELECTION" | sed -n '1p')"
PORT="$(printf '%s\n' "$SELECTION" | sed -n '2p')"
ACTION="$(printf '%s\n' "$SELECTION" | sed -n '3p')"
URL="http://127.0.0.1:$PORT"

if [[ "$ACTION" == reuse ]]; then
  echo "The matching model explorer is already responding at $URL. No process was stopped."
  echo "Case config: $CONFIG"
  exit 0
fi

echo "Model explorer: $URL"
echo "Case config: $CONFIG"
echo "Saved solutions only; this server performs no model solves. Press Ctrl-C to stop."
cd "$ROOT"
exec env NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$PYTHON" "$ROOT/code/model/tools/economics_explorer.py" --config "$CONFIG" --port "$PORT"
