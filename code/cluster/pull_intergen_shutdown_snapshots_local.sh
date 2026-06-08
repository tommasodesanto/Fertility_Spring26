#!/bin/bash
# Pull Torch shutdown snapshots back to the local repo on a timed window.
# Run locally from any directory:
#   code/cluster/pull_intergen_shutdown_snapshots_local.sh 2026-06-09T06:42:00 25 45

set -euo pipefail

TARGET_LOCAL="${1:?Usage: $0 YYYY-MM-DDTHH:MM:SS [attempts] [sleep_seconds]}"
ATTEMPTS="${2:-25}"
SLEEP_SECONDS="${3:-45}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTDIR="${REPO_DIR}/output/model/intergen_shutdown_snapshots"
LOGDIR="${OUTDIR}/logs"
REMOTE_DIR="/scratch/td2248/projects/Fertility_Spring26/code/cluster/shutdown_snapshots"

mkdir -p "${LOGDIR}" "${OUTDIR}"

DELAY_SECONDS="$(
    python3 - "${TARGET_LOCAL}" <<'PY'
from datetime import datetime
import sys
from zoneinfo import ZoneInfo

target = datetime.fromisoformat(sys.argv[1]).replace(tzinfo=ZoneInfo("Europe/Rome"))
now = datetime.now(ZoneInfo("Europe/Rome"))
print(max(0, int((target - now).total_seconds())))
PY
)"

echo "Target local time: ${TARGET_LOCAL} Europe/Rome"
echo "Initial delay seconds: ${DELAY_SECONDS}"
echo "Attempts: ${ATTEMPTS}; sleep seconds: ${SLEEP_SECONDS}"
echo "Output directory: ${OUTDIR}"
sleep "${DELAY_SECONDS}"

for attempt in $(seq 1 "${ATTEMPTS}"); do
    echo "Attempt ${attempt}/${ATTEMPTS}: $(date)"
    rsync -az --partial "torch:${REMOTE_DIR}/" "${OUTDIR}/" || true
    sleep "${SLEEP_SECONDS}"
done

echo "Finished: $(date)"
