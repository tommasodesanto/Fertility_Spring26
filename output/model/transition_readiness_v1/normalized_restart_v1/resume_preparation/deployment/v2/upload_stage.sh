#!/usr/bin/env bash
set -euo pipefail
folder=$(cd "$(dirname "$0")" && pwd)
rsync -az --exclude=results -e 'ssh -4 -oConnectTimeout=20 -oBatchMode=yes' "$folder/stage/" torch:/scratch/td2248/projects/transition_readiness_v1/normalized_resumed_fit_v2/
