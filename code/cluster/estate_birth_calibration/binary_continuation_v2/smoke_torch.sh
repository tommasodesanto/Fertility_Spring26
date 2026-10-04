#!/usr/bin/env bash
# Submit one exact two-call full-native smoke with fresh selected-point native repeat.
set -euo pipefail
remote=/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2
chain=${1:-0}
[[ "$chain" == 0 ]] || { echo 'The smoke gate is pinned to chain 0'; exit 2; }
ssh -o BatchMode=yes torch "/share/apps/anaconda3/2025.06/bin/python '$remote/submit_smoke.py' '$chain'"
