#!/usr/bin/env bash
set -euo pipefail
stage=$(cd "$(dirname "$0")/stage" && pwd)
remote=/scratch/td2248/projects/transition_readiness_v1/current_floor
# Results are immutable job evidence, never a mirror-deletion target.
rsync -az --exclude=results --exclude='*_slurm.out' -e 'ssh -4 -oConnectTimeout=20 -oBatchMode=yes' "$stage/" "torch:$remote/"
