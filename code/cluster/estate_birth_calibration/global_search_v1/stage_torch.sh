#!/usr/bin/env bash
# Stage the reviewed frozen v2 source overlay plus new exploration entrypoints.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
remote=/scratch/td2248/projects/estate_birth_global_search_20261004_v1
parent=/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2
local_deploy="$repo/output/model/experiments/birth_count_choice/estate_a_global_search_20261004_v1/deployment"
test -f "$local_deploy/stage_manifest.json" || python3 "$here/build_stage.py"
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir -p '$remote/control' '$remote/logs' '$remote/results' && cp -a '$parent/source' '$remote/source' && cp '$parent/inventory.json' '$remote/inventory.json'"
scp "$local_deploy/stage_manifest.json" "$local_deploy/incumbent.json" "$local_deploy/parent_starts.sha256" torch:"$remote/"
scp "$here/prepare_plan.py" "$here/explore.py" "$here/verify_stage.py" "$here/launch_torch.sh" "$here/submit_once.py" torch:"$remote/"
ssh -o BatchMode=yes torch "chmod 755 '$remote/launch_torch.sh' && /share/apps/anaconda3/2025.06/bin/python '$remote/verify_stage.py' --host && /share/apps/anaconda3/2025.06/bin/python '$remote/prepare_plan.py'"
scp "torch:$remote/control/plan.json" "$local_deploy/plan.json"
