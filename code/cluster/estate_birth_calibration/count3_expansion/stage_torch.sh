#!/usr/bin/env bash
# Build and stage immutable derivative, with zero model solves and no submission.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
remote=/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1
python3 "$here/build_stage.py"
archive="$repo/output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/stage.tar.gz"
ssh -o BatchMode=yes torch "test ! -e '$remote' && mkdir -p '$remote/logs'"
scp "$archive" "torch:$remote/stage.tar.gz"
ssh -o BatchMode=yes torch "cd '$remote' && tar -xzf stage.tar.gz && chmod 755 ./*.sh && /share/apps/anaconda3/2025.06/bin/python verify_stage.py --host"
ssh -o BatchMode=yes torch "cd '$remote' && /share/apps/anaconda3/2025.06/bin/python storage_check.py"
scp "torch:$remote/staging_storage.json" "$repo/output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/staging_storage.json"
