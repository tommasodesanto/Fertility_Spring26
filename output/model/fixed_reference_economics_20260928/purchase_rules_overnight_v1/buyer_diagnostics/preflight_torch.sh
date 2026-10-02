#!/usr/bin/env bash
# Exact staged-source authentication and native financial-map smoke; zero solves.
set -euo pipefail
remote=/scratch/td2248/projects/purchase_buyer_diagnostics_v3
mechanism=/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a
floor=/scratch/td2248/projects/normalized_floor_calibration_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
"$python" - "$remote" "$mechanism" "$floor" "$base" <<'PY'
import hashlib,json,sys
from pathlib import Path
for folder in map(Path,sys.argv[1:]):
    inventory=json.loads((folder/'inventory.json').read_text())
    for rel,digest in inventory['files'].items():
        source=folder/'source'/rel
        assert source.is_file() and hashlib.sha256(source.read_bytes()).hexdigest()==digest,(folder,rel)
print('all staged buyer/mechanism/floor/base source hashes verified')
PY
mkdir -p "$remote/preflight/numba_cache" "$remote/preflight/matplotlib"
export NUMBA_CACHE_DIR=/work/preflight/numba_cache MPLCONFIGDIR=/work/preflight/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$remote/preflight:/work/preflight:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']): print(rel)
PY
)
binds+=(--bind "$mechanism/source/$packet:$repo/$packet:ro")
binds+=(--bind "$mechanism/source/code/model/experiments/transition_readiness:$repo/code/model/experiments/transition_readiness:ro")
binds+=(--bind "$mechanism/source/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json:$repo/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json:ro")
binds+=(--bind "$remote/source/$packet/buyer_diagnostics:$repo/$packet/buyer_diagnostics:ro")
for arm in hard quarter; do
    apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/buyer_diagnostics/preflight_native.py" --arm "$arm" | tee "$remote/preflight/${arm}.json"
done
