#!/usr/bin/env bash
#SBATCH --job-name=purchase_buyer_readout
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/purchase_buyer_diagnostics_v6/logs/%x-%j.out
# Saved-policy postprocessing only. Submit with BUYER_MODE/BUYER_ARM etc.
set -euo pipefail
mode=${BUYER_MODE:?Set BUYER_MODE=selected or dated}
arm=${BUYER_ARM:?Set BUYER_ARM=hard or quarter}
[[ "$mode" == selected || "$mode" == dated ]] || exit 2
[[ "$arm" == hard || "$arm" == quarter ]] || exit 2
remote=/scratch/td2248/projects/purchase_buyer_diagnostics_v6
mechanism=/scratch/td2248/projects/purchase_mechanism_reviewed_93831f5a
if [[ "$mode" == dated ]]; then
  mechanism=/scratch/td2248/projects/purchase_mechanism_horizon_extension_v1
fi
selection_store=/scratch/td2248/projects/purchase_mechanism_v1
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
PY
case=${BUYER_CASE:-}
date=${BUYER_DATE:-date_000}
if [[ "$mode" == selected ]]; then
  read -r chain < <("$python" - "$selection_store/selection/manifest.json" "$selection_store/selection/selected_${arm}.json" "$selection_store" "$arm" <<'PY'
import hashlib,json,sys
from pathlib import Path
m=json.load(open(sys.argv[1]));s=json.load(open(sys.argv[2]));c=Path(sys.argv[3]);arm=sys.argv[4];e=m['arms'][arm]
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''): h.update(block)
 return h.hexdigest()
p=c/'selected_postchecks'/('chain_'+str(e['chain']))/'postcheck/completed.json'
assert s['status']=='postchecked' and s['arm']==arm and int(s['chain'])==int(e['chain'])
assert e['snapshot_remote_root']==str(p.parent.parent) and s['snapshot_remote_root']==str(p.parent.parent)
assert s['origin']==e['origin'] and s['remote_root']==e['physical_remote_root']
assert s.get('parent_remote_root')==e['parent_remote_root']
assert hashlib.sha256(Path(sys.argv[2]).read_bytes()).hexdigest()==e['selected_json_sha256']
assert hashlib.sha256(p.read_bytes()).hexdigest()==e['completed_sha256']
assert json.load(open(p))['status']=='selected_numerically_verified'
report=p.parent/'selected_postcheck/phase_b_ge/selected_root'
for name,digest in e['report_sha256'].items():
 assert hashlib.sha256((report/name).read_bytes()).hexdigest()==digest
arrays=p.parent/'selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz'
assert sha(arrays)==e['native_arrays_sha256']
print(e['chain'])
PY
)
  name="selected_${arm}"
else
  [[ "$case" == case_06_quarter_control_h48 || "$case" == case_07_quarter_temporary_h48 ]] || exit 2
  [[ "$date" == date_000 ]] || exit 2
  [[ "$case" == *"_${arm}_"* ]] || exit 2
  dated_receipt=$("$python" "$remote/source/$packet/buyer_diagnostics/verify_dated_extension.py" \
    --root "$mechanism" --store "$selection_store" --case "$case" --date "$date")
  read -r relative observed_phi < <("$python" -c \
    'import json,sys; d=json.load(sys.stdin); print(d["relative_packet"],d["observed_phi"])' \
    <<< "$dated_receipt")
  name="${case}_${date}"
fi
out="$remote/results/$name"
mkdir "$out" || { echo 'Refusing existing buyer readout'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
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
binds+=(--bind "$selection_store/selected_postchecks:$repo/$packet/results:ro")
binds+=(--bind "$selection_store/selection:$repo/$packet/collection/readout:ro")
binds+=(--bind "$mechanism/results:$repo/$packet/mechanism/results:ro")
if [[ "$mode" == selected ]]; then
  script="$repo/$packet/buyer_diagnostics/run_selected.py"
  args=(--arm "$arm" --completed "$repo/$packet/results/chain_${chain}/postcheck/completed.json"
        --selection-json "$repo/$packet/collection/readout/selected_${arm}.json"
        --selection-manifest "$repo/$packet/collection/readout/manifest.json" --out /work/results/run)
else
  script="$repo/$packet/buyer_diagnostics/run_dated.py"
  args=(--engine-root "$repo/$packet/engines/$arm"
        --packet "$repo/$packet/mechanism/results/$case/run/$relative"
        --rule "$arm" --observed-phi "$observed_phi" --out /work/results/dated_access.json)
fi
timeout --signal=TERM --kill-after=10s 1100s apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$script" "${args[@]}" > "$out/run.log" 2>&1
