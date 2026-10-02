#!/usr/bin/env bash
#SBATCH --job-name=purchase_phi
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/purchase_mechanism_horizon_extension_v1/logs/%x-%A_%a.out
set -euo pipefail

task=${SLURM_ARRAY_TASK_ID:?Array task required}
mode=${MECHANISM_SMOKE:-0}
if [[ "$mode" == 1 ]]; then
  [[ "$task" =~ ^[01]$ ]] || { echo 'Smoke task must be 0 or 1'; exit 2; }
else
  [[ "$task" =~ ^([0-9]|1[01])$ ]] || { echo 'Production task must be 0–11'; exit 2; }
fi
remote=/scratch/td2248/projects/purchase_mechanism_horizon_extension_v1
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

arms=(hard hard hard hard hard hard quarter quarter quarter quarter quarter quarter)
kinds=(control temporary permanent control temporary permanent control temporary permanent control temporary permanent)
horizons=(48 48 48 64 64 64 48 48 48 64 64 64)
if [[ "$mode" == 1 ]]; then
  arm=hard; [[ "$task" == 1 ]] && arm=quarter
  kind=control; horizon=1; name="smoke_${arm}"; cap=64
else
  arm=${arms[$task]}; kind=${kinds[$task]}; horizon=${horizons[$task]}
  printf -v label '%02d' "$task"
  name="case_${label}_${arm}_${kind}_h${horizon}"
  cap=1024
fi
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+14400))
(( deadline_epoch > 1790949600 )) && deadline_epoch=1790949600
(( deadline_epoch-start_epoch >= 600 )) || { echo 'Less than ten minutes remain before 10 AM deadline'; exit 124; }
out="$remote/results/$name"
mkdir "$out" || { echo 'Refusing existing mechanism result'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal_receipt() {
  local code=$?
  trap - EXIT
  "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" "$task" "$mode" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),start_epoch=int(sys.argv[3]),deadline_epoch=int(sys.argv[4]),task=int(sys.argv[5]),smoke=sys.argv[6]=='1',finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$code"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM

# Check every mounted source and the chosen fitted checkpoint before a model call.
"$python" - "$remote" "$selection_store" "$base" "$floor" "$arm" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,selection_store,base,floor=map(Path,sys.argv[1:5]);arm=sys.argv[5]
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''): h.update(block)
 return h.hexdigest()
for folder in (base,floor,remote):
 inv=json.loads((folder/'inventory.json').read_text())
 if folder==remote:
  assert hashlib.sha256((folder/'inventory.json').read_bytes()).hexdigest()=='a4559bf682de6e15461e6e35df595769e515d6b1c5c236900ff8fc46b9dbc24e'
 for rel,digest in inv['files'].items():
  path=folder/'source'/rel
  assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==digest,rel
manifest_path=selection_store/'selection/manifest.json'
assert hashlib.sha256(manifest_path.read_bytes()).hexdigest()=='19db8dfcc928c4cdea5810f70c4e62a9de65372bd038646d0bc2d49532913879'
selection=json.loads(manifest_path.read_text())
assert selection['target_fingerprint']=='db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1'
assert selection['weight_fingerprint']=='2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0'
entry=selection['arms'][arm]
chosen=selection_store/'selection'/('selected_'+arm+'.json')
assert hashlib.sha256(chosen.read_bytes()).hexdigest()==entry['selected_json_sha256']
candidate=json.loads(chosen.read_text())
assert candidate['status']=='postchecked' and candidate['arm']==arm and int(candidate['chain'])==entry['chain']
snapshot=selection_store/'selected_postchecks'/('chain_'+str(entry['chain']))
assert entry['snapshot_remote_root']==str(snapshot) and candidate['snapshot_remote_root']==str(snapshot)
assert candidate['origin']==entry['origin'] and candidate['remote_root']==entry['physical_remote_root']
assert candidate.get('parent_remote_root')==entry['parent_remote_root']
completed=snapshot/'postcheck/completed.json'
assert completed.is_file() and hashlib.sha256(completed.read_bytes()).hexdigest()==entry['completed_sha256']
assert json.loads(completed.read_text())['status']=='selected_numerically_verified'
report=completed.parent/'selected_postcheck/phase_b_ge/selected_root'
for name,digest in entry['report_sha256'].items():
 assert hashlib.sha256((report/name).read_bytes()).hexdigest()==digest
arrays=completed.parent/'selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz'
assert sha(arrays)==entry['native_arrays_sha256']
print(entry['chain'])
PY

chain=$("$python" - "$selection_store/selection/selected_${arm}.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]));print(int(d['chain']))
PY
)
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']): print(rel)
PY
)
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
binds+=(--bind "$remote/source/code/model/experiments/transition_readiness:$repo/code/model/experiments/transition_readiness:ro")
binds+=(--bind "$remote/source/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json:$repo/output/model/transition_readiness_v1/normalized_restart_v1/deployment/fit_plan.json:ro")
binds+=(--bind "$selection_store/selected_postchecks:$repo/$packet/results:ro")
binds+=(--bind "$selection_store/selection:$repo/$packet/collection/readout:ro")
remaining=$((deadline_epoch-$(date +%s)-15))
(( remaining > 0 )) || exit 124
arguments=(--arm "$arm" --selected-json "$repo/$packet/collection/readout/selected_${arm}.json"
 --selected-completed "$repo/$packet/results/chain_${chain}/postcheck/completed.json"
 --kind "$kind" --horizon "$horizon" --out /work/results/run
 --deadline-epoch "$deadline_epoch" --maximum-policy-calls "$cap")
[[ "$mode" == 1 ]] && arguments+=(--smoke-one-date)
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/mechanism/run_case_extended.py" "${arguments[@]}" > "$out/run.log" 2>&1
