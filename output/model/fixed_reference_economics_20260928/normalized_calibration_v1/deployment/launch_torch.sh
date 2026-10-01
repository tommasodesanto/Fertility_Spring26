#!/usr/bin/env bash
#SBATCH --job-name=norm_floor
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/normalized_calibration_v1/logs/%x-%A_%a.out
set -euo pipefail

start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+7200))
task=${SLURM_ARRAY_TASK_ID:?Requires 24-chain array}
[[ "$task" =~ ^([0-9]|1[0-9]|2[0-3])$ ]] || { echo 'Invalid chain'; exit 2; }
remote=/scratch/td2248/projects/normalized_calibration_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/normalized_calibration_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
mkdir -p "$remote/results"
out="$remote/results/chain_$task"
mkdir "$out" || { echo 'Refusing existing chain result'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal_receipt() {
  local code=$?
  trap - EXIT
  "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" "$task" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),start_epoch=int(sys.argv[3]),deadline_epoch=int(sys.argv[4]),chain=int(sys.argv[5]),finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$code"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM

"$python" - "$remote" "$base" "$inputs" "$packet" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,base,inputs=map(Path,sys.argv[1:4]);packet=sys.argv[4]
for folder in (base,remote):
 inv=json.loads((folder/'inventory.json').read_text())
 for rel,digest in inv['files'].items():
  p=folder/'source'/rel
  assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
inv=json.loads((remote/'inventory.json').read_text())
launcher_rel=packet+'/deployment/launch_torch.sh'
assert hashlib.sha256((remote/'launch_torch.sh').read_bytes()).hexdigest()==inv['files'][launcher_rel],'Entrypoint launcher drift'
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7','arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
 assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
stage=json.loads((remote/'source'/packet/'source_pins.json').read_text())
assert set(stage).issubset(inv['files']),'Pinned source missing from staged inventory'
for rel,digest in stage.items():
 assert hashlib.sha256((remote/'source'/rel).read_bytes()).hexdigest()==digest,rel
PY

# Searches cannot begin unless the one matched, full-native incumbent ROOT/REPEAT
# gate passed and its target identity matches this staged package.
"$python" - "$remote/results/normalized_incumbent_smoke/completed.json" "$remote/source/$packet/plan.json" <<'PY'
import json,sys
import hashlib
from pathlib import Path
gate=json.loads(Path(sys.argv[1]).read_text());plan=json.loads(Path(sys.argv[2]).read_text())
assert gate.get('status')=='normalized_incumbent_native_smoke_passed','Native incumbent gate did not pass'
assert gate.get('lifecycle_solves',0)>0,'Smoke did not perform native lifecycle solves'
assert gate.get('comparison',{}).get('target_fit.csv',{}).get('rows')==14,'Target row identity missing'
assert gate.get('comparison',{}).get('parameters.csv',{}).get('rows')==31,'Parameter row identity missing'
assert gate.get('plan_sha256')==hashlib.sha256(Path(sys.argv[2]).read_bytes()).hexdigest(),'Smoke plan identity mismatch'
pins=Path(sys.argv[2]).with_name('source_pins.json')
assert gate.get('source_pins_sha256')==hashlib.sha256(pins.read_bytes()).hexdigest(),'Smoke source identity mismatch'
assert plan.get('normalization',{}).get('population')==1.0,'N0=1 target contract missing'
PY

"$python" - "$out" "$start_epoch" "$deadline_epoch" "$task" "$$" <<'PY'
import json,os,sys
from pathlib import Path
Path(sys.argv[1],'launcher_start.json').write_text(json.dumps(dict(start_epoch=int(sys.argv[2]),deadline_epoch=int(sys.argv[3]),chain=int(sys.argv[4]),launcher_pid=int(sys.argv[5]),slurm_job_id=os.getenv('SLURM_JOB_ID'),threads=1,cpus=1,memory_GiB=24,maximum_objective_calls=100,final_reserve_seconds=900),indent=2)+'\n')
PY
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$remote/source/$packet:$repo/$packet:ro")
while IFS= read -r rel; do binds+=(--bind "$remote/source/$rel:$repo/$rel:ro"); done < <("$python" - "$remote/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if not rel.startswith(sys.argv[2]+'/'): print(rel)
PY
)
run_stage() {
 local stage=$1; shift
 local remaining=$((deadline_epoch-$(date +%s)-15))
 [[ "$remaining" -gt 0 ]] || { echo 'Two-hour actual-start deadline reached'; exit 124; }
 timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/run_psi.py" --chain "$task" --out "/work/results/$stage" --deadline-epoch "$deadline_epoch" "$@" > "$out/$stage.log" 2>&1
}
run_stage init --fast-objective --initialize-only
run_stage search --fast-objective
if "$python" - "$out/search/search_completed.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1])).get('selected') is not None else 1)
PY
then
 run_stage postcheck --verify-only /work/results/search/search_completed.json
fi
