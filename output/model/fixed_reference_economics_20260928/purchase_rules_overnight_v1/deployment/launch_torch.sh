#!/usr/bin/env bash
#SBATCH --job-name=purchase80
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/purchase_rules_overnight_v1/logs/%x-%A_%a.out
set -euo pipefail

start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+14400))
task=${SLURM_ARRAY_TASK_ID:?Requires 48-chain array}
[[ "$task" =~ ^([0-9]|[1-3][0-9]|4[0-7])$ ]] || { echo 'Invalid chain'; exit 2; }
remote=/scratch/td2248/projects/purchase_rules_overnight_v1
floor_remote=/scratch/td2248/projects/normalized_floor_calibration_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
mkdir -p "$remote/results"
if [[ "${CALIBRATION_SMOKE:-0}" == 1 ]]; then
 out="$remote/results/smoke_chain_$task"
else
 out="$remote/results/chain_$task"
fi
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

"$python" - "$remote" "$base" "$floor_remote" "$inputs" "$packet" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,base,floor_remote,inputs=map(Path,sys.argv[1:5]);packet=sys.argv[5]
for folder in (base,floor_remote,remote):
 inv=json.loads((folder/'inventory.json').read_text())
 for rel,digest in inv['files'].items():
  p=folder/'source'/rel
  assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
inv=json.loads((remote/'inventory.json').read_text())
launcher_rel=packet+'/deployment/launch_torch.sh'
assert hashlib.sha256((remote/'launch_torch.sh').read_bytes()).hexdigest()==inv['files'][launcher_rel],'Entrypoint launcher drift'
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7','arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
 assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
stage=json.loads((remote/'source'/packet/'engine_pins.json').read_text())
for rel,digest in stage.items():
 p=remote/'source'/packet/rel
 assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
PY

"$python" - "$out" "$start_epoch" "$deadline_epoch" "$task" "$$" <<'PY'
import json,os,sys
from pathlib import Path
Path(sys.argv[1],'launcher_start.json').write_text(json.dumps(dict(start_epoch=int(sys.argv[2]),deadline_epoch=int(sys.argv[3]),chain=int(sys.argv[4]),launcher_pid=int(sys.argv[5]),slurm_job_id=os.getenv('SLURM_JOB_ID'),threads=1,cpus=1,memory_GiB=24,maximum_objective_calls=250,final_reserve_seconds=900),indent=2)+'\n')
PY
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
while IFS= read -r rel; do binds+=(--bind "$floor_remote/source/$rel:$repo/$rel:ro"); done < <("$python" - "$floor_remote/inventory.json" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']): print(rel)
PY
)
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
 [[ "$remaining" -gt 0 ]] || { echo 'Four-hour actual-start deadline reached'; exit 124; }
 timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "$repo/$packet/run_psi.py" --chain "$task" --out "/work/results/$stage" --deadline-epoch "$deadline_epoch" "$@" > "$out/$stage.log" 2>&1
}
run_stage init --fast-objective --initialize-only
if [[ "${CALIBRATION_SMOKE:-0}" == 1 ]]; then
 run_stage smoke --fast-objective --smoke-only
 exit 0
fi
run_stage search --fast-objective
if "$python" - "$out/search/search_completed.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1])).get('selected') is not None else 1)
PY
then
 run_stage postcheck --verify-only /work/results/search/search_completed.json
fi
