#!/usr/bin/env bash
#SBATCH --job-name=floor_one_shock
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
set -euo pipefail
mode=; seconds=; plan=; label=
while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) mode=$2; shift 2;;
    --seconds) seconds=$2; shift 2;;
    --plan) plan=$2; shift 2;;
    --label) label=$2; shift 2;;
    *) echo "Unknown launcher argument: $1" >&2; exit 2;;
  esac
done
[[ "$mode" =~ ^(preflight|smoke|fit)$ && "$seconds" =~ ^[0-9]+$ && "$seconds" -ge 30 && "$label" =~ ^[a-zA-Z0-9_-]+$ ]] || { echo 'Explicit mode, bounded seconds and new label required'; exit 2; }
[[ "$mode" == preflight || -n "$plan" ]] || { echo 'Smoke/fit requires a pinned plan'; exit 2; }
remote=/scratch/td2248/projects/transition_readiness_v1/current_floor
calibration=/scratch/td2248/projects/utility_floor_psi_continuation_v1
base=/scratch/td2248/projects/grid_resolution_credit053_v2
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
packet=output/model/fixed_reference_economics_20260928/utility_floor_psi_continuation_v1
code=code/model/experiments/transition_readiness
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
start_epoch=$(date +%s); deadline_epoch=$((start_epoch+seconds))
mkdir -p "$remote/results"; out="$remote/results/$label"
mkdir "$out" || { echo 'Refusing existing output'; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal_receipt() {
 local code=$?; trap - EXIT
 [[ -z "${heartbeat_pid:-}" ]] || kill "$heartbeat_pid" 2>/dev/null || true
 "$python" - "$out" "$code" "$start_epoch" "$deadline_epoch" "$mode" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1],'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(sys.argv[2]),start_epoch=int(sys.argv[3]),deadline_epoch=int(sys.argv[4]),mode=sys.argv[5],finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
 exit "$code"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM
(while true; do date -u +%Y-%m-%dT%H:%M:%SZ > "$out/launcher_heartbeat.txt"; sleep 300; done) >/dev/null 2>&1 < /dev/null &
heartbeat_pid=$!
"$python" - "$remote" "$base" "$calibration" "$inputs" "$repo" "$plan" <<'PY'
import hashlib,json,sys
from pathlib import Path
remote,base,calibration,inputs,repo=map(Path,sys.argv[1:6])
for folder in (base,calibration,remote):
 for rel,digest in json.loads((folder/'inventory.json').read_text())['files'].items():
  assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest,rel
own=json.loads((remote/'inventory.json').read_text())['files']
launcher='code/model/experiments/transition_readiness/floor_launch.sh'
assert hashlib.sha256((remote/'floor_launch.sh').read_bytes()).hexdigest()==own[launcher], 'Launcher drift'
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7','arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
 assert hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest,name
if sys.argv[6]:
 plan=Path(sys.argv[6]);relative=str(plan.relative_to(repo))
 assert relative in own,'Plan is not inventoried'
PY
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
binds=(--bind "$frozen:$repo:ro" --bind "$out:/work/results:rw" --bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
while IFS= read -r rel; do binds+=(--bind "$base/source/$rel:$repo/$rel:ro"); done < "$base/mounts.txt"
binds+=(--bind "$calibration/source/$packet:$repo/$packet:ro")
while IFS= read -r rel; do binds+=(--bind "$calibration/source/$rel:$repo/$rel:ro"); done < <("$python" - "$calibration/inventory.json" "$packet" <<'PY'
import json,sys
for rel in sorted(json.load(open(sys.argv[1]))['files']):
 if not rel.startswith(sys.argv[2]+'/'):print(rel)
PY
)
# Own overlay is last: explicit inventories preserve all calibrated native mounts.
while IFS= read -r rel; do binds+=(--bind "$remote/source/$rel:$repo/$rel:ro"); done < "$remote/mounts.txt"
"$python" - "$out" "$start_epoch" "$deadline_epoch" "$mode" <<'PY'
import json,os,sys
from pathlib import Path
Path(sys.argv[1],'launcher_start.json').write_text(json.dumps(dict(start_epoch=int(sys.argv[2]),deadline_epoch=int(sys.argv[3]),mode=sys.argv[4],slurm_job_id=os.getenv('SLURM_JOB_ID'),threads=1,cpus=4,memory_GiB=96),indent=2)+'\n')
PY
if [[ "$mode" == preflight ]]; then
 handoff=output/model/transition_readiness_v1/current_floor_handoff/handoff.json
 digest=$("$python" - "$remote/inventory.json" "$handoff" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['files'][sys.argv[2]])
PY
)
 args=("$repo/$code/floor_runtime.py" --handoff "$repo/$handoff" --handoff-sha256 "$digest" --output /work/results --native-import)
else
 args=("$repo/$code/one_shock_floor.py" --plan "$plan" --output /work/results)
 [[ "$mode" == smoke ]] && args+=(--native-smoke) || args+=(--execute)
fi
remaining=$((deadline_epoch-$(date +%s)-10)); [[ "$remaining" -gt 0 ]] || exit 124
# Finite external wall limit complements controller per-stage and policy-call caps.
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" "${args[@]}" > "$out/run.log" 2>&1
