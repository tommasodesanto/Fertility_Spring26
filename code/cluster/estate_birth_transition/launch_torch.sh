#!/usr/bin/env bash
# Separate zero-solve preflight, actual exact-loop smoke, and lead-reviewed fit.
#SBATCH --job-name=estate_transition
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cl
#SBATCH --output=/scratch/td2248/projects/current_estate_transition_20261003_v5/logs/%x-%j.out
set -euo pipefail
remote=/scratch/td2248/projects/current_estate_transition_20261003_v5
base=/scratch/td2248/projects/grid_resolution_credit053_v2
floor_remote=/scratch/td2248/projects/normalized_floor_calibration_v1
frozen=/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project
inputs=/scratch/td2248/projects/publication_refactor_20260929/export_v1/inputs
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
mode=${1:-${TRANSITION_MODE:-preflight}}
[[ "$mode" =~ ^(preflight|smoke|fit)$ ]] || { echo 'Invalid mode'; exit 2; }
[[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'Array execution forbidden'; exit 2; }
if [[ "$mode" == preflight ]]; then seconds=${TRANSITION_WALL_SECONDS:-600}; else
  [[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Smoke/fit require one Slurm job'; exit 2; }
  seconds=${TRANSITION_WALL_SECONDS:?Explicit wall budget required}
fi
[[ "$seconds" =~ ^[0-9]+$ && "$seconds" -gt 0 && "$seconds" -le 21600 ]] || { echo 'Invalid bounded wall budget'; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-8}" == 8 ]] || { echo 'Exactly eight allocated CPUs required'; exit 2; }
[[ "${SLURM_MEM_PER_NODE:-98304}" -le 98304 ]] || { echo 'Memory exceeds 96 GiB'; exit 2; }
module load anaconda3/2025.06
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
start_epoch=$(date +%s); deadline_epoch=$((start_epoch+seconds))
mkdir -p "$remote/results"
out="$remote/results/$mode"
mkdir "$out" || { echo "Refusing existing output: $out"; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
terminal_receipt() {
  local status=$?; trap - EXIT
  [[ -z "${heartbeat_pid:-}" ]] || kill "$heartbeat_pid" 2>/dev/null || true
  "$python" - "$out" "$status" "$mode" "$start_epoch" "$deadline_epoch" <<'PY'
import json,os,sys,time
from pathlib import Path
out,status,mode,start,deadline=sys.argv[1:]
Path(out,'launcher_terminal.json').write_text(json.dumps(dict(exit_code=int(status),mode=mode,
    start_epoch=int(start),deadline_epoch=int(deadline),finished_epoch=time.time(),
    slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),indent=2)+'\n')
PY
  exit "$status"
}
trap terminal_receipt EXIT
trap 'exit 143' TERM
"$python" - "$out" "$remote" "$mode" "$start_epoch" "$deadline_epoch" "$seconds" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
out,remote,mode,start,deadline,seconds=sys.argv[1:]
Path(out,'launcher_start.json').write_text(json.dumps(dict(mode=mode,start_epoch=int(start),deadline_epoch=int(deadline),
    slurm_job_id=os.getenv('SLURM_JOB_ID'),cpus=8,numba_threads=8,blas_threads=1,memory_GiB=96,wall_seconds=int(seconds),
    inventory_sha256=hashlib.sha256((Path(remote)/'inventory.json').read_bytes()).hexdigest()),indent=2)+'\n')
PY
(while true; do date -u +%Y-%m-%dT%H:%M:%SZ > "$out/launcher_heartbeat.txt"; sleep 300; done) >/dev/null 2>&1 < /dev/null &
heartbeat_pid=$!
"$python" "$remote/deploy.py" verify --stage "$remote" > "$out/host_verification.json"
"$python" - "$base" "$floor_remote" "$inputs" <<'PY'
import hashlib,json,sys
from pathlib import Path
for folder in map(Path,sys.argv[1:3]):
    for rel,digest in json.loads((folder/'inventory.json').read_text())['files'].items():
        assert hashlib.sha256((folder/'source'/rel).read_bytes()).hexdigest()==digest,rel
root=Path(sys.argv[3])
for name,digest in {'bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7',
 'arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0'}.items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
PY
if [[ "$mode" == fit ]]; then
 "$python" "$remote/deploy.py" verify-gate --stage "$remote" --gate "$remote/fit_review_gate.json" > "$out/review_gate_verification.json"
fi
binds=(--bind "$frozen:$repo:ro" --bind "$remote:/work/deployment:ro")
"$python" "$remote/deploy.py" mounts --stage "$remote" --base "$base" --floor "$floor_remote" --repo "$repo" > "$out/bind_manifest.json"
while IFS= read -r binding; do binds+=(--bind "$binding"); done < <("$python" - "$out/bind_manifest.json" <<'PYBIND'
import json,sys
for binding in json.load(open(sys.argv[1]))['bindings']:print(binding)
PYBIND
)
# Avoid exporting inherited duplicate bind lists into Apptainer's child exec.
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH
# Fixed input directory wins after any broader staged output-directory mount.
"$python" - "$remote/inventory.json" "$inputs" <<'PYINPUT'
import hashlib,json,sys
from pathlib import Path
own=json.load(open(sys.argv[1]))['files'];prefix='output/model/publication_refactor_20260929/local_export_v1/inputs/'
for rel,digest in own.items():
 if rel.startswith(prefix):assert hashlib.sha256((Path(sys.argv[2])/rel[len(prefix):]).read_bytes()).hexdigest()==digest,'Fixed-input pin conflict: '+rel
PYINPUT
binds+=(--bind "$inputs:$repo/output/model/publication_refactor_20260929/local_export_v1/inputs:ro")
binds+=(--bind "$out:/work/results:rw")
export NUMBA_CACHE_DIR=/work/results/numba_cache MPLCONFIGDIR=/work/results/matplotlib
apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" /work/deployment/deploy.py verify \
 --stage /work/deployment --mounted-root "$repo" > "$out/container_verification.json"
plan_mode=$mode; [[ "$mode" == preflight ]] && plan_mode=smoke
plan_rel=$("$python" - "$remote/inventory.json" "$plan_mode" "$seconds" "$mode" <<'PY'
import json,sys
from pathlib import Path
inv=json.load(open(sys.argv[1]));rel=inv['plans'][sys.argv[2]]['path']
plan=json.load(open(Path(sys.argv[1]).parent/'source'/rel))
if sys.argv[4]!='preflight':assert plan['budget']['total_seconds']+120<=int(sys.argv[3]),'Plan exceeds external wall budget'
print(rel)
PY
)
args=(--plan "$repo/$plan_rel" --output /work/results/run)
case "$mode" in preflight) args+=(--preflight);; smoke) args+=(--native-smoke);; fit) args+=(--execute);; esac
remaining=$((deadline_epoch-$(date +%s)-15)); [[ "$remaining" -gt 0 ]] || exit 124
timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$repo/code/model/experiments/birth_count_choice/transition.py" "${args[@]}" > "$out/driver.log" 2>&1
if [[ "$mode" == preflight ]]; then
 remaining=$((deadline_epoch-$(date +%s)-15)); [[ "$remaining" -gt 0 ]] || exit 124
 timeout --signal=TERM --kill-after=10s "${remaining}s" apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
  - "$repo/$plan_rel" "$repo/code/model" /work/results/native_constructor <<'PYNATIVE' > "$out/native_constructor.log" 2>&1
import json,sys
from pathlib import Path
sys.path.insert(0,sys.argv[2])
from experiments.birth_count_choice.transition_runtime import CurrentEstateARuntime
from numba import get_num_threads
assert get_num_threads()==8,'Native preflight requires eight Numba threads'
plan=json.load(open(sys.argv[1]));folder=Path(sys.argv[3])
rt=CurrentEstateARuntime.from_handoff(plan['handoff'],folder)
assert rt.total_native_calls==0,'Zero-solve constructor unexpectedly made native calls'
assert rt.identity()==plan['identity'],'Constructor runtime identity differs from pinned plan'
receipt=dict(status='PASS_ZERO_SOLVES_NATIVE_CONSTRUCTOR',policy_calls=0,bellman_calls=0,kfe_calls=0,numba_threads=get_num_threads(),
    identity=rt.identity(),scientific_validation=False,reference_reconstruction_verified=False)
(folder/'native_status.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(json.dumps(receipt,indent=2))
PYNATIVE
fi
