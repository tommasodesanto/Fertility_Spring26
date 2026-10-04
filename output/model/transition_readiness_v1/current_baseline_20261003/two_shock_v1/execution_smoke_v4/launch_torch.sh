#!/usr/bin/env bash
#SBATCH --job-name=estate_two_shock
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cl
#SBATCH --output=/scratch/td2248/projects/current_estate_two_shock_20261004_v4/logs/%x-%j.out
set -euo pipefail
local_root=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v4
remote=/scratch/td2248/projects/current_estate_two_shock_20261004_v4
repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
mode=${1:-preflight}
name=${2:-${SLURM_JOB_ID:-}}
[[ "$mode" =~ ^(preflight|smoke|fit)$ ]] || { echo 'Invalid mode'; exit 2; }
[[ "$name" =~ ^[A-Za-z0-9_-]+$ ]] || { echo 'Explicit unique run name required'; exit 2; }
[[ -z "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'Arrays forbidden'; exit 2; }
[[ "${SLURM_CPUS_PER_TASK:-1}" == 1 ]] || { echo 'Exactly one CPU required'; exit 2; }
[[ "${SLURM_MEM_PER_NODE:-24576}" -le 24576 ]] || { echo '24GiB memory cap exceeded'; exit 2; }
case "$mode" in preflight) seconds=600;; smoke) seconds=1800;; fit) seconds=21600;; esac
if [[ "$mode" != preflight ]]; then [[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Smoke/fit require one Slurm job'; exit 2; }; fi
smoke_pin=${3:-}
if [[ "$mode" == fit ]]; then [[ -n "$smoke_pin" ]] || { echo 'Fit requires exact passed smoke pin JSON'; exit 2; }; fi
module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH
mkdir -p "$remote/jobs" "$remote/results"
if [[ "$mode" != preflight ]]; then
 claim="$remote/jobs/${mode}.claim"
 mkdir "$claim" || { echo "Refusing duplicate or unknown $mode outcome: $claim"; exit 2; }
 printf '%s\n' "$name" > "$claim/run_name"
 printf '%s\n' "${SLURM_JOB_ID:-}" > "$claim/slurm_job_id"
fi
job="$remote/jobs/${mode}_${name}"
mkdir "$job" || { echo "Refusing existing or unknown job: $job"; exit 2; }
out="$remote/results/${mode}_${name}"
mkdir "$out" || { echo "Refusing existing or unknown output: $out"; exit 2; }
mkdir "$out/numba_cache" "$out/matplotlib"
visible_out="$local_root/results/${mode}_${name}"
start_epoch=$(date +%s); deadline_epoch=$((start_epoch+seconds))
child_pid=''; heartbeat_pid=''
stop_child() {
 if [[ -n "$child_pid" ]] && kill -0 "$child_pid" 2>/dev/null; then
  kill -TERM -- "-$child_pid" 2>/dev/null || true
  for _ in {1..10}; do kill -0 "$child_pid" 2>/dev/null || break; sleep 1; done
  kill -KILL -- "-$child_pid" 2>/dev/null || true
  wait "$child_pid" 2>/dev/null || true
 fi
}
terminal() {
 status=$?; trap - EXIT TERM INT
 stop_child
 [[ -z "$heartbeat_pid" ]] || kill "$heartbeat_pid" 2>/dev/null || true
 "$python" - "$job/launcher_terminal.json" "$status" "$mode" "$start_epoch" "$deadline_epoch" <<'PYRECEIPT'
import json,os,sys,time
from pathlib import Path
path,status,mode,start,deadline=sys.argv[1:]
Path(path).write_text(json.dumps(dict(exit_code=int(status),mode=mode,start_epoch=int(start),deadline_epoch=int(deadline),finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),no_auto_retry=True),sort_keys=True,indent=2)+'\n')
PYRECEIPT
 cp "$job/launcher_terminal.json" "$out/launcher_terminal.json"
 exit "$status"
}
trap terminal EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
"$python" - "$job/launcher_start.json" "$remote/inventory.json" "$mode" "$seconds" "$start_epoch" "$deadline_epoch" "$visible_out" <<'PYRECEIPT'
import hashlib,json,os,sys
from pathlib import Path
path,inventory,mode,seconds,start,deadline,output=sys.argv[1:]
Path(path).write_text(json.dumps(dict(mode=mode,wall_seconds=int(seconds),start_epoch=int(start),deadline_epoch=int(deadline),output=output,cpus=1,numba_threads=1,blas_threads=1,memory_gib=24,slurm_job_id=os.getenv('SLURM_JOB_ID'),inventory_sha256=hashlib.sha256(Path(inventory).read_bytes()).hexdigest(),no_auto_retry=True),sort_keys=True,indent=2)+'\n')
PYRECEIPT
cp "$job/launcher_start.json" "$out/launcher_start.json"
(while true; do date -u +%Y-%m-%dT%H:%M:%SZ > "$job/launcher_heartbeat.txt"; sleep 60; done) >/dev/null 2>&1 < /dev/null &
heartbeat_pid=$!
run_bounded() {
 local remaining=$((deadline_epoch-$(date +%s)-15))
 [[ "$remaining" -gt 0 ]] || return 124
 setsid timeout --signal=TERM --kill-after=10s "${remaining}s" "$@" <&0 & child_pid=$!
 wait "$child_pid"; local status=$?; child_pid=''; return "$status"
}
run_bounded "$python" "$remote/prepare_two_shock_torch.py" verify --stage "$remote" > "$out/host_verification.json"
binds=(--bind "$remote/frozen/source:$repo:ro" --bind "$remote:$local_root:ro"
 --bind "$remote/jobs:$local_root/jobs:rw" --bind "$remote/results:$local_root/results:rw")
# Remote stage creation must include empty jobs/results mount points.
export NUMBA_CACHE_DIR="$visible_out/numba_cache" MPLCONFIGDIR="$visible_out/matplotlib"
if [[ "$mode" == fit ]]; then
 # Only evidence inside this package's stable results namespace is exposed.
 "$python" - "$smoke_pin" "$local_root/results" <<'PYPIN'
import json,sys
from pathlib import Path
p=json.loads(sys.argv[1]);assert set(p)=={'path','sha256'}
assert Path(p['path']).is_absolute() and Path(p['path']).is_relative_to(Path(sys.argv[2])),'Smoke receipt must be collected into stable package results namespace'
PYPIN
fi

run_bounded apptainer exec "${binds[@]}" --pwd "$local_root" "$image" "$python" \
 "$local_root/prepare_two_shock_torch.py" verify --stage "$local_root" > "$out/container_verification.json"
manifest="$local_root/inputs/${mode}_manifest.json"
[[ "$mode" != preflight ]] || manifest="$local_root/inputs/smoke_manifest.json"
driver="$local_root/frozen/source/code/model/experiments/birth_count_choice/two_shock.py"
args=(--manifest "$manifest" --output "$visible_out/run")
case "$mode" in preflight) args+=(--preflight);; smoke) args+=(--smoke);; fit) args+=(--run --smoke-receipt-pin "$smoke_pin");; esac
run_bounded apptainer exec "${binds[@]}" --pwd "$local_root" "$image" "$python" "$driver" "${args[@]}" > "$out/driver.log" 2>&1
if [[ "$mode" == preflight ]]; then
 run_bounded apptainer exec "${binds[@]}" --pwd "$local_root" "$image" "$python" - "$driver" "$manifest" "$visible_out/native_constructor" <<'PYNATIVE' > "$out/native_constructor.log" 2>&1
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(sys.argv[1]).parent))
import two_shock as d
import two_shock_runtime as r
from numba import get_num_threads
assert get_num_threads()==1
plan=json.load(open(sys.argv[2]));d.preflight(plan)
out=Path(sys.argv[3]);runtime=r.build_runtime(plan=plan,output=out,smoke=True)
assert runtime.rt.total_native_calls==0
assert runtime.rt.identity()==plan['identity']
d.write(out/'native_status.json',dict(status='PASS_ZERO_SOLVES_NATIVE_CONSTRUCTOR',native_calls=0,numba_threads=1,identity=runtime.rt.identity(),scientific_validation=False))
PYNATIVE
fi
