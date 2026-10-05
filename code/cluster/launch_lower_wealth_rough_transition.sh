#!/usr/bin/env bash
#SBATCH --job-name=estate_lower_rough
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cs
#SBATCH --output=/scratch/td2248/projects/estate_a_lower_wealth_rough_transition_20261005_v1/logs/%x-%j.out
set -euo pipefail

repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
local_v5="$repo/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5"
local_addon="$repo/output/model/transition_readiness_v1/estate_a_lower_wealth_20261005/rough_two_shock_v1/deployment"
remote_v5=/scratch/td2248/projects/current_estate_two_shock_20261004_v5
remote_overnight=/scratch/td2248/projects/estate_birth_overnight_20261004_v1
remote_addon=/scratch/td2248/projects/estate_a_lower_wealth_rough_transition_20261005_v1
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
task=estate_a_lower_wealth_rough_h16_v1

[[ "$#" == 1 && "$1" =~ ^[[:xdigit:]]{64}$ ]] || { echo 'Pass exactly the pinned package.json SHA256' >&2; exit 2; }
expected_package_sha=$1
[[ -n "${SLURM_JOB_ID:-}" && -z "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'One non-array Slurm job required' >&2; exit 2; }
[[ "${SLURM_NTASKS:-1}" == 1 && "${SLURM_CPUS_PER_TASK:-1}" == 1 ]] || { echo 'Exactly one task and CPU required' >&2; exit 2; }
[[ "${SLURM_MEM_PER_NODE:-24576}" =~ ^[0-9]+$ && "${SLURM_MEM_PER_NODE:-24576}" -le 24576 ]] || { echo '24 GiB memory allocation required' >&2; exit 2; }
[[ -d "$remote_v5/frozen/source" && -f "$remote_overnight/inventory.json" && -d "$remote_addon/files" ]] || { echo 'Pinned source or add-on missing' >&2; exit 2; }

module load anaconda3/2025.06
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
export PYTHONPATH="$repo/code/model/experiments/birth_count_choice:$repo/code/model:$local_addon/files"
export APPTAINERENV_PYTHONPATH="$PYTHONPATH"
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH

# This persistent claim is never removed, including on a failed run.
claim="$remote_addon/jobs/$task.claim"
mkdir "$claim" || { echo "Refusing duplicate or unknown submission: $claim" >&2; exit 2; }
printf '%s\n' "$SLURM_JOB_ID" > "$claim/slurm_job_id"
printf '%s\n' "$expected_package_sha" > "$claim/package_sha256"
job="$remote_addon/jobs/${task}_${SLURM_JOB_ID}"
out="$remote_addon/results/$task"
mkdir "$job" "$out" "$out/numba_cache" "$out/matplotlib" || exit 2
visible="$local_addon/results/$task"
export NUMBA_CACHE_DIR="$visible/numba_cache" MPLCONFIGDIR="$visible/matplotlib"
start_epoch=$(date +%s)
deadline_epoch=$((start_epoch+14280))
child_pid=''
heartbeat_pid=''
disk_guard_pid=''

stop_child() {
 if [[ -n "$child_pid" ]] && kill -0 "$child_pid" 2>/dev/null; then
  kill -TERM -- "-$child_pid" 2>/dev/null || true
  for _ in {1..10}; do kill -0 "$child_pid" 2>/dev/null || break; sleep 1; done
  kill -KILL -- "-$child_pid" 2>/dev/null || true
  wait "$child_pid" 2>/dev/null || true
 fi
}
terminal() {
 local status=$?
 trap - EXIT TERM INT
 stop_child
 if [[ -n "$heartbeat_pid" ]]; then kill "$heartbeat_pid" 2>/dev/null || true; wait "$heartbeat_pid" 2>/dev/null || true; fi
 if [[ -n "$disk_guard_pid" ]]; then kill "$disk_guard_pid" 2>/dev/null || true; wait "$disk_guard_pid" 2>/dev/null || true; fi
 "$python" - "$job/launcher_terminal.json" "$status" "$start_epoch" "$deadline_epoch" <<'PY'
import json,os,sys,time
from pathlib import Path
path,status,start,deadline=sys.argv[1:]
Path(path).write_text(json.dumps(dict(status='complete' if int(status)==0 else 'failed',exit_code=int(status),
 start_epoch=int(start),deadline_epoch=int(deadline),finished_epoch=time.time(),slurm_job_id=os.getenv('SLURM_JOB_ID'),
 no_auto_retry=True,task_id='estate_a_lower_wealth_rough_h16_v1',
 failure_reason='own_output_cap_exceeded' if Path(path).parent.joinpath('disk_cap_exceeded.json').exists() else None),
 sort_keys=True,indent=2)+'\n')
PY
 cp "$job/launcher_terminal.json" "$out/launcher_terminal.json"
 exit "$status"
}
trap terminal EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

# A stage cap always consumes the same global clock; none resets the deadline.
run_bounded() {
 local cap=$1; shift
 local remaining=$((deadline_epoch-$(date +%s)-60))
 [[ "$remaining" -gt 0 ]] || return 124
 if [[ "$cap" -gt "$remaining" ]]; then cap=$remaining; fi
 setsid timeout --signal=TERM --kill-after=10s "${cap}s" "$@" <&0 & child_pid=$!
 local status=0
 wait "$child_pid" || status=$?
 child_pid=''
 return "$status"
}

"$python" - "$remote_addon/package.json" "$expected_package_sha" "$remote_addon/files" "$remote_v5" "$remote_overnight" > "$out/package_verification.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
manifest,expected,root,v5,overnight=Path(sys.argv[1]),sys.argv[2],Path(sys.argv[3]),Path(sys.argv[4]),Path(sys.argv[5])
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
if sha(manifest)!=expected:raise SystemExit('Submitted package pin differs')
package=json.loads(manifest.read_text())
required={'lower_wealth_transition_reference.py','rough_two_shock.py','rough_two_shock_runtime.py'}
if package.get('schema')!='estate_a_lower_wealth_rough_package_v1' or set(package['files'])!=required:
 raise SystemExit('Wrong package file set')
for name,digest in package['files'].items():
 if sha(root/name)!=digest:raise SystemExit('Add-on file bytes differ: '+name)
inputs={'v5_inventory':v5/'inventory.json','v5_fit_manifest':v5/'inputs/fit_manifest.json',
 'overnight_inventory':overnight/'inventory.json',
 'selected_json':overnight/'results/production_binary_chain_11/run/0149_nm/phase_b_ge/selected.json',
 'selected_repeat_arrays':overnight/'results/production_binary_chain_11/run/0149_nm/phase_b_ge/selected_repeat/stage/solution_arrays.npz'}
if set(package['inputs'])!=set(inputs):raise SystemExit('Package input pin set differs')
for name,path in inputs.items():
 if sha(path)!=package['inputs'][name]:raise SystemExit('Pinned input bytes differ: '+name)
print(json.dumps(dict(verified=True,package_sha256=expected,files=package['files']),sort_keys=True))
PY
"$python" - "$job/launcher_start.json" "$start_epoch" "$deadline_epoch" "$expected_package_sha" <<'PY'
import json,os,sys
from pathlib import Path
path,start,deadline,pin=sys.argv[1:]
Path(path).write_text(json.dumps(dict(task_id='estate_a_lower_wealth_rough_h16_v1',start_epoch=int(start),
 deadline_epoch=int(deadline),external_wall_seconds=14400,internal_wall_seconds=14280,cpus=1,memory_gib=24,
 numba_threads=1,blas_threads=1,package_sha256=pin,slurm_job_id=os.getenv('SLURM_JOB_ID'),
 no_auto_retry=True,rough_diagnostic=True),sort_keys=True,indent=2)+'\n')
PY
cp "$job/launcher_start.json" "$out/launcher_start.json"
timeout 20s myquota > "$job/initial_quota.txt"
"$python" - "$job/initial_quota.txt" > "$job/quota_gate.json" <<'PY'
import json,re,sys
from pathlib import Path
line=next((x for x in re.sub(r'\x1b\[[0-9;]*m','',Path(sys.argv[1]).read_text()).splitlines()
           if x.startswith('/scratch ')),None)
if line is None:raise SystemExit('Scratch quota line missing')
limit=re.search(r'(\d+(?:\.\d+)?)TB/',line)
usage=re.search(r'\d+(?:\.\d+)?TB\((\d+(?:\.\d+)?)%\)',line)
if not limit or not usage:raise SystemExit('Scratch quota format changed')
headroom_gb=float(limit.group(1))*1000*(1-float(usage.group(1))/100)
if headroom_gb<34.4:raise SystemExit('Less than 32 GiB scratch quota headroom')
print(json.dumps(dict(status='PASS',headroom_gb_from_report=headroom_gb,minimum_gb=34.4),sort_keys=True))
PY
(while true; do date -u +%Y-%m-%dT%H:%M:%SZ > "$job/launcher_heartbeat.txt"; sleep 60; done) >/dev/null 2>&1 < /dev/null &
heartbeat_pid=$!
launcher_pid=$$
(while true; do
 used_kib=$(du -sk "$out" 2>/dev/null | awk '{print $1}') || used_kib=''
 [[ "$used_kib" =~ ^[0-9]+$ ]] || { sleep 60; continue; }
 if [[ "$used_kib" -gt 16777216 ]]; then
  printf '{"reason":"own_output_cap_exceeded","limit_gib":16,"observed_kib":%s}\n' "$used_kib" > "$job/disk_cap_exceeded.json"
  kill -TERM "$launcher_pid"
  exit 0
 fi
 sleep 60
done) >/dev/null 2>&1 < /dev/null &
disk_guard_pid=$!

# The original v5 package and overnight case remain read-only. The isolated
# add-on holds only new programs, control receipts and outputs.
binds=(--bind "$remote_v5/frozen/source:$repo:ro" --bind "$remote_v5:$local_v5:ro"
 --bind "$remote_overnight:$remote_overnight:ro" --bind "$remote_addon:$local_addon:ro"
 --bind "$remote_addon/results:$local_addon/results:rw")
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" - \
 "$local_addon/package.json" "$expected_package_sha" "$local_addon/files" <<'PY' > "$out/package_container_verification.json"
import hashlib,json,sys
from pathlib import Path
manifest,expected,root=Path(sys.argv[1]),sys.argv[2],Path(sys.argv[3])
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
if sha(manifest)!=expected:raise SystemExit('Container package pin differs')
package=json.loads(manifest.read_text())
for name,digest in package['files'].items():
 if sha(root/name)!=digest:raise SystemExit('Container add-on bytes differ: '+name)
print(json.dumps(dict(verified=True,package_sha256=expected),sort_keys=True))
PY
# The retained packages were already authenticated. Recheck only the exact
# inputs above and let the bridge/rough native preflights authenticate sources
# actually used by this new run.
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/lower_wealth_transition_reference.py" --run-root "$remote_overnight/results/production_binary_chain_11/run" \
 --inventory "$remote_overnight/inventory.json" --source-root "$repo" --output-case "$visible/bridge_case" \
 --budget-seconds 1800 --preflight-only > "$out/bridge_preflight.json"
run_bounded 1800 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/lower_wealth_transition_reference.py" --run-root "$remote_overnight/results/production_binary_chain_11/run" \
 --inventory "$remote_overnight/inventory.json" --source-root "$repo" --output-case "$visible/bridge_case" \
 --budget-seconds 1800 > "$out/bridge.log" 2>&1

run_bounded 300 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --prepare-base-plan --bridge-case "$visible/bridge_case" \
 --source-template-manifest "$local_v5/inputs/fit_manifest.json" --output "$visible/base_plan.json" \
 > "$out/base_plan.log" 2>&1
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --prepare-smoke --base-plan "$visible/base_plan.json" \
 --deadline-epoch "$deadline_epoch" --output "$visible/smoke_manifest.json" > "$out/smoke_prepare.log" 2>&1
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --manifest "$visible/smoke_manifest.json" \
 --output "$visible/smoke_preflight" --preflight > "$out/smoke_preflight.json"
run_bounded 900 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --manifest "$visible/smoke_manifest.json" \
 --output "$visible/smoke" --smoke > "$out/smoke.log" 2>&1

"$python" - "$out/bridge_case/bridge_receipt.json" "$out/smoke/smoke_receipt.json" "$out/base_plan.json" "$out/smoke_manifest.json" > "$out/shared_call_check.json" <<'PY'
import json,sys
from pathlib import Path
bridge=json.loads(Path(sys.argv[1]).read_text());smoke=json.loads(Path(sys.argv[2]).read_text())
base=json.loads(Path(sys.argv[3]).read_text());manifest=json.loads(Path(sys.argv[4]).read_text())
if smoke.get('status')!='execution_passed' or smoke.get('schema')!=manifest.get('schema') or \
   smoke.get('source_identity')!=base['identity'] or manifest.get('identity')!=base['identity']:
 raise SystemExit('Exact same-source execution smoke did not pass')
b=int(bridge['lifecycle_solves']);s=int(smoke['actual_policy_calls'])
if b!=1 or not 0<=s<=100 or b+s+19899>20000:raise SystemExit('Common native-call cap exceeded')
print(json.dumps(dict(bridge=b,smoke=s,fit_cap=19899,aggregate_cap=20000,checked=True),sort_keys=True))
PY
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --prepare --base-plan "$visible/base_plan.json" \
 --deadline-epoch "$deadline_epoch" --output "$visible/fit_manifest.json" > "$out/fit_prepare.log" 2>&1
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --manifest "$visible/fit_manifest.json" \
 --output "$visible/fit_preflight" --preflight > "$out/fit_preflight.json"
run_bounded 14280 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --manifest "$visible/fit_manifest.json" \
 --output "$visible/fit" --run > "$out/fit.log" 2>&1
"$python" - "$out/bridge_case/bridge_receipt.json" "$out/smoke/smoke_receipt.json" \
 "$out/fit/complete.json" "$out/fit_manifest.json" "$remote_addon/package.json" > "$out/final_common_ledger.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
bridge,smoke,fit,manifest,package=[json.loads(Path(p).read_text()) for p in sys.argv[1:]]
if fit.get('status')!='rough_matched' or fit.get('source_identity')!=manifest.get('identity'):
 raise SystemExit('Rough fit completion/source identity differs')
if fit.get('actual_policy_calls')!=fit.get('native_actual_policy_calls'):
 raise SystemExit('Fit controller/native calls differ')
for key,file in (('rough_driver','rough_two_shock.py'),('rough_runtime','rough_two_shock_runtime.py')):
 if manifest['source_pins'][key]['sha256']!=package['files'][file]:
  raise SystemExit('Final rough source pin differs: '+key)
calls=int(bridge['lifecycle_solves'])+int(smoke['actual_policy_calls'])+int(fit['native_actual_policy_calls'])
if calls>20000:raise SystemExit('Shared actual native-call cap exceeded')
print(json.dumps(dict(status='PASS',bridge_calls=bridge['lifecycle_solves'],
 smoke_calls=smoke['actual_policy_calls'],fit_calls=fit['native_actual_policy_calls'],
 aggregate_actual_calls=calls,aggregate_cap=20000,rough_source_pins=manifest['source_pins']),sort_keys=True))
PY
