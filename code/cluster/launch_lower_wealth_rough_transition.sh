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
#SBATCH --output=/scratch/td2248/projects/estate_a_lower_wealth_rough_transition_20261005_v4/logs/%x-%j.out
set -euo pipefail

repo=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
local_v5="$repo/output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/execution_smoke_v5"
local_addon="$repo/output/model/transition_readiness_v1/estate_a_lower_wealth_20261005/rough_two_shock_v1/deployment_v4"
local_v3="$repo/output/model/transition_readiness_v1/estate_a_lower_wealth_20261005/rough_two_shock_v1/deployment_v3"
remote_v5=/scratch/td2248/projects/current_estate_two_shock_20261004_v5
remote_v3=/scratch/td2248/projects/estate_a_lower_wealth_rough_transition_20261005_v3
remote_addon=/scratch/td2248/projects/estate_a_lower_wealth_rough_transition_20261005_v4
python=/share/apps/anaconda3/2025.06/bin/python
image=/share/apps/images/ubuntu-24.04.4.sif
task=estate_a_lower_wealth_rough_h16_v4
prior_budget_charge_upper_bound=103
original_deadline_epoch=1791238908

[[ "$#" == 1 && "$1" =~ ^[[:xdigit:]]{64}$ ]] || { echo 'Pass exactly the pinned package.json SHA256' >&2; exit 2; }
expected_package_sha=$1
[[ -n "${SLURM_JOB_ID:-}" && -z "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo 'One non-array Slurm job required' >&2; exit 2; }
[[ "${SLURM_NTASKS:-1}" == 1 && "${SLURM_CPUS_PER_TASK:-1}" == 1 ]] || { echo 'Exactly one task and CPU required' >&2; exit 2; }
[[ "${SLURM_MEM_PER_NODE:-24576}" =~ ^[0-9]+$ && "${SLURM_MEM_PER_NODE:-24576}" -le 24576 ]] || { echo '24 GiB memory allocation required' >&2; exit 2; }
[[ -d "$remote_v5/frozen/source" && -f "$remote_v3/results/estate_a_lower_wealth_rough_h16_v3/base_plan.json" && -d "$remote_addon/files" ]] || { echo 'Pinned source, v3 base plan or add-on missing' >&2; exit 2; }

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
if [[ "$deadline_epoch" -gt "$original_deadline_epoch" ]]; then deadline_epoch=$original_deadline_epoch; fi
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
 no_auto_retry=True,task_id='estate_a_lower_wealth_rough_h16_v4',
 prior_attempt_job_ids=['19238662','19239226','19239878'],prior_budget_charge_upper_bound=103,
 prior_actual_native_calls_lower_bound=32,
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

"$python" - "$remote_addon/package.json" "$expected_package_sha" "$remote_addon/files" "$remote_v3" > "$out/package_verification.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
manifest,expected,root,v3=Path(sys.argv[1]),sys.argv[2],Path(sys.argv[3]),Path(sys.argv[4])
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
if sha(manifest)!=expected:raise SystemExit('Submitted package pin differs')
package=json.loads(manifest.read_text())
required={'rough_two_shock.py','rough_two_shock_runtime.py'}
if package.get('schema')!='estate_a_lower_wealth_rough_fit_only_package_v4' or set(package['files'])!=required:
 raise SystemExit('Wrong package file set')
for name,digest in package['files'].items():
 if sha(root/name)!=digest:raise SystemExit('Add-on file bytes differ: '+name)
case=v3/'results/estate_a_lower_wealth_rough_h16_v3'
inputs={'v3_base_plan':case/'base_plan.json','v3_handoff':case/'base_plan_handoff.json',
 'v3_bridge_receipt':case/'bridge_case/bridge_receipt.json',
 'v3_input_contract':case/'bridge_case/input_contract.json','v3_metadata':case/'bridge_case/metadata.json'}
if set(package['inputs'])!=set(inputs):raise SystemExit('Package input pin set differs')
for name,path in inputs.items():
 if sha(path)!=package['inputs'][name]:raise SystemExit('Pinned input bytes differ: '+name)
base=json.loads(inputs['v3_base_plan'].read_text())
if base['handoff']['sha256']!=package['inputs']['v3_handoff'] or base['initial_psi']!=0.176948250201189:
 raise SystemExit('Retained v3 reference identity differs')
print(json.dumps(dict(verified=True,package_sha256=expected,files=package['files']),sort_keys=True))
PY
"$python" - "$job/launcher_start.json" "$start_epoch" "$deadline_epoch" "$expected_package_sha" <<'PY'
import json,os,sys
from pathlib import Path
path,start,deadline,pin=sys.argv[1:]
Path(path).write_text(json.dumps(dict(task_id='estate_a_lower_wealth_rough_h16_v4',start_epoch=int(start),
 deadline_epoch=int(deadline),external_wall_seconds=14400,internal_wall_seconds=14280,cpus=1,memory_gib=24,
 numba_threads=1,blas_threads=1,package_sha256=pin,slurm_job_id=os.getenv('SLURM_JOB_ID'),
 no_auto_retry=True,rough_diagnostic=True,fit_only=True,
 prior_attempt_job_ids=['19238662','19239226','19239878'],
 prior_budget_charge_upper_bound=103,prior_actual_native_calls_lower_bound=32,
 original_shared_deadline_epoch=1791238908),sort_keys=True,indent=2)+'\n')
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

# Reuse the authenticated v3 bridge case and base plan read-only. The fit
# rebuilds its own reference, seeds, and candidate evaluations from that case.
binds=(--bind "$remote_v5/frozen/source:$repo:ro" --bind "$remote_v5:$local_v5:ro"
 --bind "$remote_v3:$local_v3:ro" --bind "$remote_addon:$local_addon:ro"
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
v3_base_plan="$local_v3/results/estate_a_lower_wealth_rough_h16_v3/base_plan.json"
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --prepare --base-plan "$v3_base_plan" \
 --deadline-epoch "$deadline_epoch" --output "$visible/fit_manifest.json" > "$out/fit_prepare.log" 2>&1
run_bounded 120 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --manifest "$visible/fit_manifest.json" \
 --output "$visible/fit_preflight" --preflight > "$out/fit_preflight.json"
run_bounded 14280 apptainer exec "${binds[@]}" --pwd "$repo" "$image" "$python" \
 "$local_addon/files/rough_two_shock.py" --manifest "$visible/fit_manifest.json" \
 --output "$visible/fit" --run > "$out/fit.log" 2>&1
"$python" - "$out/fit/complete.json" "$out/fit_manifest.json" \
 "$remote_addon/package.json" > "$out/final_common_ledger.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
fit,manifest,package=[json.loads(Path(p).read_text()) for p in sys.argv[1:]]
if fit.get('status')!='rough_matched' or fit.get('source_identity')!=manifest.get('identity'):
 raise SystemExit('Rough fit completion/source identity differs')
if fit.get('actual_policy_calls')!=fit.get('native_actual_policy_calls'):
 raise SystemExit('Fit controller/native calls differ')
for key,file in (('rough_driver','rough_two_shock.py'),('rough_runtime','rough_two_shock_runtime.py')):
 if manifest['source_pins'][key]['sha256']!=package['files'][file]:
  raise SystemExit('Final rough source pin differs: '+key)
fit_calls=int(fit['native_actual_policy_calls'])
charged=103+fit_calls
if fit_calls>19897 or charged>20000:raise SystemExit('Shared native-call budget charge exceeded')
print(json.dumps(dict(status='PASS',prior_job_ids=['19238662','19239226','19239878'],
 prior_budget_charge_upper_bound=103,prior_actual_native_calls_lower_bound=32,
 prior_actual_final_unknown=True,fit_actual_calls=fit_calls,aggregate_budget_charge=charged,
 aggregate_cap=20000,rough_source_pins=manifest['source_pins']),sort_keys=True))
PY
