#!/usr/bin/env python3
"""Build a fresh authenticated two-shock package; never submit or overwrite a stage.

Torch mounts the package at its original absolute local path. Control pins and
result receipts therefore have one namespace on both machines.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE_REL = 'output/model/transition_readiness_v1/current_baseline_20261003'
EXPECTED_INVENTORY = '15cc5f6036d57b985ea2ea9098f015dd527d2ce11ec94bd458a49913f017d977'
EXPECTED_MANIFEST = '07d84336a3112b251afe505908113d9c00585b91c34bd50f0dee108435db496d'
EXPECTED_HELPER = 'ad4d8407b8a192f21a39c28c916f5ca36b565733ab5c85e6d31682d540d23ecd'
MANIFEST_REL = 'output/model/overnight_calibration_20260928/contract_v1/source_manifest.json'
EXPECTED_INPUTS = {
 'output/model/publication_refactor_20260929/local_export_v1/inputs/bundle.json':'427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7',
 'output/model/publication_refactor_20260929/local_export_v1/inputs/arrays.npz':'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0',
 'output/model/fertility_identification_20260928/fixed_reference_manifest.json':'147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4',
 'output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py':'96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44',
}
DRIVER_REL = 'code/model/experiments/birth_count_choice/two_shock.py'
RUNTIME_REL = 'code/model/experiments/birth_count_choice/two_shock_runtime.py'
HELPER_REL = 'code/cluster/estate_birth_transition/two_shock_source_overlay.py'
LOCAL_LAUNCHER_REL = 'code/cluster/estate_birth_transition/launch_continuation.py'
REMOTE_DEFAULT = '/scratch/td2248/projects/current_estate_two_shock_20261004_v1'
SCHEMA = 'authenticated_two_shock_torch_package_v1'


def require(ok, message):
    if not ok: raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def read(path): return json.loads(Path(path).read_text())
def dump(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
def pin(path): return dict(path=str(Path(path).resolve()),sha256=sha(path))

def safe_rel(value):
    p=Path(value)
    require(not p.is_absolute() and '..' not in p.parts and str(p) not in ('','.'),'Unsafe package path: '+str(p))
    return p


def authenticate_files(root, files):
    for rel,digest in sorted(files.items()):
        p=Path(root)/safe_rel(rel)
        require(p.is_file() and not p.is_symlink() and sha(p)==digest,'Authenticated file mismatch: '+str(p))


def copy_exact(root, destination, files):
    for rel,digest in sorted(files.items()):
        source=Path(root)/safe_rel(rel);target=Path(destination)/rel
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        require(sha(target)==digest,'Copied bytes mismatch: '+rel)
        target.chmod(0o444)


def relocate_plan(plan, canonical_root, source):
    """Relocate only declared path pins and authenticate every underlying content."""
    result=copy.deepcopy(plan)
    items=list(result['source_files'].values())+[result['handoff']]
    items += [result['target_contract'][key] for key in ('annual','blocks')]
    for item in items:
        rel=Path(item['path']).relative_to(canonical_root);safe_rel(rel)
        item['path']=str(Path(source)/rel)
        require(sha(item['path'])==item['sha256'],'Relocated plan content mismatch: '+str(rel))
    return result


def relocate_overlay(meta, manifest, destination, canonical_root):
    require(meta.get('schema')=='exact_read_only_source_overlay_v1','Wrong overlay schema')
    require(meta.get('manifest_sha256')==EXPECTED_MANIFEST and meta.get('all_source_hashes_verified') is True,'Wrong overlay authentication')
    require(meta.get('original_root')==str(canonical_root),'Wrong overlay original root')
    files=manifest['files'];require(len(files)==1241,'Wrong overlay file count')
    old=Path(meta['snapshot_root'])
    # The original helper recorded resolved paths: legacy workspace aliases
    # can name the exact same frozen file through an archive symlink. Check
    # that original namespace exactly; the fresh namespace is built separately.
    expected={str((Path(canonical_root)/safe_rel(rel)).resolve()):str((old/rel).resolve()) for rel in files}
    require(meta.get('mapping')==expected and meta.get('manifest_file_count')==len(files),'Overlay path mapping differs')
    result=copy.deepcopy(meta);snapshot=Path(destination)/'frozen/source_overlay/files'
    result.update(snapshot_root=str(snapshot),manifest_path=str(Path(destination)/'frozen/source'/MANIFEST_REL),
        mapping={str(Path(canonical_root)/rel):str(snapshot/rel) for rel in files})
    return result


def prepare_plans(destination, python):
    destination=Path(destination);source=destination/'frozen/source'
    # A fresh process imports the staged driver, so active concurrent edits can
    # never be hidden by a cached module in the packaging process.
    program='''import json,sys\nfrom pathlib import Path\nsys.path.insert(0,sys.argv[1])\nimport two_shock as d\nb=Path(sys.argv[2]);base=d.pin(b/'inputs/base_fit_plan.json')\noverlay=dict(helper=d.pin(b/'frozen/source/code/cluster/estate_birth_transition/two_shock_source_overlay.py'),manifest=d.pin(b/'frozen/source_overlay/overlay.json'))\nfor mode in ('smoke','fit'):\n p=d.prepare_manifest(base,smoke=mode=='smoke',legacy_source_overlay=overlay)\n d.validate_mode_controls(p)\n d.write(b/'inputs'/f'{mode}_manifest.json',p)\n'''
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',NUMBA_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    subprocess.run([str(python),'-c',program,str(source/Path(DRIVER_REL).parent),str(destination)],env=env,check=True,timeout=120)


def build(args):
    repo=Path(args.repo).resolve();base=repo/BASE_REL;dest=Path(args.destination).resolve()
    require(not dest.exists(),'Fresh destination required; refusing existing stage: '+str(dest))
    remote=Path(args.remote_root)
    require(remote.is_absolute() and ':' not in str(remote) and '\n' not in str(remote),'Absolute safe remote root required')
    require(str(remote)!=str(dest),'Remote storage must differ from container-visible local namespace')
    inventory_path=base/'deployment_v9/inventory.json'
    require(sha(inventory_path)==EXPECTED_INVENTORY,'Authoritative v9 inventory SHA mismatch')
    inv=read(inventory_path);canonical=Path(inv['root'])
    require(canonical==repo,'Canonical repository root differs; explicit original namespace required')
    require(len(inv['files'])==709,'Wrong native source count')
    frozen=base/'two_shock_v1/attempt3/frozen/source'
    authenticate_files(frozen,inv['files']);authenticate_files(frozen,EXPECTED_INPUTS)
    # Authenticate v9 metadata/launch ancestry without consulting mutable source.
    authenticate_files(base/'deployment_v9',inv['entrypoints'])
    manifest_path=repo/MANIFEST_REL
    require(sha(manifest_path)==EXPECTED_MANIFEST,'Nested original manifest SHA mismatch')
    manifest=read(manifest_path)
    overlay_base=base/'two_shock_v1/frozen/source_overlay';meta=read(overlay_base/'overlay.json')
    newmeta=relocate_overlay(meta,manifest,dest,canonical)
    require(Path(meta['snapshot_root']).resolve()==(overlay_base/'files').resolve(),'Overlay snapshot escapes frozen package')
    authenticate_files(overlay_base/'files',manifest['files'])
    require(sha(frozen/HELPER_REL)==EXPECTED_HELPER and sha(repo/HELPER_REL)==EXPECTED_HELPER,'Read-only overlay helper changed')
    additions={rel:sha(repo/rel) for rel in (DRIVER_REL,RUNTIME_REL,LOCAL_LAUNCHER_REL)}
    require(not (set(additions)&set(inv['files'])),'Reviewed additions would overwrite original native source')
    # All immutable inputs authenticate before creating any destination content.
    dest.mkdir(parents=True);source=dest/'frozen/source'
    copy_exact(frozen,source,inv['files']);copy_exact(frozen,source,EXPECTED_INPUTS)
    copy_exact(frozen,source,{HELPER_REL:EXPECTED_HELPER})
    copy_exact(repo,source,additions)
    copy_exact(repo,source,{MANIFEST_REL:EXPECTED_MANIFEST})
    copy_exact(overlay_base/'files',dest/'frozen/source_overlay/files',manifest['files'])
    dump(dest/'frozen/source_overlay/overlay.json',newmeta)
    copy_exact(base/'deployment_v9',dest/'provenance',{'inventory.json':EXPECTED_INVENTORY})
    (dest/'provenance/inventory.json').rename(dest/'provenance/v9_inventory.json')
    copy_exact(base/'deployment_v9',dest/'provenance/v9_entrypoints',inv['entrypoints'])
    for mode in ('fit','smoke'):
        rel=inv['plans'][mode]['path'];require(inv['files'][rel]==inv['plans'][mode]['sha256'],'Original plan inventory pin differs')
        dump(dest/'inputs'/f'base_{mode}_plan.json',relocate_plan(read(frozen/rel),canonical,source))
    prepare_plans(dest,args.python)
    # Verify that plans constructed by the fresh driver pin precisely its bytes.
    for mode in ('smoke','fit'):
        plan=read(dest/'inputs'/f'{mode}_manifest.json')
        require(plan['source_pins']['two_shock_driver']==pin(source/DRIVER_REL),'Fresh driver source pin differs')
        require(plan['source_pins']['two_shock_runtime']==pin(source/RUNTIME_REL),'Fresh runtime source pin differs')
        require(plan['identity']==inv['identity'],'Native identity changed')
    shutil.copyfile(Path(__file__),dest/'prepare_two_shock_torch.py')
    (dest/'jobs').mkdir();(dest/'results').mkdir()
    launcher=dest/'launch_torch.sh';launcher.write_text(launcher_text(dest,remote,canonical));launcher.chmod(0o555)
    files={str(p.relative_to(dest)):sha(p) for p in sorted(dest.rglob('*')) if p.is_file()}
    inventory=dict(schema=SCHEMA,local_root=str(dest),remote_root=str(remote),canonical_root=str(canonical),files=files,
        native_inventory_sha256=EXPECTED_INVENTORY,native_file_count=709,auxiliary_file_count=4,overlay_manifest_sha256=EXPECTED_MANIFEST,
        overlay_file_count=1241,identity=inv['identity'],plans={m:pin(dest/'inputs'/f'{m}_manifest.json') for m in ('smoke','fit')},
        resources=dict(cpus=1,numba_threads=1,blas_threads=1,memory_gib=24,maximum_wall_seconds=21600),
        no_model_solves=True,no_submission=True)
    dump(dest/'inventory.json',inventory)
    receipt=verify(dest)
    dump(dest/'stage_receipt.json',dict(status='prepared_no_submission',verification=receipt,inventory_sha256=sha(dest/'inventory.json')))
    for p in dest.rglob('*'):
        if p.is_file() and p.name not in ('launch_torch.sh',):p.chmod(0o444)
    return dict(destination=str(dest),remote_root=str(remote),inventory_sha256=sha(dest/'inventory.json'),**receipt)


def verify(stage):
    stage=Path(stage);inv=read(stage/'inventory.json')
    require(inv.get('schema')==SCHEMA,'Wrong package schema')
    authenticate_files(stage,inv['files'])
    original=read(stage/'provenance/v9_inventory.json')
    require(sha(stage/'provenance/v9_inventory.json')==EXPECTED_INVENTORY,'Original v9 inventory changed')
    authenticate_files(stage/'frozen/source',original['files'])
    authenticate_files(stage/'frozen/source',EXPECTED_INPUTS)
    require(sha(stage/'frozen/source'/MANIFEST_REL)==EXPECTED_MANIFEST,'Original overlay manifest changed')
    authenticate_files(stage/'frozen/source_overlay/files',read(stage/'frozen/source'/MANIFEST_REL)['files'])
    require(read(stage/'inputs/fit_manifest.json')['identity']==inv['identity'],'Fit native identity differs')
    require(read(stage/'inputs/smoke_manifest.json')['identity']==inv['identity'],'Smoke native identity differs')
    return dict(status='PASS_ZERO_SOLVES_PACKAGE',files=len(inv['files']),native_calls=0,model_solves=0,scientific_validation=False)


def launcher_text(local,remote,canonical):
    # Constants are shell-quoted; receipt JSON and plan pins are structured args.
    values={'LOCAL':shlex.quote(str(local)),'REMOTE':shlex.quote(str(remote)),'CANONICAL':shlex.quote(str(canonical)),
            'LOGROOT':str(remote/'logs')}
    return LAUNCHER.replace('@LOGROOT@',values['LOGROOT']).replace('@LOCAL@',values['LOCAL']).replace('@REMOTE@',values['REMOTE']).replace('@CANONICAL@',values['CANONICAL'])


LAUNCHER=r'''#!/usr/bin/env bash
#SBATCH --job-name=estate_two_shock
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --signal=B:TERM@30
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cl
#SBATCH --output=@LOGROOT@/%x-%j.out
set -euo pipefail
local_root=@LOCAL@
remote=@REMOTE@
repo=@CANONICAL@
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
'''


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build');b.add_argument('--repo',default=str(REPO));b.add_argument('--destination',default=str(REPO/BASE_REL/'two_shock_v1/execution_smoke_v1'))
    b.add_argument('--remote-root',default=REMOTE_DEFAULT);b.add_argument('--python',default=sys.executable)
    v=sub.add_parser('verify');v.add_argument('--stage',required=True)
    args=p.parse_args(argv);result=build(args) if args.command=='build' else verify(args.stage)
    print(json.dumps(result,sort_keys=True,indent=2));return result
if __name__=='__main__':main()
