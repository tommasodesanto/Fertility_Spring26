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
import re
import signal
import time
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


REFERENCE_REL = 'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
PAIR_LOCK_SHA = '6c5e14b40eba0ca63911f5a4514d3a63832259f2f36458b05c57d801eb0dadad'
CASE_CHECKPOINT_SHA = '090c9ebda662bf7837c4f4cf1d816159bc9d203a9babe7c70d00d0c4be1e575e'
REVIEW_ARCHIVE_SHA = '7bc2ddf44e15751aa1b3790783ec4e9c56a96d431ee3e190cd6336b68caf5e5f'
REVIEW_MANIFEST_SHA = 'c6a2aa6909e4ef7afe8610a2ee3a70e92a4ddf6b1c5333bdba1eb5d323fe319d'
CASE_LEAVES = {'parameters.csv':'48e01f453c90e10f0527afe05b1b751a67fb6b056e69d1c66ed57696f93fa92c',
               'receipt.json':'c15cdd5862106cf95c392d4e2d75adfc786625b5108f5cc3c660b4c3cc1d76de'}


def reference_input_pins(manifest, canonical_root, read_pin):
    """The explicit zero-constructor authentication graph, without imports/solves.

    Follow only pins consumed by fixed_price.authenticate, evening.setup,
    current_transition.setup and pair_runtime/read_contract. In particular the
    base contract's obsolete /scratch pins are not used by pair_runtime: its
    reviewed relocated reference_root and PARENT_LOCK define that namespace.
    """
    root=Path(canonical_root);files={};provenance={}
    def add(item,label):
        require(set(item)=={'path','sha256'},'Exact reference dependency pin required: '+label)
        rel=str(safe_rel(Path(item['path']).relative_to(root)))
        require(rel not in files or files[rel]==item['sha256'],'Conflicting reference input hashes: '+rel)
        files[rel]=item['sha256'];provenance.setdefault(rel,[]).append(label)
        return item
    def load(item,label):return read_pin(add(item,label))
    controls={k:load(manifest[k],'fixed_reference_manifest.'+k) for k in
        ('contract','objective','source_manifest','source_contract','native_ancestry_contract')}
    contract=controls['contract'];native=controls['native_ancestry_contract']
    for name,item in contract['files'].items():add(item,'reference contract.files.'+name)
    for name,digest in manifest['artifact_hashes'].items():
        add(dict(path=str(Path(manifest['local_export'])/safe_rel(name)),sha256=digest),'reference export artifact '+name)
    add(dict(path=str(Path(manifest['local_export'])/'initial_state.pkl.gz'),sha256=manifest['checkpoint']['sha256']),'reference checkpoint')
    for name,item in native['files'].items():add(item,'native ancestry.files.'+name)
    base=load(native['base_contract'],'native ancestry.base_contract')
    load(native['objective'],'native ancestry.objective')
    native_inventory=load(native['source_manifest'],'native ancestry.source_manifest')
    for rel,digest in native_inventory['files'].items():
        add(dict(path=str(Path(native['source_root'])/safe_rel(rel)),sha256=digest),'native source inventory '+rel)
    case=Path(contract['reference_case'])
    add(dict(path=str(case/'initial_state.pkl.gz'),sha256=CASE_CHECKPOINT_SHA),'pinned current_transition_runtime.CHECKPOINT_SHA')
    for name,digest in CASE_LEAVES.items():
        add(dict(path=str(case/name),sha256=digest),'reviewed ZIP '+REVIEW_ARCHIVE_SHA+' / SOURCE_MANIFEST '+REVIEW_MANIFEST_SHA)
    pair=Path(base['reference_root'])
    lock=load(dict(path=str(pair/'inputs/launch_lock.json'),sha256=PAIR_LOCK_SHA),'pinned recovery_runner.PARENT_LOCK')
    for rel,digest in lock['runtime_file_sha256'].items():
        add(dict(path=str(pair/safe_rel(rel)),sha256=digest),'pair lock.runtime_file_sha256 '+rel)
    for name,key in [('objective.json','objective_sha256'),('proposal_bank.json','proposal_bank_sha256'),('source_manifest.json','source_manifest_sha256')]:
        value=load(dict(path=str(pair/'inputs'/name),sha256=lock[key]),'pair lock.'+key)
        if name=='source_manifest.json':
            require(value['source_root']==str(pair/'source'),'Pair source inventory root differs')
            for rel,digest in value['files'].items():
                add(dict(path=str(pair/'source'/safe_rel(rel)),sha256=digest),'pair source inventory '+rel)
    add(dict(path=str(pair/'ancestor_commute.py'),sha256=lock['ancestor_sha256']),'pair lock.ancestor_sha256')
    tax=pair.parent/'paygo_tax_comparison_20260924/run_paygo_two_rate.py'
    add(dict(path=str(tax),sha256=lock['tax_driver_sha256']),'pair lock.tax_driver_sha256')
    return files,provenance


def authenticated_reference_sources(repo, frozen, overlay, original_inventory, overlay_inventory):
    """Resolve every required byte using its independent pin; never default it."""
    repo=Path(repo);candidates={}
    for directory,inventory in ((Path(frozen),original_inventory),(Path(overlay),overlay_inventory)):
        for rel,digest in inventory.items():candidates.setdefault(digest,[]).append(directory/safe_rel(rel))
    chosen={}
    def source(item):
        rel=str(Path(item['path']).relative_to(repo));digest=item['sha256']
        for path in [Path(frozen)/rel,Path(overlay)/rel,*candidates.get(digest,[]),repo/rel]:
            if path.is_file() and sha(path)==digest:
                chosen[rel]=path;return path
        raise ValueError('No exact authenticated reference input bytes: '+rel+' SHA256='+digest)
    manifest=read(source(dict(path=str(repo/REFERENCE_REL),sha256=EXPECTED_INPUTS[REFERENCE_REL])))
    files,provenance=reference_input_pins(manifest,repo,lambda item:read(source(item)))
    for rel,digest in files.items():source(dict(path=str(repo/rel),sha256=digest))
    return files,provenance,chosen


def materialize_overlay_sources(source_root, snapshot_root, files):
    """Make original filenames discoverable without replacing existing bytes.

    The installed read-only overlay remains authoritative for mapped reads and
    SourceFileLoader.get_code. Physical files only supply the missing filesystem
    entries needed by Python's module finder and source-presence checks.
    """
    source=Path(source_root);snapshot=Path(snapshot_root)
    authenticate_files(snapshot,files)
    records={}
    for rel,digest in sorted(files.items()):
        path=source/safe_rel(rel)
        if path.exists():
            require(path.is_file() and not path.is_symlink(),'Overlay discovery path is not an ordinary file: '+rel)
            physical=sha(path)
            classification='existing_identical' if physical==digest else 'existing_different_overlay_precedence'
        else:
            path.parent.mkdir(parents=True,exist_ok=True)
            # Exclusive creation prevents a concurrent addition from being
            # overwritten between the presence check and the physical copy.
            with path.open('xb') as output, (snapshot/rel).open('rb') as input_stream:
                shutil.copyfileobj(input_stream,output,1<<20)
            path.chmod(0o444);physical=sha(path)
            require(physical==digest,'Overlay discovery copy changed: '+rel)
            classification='materialized_missing'
        records[rel]=dict(classification=classification,physical_sha256=physical,overlay_sha256=digest,
            read_and_loader_precedence='unchanged exact read-only overlay')
    return dict(schema='exact_overlay_discovery_materialization_v1',files=records,
        counts={name:sum(item['classification']==name for item in records.values()) for name in
          ('materialized_missing','existing_identical','existing_different_overlay_precedence')},
        existing_files_overwritten=0)


def verify_overlay_materialization(source, files, receipt):
    require(receipt.get('schema')=='exact_overlay_discovery_materialization_v1' and
        receipt.get('existing_files_overwritten')==0,'Overlay discovery provenance differs')
    require(set(receipt['files'])==set(files),'Overlay discovery filename set differs')
    for rel,digest in files.items():
        record=receipt['files'][rel];classification=record['classification']
        require(classification in ('materialized_missing','existing_identical','existing_different_overlay_precedence'),'Unknown overlay namespace priority')
        require(record['overlay_sha256']==digest and sha(Path(source)/rel)==record['physical_sha256'],'Overlay discovery physical/source pin differs: '+rel)
        require(classification=='existing_different_overlay_precedence' or record['physical_sha256']==digest,'Overlay discovery bytes differ: '+rel)
        if classification=='existing_different_overlay_precedence':require(record['physical_sha256']!=digest,'Different-source precedence classified incorrectly')
        require(not ((Path(source)/rel).stat().st_mode & 0o222),'Overlay discovery source is writable: '+rel)
    expected={name:sum(item['classification']==name for item in receipt['files'].values()) for name in
        ('materialized_missing','existing_identical','existing_different_overlay_precedence')}
    require(receipt['counts']==expected,'Overlay discovery counts differ')


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
    reference_files,reference_provenance,reference_sources=authenticated_reference_sources(repo,frozen,overlay_base/'files',inv['files'],manifest['files'])
    original_pins={**inv['files'],**EXPECTED_INPUTS,HELPER_REL:EXPECTED_HELPER,MANIFEST_REL:EXPECTED_MANIFEST}
    for rel in set(reference_files)&set(original_pins):
        require(reference_files[rel]==original_pins[rel],'Reference input conflicts with immutable original pin: '+rel)
    review_manifest=repo/'output/model/review_bundle_20261003/Fertility_Model_Review_20261003/SOURCE_MANIFEST.sha256'
    require(sha(review_manifest)==REVIEW_MANIFEST_SHA,'Reviewed case-leaf authority manifest changed')
    additions={rel:sha(repo/rel) for rel in (DRIVER_REL,RUNTIME_REL,LOCAL_LAUNCHER_REL)}
    require(not (set(additions)&set(inv['files'])),'Reviewed additions would overwrite original native source')
    # All immutable inputs authenticate before creating any destination content.
    dest.mkdir(parents=True);source=dest/'frozen/source'
    copy_exact(frozen,source,inv['files']);copy_exact(frozen,source,EXPECTED_INPUTS)
    copy_exact(frozen,source,{HELPER_REL:EXPECTED_HELPER})
    copy_exact(repo,source,additions)
    copy_exact(repo,source,{MANIFEST_REL:EXPECTED_MANIFEST})
    for rel,digest in reference_files.items():
        existing=source/rel
        if existing.exists():require(sha(existing)==digest,'Reference input conflicts with native source: '+rel)
        if not existing.exists():
            existing.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(reference_sources[rel],existing);existing.chmod(0o444)
        require(sha(existing)==digest,'Reference input copy mismatch: '+rel)
    dump(dest/'provenance/reference_inputs.json',dict(files=reference_files,authority=reference_provenance,reviewed_case_leaf_authority=dict(archive_sha256=REVIEW_ARCHIVE_SHA,source_manifest_sha256=REVIEW_MANIFEST_SHA,no_baseline_adoption=True)))
    shutil.copyfile(review_manifest,dest/'provenance/reviewed_case_leaf_manifest.sha256')
    require(sha(dest/'provenance/reviewed_case_leaf_manifest.sha256')==REVIEW_MANIFEST_SHA,'Copied reviewed leaf authority changed')
    copy_exact(overlay_base/'files',dest/'frozen/source_overlay/files',manifest['files'])
    discovery=materialize_overlay_sources(source,dest/'frozen/source_overlay/files',manifest['files'])
    dump(dest/'provenance/overlay_discovery.json',discovery)
    # Confirm every protected native/input/reviewed file retained its exact byte pin.
    authenticate_files(source,{**original_pins,**additions,**reference_files})
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
        overlay_file_count=1241,overlay_discovery_counts=discovery['counts'],reference_input_count=len(reference_files),identity=inv['identity'],plans={m:pin(dest/'inputs'/f'{m}_manifest.json') for m in ('smoke','fit')},
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
    verify_overlay_materialization(stage/'frozen/source',read(stage/'frozen/source'/MANIFEST_REL)['files'],read(stage/'provenance/overlay_discovery.json'))
    require(sha(stage/'provenance/reviewed_case_leaf_manifest.sha256')==REVIEW_MANIFEST_SHA,'Reviewed case-leaf authority differs')
    reference_receipt=read(stage/'provenance/reference_inputs.json')
    reference_manifest=read(stage/'frozen/source'/REFERENCE_REL)
    expected_reference,_=reference_input_pins(reference_manifest,inv['canonical_root'],lambda item:read(stage/'frozen/source'/Path(item['path']).relative_to(inv['canonical_root'])))
    require(reference_receipt['files']==expected_reference,'Reference input graph differs')
    authenticate_files(stage/'frozen/source',expected_reference)
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


COMPLETION_SCHEMA = 'two_shock_constructor_completion_v1'
COMPLETION_INVENTORY = 'c2618b8b263ed90c82ca60f40b1f02c8836bc6ebeede1db3e89087ea7bbc9ce3'
COMPLETION_STAGE = REPO/BASE_REL/'two_shock_v1/execution_smoke_v4'
COMPLETION_REMOTE = '/scratch/td2248/projects/current_estate_two_shock_20261004_v4'
COMPLETION_PROOFS = ('launcher_start.json','launcher_terminal.json','host_verification.json',
                     'container_verification.json','run/preflight.json')


def constructor_body(launcher):
    """Extract the unique original heredoc without changing a single byte."""
    lines=launcher.splitlines(keepends=True)
    starts=[i for i,line in enumerate(lines) if "<<'PYNATIVE'" in line]
    ends=[i for i,line in enumerate(lines) if line.rstrip('\r\n')=='PYNATIVE']
    require(len(starts)==len(ends)==1 and starts[0]<ends[0], 'Unique bounded PYNATIVE markers required')
    body=''.join(lines[starts[0]+1:ends[0]])
    require(len(body.encode())<16384, 'Constructor body exceeds bounded extraction')
    for strict in ('plan=json.load(open(sys.argv[2]));d.preflight(plan)',
                   'runtime=r.build_runtime(plan=plan,output=out,smoke=True)',
                   'assert get_num_threads()==1', 'assert runtime.rt.total_native_calls==0',
                   "assert runtime.rt.identity()==plan['identity']", "out/'native_status.json'"):
        require(strict in body, 'Original strict constructor check missing: '+strict)
    return body


def completion_fingerprints(plan):
    """Pure JSON equivalent of the authenticated driver's smoke fingerprints."""
    require(plan.get('smoke') is True, 'Pinned smoke plan required')
    digest=lambda value: hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    controls={k:plan[k] for k in ('gates','seed','fit','endpoint','path','budget','horizons','smoke_seed_endpoint_padding')}
    return dict(source=digest(plan['source_pins']),contract=digest(dict({k:plan[k] for k in
        ('schema','kind','identity','baseline_psi','psi_bound_ratios','stages','rows','weights','target_contract')},
        legacy_source_overlay=plan.get('legacy_source_overlay'))),controls=digest(controls),
        empirical_controls=digest(plan['empirical_controls']),smoke_controls=digest(dict(
        protocol=plan['smoke_protocol'],numerical_controls=controls,fixed_levels=[plan['baseline_psi']]*2,
        scalar_optimizer=False,derivative_probes=False)))


def validate_completion_proofs(stage,proof_rel,inventory,plan):
    proof=Path(stage)/safe_rel(proof_rel)
    values={name:read(proof/name) for name in COMPLETION_PROOFS}
    start=values['launcher_start.json'];terminal=values['launcher_terminal.json']
    require(start.get('mode')=='preflight' and start.get('wall_seconds')==600 and
        start.get('inventory_sha256')==COMPLETION_INVENTORY and start.get('cpus')==1 and
        start.get('memory_gib')==24 and start.get('numba_threads')==start.get('blas_threads')==1 and
        start.get('no_auto_retry') is True, 'Original preflight start/resources differ')
    require(type(start.get('start_epoch')) is int and start.get('deadline_epoch')==start['start_epoch']+600,
        'Original preflight deadline differs')
    require(start.get('output')==str(Path(inventory['local_root'])/proof_rel), 'Original preflight output differs')
    require(terminal.get('mode')=='preflight' and type(terminal.get('exit_code')) is int and
        terminal['exit_code']==124 and terminal.get('no_auto_retry') is True and
        terminal.get('start_epoch')==start['start_epoch'] and terminal.get('deadline_epoch')==start['deadline_epoch'] and
        terminal.get('slurm_job_id')==start.get('slurm_job_id') and
        type(terminal.get('finished_epoch')) in (int,float) and
        start['start_epoch']<=terminal['finished_epoch']<=start['deadline_epoch'], 'Known failed terminal124 required')
    for name in ('host_verification.json','container_verification.json'):
        v=values[name]
        require(v.get('status')=='PASS_ZERO_SOLVES_PACKAGE' and v.get('files')==4599 and
            type(v.get('native_calls')) is int and v['native_calls']==0 and
            type(v.get('model_solves')) is int and v['model_solves']==0 and
            v.get('scientific_validation') is False, 'Actual passed 4599-file zero-call proof required: '+name)
    v=values['run/preflight.json']
    require(v.get('status')=='PASS' and type(v.get('native_calls')) is int and v['native_calls']==0 and
        v.get('scientific_validation') is False and v.get('production_ready') is False and
        v.get('schema')==plan['schema'] and v.get('fingerprints')==completion_fingerprints(plan) and
        v.get('horizons')==plan['horizons'] and v.get('total_seconds')==plan['budget']['total_seconds'] and
        v.get('policy_call_stop_cap')==plan['budget']['maximum_policy_calls'], 'Driver PASS fingerprints differ')
    return values


def completion_stage_pins(stage,inventory):
    names=['launch_torch.sh','prepare_two_shock_torch.py','inputs/smoke_manifest.json',
        'inputs/fit_manifest.json','inputs/base_fit_plan.json','frozen/source/'+DRIVER_REL,
        'frozen/source/'+RUNTIME_REL,'frozen/source/'+HELPER_REL,'frozen/source_overlay/overlay.json']
    require(all(name in inventory['files'] for name in names), 'Required package control pin missing')
    pins={name:inventory['files'][name] for name in names}
    authenticate_files(stage,pins)
    return pins


def prepare_constructor_completion(args):
    stage=Path(args.stage).resolve();dest=Path(args.destination).resolve()
    require(not dest.exists(), 'Fresh destination required: '+str(dest))
    require(stage!=dest and not dest.is_relative_to(stage), 'Completion packet must be outside immutable stage')
    require(sha(stage/'inventory.json')==COMPLETION_INVENTORY, 'Original v4 inventory SHA differs')
    inv=read(stage/'inventory.json')
    require(inv.get('schema')==SCHEMA and inv.get('local_root')==str(stage) and
        inv.get('remote_root')==COMPLETION_REMOTE and len(inv['files'])==4599, 'Original v4 package identity differs')
    pins=completion_stage_pins(stage,inv);plan=read(stage/'inputs/smoke_manifest.json')
    require(plan['identity']==inv['identity'], 'Original native identity differs')
    for name,rel in [('two_shock_driver',DRIVER_REL),('two_shock_runtime',RUNTIME_REL)]:
        require(plan['source_pins'][name]==dict(path=str(stage/'frozen/source'/rel),sha256=pins['frozen/source/'+rel]),
            'Original driver/runtime plan source pin differs')
    body=constructor_body((stage/'launch_torch.sh').read_text())
    proof_rel=str(safe_rel(args.proof_relative));validate_completion_proofs(stage,proof_rel,inv,plan)
    require(re.fullmatch(r'[A-Za-z0-9_-]+',args.name) is not None, 'Explicit safe completion name required')
    remote_packet=Path(args.remote_packet)
    require(remote_packet.is_absolute() and remote_packet!=Path(inv['remote_root']) and
        not remote_packet.is_relative_to(Path(inv['remote_root'])), 'Fresh separate remote packet required')
    proof_pins={str(Path(proof_rel)/name):sha(stage/proof_rel/name) for name in COMPLETION_PROOFS}
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'constructor.py').write_text(body)
    shutil.copyfile(Path(__file__),dest/'prepare_two_shock_torch.py')
    for name in COMPLETION_PROOFS:
        target=dest/'proofs'/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(stage/proof_rel/name,target)
    manifest=dict(schema=COMPLETION_SCHEMA,stage_local=str(stage),stage_remote=inv['remote_root'],
        packet_local=str(dest),packet_remote=str(remote_packet),canonical_root=inv['canonical_root'],
        inventory_sha256=COMPLETION_INVENTORY,stage_pins=pins,proof_relative=proof_rel,proof_pins=proof_pins,
        name=args.name,identity=plan['identity'],fingerprints=completion_fingerprints(plan),
        constructor_sha256=sha(dest/'constructor.py'),original_overall_preflight='FAILED_EXIT_124',
        resources=dict(cpus=1,memory_gib=24,numba_threads=1,blas_threads=1,wall_seconds=600),
        no_model_solves=True,no_smoke_or_fit=True,no_auto_retry=True)
    dump(dest/'manifest.json',manifest)
    files={str(path.relative_to(dest)):sha(path) for path in dest.rglob('*') if path.is_file()}
    dump(dest/'inventory.json',dict(schema=COMPLETION_SCHEMA,files=files))
    inventory_sha=sha(dest/'inventory.json')
    script='''#!/usr/bin/env bash
#SBATCH --job-name=constructor_completion
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --account=torch_pr_570_general
#SBATCH --partition=cl
set -euo pipefail
module load anaconda3/2025.06
exec /share/apps/anaconda3/2025.06/bin/python {helper} run-constructor-completion --packet {packet} --inventory-sha256 {digest}
'''.format(helper=shlex.quote(str(remote_packet/'prepare_two_shock_torch.py')),packet=shlex.quote(str(remote_packet)),digest=inventory_sha)
    (dest/'launch_constructor.sh').write_text(script)
    # The launcher embeds the inventory digest; its independent pin is returned.
    for path in dest.rglob('*'):
        if path.is_file():path.chmod(0o555 if path.name=='launch_constructor.sh' else 0o444)
    return dict(status='PREPARED_ZERO_SOLVE_COMPLETION_NOT_RUN',destination=str(dest),remote_packet=str(remote_packet),
        inventory_sha256=inventory_sha,launcher_sha256=sha(dest/'launch_constructor.sh'),files=len(files),
        stage_inventory_sha256=COMPLETION_INVENTORY,original_overall_preflight='FAILED_EXIT_124')


def verify_constructor_completion(packet,inventory_sha256):
    """Bounded proof/control checks only; never scan the old 4599-file tree."""
    packet=Path(packet)
    require(sha(packet/'inventory.json')==inventory_sha256, 'Completion inventory SHA differs')
    inv=read(packet/'inventory.json');require(inv.get('schema')==COMPLETION_SCHEMA, 'Completion inventory schema differs')
    require(1<=len(inv['files'])<=16, 'Completion packet is not bounded')
    authenticate_files(packet,inv['files']);m=read(packet/'manifest.json')
    require(m.get('schema')==COMPLETION_SCHEMA and m.get('inventory_sha256')==COMPLETION_INVENTORY and
        m.get('resources')==dict(cpus=1,memory_gib=24,numba_threads=1,blas_threads=1,wall_seconds=600) and
        m.get('original_overall_preflight')=='FAILED_EXIT_124' and m.get('no_smoke_or_fit') is True and
        m.get('no_auto_retry') is True, 'Completion contract/resources differ')
    require(str(packet) in (m['packet_local'],m['packet_remote']), 'Completion packet namespace differs')
    stage=Path(m['stage_remote'] if str(packet)==m['packet_remote'] else m['stage_local'])
    require(sha(stage/'inventory.json')==COMPLETION_INVENTORY, 'Original v4 inventory SHA differs')
    original=read(stage/'inventory.json');pins=completion_stage_pins(stage,original)
    require(pins==m['stage_pins'], 'Original package control pins differ')
    authenticate_files(stage,m['proof_pins'])
    for name in COMPLETION_PROOFS:
        require(sha(packet/'proofs'/name)==m['proof_pins'][str(Path(m['proof_relative'])/name)], 'Copied proof differs')
    plan=read(stage/'inputs/smoke_manifest.json')
    validate_completion_proofs(stage,m['proof_relative'],original,plan)
    require(m['identity']==original['identity']==plan['identity'] and m['fingerprints']==completion_fingerprints(plan),
        'Completion native/target identity differs')
    require((packet/'constructor.py').read_text()==constructor_body((stage/'launch_torch.sh').read_text()) and
        sha(packet/'constructor.py')==m['constructor_sha256'], 'Exact original constructor body differs')
    return m


def bounded_completion_process(command,*,environment,stdin,stdout,seconds):
    require(0<seconds<=585, 'Bounded constructor process budget required')
    child=subprocess.Popen(command,env=environment,stdin=stdin,stdout=stdout,stderr=subprocess.STDOUT,start_new_session=True)
    try:
        return child.wait(timeout=seconds)
    except (subprocess.TimeoutExpired,KeyboardInterrupt):
        status=124
    finally:
        # Clean the entire group even if its leader exited leaving descendants.
        try:os.killpg(child.pid,signal.SIGTERM)
        except ProcessLookupError:pass
        try:child.wait(timeout=5)
        except subprocess.TimeoutExpired:pass
        try:os.killpg(child.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        child.wait(timeout=5)
    return status


def run_constructor_completion(args):
    import resource
    start=time.monotonic();started=time.time();packet=Path(args.packet).resolve()
    m=verify_constructor_completion(packet,args.inventory_sha256)
    require(str(packet)==m['packet_remote'], 'Constructor completion runs only in its remote namespace')
    require(not os.getenv('SLURM_ARRAY_TASK_ID') and os.getenv('SLURM_CPUS_PER_TASK','1')=='1' and
        int(os.getenv('SLURM_MEM_PER_NODE','24576'))<=24576, 'One CPU / 24GiB non-array context required')
    stage=Path(m['stage_remote']);jobs=stage/'jobs';results=stage/'results'
    claim=jobs/'constructor_completion.claim';claim.mkdir(exist_ok=False)
    dump(claim/'claim.json',dict(name=m['name'],inventory_sha256=args.inventory_sha256,started_epoch=started,no_auto_retry=True))
    out=results/m['name'];out.mkdir(exist_ok=False)
    job=jobs/m['name'];job.mkdir(exist_ok=False)
    receipt=dict(mode='constructor_completion',start_epoch=started,deadline_epoch=started+600,
        inventory_sha256=args.inventory_sha256,stage_inventory_sha256=COMPLETION_INVENTORY,
        original_overall_preflight='FAILED_EXIT_124',output=str(Path(m['stage_local'])/'results'/m['name']),
        resources=m['resources'],no_auto_retry=True,scientific_validation=False)
    dump(out/'launcher_start.json',receipt);dump(job/'launcher_start.json',receipt)
    status=1;native=None
    def terminated(signum,frame):raise KeyboardInterrupt
    previous={sig:signal.signal(sig,terminated) for sig in (signal.SIGTERM,signal.SIGINT)}
    try:
        cap=24*1024**3;soft,hard=resource.getrlimit(resource.RLIMIT_AS)
        resource.setrlimit(resource.RLIMIT_AS,(min(cap,hard) if hard!=resource.RLIM_INFINITY else cap,hard))
        environment=dict(os.environ,NUMBA_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
            VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
        for key in ('APPTAINER_BIND','APPTAINER_BINDPATH','SINGULARITY_BIND','SINGULARITY_BINDPATH'):environment.pop(key,None)
        for name in ('numba_cache','matplotlib'):(out/name).mkdir()
        visible=Path(m['stage_local'])/'results'/m['name']
        environment.update(NUMBA_CACHE_DIR=str(visible/'numba_cache'),MPLCONFIGDIR=str(visible/'matplotlib'))
        command=['apptainer','exec','--bind',f"{stage}/frozen/source:{m['canonical_root']}:ro",'--bind',f"{stage}:{m['stage_local']}:ro",
            '--bind',f"{jobs}:{m['stage_local']}/jobs:rw",'--bind',f"{results}:{m['stage_local']}/results:rw",'--pwd',m['stage_local'],
            '/share/apps/images/ubuntu-24.04.4.sif','/share/apps/anaconda3/2025.06/bin/python','-',
            str(Path(m['stage_local'])/'frozen/source'/DRIVER_REL),str(Path(m['stage_local'])/'inputs/smoke_manifest.json'),str(visible/'native_constructor')]
        remaining=585-(time.monotonic()-start)
        with (packet/'constructor.py').open('rb') as body,(out/'native_constructor.log').open('wb') as log:
            status=bounded_completion_process(command,environment=environment,stdin=body,stdout=log,seconds=remaining)
        if status==0:
            native=read(out/'native_constructor/native_status.json')
            require(native.get('status')=='PASS_ZERO_SOLVES_NATIVE_CONSTRUCTOR' and type(native.get('native_calls')) is int and
                native['native_calls']==0 and native.get('numba_threads')==1 and native.get('identity')==m['identity'] and
                native.get('scientific_validation') is False, 'Actual zero-call constructor native_status required')
    except KeyboardInterrupt:status=143
    except Exception:
        status=1;raise
    finally:
        for sig,handler in previous.items():signal.signal(sig,handler)
        terminal=dict(receipt,exit_code=status,finished_epoch=time.time(),native_status_pin=pin(out/'native_constructor/native_status.json') if native else None,
            status='PASS_ZERO_SOLVES_CONSTRUCTOR_COMPLETION' if status==0 and native else 'FAILED_CONSTRUCTOR_COMPLETION',
            composite_ready=status==0 and native is not None)
        dump(out/'launcher_terminal.json',terminal);dump(job/'launcher_terminal.json',terminal)
    require(status==0 and native is not None, 'Constructor completion failed: exit '+str(status))
    return terminal


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build');b.add_argument('--repo',default=str(REPO));b.add_argument('--destination',default=str(REPO/BASE_REL/'two_shock_v1/execution_smoke_v1'))
    b.add_argument('--remote-root',default=REMOTE_DEFAULT);b.add_argument('--python',default=sys.executable)
    v=sub.add_parser('verify');v.add_argument('--stage',required=True)
    c=sub.add_parser('prepare-constructor-completion')
    c.add_argument('--stage',default=str(COMPLETION_STAGE));c.add_argument('--destination',required=True)
    c.add_argument('--remote-packet',required=True);c.add_argument('--proof-relative',default='results/preflight_preflight_v1')
    c.add_argument('--name',default='constructor_completion_v1')
    c=sub.add_parser('verify-constructor-completion');c.add_argument('--packet',required=True);c.add_argument('--inventory-sha256',required=True)
    c=sub.add_parser('run-constructor-completion');c.add_argument('--packet',required=True);c.add_argument('--inventory-sha256',required=True)
    args=p.parse_args(argv)
    if args.command=='build':result=build(args)
    elif args.command=='verify':result=verify(args.stage)
    elif args.command=='prepare-constructor-completion':result=prepare_constructor_completion(args)
    elif args.command=='verify-constructor-completion':
        m=verify_constructor_completion(args.packet,args.inventory_sha256)
        result=dict(status='PASS_ZERO_SOLVES_COMPLETION_PACKET',stage_inventory_sha256=m['inventory_sha256'],native_calls=0)
    else:result=run_constructor_completion(args)
    print(json.dumps(result,sort_keys=True,indent=2));return result
if __name__=='__main__':main()
