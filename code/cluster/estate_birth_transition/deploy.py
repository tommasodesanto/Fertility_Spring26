#!/usr/bin/env python3
"""Immutable Estate-A transition packaging and explicit, separately gated submission."""
from __future__ import annotations
import argparse, copy, gzip, hashlib, io, json, math, os, shlex, shutil, subprocess, tarfile, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Copied deployers run at /work/deployment (or a remote stage), where a
# repository-parent index is undefined. The inventoried canonical root owns
# the container's mounted namespace; source-checkout use retains its own root.
if (HERE/'inventory.json').is_file():
    ROOT = Path(json.loads((HERE/'inventory.json').read_text())['root'])
    if not ROOT.is_absolute():raise RuntimeError('Staged inventory root must be absolute')
elif HERE.parts[-3:] == ('code','cluster','estate_birth_transition'):
    ROOT = HERE.parents[2]
else:
    raise RuntimeError('Copied deployer requires its adjacent inventory.json')
REMOTE = '/scratch/td2248/projects/current_estate_transition_20261003_v6'
OUT = ROOT/'output/model/transition_readiness_v1/current_baseline_20261003/deployment_v6'
PARENT = ROOT/'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/stage.tar.gz'
PARENT_SHA = '5cb99a84f1aa51dab49462292ab31756f9183776e30d2d327844ae2af2c5fc44'
CASE = 'output/model/experiments/birth_count_choice/estate_a_v1/single/cases/20261003T212605039706Z_a739edc3'
EXCLUDED = {'.pyc', '.nbc', '.nbi'}

def sha(blob): return hashlib.sha256(blob).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def write(path, data): Path(path).write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
def require(condition, message):
    if not condition: raise RuntimeError(message)
def run(command): return subprocess.run(command, check=True, text=True, capture_output=True).stdout.strip()
def ssh(command): return run(['ssh', '-o', 'BatchMode=yes', 'torch', command])
def quote(value): return shlex.quote(str(value))
def relative(path): return str(Path(path).resolve().relative_to(ROOT))

def build(args):
    out = Path(args.stage).resolve()
    require(getattr(args,'panel_config',None), 'v6 build requires an explicit pinned panel config')
    require(not out.exists() or not any(out.iterdir()), 'Refusing nonempty deployment directory')
    previous=Path(args.from_stage).resolve() if getattr(args,'from_stage',None) else None
    if previous:
        verify(argparse.Namespace(stage=previous,mounted_root=None))
        prior=read(previous/'inventory.json')
        require(prior['remote_root']!=REMOTE, 'Source delta requires a distinct previous remote root')
        source={rel:(previous/'source'/rel).read_bytes() for rel in prior['files']}
        parent_archive_sha256=prior['parent_archive_sha256']
    else:
        require(sha(PARENT.read_bytes()) == PARENT_SHA, 'Tested immutable parent archive drift')
        with tarfile.open(PARENT) as archive:
            prior = json.load(archive.extractfile('inventory.json'))
            source = {n.removeprefix('source/'):archive.extractfile(n).read() for n in archive.getnames()
                      if n.startswith('source/') and archive.getmember(n).isfile()}
        require({k:sha(v) for k,v in source.items()} == prior['files'], 'Immutable parent inventory drift')
        parent_archive_sha256=PARENT_SHA
    def add_file(path):
        path = Path(path).resolve(); source[relative(path)] = path.read_bytes()
    def add_directory(path):
        for file in sorted((ROOT/path).rglob('*')):
            if file.is_file() and '__pycache__' not in file.parts and file.suffix not in EXCLUDED: add_file(file)
    for directory in ('code/model/experiments/birth_count_choice', 'code/model/experiments/transition_readiness',
                      'code/model/production/reference_inputs', CASE): add_directory(directory)
    plan_paths = {'smoke':relative(args.smoke_plan), 'fit':relative(args.fit_plan)}
    # Traverse declared pins only: never copy historical controller outputs or Jacobians.
    def collect(value):
        if isinstance(value, dict):
            if isinstance(value.get('path'), str) and isinstance(value.get('sha256'), str):
                path = Path(value['path']); require(path.is_file(), 'Missing plan pin: '+str(path))
                require(sha(path.read_bytes()) == value['sha256'], 'Declared plan pin drift: '+str(path)); add_file(path)
            for key,item in value.items():
                if key == 'source_pins':
                    for rel,digest in item.items():
                        path=ROOT/rel; require(sha(path.read_bytes())==digest, 'Source pin drift: '+rel); add_file(path)
                else: collect(item)
        elif isinstance(value, list):
            for item in value: collect(item)
    plans={}
    for mode,rel in plan_paths.items():
        add_file(ROOT/rel); plans[mode]=read(ROOT/rel); collect(plans[mode])
        handoff=read(plans[mode]['handoff']['path']); collect(handoff)
        require(handoff['saved_case']==CASE, 'Wrong saved Estate-A baseline')
        for name,digest in {**handoff.get('saved_files',{}),**handoff.get('standard_plot_pins',{})}.items():
            require(sha(source[CASE+'/'+name])==digest, 'Saved case pin drift: '+name)
        require(plans[mode]['current_context_schema']=='current_estate_a_one_permanent_v1', 'Wrong current contract')
        require(plans[mode]['execution_enabled'] is True, 'Execution must be explicit in pinned plan')
    require(plans['smoke']['identity']==plans['fit']['identity'], 'Smoke/fit identities differ')
    require(plans['smoke']['horizons']==[6,8] and plans['smoke']['path']['max_evaluations']==6, 'Wrong exact-loop smoke')
    panel_path=relative(args.panel_config) if getattr(args,'panel_config',None) else None
    if panel_path:
        add_file(ROOT/panel_path);collect(read(ROOT/panel_path))
    names=('deploy.py','launch_torch.sh','panel_launch_torch.sh') if panel_path else ('deploy.py','launch_torch.sh')
    entrypoints={name:sha((HERE/name).read_bytes()) for name in names}
    inventory=dict(schema='estate_a_transition_stage_v1', remote_root=REMOTE, root=str(ROOT),
        parent_archive_sha256=parent_archive_sha256, files={k:sha(v) for k,v in sorted(source.items())},
        entrypoints=entrypoints, plans={k:dict(path=v,sha256=sha(source[v])) for k,v in plan_paths.items()},
        identity=plans['fit']['identity'], resources=dict(cpus=8,numba_threads=8,blas_threads=1,memory_GiB=96,maximum_wall_seconds=21600),
        no_old_transition_state_or_jacobian_reuse=True)
    if panel_path:
        inventory['panel_config']=dict(path=panel_path,sha256=sha(source[panel_path]))
        validate_panel(inventory,read(ROOT/panel_path),plans['fit'])
        inventory['panel_resources']=dict(nodes_per_task=1,ntasks=1,cpus_per_task=8,memory_GiB=96,guesses=12,maximum_concurrent_nodes=12,array='0-11%12')
    out.mkdir(parents=True, exist_ok=True)
    if previous:shutil.copytree(previous/'source',out/'source',copy_function=os.link)
    write(out/'inventory.json', inventory)
    changed=[k for k,v in source.items() if not previous or prior['files'].get(k)!=sha(v)]
    entries={'source/'+k:source[k] for k in changed}
    entries.update({name:(HERE/name).read_bytes() for name in entrypoints})
    entries['inventory.json']=(out/'inventory.json').read_bytes()
    # Keep the local payload identical to the archive: later verification and
    # submission read plans/source directly from this immutable stage.
    for name,blob in sorted(entries.items()):
        target=out/name
        require(not Path(name).is_absolute() and '..' not in Path(name).parts, 'Unsafe payload path: '+name)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():target.unlink()  # Detach delta source hardlinks.
        target.write_bytes(blob)
        target.chmod(0o755 if name.endswith('.sh') else 0o644)
    with (out/'stage.tar.gz').open('wb') as raw, gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
        with tarfile.open(fileobj=gz,mode='w') as archive:
            for name,blob in sorted(entries.items()):
                info=tarfile.TarInfo(name);info.size=len(blob);info.mtime=0;info.mode=0o755 if name.endswith('.sh') else 0o644
                archive.addfile(info,io.BytesIO(blob))
    verification=verify(argparse.Namespace(stage=out,mounted_root=None))
    receipt=dict(status='prepared_no_submission',local_verification=verification,remote_root=REMOTE,source_files=len(source),
        archive_sha256=sha((out/'stage.tar.gz').read_bytes()),inventory_sha256=sha((out/'inventory.json').read_bytes()))
    if previous:receipt.update(remote_parent=prior['remote_root'],remote_parent_inventory_sha256=sha((previous/'inventory.json').read_bytes()),
        changed_source_files=changed,source_delta=True)
    write(out/'stage_receipt.json', receipt); return receipt

def validate_panel(inventory, config, plan):
    require(config.get('schema')=='estate_a_transition_panel_v1', 'Wrong panel config schema')
    require(config.get('identity')==inventory['identity'], 'Panel identity differs from current runtime')
    pin=config.get('plan',{})
    require(set(pin)=={'path','sha256'} and relative(pin['path'])==inventory['plans']['fit']['path'] and
            pin['sha256']==inventory['plans']['fit']['sha256'], 'Panel must pin the exact fit plan')
    require(plan['mode']=='diagnostic' and plan['horizons']==[24,32] and plan['budget']['total_seconds']<=21480, 'Panel requires retained bounded 24/32 diagnostic fit plan')
    guesses=config.get('guesses',[])
    require(len(guesses)==12 and [row.get('index') for row in guesses]==list(range(12)), 'Exactly 12 ordered independent guesses required')
    for row in guesses:
        psi=row.get('psi')
        require(type(psi) in (float,int) and math.isfinite(psi) and
                plan['initial_psi']*plan['psi_bound_ratios'][0]<psi<plan['initial_psi']*plan['psi_bound_ratios'][1], 'Panel guess outside retained preference domain')
    require(len({float(row['psi']) for row in guesses})==12,'Panel guesses must be distinct')
    panel=config.get('panel_source',{})
    require(set(panel)=={'path','sha256'} and relative(panel['path'])=='code/model/experiments/birth_count_choice/transition_panel.py' and inventory['files'].get(relative(panel['path']))==panel['sha256'], 'Separate panel-driver source pin required')
    return dict(status='passed_zero_solves',guesses=12,maximum_concurrent_nodes=12,policy_parameter='psi_child')

def panel_verify(args):
    stage=Path(args.stage);inv=read(stage/'inventory.json');root=Path(args.mounted_root) if args.mounted_root else stage/'source'
    panel=inv['panel_config'];require(sha((root/panel['path']).read_bytes())==panel['sha256'],'Panel config drift')
    return validate_panel(inv,read(root/panel['path']),read(root/inv['plans']['fit']['path']))

def verify(args):
    stage=Path(args.stage); inventory=read(stage/'inventory.json')
    root=Path(args.mounted_root) if args.mounted_root else stage/'source'
    for rel,digest in inventory['files'].items():
        require((root/rel).is_file() and sha((root/rel).read_bytes())==digest, 'Source/data drift: '+rel)
    for name,digest in inventory['entrypoints'].items():
        require(sha((stage/name).read_bytes())==digest, 'Entrypoint drift: '+name)
    return dict(status='passed_zero_solves',source_files=len(inventory['files']),inventory_sha256=sha((stage/'inventory.json').read_bytes()))

def mount_bindings(stage, base, floor, repo):
    """Consolidate complete pinned packages; preserve exact final file precedence."""
    stage,base,floor=map(Path,(stage,base,floor)); repo=str(repo)
    own=read(stage/'inventory.json')['files']; final={}
    for layer in (base,floor,stage):
        for rel in read(layer/'inventory.json')['files']:final[rel]=layer/'source'/rel
    # Only complete, explicitly shipped packages and current input subtrees may
    # replace directories. Historical output inventories are sparse overlays:
    # mounting their parent hides authenticated siblings in the frozen root.
    approved_groups=(
        'code/model/experiments/birth_count_choice',
        'code/model/experiments/purchase_timing_sandbox',
        'code/model/experiments/transition_readiness',
        'code/model/refactor_lab',
        'code/model/production/reference_inputs',
        CASE,
        'output/model/transition_readiness_v1/current_baseline_20261003/plans',
        'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source/small_credit_lab',
    )
    groups={group for group in approved_groups if any(rel.startswith(group+'/') for rel in own)}
    # Avoid nested directory mounts: each staged file has a unique directory owner.
    groups=sorted(g for g in groups if not any(g.startswith(other+'/') for other in groups if other!=g))
    pairs=[(str(stage/'source'/g),repo+'/'+g) for g in groups]
    for rel,source in sorted(final.items()):
        covered=any(rel.startswith(g+'/') for g in groups)
        if source==stage/'source'/rel and covered:continue
        pairs.append((str(source),repo+'/'+rel))
    require(len({target for _,target in pairs})==len(pairs),'Duplicate bind target')
    # Prove the effective origin of every original pin after consolidation.
    for rel,source in final.items():
        target=repo+'/'+rel;resolved=None
        for origin,destination in pairs:
            if target==destination:resolved=Path(origin)
            elif target.startswith(destination+'/'):resolved=Path(origin)/target[len(destination)+1:]
        require(resolved==source,'Consolidated bind changed source: '+rel)
    bindings=[origin+':'+target+':ro' for origin,target in pairs]
    byte_count=sum(len(b.encode())+1 for b in bindings)
    require(byte_count<=100000,'Source bind arguments exceed 100 KB; fixed launcher binds reserve 10 KB below 128 KiB')
    return dict(bindings=bindings,binding_count=len(bindings),source_bind_bytes=byte_count,
                verified_staged_files=len(own),effective_pinned_files=len(final),duplicate_targets=0)

def revise_launcher(args):
    require(False, 'v6 panel deployment requires build --from-stage with both unchanged plans and panel config')
    previous=Path(args.from_stage).resolve();out=Path(args.stage).resolve()
    verify(argparse.Namespace(stage=previous,mounted_root=None))
    require(not out.exists() or not any(out.iterdir()),'Refusing nonempty revised stage')
    inventory=read(previous/'inventory.json');parent_remote=inventory['remote_root']
    replacement=None
    if getattr(args,'fit_plan',None):
        rel=inventory['plans']['fit']['path'];old_plan=read(previous/'source'/rel)
        blob=Path(args.fit_plan).read_bytes();new_plan=json.loads(blob)
        require(new_plan['identity']==inventory['identity'], 'Fit override changes current runtime identity')
        left,right=copy.deepcopy(old_plan),copy.deepcopy(new_plan)
        for item in (left,right):
            item['fit'].pop('max_evaluations');item['budget'].pop('maximum_policy_calls')
        require(left==right, 'Fit override may change only fit evaluation and native-call caps')
        require(type(new_plan['fit']['max_evaluations']) is int and new_plan['fit']['max_evaluations']>=5, 'Invalid fit evaluation cap')
        require(type(new_plan['budget']['maximum_policy_calls']) is int and new_plan['budget']['maximum_policy_calls']>0, 'Invalid native-call cap')
        replacement=(rel,blob)
    require(parent_remote!=REMOTE,'Revision requires a distinct previous remote root')
    out.mkdir(parents=True,exist_ok=True)
    # Source hardlinks are read-only inputs; this command never modifies them.
    shutil.copytree(previous/'source',out/'source',copy_function=os.link)
    if replacement:
        rel,blob=replacement;target=out/'source'/rel
        # Detach this file before writing: the unchanged payload uses hardlinks.
        target.unlink();target.write_bytes(blob)
        inventory['files'][rel]=sha(blob);inventory['plans']['fit']['sha256']=sha(blob)
    inventory['remote_root']=REMOTE
    inventory['resources']=dict(cpus=8,numba_threads=8,blas_threads=1,memory_GiB=96,maximum_wall_seconds=21600)
    inventory['entrypoints']={name:sha((HERE/name).read_bytes()) for name in ('deploy.py','launch_torch.sh')}
    write(out/'inventory.json',inventory)
    entries={'inventory.json':(out/'inventory.json').read_bytes()}
    if replacement:entries['source/'+replacement[0]]=replacement[1]
    for name in inventory['entrypoints']:
        entries[name]=(HERE/name).read_bytes();(out/name).write_bytes(entries[name])
        (out/name).chmod(0o755 if name.endswith('.sh') else 0o644)
    with (out/'stage.tar.gz').open('wb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
        with tarfile.open(fileobj=gz,mode='w') as archive:
            for name,blob in sorted(entries.items()):
                info=tarfile.TarInfo(name);info.size=len(blob);info.mtime=0;info.mode=0o755 if name.endswith('.sh') else 0o644
                archive.addfile(info,io.BytesIO(blob))
    verification=verify(argparse.Namespace(stage=out,mounted_root=None))
    receipt=dict(status='prepared_launcher_delta_no_submission',remote_root=REMOTE,remote_parent=parent_remote,
       remote_parent_inventory_sha256=sha((previous/'inventory.json').read_bytes()),
       source_files=len(inventory['files']),source_inventory_unchanged=replacement is None,
       changed_source_files=[] if replacement is None else [replacement[0]],
       runtime_sources_and_smoke_plan_unchanged=True,local_verification=verification,
       archive_sha256=sha((out/'stage.tar.gz').read_bytes()),inventory_sha256=sha((out/'inventory.json').read_bytes()))
    write(out/'stage_receipt.json',receipt);return receipt

def stage_remote(args):
    stage=Path(args.stage); verify(argparse.Namespace(stage=stage,mounted_root=None))
    require(read(stage/'inventory.json')['remote_root']==REMOTE, 'Stage belongs to a different remote root')
    archive=stage/'stage.tar.gz';receipt=read(stage/'stage_receipt.json')
    require(sha(archive.read_bytes())==receipt['archive_sha256'], 'Archive drift')
    parent=receipt.get('remote_parent')
    if parent:
        require(parent!=REMOTE and parent.startswith('/scratch/td2248/projects/'), 'Invalid delta parent')
        check='import hashlib,sys;assert hashlib.sha256(open(sys.argv[1],"rb").read()).hexdigest()==sys.argv[2]'
        ssh('/share/apps/anaconda3/2025.06/bin/python -c '+quote(check)+' '+quote(parent+'/inventory.json')+' '+quote(receipt['remote_parent_inventory_sha256']))
        ssh('/share/apps/anaconda3/2025.06/bin/python '+quote(parent+'/deploy.py')+' verify --stage '+quote(parent))
    ssh('test ! -e '+quote(REMOTE)+' && mkdir -p '+quote(REMOTE+'/logs'))
    if parent:ssh('cp -a --reflink=auto '+quote(parent+'/source')+' '+quote(REMOTE+'/source'))
    run(['scp',str(archive),'torch:'+REMOTE+'/stage.tar.gz'])
    result=ssh('cd '+quote(REMOTE)+' && tar -xzf stage.tar.gz && /share/apps/anaconda3/2025.06/bin/python deploy.py verify --stage '+quote(REMOTE))
    write(stage/'remote_stage_receipt.json',{**receipt,'status':'staged_no_submission','verification':result})
    return dict(status='staged_no_submission',verification=result)

def gate_check(stage, gate_path):
    inv=read(stage/'inventory.json'); gate=read(gate_path)
    require(gate.get('reviewed_by_lead') is True and bool(gate.get('reviewer')), 'Explicit lead smoke review required')
    require(gate.get('inventory_sha256')==sha((stage/'inventory.json').read_bytes()), 'Review binds a different stage')
    require(gate.get('plans')==inv['plans'] and gate.get('identity')==inv['identity'], 'Review source/plan identity differs')
    if str(stage) == inv['remote_root']:
        native=Path(gate.get('native_remote_root',str(stage)))
        if native!=stage:
            require(gate.get('bridge_kind')=='identical_native_v5_panel_v6', 'Explicit unchanged-native bridge required')
            require(sha((native/'inventory.json').read_bytes())==gate['native_inventory_sha256'],'Native smoke inventory changed')
        smoke_path=native/'results/smoke/run/smoke_receipt.json'
        require(smoke_path.is_file() and sha(smoke_path.read_bytes())==gate['smoke_receipt_sha256'], 'Actual staged smoke receipt differs')
    return gate

def review(args):
    stage=Path(args.stage); inv=read(stage/'inventory.json'); smoke=read(args.smoke_receipt)
    bridge={}
    if inv.get('panel_config'):
        require(args.native_stage and args.panel_test_receipt, 'Panel review requires authenticated native stage and mocked panel test receipt')
        native=Path(args.native_stage);verify(argparse.Namespace(stage=native,mounted_root=None));old=read(native/'inventory.json')
        require(old['identity']==inv['identity'] and old['plans']==inv['plans'], 'Native identity or plans changed; a new native smoke is required')
        require(all(inv['files'].get(rel)==digest for rel,digest in old['files'].items()), 'Existing native/source/data pins changed; cross-stage smoke forbidden')
        test=read(args.panel_test_receipt);config=read(stage/'source'/inv['panel_config']['path'])
        require(test.get('status')=='PASS' and test.get('panel_source_sha256')==config['panel_source']['sha256'] and test.get('controller_run_reused') is True and test.get('effective_plan_only_added_fit_start_psi') is True, 'Panel mock test must pin current driver and prove unchanged controller/effective plan')
        bridge=dict(bridge_kind='identical_native_v5_panel_v6',native_remote_root=old['remote_root'],
            native_inventory_sha256=sha((native/'inventory.json').read_bytes()),panel_test_receipt_sha256=sha(Path(args.panel_test_receipt).read_bytes()))
    require(smoke.get('status')=='PASS' and smoke.get('identity')==inv['identity'], 'Actual smoke did not pass matching identity')
    require(smoke.get('native_setup_verified') is True and bool(smoke.get('native_root_gates')) and all(smoke['native_root_gates'].values()), 'Actual native root/replay/accounting gates required')
    require(smoke.get('fresh_twelve_date_seed') is True and smoke.get('maximum_maps')==6, 'Fresh measured seed/exact smoke required')
    require(not Path(args.gate).exists(), 'Refusing existing review gate')
    gate=dict(reviewed_by_lead=True,reviewer=args.reviewer,reviewed_epoch=time.time(),identity=inv['identity'],plans=inv['plans'],
        inventory_sha256=sha((stage/'inventory.json').read_bytes()),smoke_receipt_sha256=sha(Path(args.smoke_receipt).read_bytes()))
    gate.update(bridge)
    write(args.gate,gate); return gate

def submit(args):
    stage=Path(args.stage); inv=read(stage/'inventory.json');require((stage/'remote_stage_receipt.json').is_file(), 'Stage first')
    require(inv['remote_root']==REMOTE, 'Submission stage belongs to a different remote root')
    require(0 < args.wall_seconds <= 21600, 'Explicit wall budget must be within six hours')
    plan_mode='fit' if args.mode=='panel' else args.mode
    plan=read(stage/'source'/inv['plans'][plan_mode]['path'])
    require(plan['budget']['total_seconds']+120 <= args.wall_seconds, 'Allow 120 seconds for verification/launcher overhead')
    if args.mode in ('fit','panel'): gate_check(stage,args.gate)
    if inv.get('panel_config'):
        require(args.mode in ('smoke','panel'), 'Panel stage supports only smoke or independent panel')
        panel_verify(argparse.Namespace(stage=stage,mounted_root=None))
    receipt_path=stage/(args.mode+'_submission_receipt.json')
    fd=os.open(receipt_path,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o644)
    with os.fdopen(fd,'w') as handle: json.dump(dict(status='submission_guard_created',mode=args.mode,time=time.time()),handle)
    try:
        if args.mode in ('fit','panel'): run(['scp',args.gate,'torch:'+REMOTE+'/fit_review_gate.json'])
        command='cd '+quote(REMOTE)+' && sbatch --parsable --partition=cl --time='+f'{args.wall_seconds//3600:02d}:{args.wall_seconds%3600//60:02d}:{args.wall_seconds%60:02d}'+' --export=ALL,TRANSITION_MODE='+args.mode+',TRANSITION_WALL_SECONDS='+str(args.wall_seconds)+(' --nodes=1 --ntasks=1 --array=0-11%12' if args.mode=='panel' else ' --nodes=1 --ntasks=1')+(' panel_launch_torch.sh' if inv.get('panel_config') else ' launch_torch.sh')
        job=ssh(command).split(';')[0]; require(job.isdigit(), 'Unrecognized sbatch reply')
        receipt=dict(status='submitted',mode=args.mode,job_id=job,wall_seconds=args.wall_seconds,submitted_epoch=time.time(),
                     inventory_sha256=sha((stage/'inventory.json').read_bytes()),remote_root=REMOTE,partition='cl',cpus=8,numba_threads=8,memory_GiB=96,nodes_per_task=1,ntasks=1,array='0-11%12' if args.mode=='panel' else None,maximum_concurrent_nodes=12 if args.mode=='panel' else 1,no_auto_retry=True)
    except Exception as exc:
        write(receipt_path,dict(status='submission_failed_or_unknown_no_retry',error=str(exc),mode=args.mode));raise
    write(receipt_path,receipt); return receipt

def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='action',required=True)
    for action in ('build','revise-launcher','mounts','verify','verify-panel','verify-gate','stage','review-smoke','submit'):
        item=sub.add_parser(action);item.add_argument('--stage',default=str(OUT))
        if action=='revise-launcher':
            item.add_argument('--from-stage',required=True);item.add_argument('--fit-plan')
        if action=='mounts':
            item.add_argument('--base',required=True);item.add_argument('--floor',required=True);item.add_argument('--repo',required=True)
        if action=='build':
            item.add_argument('--from-stage');item.add_argument('--panel-config')
            item.add_argument('--smoke-plan',required=True);item.add_argument('--fit-plan',required=True)
        if action in ('verify','verify-panel'):item.add_argument('--mounted-root')
        if action=='verify-gate':item.add_argument('--gate',required=True)
        if action=='review-smoke':
            item.add_argument('--smoke-receipt',required=True);item.add_argument('--gate',required=True);item.add_argument('--reviewer',required=True);item.add_argument('--native-stage');item.add_argument('--panel-test-receipt')
        if action=='submit':
            item.add_argument('--mode',choices=('smoke','fit','panel'),required=True);item.add_argument('--wall-seconds',type=int,required=True);item.add_argument('--gate')
    args=parser.parse_args(); actions={'build':build,'revise-launcher':revise_launcher,'mounts':lambda a: mount_bindings(a.stage,a.base,a.floor,a.repo),'verify':verify,'verify-panel':panel_verify,'verify-gate':lambda a: gate_check(Path(a.stage),a.gate),'stage':stage_remote,'review-smoke':review,'submit':submit}
    print(json.dumps(actions[args.action](args),indent=2))
if __name__=='__main__':main()
