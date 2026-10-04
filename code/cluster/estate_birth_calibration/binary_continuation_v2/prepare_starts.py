"""Build a 20-start one-birth continuation plan from authenticated saved points."""
import argparse, hashlib, json, math, random
from pathlib import Path

RECOVERY_SHA='974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22'
PARENT_SHA='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
RECOVERY_JOB='19141024'
PARENT_JOB='19127370'
PARENT_TASK_JOB_IDS=('19127379','19127380','19127381','19127382','19127383')
PLAN_REL=Path('source/output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/plan_binary.json')
PARENT_PLAN_REL=Path('source/output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json')

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def canonical(x): return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def normalized_distance(a,b,bounds):
    return math.sqrt(sum(((a[k]-b[k])/(bounds[k][1]-bounds[k][0]))**2 for k in bounds))

def recovery_points(root,inv,base_plan):
    points=[]; records=[]
    for chain in range(5):
        folder=root/'results'/f'production_binary_chain_{chain}'
        launch=read(folder/'launcher_start.json');terminal=read(folder/'launcher_terminal.json')
        contract=read(folder/'run/start_contract.json');done=read(folder/'run/completed.json')
        search=read(folder/'run/search_completed.json');native=read(folder/'run/native_postcheck/completed.json')
        require(launch['mode']=='production' and launch['arm']=='binary' and launch['chain']==chain and launch['stage_inventory_sha256']==RECOVERY_SHA,'Recovery launch identity drift')
        require(terminal['exit_code']==0 and terminal['arm']=='binary' and terminal['chain']==chain and terminal['slurm_job_id']==launch['slurm_job_id'],'Recovery task failed')
        require(done['status']=='selected_numerically_verified' and done['arm']=='binary' and done['birth_cap']==1 and done['chain']==chain,'Recovery selected gate failed')
        require(native['status']=='full_native_postcheck_passed' and native['search_receipt_sha256']==sha(folder/'run/search_completed.json'),'Recovery selected native receipt drift')
        require(done['repeat']['status']=='exact_full_ge_repeat_passed' and len(done['repeat']['standard_plot_hashes'])==17,'Recovery repeat gate failed')
        require(len(done['target_fit'])==14 and len(done['parameters'])==31 and done['selected_postcheck']['status']=='passed','Recovery report contract drift')
        require(all(done[k]==native[k]==inv[k] for k in ('target_fingerprint','weight_fingerprint')),'Recovery target/weight fingerprints drift')
        require(contract['arm']=='binary' and contract['birth_cap']==1 and contract['chain']==chain and contract['target_fingerprint']==inv['target_fingerprint'] and contract['weight_fingerprint']==inv['weight_fingerprint'],'Recovery start contract drift')
        require(contract['starts_file_sha256']==sha(root/PLAN_REL) and search['starts_file_sha256']==contract['starts_file_sha256'],'Recovery source plan drift')
        point=done['selected']['parameters'];reported={r['parameter']:float(r['estimate']) for r in done['parameters'] if r['parameter'] in point}
        require(point==reported,'Recovery selected parameter table drift')
        points.append(point);records.append(dict(kind='verified_recovery_selected',source_job_id=launch['slurm_job_id'],chain=chain,
            source_inventory_sha256=RECOVERY_SHA,completed_sha256=sha(folder/'run/completed.json'),
            native_postcheck_sha256=sha(folder/'run/native_postcheck/completed.json'),search_sha256=sha(folder/'run/search_completed.json'),loss=done['native_loss']))
    return points,records

def parent_points(root,inv,plan):
    require(sha(root/'inventory.json')==PARENT_SHA,'Parent inventory bytes drift')
    old=read(root/PARENT_PLAN_REL)
    require(sha(root/PARENT_PLAN_REL)==inv['start_plan_sha256'] and old['target_fingerprint']==inv['target_fingerprint'] and old['weight_fingerprint']==inv['weight_fingerprint'],'Parent plan identity drift')
    require(old['bounds']==plan['bounds'] and old['target_contract']==plan['target_contract'],'Parent/continuation contract drift')
    points=[];records=[]
    for chain in range(5):
        folder=root/'results'/f'production_binary_chain_{chain}'
        launch=read(folder/'launcher_start.json');terminal=read(folder/'launcher_terminal.json')
        contract=read(folder/'run/start_contract.json');best=read(folder/'run/best_so_far.json')['best'];cases=read(folder/'run/cases.json')
        require(launch['mode']=='production' and launch['arm']=='binary' and launch['chain']==chain and launch['stage_inventory_sha256']==PARENT_SHA,'Parent launch identity drift')
        require(str(terminal['slurm_job_id'])==str(launch['slurm_job_id'])==PARENT_TASK_JOB_IDS[chain] and terminal['exit_code']!=0,'Parent array task accounting identity/status drift')
        require(contract['arm']=='binary' and contract['birth_cap']==1 and contract['chain']==chain and contract['starts_file_sha256']==inv['start_plan_sha256'],'Parent start contract drift')
        require(contract['target_fingerprint']==inv['target_fingerprint']==plan['target_fingerprint'] and contract['weight_fingerprint']==inv['weight_fingerprint']==plan['weight_fingerprint'] and contract['bounds']==plan['bounds'],'Parent target/bounds drift')
        require(best['status']=='passed' and best in cases and best==min((c for c in cases if c.get('status')=='passed'),key=lambda c:c['loss']),'Parent saved feasible checkpoint drift')
        point=best['parameters'];require(set(point)==set(plan['bounds']),'Parent checkpoint coordinates drift')
        points.append(point);records.append(dict(kind='parent_saved_feasible_checkpoint',parent_array_job_id=PARENT_JOB,
            parent_array_task_job_id=launch['slurm_job_id'],chain=chain,
            source_inventory_sha256=PARENT_SHA,best_sha256=sha(folder/'run/best_so_far.json'),cases_sha256=sha(folder/'run/cases.json'),
            start_contract_sha256=sha(folder/'run/start_contract.json'),provisional_loss=best['loss'],verified_native=False))
    return points,records

def perturbations(center,bounds,occupied):
    rng=random.Random(20261004);out=[];metadata=[]
    for j in range(10):
        for attempt in range(1000):
            radius=(0.035,0.05,0.065,0.08,0.095)[j%5]
            z=[rng.gauss(0,1) for _ in bounds]
            norm=math.sqrt(sum(v*v for v in z)) or 1.
            point={}
            for (k,(lo,hi)),v in zip(bounds.items(),z):
                x=center[k]+radius*(hi-lo)*v/norm
                point[k]=float(min(hi,max(lo,x)))
            if canonical(point) in {canonical(p) for p in occupied+out}: continue
            if occupied+out and min(normalized_distance(point,p,bounds) for p in occupied+out)<0.005: continue
            out.append(point);metadata.append(dict(kind='bounded_perturbation',center='verified_recovery_chain_1',radius_fraction=radius,draw=j,attempt=attempt,seed=20261004));break
        else: raise RuntimeError('Could not produce a distinct bounded perturbation')
    return out,metadata

def prepare(recovery,parent,out):
    rinv=read(recovery/'inventory.json');pinv=read(parent/'inventory.json')
    require(sha(recovery/'inventory.json')==RECOVERY_SHA,'Recovery stage inventory drift')
    require(sha(parent/'inventory.json')==PARENT_SHA,'Parent stage inventory drift')
    base=read(recovery/PLAN_REL);old=read(parent/PARENT_PLAN_REL)
    require(base['recovery_arm']=='binary' and len(base['starts'])==5,'Recovery binary plan drift')
    require(base['target_fingerprint']==rinv['target_fingerprint']==old['target_fingerprint'] and base['weight_fingerprint']==rinv['weight_fingerprint']==old['weight_fingerprint'],'Target/weight mismatch across authenticated stages')
    require(base['bounds']==old['bounds'] and base['target_contract']==old['target_contract'],'Parameter or target contract mismatch')
    rpoints,rmeta=recovery_points(recovery,rinv,base)
    ppoints,pmeta=parent_points(parent,pinv,base)
    center=rpoints[1]
    perturb,pmeta2=perturbations(center,base['bounds'],rpoints+ppoints)
    starts=rpoints+ppoints+perturb;provenance=rmeta+pmeta+pmeta2
    require(len(starts)==20 and len({canonical(x) for x in starts})==20,'Start count/uniqueness failure')
    distances=[normalized_distance(starts[i],starts[j],base['bounds']) for i in range(20) for j in range(i)]
    require(min(distances)>0.005,'Start points are not sufficiently distinct')
    plan=dict(base)
    plan.update(starts=starts,arms={'binary':1},start_provenance=[dict(index=i,**x) for i,x in enumerate(provenance)],
        continuation_arm='binary',recovery_job_id=RECOVERY_JOB,recovery_inventory_sha256=RECOVERY_SHA,
        parent_job_id=PARENT_JOB,parent_array_task_job_ids=list(PARENT_TASK_JOB_IDS),parent_inventory_sha256=PARENT_SHA,recovery_root='/work/recovery',parent_root='/work/parent',
        source_provenance=dict(recovery_inventory_sha256=RECOVERY_SHA,recovery_stage_inventory_sha256=sha(recovery/'inventory.json'),
          parent_inventory_sha256=PARENT_SHA,parent_plan_sha256=sha(parent/PARENT_PLAN_REL),recovery_plan_sha256=sha(recovery/PLAN_REL)),
        pairwise_normalized_distances=distances,deterministic_seed=20261004,no_auto_retry=True,no_scientific_adoption=True,
        continuation_design='5 verified recovery selected vectors; 5 authenticated provisional parent feasible checkpoints; 10 bounded perturbations around verified recovery chain 1')
    out.mkdir(parents=True,exist_ok=True)
    require(not (out/'starts.json').exists() and not (out/'starts_receipt.json').exists(),'Refusing preexisting start plan')
    b=(json.dumps(plan,sort_keys=True,indent=2,allow_nan=False)+'\n').encode();(out/'starts.json').write_bytes(b)
    receipt=dict(status='20_authenticated_distinct_starts_prepared',starts_sha256=hashlib.sha256(b).hexdigest(),count=20,
        categories=dict(verified_recovery_selected=5,parent_saved_feasible_checkpoint=5,bounded_perturbation=10),
        minimum_normalized_pairwise_distance=min(distances),recovery_inventory_sha256=RECOVERY_SHA,parent_inventory_sha256=PARENT_SHA,
        target_fingerprint=rinv['target_fingerprint'],weight_fingerprint=rinv['weight_fingerprint'],no_model_solves=True)
    (out/'starts_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

def audit(plan,recovery,parent):
    rinv=read(recovery/'inventory.json');pinv=read(parent/'inventory.json')
    require(sha(recovery/'inventory.json')==RECOVERY_SHA and sha(parent/'inventory.json')==PARENT_SHA,'Source inventory drift')
    base=read(recovery/PLAN_REL)
    rp,rm=recovery_points(recovery,rinv,base);pp,pm=parent_points(parent,pinv,base)
    perturb,pm2=perturbations(rp[1],base['bounds'],rp+pp)
    points=rp+pp+perturb;meta=rm+pm+pm2
    require(plan['starts']==points and plan['start_provenance']==[dict(index=i,**x) for i,x in enumerate(meta)],'Plan differs from authenticated checkpoint selection')
    require(plan['target_fingerprint']==rinv['target_fingerprint'] and plan['weight_fingerprint']==rinv['weight_fingerprint'],'Plan objective fingerprint drift')
    return dict(status='all_20_starts_reauthenticated',recovery_selected=5,parent_feasible=5,perturbations=10)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--recovery',type=Path,required=True);p.add_argument('--parent',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();print(json.dumps(prepare(a.recovery,a.parent,a.out),indent=2))
