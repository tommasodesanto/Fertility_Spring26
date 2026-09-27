#!/usr/bin/env python3
"""Bounded local continuation using the unchanged frozen objective and supervisor.

The separate plan selects candidate points and wall-clock limits only. Economic
parameters/targets/gates remain in the authenticated scientific contract.
"""
from __future__ import annotations
import argparse, collections, hashlib, importlib.util, json, math, os, signal, sys, time
from pathlib import Path

THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_name(p.name+'.new')
    q.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n');q.replace(p)
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def validate_plan(plan,objective):
    assert plan['workers']==2
    assert isinstance(plan['authorization'],str) and plan['authorization'].strip()
    assert plan['supersedes_overnight_deadline'] is True
    assert 0<plan['search_seconds']<plan['total_seconds']<=1740
    assert plan['total_seconds']-plan['search_seconds']>=600
    assert 0<plan['case_seconds']<=480
    assert 1<=len(plan['cases'])<=12
    bounds={r['parameter']:(r['lower'],r['upper']) for r in objective['parameter_restrictions']}
    names=[]
    for row in [dict(id='anchor',point=plan['point'])]+plan['cases']:
        name=row['id'];assert name and all(c.isalnum() or c in '_-' for c in name)
        assert not name.startswith(('smoke_','repeat_'));names.append(name)
        assert set(row['point'])==set(bounds)
        for key,value in row['point'].items(): assert math.isfinite(value) and bounds[key][0]<=value<=bounds[key][1],key
    assert len(names)==len(set(names))
    if 'absolute_end_epoch' in plan: assert math.isfinite(plan['absolute_end_epoch'])

def request(driver,c,contract_sha,row,stage,end,graphs):
    context=dict(candidate_id=row['id'],stage=stage,contract_sha256=contract_sha,source_sha256=c['files']['driver']['sha256'],target_sha256=c['objective']['sha256'],point_sha256=driver.canon(row['point']))
    return dict(id=row['id'],point=row['point'],context=context,stage=stage,point_sha256=driver.canon(row['point']),scientific_candidate_id=driver.candidate_id(c,row['point']),normalization_inputs=driver.norm_inputs(c),contract_sha256=contract_sha,controller_pid=os.getpid(),deadline_epoch=end,graphs=graphs)

def setup(plan_path):
    plan=read(plan_path);contract=Path(plan['contract']).resolve()
    assert sha(contract)==plan['contract_sha256'];assert sha(__file__)==plan['runner_sha256']
    c=read(contract);assert c['execution']['kind']=='local'
    os.environ['EXPECTED_UTILITY_OVERNIGHT_SHA256']=sha(contract)
    os.environ['E5F_LOCAL_EXECUTION_AUTHORIZATION']=c['execution']['authorization_id']
    os.environ.update({k:'1' for k in THREADS});os.environ['MPLBACKEND']='Agg'
    driver=load('frozen_local_continuation_driver',c['files']['driver']['path'])
    _,objective=driver.verify(contract);driver.verify_execution(c);validate_plan(plan,objective)
    anchor=Path(plan['anchor_case']).resolve()
    assert set(plan['anchor_files'])=={'receipt.json','initial_state.pkl.gz','target_fit.csv','parameters.csv'}
    for name,expected in plan['anchor_files'].items(): assert sha(anchor/name)==expected,name
    r=read(anchor/'receipt.json')
    assert r['point']==plan['point']
    assert r['target_system_sha256']==c['objective']['sha256']
    assert r['source_manifest_sha256']==c['source_manifest']['sha256']
    assert r['case_checkpoint_sha256']==sha(anchor/'initial_state.pkl.gz')
    assert len(driver.keyed_csv(anchor/'target_fit.csv','moment'))==14
    assert len(driver.keyed_csv(anchor/'parameters.csv','parameter'))==31
    sys.path.insert(0,c['runtime_tools'])
    supervision=load('frozen_local_continuation_supervision',c['files']['recovery_search']['path'])
    return plan,contract,c,anchor,driver,supervision

def run(plan_path,out,validate_only=False):
    plan,contract,c,anchor,driver,supervision=setup(plan_path)
    if validate_only:
        print(json.dumps(dict(status='plan_verified_no_solves',workers=2,candidates=len(plan['cases']))));return
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    start=time.time();end=min(start+plan['total_seconds'],plan.get('absolute_end_epoch',math.inf))
    cutoff=min(start+plan['search_seconds'],end-(plan['total_seconds']-plan['search_seconds']))
    assert start<cutoff
    clock=dict(start=start,search_cutoff=cutoff,end=end)
    records=[];batches=[];best=None;fatal=False;stage='setup';last_heartbeat=0
    write(out/'launch.json',dict(pid=os.getpid(),plan=plan,plan_sha256=sha(plan_path),clock=clock))
    print(json.dumps(dict(controller_pid=os.getpid(),output=str(out),clock=clock)),flush=True)
    def heartbeat(**progress):
        nonlocal last_heartbeat
        if time.time()-last_heartbeat>=5 or progress.pop('force',False):
            write(out/'heartbeat.json',dict(epoch=time.time(),stage=stage,completed=len(records),best_loss=best['loss'] if best else None,**progress));last_heartbeat=time.time()
    def batch(rows,phase,deadline,graphs=False):
        nonlocal stage,best,fatal
        stage=phase
        def launch(row,batch_deadline):
            driver.verify(contract)
            req=request(driver,c,sha(contract),row,phase,min(batch_deadline,time.time()+plan['case_seconds']),graphs)
            path=out/(row['id']+'.request.json');write(path,req)
            row.update(context=req['context'],request=req,request_path=str(path))
            cmd=[sys.executable,c['files']['driver']['path'],'--stage','evaluate','--contract',str(contract),'--output',str(out/row['id']),'--request',str(path)]
            return supervision.ManagedProcess(cmd,out/(row['id']+'.log'),req['deadline_epoch'],os.environ.copy())
        def finish(row,process,code):
            nonlocal best,fatal
            try:
                driver.verify(contract);status,data,error=driver.classify(out/row['id'],c,row,process,code)
            except Exception as exc: status,data,error='fatal',{},dict(error=str(exc),classification='controller_integrity_failure')
            rec=dict(case=row['id'],status=status,point=row['point'],request_path=row['request_path'],completed_epoch=time.time(),execution=dict(returncode=code,deadline=process.deadline,owned_timeout=process.observed_running_at_expiry and process.deadline_kill_reaped),error=error,**data)
            records.append(rec)
            if status=='fatal' or phase in ('smoke','repeat') and status!='success':fatal=True
            rec['halt_new_dispatch']=fatal
            if status=='success' and phase!='repeat' and (best is None or rec['loss']<best['loss']):best=rec;write(out/'best_so_far.json',best)
            write(out/'latest_completed.json',rec);write(out/'checkpoint.json',dict(records=records,best=best,clock=clock,fatal=fatal));heartbeat(force=True)
            return rec
        result=supervision.run_batch(rows,workers=2,deadline=deadline,launch=launch,finish=finish,heartbeat=heartbeat,poll_seconds=1,allowed_statuses={'success'} if phase in ('smoke','repeat') else {'success','inadmissible','censored_timeout','censored_late_completion'},guard=lambda:'fatal_stop' if fatal else None)
        if any(r['status']=='failed' for r in result['results']):fatal=True
        batches.append(dict(stage=phase,**result));write(out/'batches.json',batches);return result
    status='incomplete';error=None
    try:
        smoke=batch([dict(id=f'smoke_{i}',point=plan['point']) for i in (1,2)],'smoke',cutoff)
        assert not fatal and smoke['complete'] and len(smoke['results'])==2,'Anchor smoke failed'
        comparisons=[driver.compare_tables(anchor,r['case_path']) for r in smoke['results']]
        write(out/'smoke_complete.json',dict(status='exact_anchor_smokes_passed',comparisons=comparisons,records=smoke['results']))
        batch([dict(r) for r in plan['cases']],'search',cutoff)
        assert not fatal,'Fatal search result; stopped without recovery or relabeling'
        selected=dict(best);driver.verify_recovery_integrity(c,contract,selected)
        write(out/'selected.json',dict(selected=selected,frozen_before_repeats=True))
        repeat=batch([dict(id=f'repeat_{i}',point=selected['point']) for i in (1,2)],'repeat',end-30,graphs=True)
        assert not fatal and repeat['complete'] and len(repeat['results'])==2,'Final repeats incomplete'
        driver.verify_recovery_integrity(c,contract,selected)
        driver.export_selected(selected,repeat['results'],out/'selected_export',c,end)
        status='bounded_search_complete'
    except Exception as exc:
        error=dict(type=type(exc).__name__,message=str(exc));status='fatal_stop' if fatal else 'incomplete'
    finally:
        attempted={r['case'] for r in records}
        write(out/'complete.json',dict(status=status,error=error,records=records,best=best,counts=dict(collections.Counter(r['status'] for r in records)),not_run=[r['id'] for r in plan['cases'] if r['id'] not in attempted],clock=clock,end_epoch=time.time()))
        heartbeat(force=True)
    if status!='bounded_search_complete':raise RuntimeError(f'{status}: {error}')

def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--validate-only',action='store_true');a=p.parse_args()
    if a.validate_only:
        run(a.plan,a.output,True);return
    plan=read(a.plan)
    remaining=min(plan['total_seconds'],plan.get('absolute_end_epoch',math.inf)-time.time())
    if remaining<=0:raise TimeoutError('Local continuation plan expired')
    def stop(signum,frame):raise TimeoutError('Controller interrupted or absolute budget reached: signal '+str(signum))
    for signum in (signal.SIGALRM,signal.SIGTERM,signal.SIGINT):signal.signal(signum,stop)
    signal.setitimer(signal.ITIMER_REAL,remaining)
    try:run(a.plan,a.output)
    finally:signal.setitimer(signal.ITIMER_REAL,0)
if __name__=='__main__':main()
