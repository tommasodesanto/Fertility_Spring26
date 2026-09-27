"""Bounded, authenticated central differences of the normalized E5F objective.

No model changes. Local dimensionless coordinates: logarithms except beta,
housing loading and child curvature divided by fixed economic scales.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, math, os, signal, subprocess, sys, time
from pathlib import Path

SCALES={'beta_annual':.05,'delta_alpha_jump':.25,'child_benefit_curvature':.8}
def read(p): return json.loads(Path(p).read_text())
def write(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_suffix(p.suffix+'.new');q.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n');q.replace(p)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def make_points(point,restrictions):
    rows=[];bounds={r['parameter']:(r['lower'],r['upper']) for r in restrictions}
    for name,value in point.items():
        h=.01 if name in SCALES else .02
        for sign in (-1,1):
            p=dict(point);p[name]=value+sign*h*SCALES[name] if name in SCALES else value*math.exp(sign*h)
            assert bounds[name][0]<=p[name]<=bounds[name][1],name
            rows.append(dict(id=name+('_minus' if sign<0 else '_plus'),parameter=name,sign=sign,h=h,coordinate='linear/'+str(SCALES[name]) if name in SCALES else 'log',point=p))
    return rows

def collect(plan,out):
    import numpy as np
    out=Path(out);anchor=list(csv.DictReader(open(Path(plan['anchor_case'])/'target_fit.csv')))
    names=[r['moment'] for r in anchor];columns=[];details=[];psi=[]
    for name in plan['point']:
        pair=[r for r in plan['cases'] if r['parameter']==name];mm,pp=pair
        vals=[];benefits=[]
        for r in pair:
            folder=out/r['id']/'case'
            if not (folder/'receipt.json').exists(): return None
            table=list(csv.DictReader(open(folder/'target_fit.csv')))
            assert [x['moment'] for x in table]==names
            vals.append(np.array([float(x['model']) for x in table]))
            benefits.append(read(folder/'receipt.json')['normalization']['psi_child'])
        derivative=(vals[1]-vals[0])/(2*pp['h']);columns.append(derivative)
        psi.append((benefits[1]-benefits[0])/(2*pp['h']))
        curvature=(vals[1]+vals[0]-2*np.array([float(x['model']) for x in anchor]))/(pp['h']**2)
        for i,m in enumerate(names): details.append(dict(parameter=name,moment=m,coordinate=pp['coordinate'],step=pp['h'],derivative=derivative[i],second_difference=curvature[i],minus_model=vals[0][i],plus_model=vals[1][i],weight=float(anchor[i]['weight'] or 0)))
    raw=np.column_stack(columns);weights=np.sqrt([float(r['weight'] or 0) for r in anchor]);weighted=raw*weights[:,None]
    u,s,vh=np.linalg.svd(weighted,full_matrices=False)
    gaps=np.array([float(r['gap']) for r in anchor])*weights
    result=dict(status='complete',parameters=list(plan['point']),moments=names,raw_jacobian=raw.tolist(),weighted_jacobian=weighted.tolist(),singular_values=s.tolist(),right_singular_vectors=vh.tolist(),local_condition_number=float(s[0]/s[-1]) if s[-1]>0 else None,loss_gradient=(2*weighted.T@gaps).tolist(),normalized_benefit_derivatives=psi,caveat='Local finite differences on a discrete grid; no global identification or optimizer convergence claim. Singular values depend on documented coordinate scaling. Step-halving remains unperformed.')
    write(out/'jacobian.json',result)
    with (out/'jacobian.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(details[0]));w.writeheader();w.writerows(details)
    return result

def run(plan_path,out):
    plan=read(plan_path);out=Path(out);out.mkdir(parents=True,exist_ok=False)
    contract=Path(plan['contract']);c=read(contract)
    assert sha(contract)==plan['contract_sha256'];assert sha(__file__)==plan['runner_sha256']
    env=os.environ.copy();env['EXPECTED_UTILITY_OVERNIGHT_SHA256']=sha(contract);env['E5F_LOCAL_EXECUTION_AUTHORIZATION']=c['execution']['authorization_id']
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):env[k]='1'
    os.environ.update(env)
    spec=importlib.util.spec_from_file_location('frozen_jacobian_driver',c['files']['driver']['path']);driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    driver.verify(contract);driver.verify_execution(c)
    start=time.time();deadline=start+plan['total_seconds'];active={};records=[];pending=list(plan['cases']);failed=False
    write(out/'launch.json',dict(plan=plan,plan_sha256=sha(plan_path),pid=os.getpid(),start=start,deadline=deadline))
    while pending or active:
        # First pair is an exact-controller-loop smoke before remaining dispatch.
        smoke_ok=len(records)>=2 and all(r['status']=='success' for r in records[:2])
        while pending and not failed and len(active)<2 and (smoke_ok or len(records)+len(active)<2) and time.time()<deadline:
            row=pending.pop(0);name=row['id'];p=row['point'];end=min(deadline,time.time()+plan['case_seconds'])
            context=dict(candidate_id=name,stage='jacobian',contract_sha256=sha(contract),source_sha256=c['files']['driver']['sha256'],target_sha256=c['objective']['sha256'],point_sha256=driver.canon(p))
            request=dict(id=name,point=p,context=context,stage='jacobian',point_sha256=driver.canon(p),scientific_candidate_id=driver.candidate_id(c,p),normalization_inputs=driver.norm_inputs(c),contract_sha256=sha(contract),controller_pid=os.getpid(),deadline_epoch=end,graphs=False)
            path=out/(name+'.request.json');write(path,request);log=(out/(name+'.log')).open('w')
            proc=subprocess.Popen([sys.executable,c['files']['driver']['path'],'--stage','evaluate','--contract',str(contract),'--output',str(out/name),'--request',str(path)],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            active[name]=(proc,log,request,end)
        for name,(proc,log,req,end) in list(active.items()):
            timeout=time.time()>=end and proc.poll() is None
            if timeout:
                os.killpg(proc.pid,signal.SIGTERM)
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait(timeout=5)
            if proc.poll() is None:continue
            log.close();status='failure';data={}
            try:
                if proc.returncode!=0 or timeout:raise RuntimeError('Evaluator timeout' if timeout else f'Evaluator exit {proc.returncode}')
                data=driver.validate_success(out/name,c,req);status='success'
            except Exception as exc:data={'error':str(exc)};failed=True
            rec=dict(id=name,status=status,completed=time.time(),**data);records.append(rec);del active[name]
            write(out/'latest_completed.json',rec)
            successes=[r for r in records if r['status']=='success']
            if successes:write(out/'best_so_far.json',min(successes,key=lambda r:r['loss']))
            write(out/'checkpoint.json',dict(records=records,pending=[r['id'] for r in pending],halted=failed))
        write(out/'heartbeat.json',dict(time=time.time(),active={n:p[0].pid for n,p in active.items()},completed=len(records),pending=len(pending),failed=failed,deadline=deadline))
        if (failed or time.time()>=deadline) and not active:break
        if active:time.sleep(5)
    if not failed and not pending:collect(plan,out)
    write(out/'complete.json',dict(status='complete' if not failed and not pending else 'incomplete',records=records,not_run=[r['id'] for r in pending],end=time.time()))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--collect-only',action='store_true');a=ap.parse_args()
    if a.collect_only:collect(read(a.plan),a.output)
    else:run(a.plan,a.output)
if __name__=='__main__':main()
