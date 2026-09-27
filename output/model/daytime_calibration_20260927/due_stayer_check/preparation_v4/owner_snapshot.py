import json,time,hashlib,pathlib,subprocess,os,signal
root=pathlib.Path('/Users/tommasodesanto/.codex/worktrees/due-existing-owner-credit/Fertility_Spring26')
base=root/'output/model/daytime_calibration_20260927/due_stayer_check'
prep=base/'preparation_v4';plan=json.loads((prep/'plan_draft.json').read_text());plan['execution_authorized']=True
plan['status']='lead_reviewed_bounded_same_price_match';plan['absolute_end_epoch']=1790536354.173106
p=root/'code/model/tools/e5f_due_purchase_audit.py'
plan['audit_contract'].update(purchase_audit_path=str(p),purchase_audit_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
plan['lead_review']={'scope':'DUE arm after explicit saved baseline review;39 focused tests and17 baseline plots reviewed','focused_tests':40,'economic_change':'none in baseline, DUE only in later same-parameter arm','expected_seconds_per_solve':60,'max_solves':4}
review=base/'baseline_v2_review/complete.json';assert json.loads(review.read_text())['parameters_exact']==31
plan['baseline_review']={'path':str(review),'sha256':hashlib.sha256(review.read_bytes()).hexdigest()}
planpath=prep/'plan_reviewed.json';planpath.write_text(json.dumps(plan,indent=2)+'\n');planpath.chmod(0o444)
for p in (prep/'source_snapshot').rglob('*.py'):p.chmod(0o444)
out=base/'due_v2';log=base/'due_v2.log';deadline=min(time.time()+300,plan['absolute_end_epoch'])
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']});env['E5F_LOCAL_EXECUTION_AUTHORIZATION']='tommaso_authorized_20260927_local_primary_continuation_v1'
cmd=['/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python',str(root/'code/model/tools/run_e5f_due_stayer_matched_check.py'),'--plan',str(planpath),'--arm','due','--output',str(out)]
with log.open('w') as f:
 child=subprocess.Popen(cmd,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
 launch={'pid':child.pid,'controller_pid':os.getpid(),'start_epoch':time.time(),'deadline_epoch':deadline,'plan_sha256':hashlib.sha256(planpath.read_bytes()).hexdigest(),'command':cmd,'workers':1,'threads':1}
 (base/'due_v2_launch.json').write_text(json.dumps(launch,indent=2)+'\n');print(json.dumps(launch),flush=True)
 try:code=child.wait(timeout=max(0,deadline-time.time()));status='completed' if code==0 else 'failed'
 except subprocess.TimeoutExpired:
  os.killpg(child.pid,signal.SIGKILL);code=child.wait();status='owned_timeout'
 (base/'due_v2_process.json').write_text(json.dumps({'status':status,'exit_code':code,'end_epoch':time.time()},indent=2)+'\n');print(status,code,flush=True)
