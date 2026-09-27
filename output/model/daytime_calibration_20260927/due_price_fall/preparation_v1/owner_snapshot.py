import json,time,hashlib,pathlib,subprocess,os,signal
root=pathlib.Path('/Users/tommasodesanto/.codex/worktrees/due-existing-owner-credit/Fertility_Spring26')
base=root/'output/model/daytime_calibration_20260927/due_price_fall';prep=base/'preparation_v1'
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
start=time.time();end=start+900
receipt=prep/'parent_supervisor.json'
write(receipt,dict(status='supervising',pid=os.getpid(),maximum_active_children=1,seconds_per_arm=300,start_epoch=start,absolute_end_epoch=end,owner_source_sha256=sha(pathlib.Path(__file__))))
receipt.chmod(0o444)
plan=json.loads((prep/'plan_draft.json').read_text());plan.update(execution_authorized=True,absolute_end_epoch=end,parent_supervisor_receipt=str(receipt),parent_supervisor_sha256=sha(receipt),lead_review='Reviewed dated operator, fixed-pension PE scope, exact g0, original queue, origin-aware estates;31 focused tests pass; 2 Bellman calls only')
pp=prep/'plan_reviewed.json';write(pp,plan);pp.chmod(0o444)
for p in (prep/'source_snapshot').rglob('*.py'):p.chmod(0o444)
env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']});env['E5F_LOCAL_EXECUTION_AUTHORIZATION']='tommaso_authorized_20260927_local_primary_continuation_v1'
records=[]
for arm in plan['arms']:
 if time.time()>=end:break
 out=base/(arm+'_v1');deadline=min(time.time()+300,end)
 cmd=['/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/.venv/bin/python',str(root/'code/model/tools/run_e5f_due_price_fall_check.py'),'--plan',str(pp),'--arm',arm,'--output',str(out)]
 with (base/(arm+'_v1.log')).open('w') as log:
  child=subprocess.Popen(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  record=dict(arm=arm,pid=child.pid,controller_pid=os.getpid(),start_epoch=time.time(),deadline_epoch=deadline,command=cmd,workers=1,threads=1,plan_sha256=sha(pp));write(base/(arm+'_launch.json'),record);print(json.dumps(record),flush=True)
  try:
   code=child.wait(timeout=max(0,deadline-time.time()));status='completed' if code==0 and time.time()<=deadline else 'failed_or_late'
  except subprocess.TimeoutExpired:
   os.killpg(child.pid,signal.SIGKILL);code=child.wait();status='owned_timeout'
  record.update(exit_code=code,status=status,end_epoch=time.time());records.append(record);write(base/'checkpoint.json',dict(records=records,absolute_end_epoch=end))
  print(arm,status,code,flush=True)
  # Scientific infeasibility permits the independent matched arm; unknown workflow errors stop.
  if code!=0:
   failure=out/'failure.json'
   if not failure.exists() or json.loads(failure.read_text()).get('error_type')!='InheritedDistributionInfeasible':break
write(base/'process_complete.json',dict(records=records,end_epoch=time.time(),absolute_end_epoch=end,remaining_arms=[a for a in plan['arms'] if a not in [r['arm'] for r in records]]))
