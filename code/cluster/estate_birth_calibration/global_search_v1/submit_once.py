"""Exactly-once smoke, native-control, and production submissions."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
ROOT=Path('/scratch/td2248/projects/estate_birth_global_search_20261004_v1')
PY='/share/apps/anaconda3/2025.06/bin/python'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def gate(mode):
    name={'smoke':'smoke','preflight':'preflight'}[mode]
    r=ROOT/'results'/f'{name}_task_0'
    d=json.loads((r/'run/completed.json').read_text())
    t=json.loads((r/'launcher_terminal.json').read_text())
    expect={'smoke':'smoke_loop_passed_zero_solves','preflight':'native_preflight_passed'}[mode]
    assert d['status']==expect and d['completed_cases']==(2 if mode=='smoke' else 1)
    assert d['plan_sha256']==sha(ROOT/'control/plan.json') and t['exit_code']==0
    assert d['target_fingerprint']=='c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70'
    assert d['weight_fingerprint']=='f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4'
    return dict(status='passed',mode=mode,result=str(r),best_loss=d['best']['loss'] if d['best'] else None)
def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['smoke','preflight','production'],required=True);a=p.parse_args();mode=a.mode
    subprocess.run([PY,str(ROOT/'verify_stage.py'),'--host'],check=True)
    plan=json.loads((ROOT/'control/plan.json').read_text())
    assert plan['task_count']==16 and len(plan['points'])==64 and plan['task_wall_seconds']==5400 and plan['case_budget_seconds']==1200
    if mode=='preflight':gate('smoke')
    if mode=='production':gate('smoke');gate('preflight')
    control=ROOT/'control';receipt=control/f'{mode}_submission.json';lock=control/f'{mode}_submission.lock'
    if receipt.exists() or lock.exists():raise SystemExit('Submission already attempted; audit Slurm and receipt rather than resubmit')
    free=os.statvfs('/scratch/td2248').f_bavail*os.statvfs('/scratch/td2248').f_frsize
    reserve=10000*1024**3
    if free<reserve:raise SystemExit('Scratch below 10,000 GiB reserve')
    lock.mkdir()
    array='0-15%16' if mode=='production' else '0-0%1'
    wall={'smoke':'00:05:00','preflight':'00:25:00','production':'01:30:00'}[mode]
    state=dict(status='submitting',mode=mode,array=array,wall=wall,plan_sha256=sha(control/'plan.json'),
               stage_manifest_sha256=sha(ROOT/'stage_manifest.json'),free_bytes=free,reserve_bytes=reserve,
               no_auto_retry=True,no_auto_extension=True,submit_epoch=time.time())
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    args=['sbatch','--parsable',f'--time={wall}',f'--array={array}',f'--export=ALL,ESTATE_RUN_MODE={mode}']
    if mode=='production':args.append('--begin=now+5minutes')
    args.append(str(ROOT/'launch_torch.sh'))
    result=subprocess.run(args,capture_output=True,text=True,check=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit():raise RuntimeError('Ambiguous Slurm submission; receipt remains submitting')
    state.update(status='submitted',job_id=job,sbatch_stdout=result.stdout.strip(),
                 release_status='delayed_start_receipt_pinned' if mode=='production' else 'automatic')
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps(state))
if __name__=='__main__':main()
