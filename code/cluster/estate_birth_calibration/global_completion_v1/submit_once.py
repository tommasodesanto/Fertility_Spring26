"""Exactly-once dependent submission after the zero-solve execution smoke."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path('/scratch/td2248/projects/estate_birth_global_completion_20261005_v1')
PY='/share/apps/anaconda3/2025.06/bin/python'
DEPENDENCY='afterany:19194495:19194496'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    subprocess.run([PY,str(ROOT/'verify_stage.py'),'--host'],check=True,stdout=subprocess.DEVNULL)
    plan=json.loads((ROOT/'control/plan.json').read_text())
    assert plan['task_count']==9 and len(plan['points'])==64 and len(sum(plan['chunks'],[]))==36
    assert plan['task_wall_seconds']==5400 and plan['case_budget_seconds']==1200
    assert plan['target_fingerprint']=='c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70'
    assert plan['weight_fingerprint']=='f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4'
    run=ROOT/'results/smoke_task_0'
    completed=json.loads((run/'run/completed.json').read_text())
    terminal=json.loads((run/'launcher_terminal.json').read_text())
    assert completed['status']=='smoke_loop_passed_zero_solves' and completed['completed_cases']==2
    assert completed['plan_sha256']==sha(ROOT/'control/plan.json') and terminal['exit_code']==0
    cases=json.loads((run/'run/cases.json').read_text())
    assert [c['index'] for c in cases]==[1,2] and all(c['case_result']['lifecycle_solves']==0 for c in cases)
    assert time.time()<1791223200-5400,'Not enough time before global 18UTC cutoff'
    jobs=subprocess.check_output(['squeue','-u',os.environ['USER'],'-h','-o','%j'],text=True).splitlines()
    assert not any(j=='estateglobalcompletion' for j in jobs),'Possible duplicate global completion array'
    free=os.statvfs('/scratch/td2248').f_bavail*os.statvfs('/scratch/td2248').f_frsize
    assert free>=10000*1024**3,'Scratch below 10,000 GiB reserve'
    control=ROOT/'control';receipt=control/'production_submission.json';lock=control/'production_submission.lock'
    if receipt.exists() or lock.exists():raise SystemExit('Submission already attempted; audit receipt/Slurm')
    lock.mkdir()
    state=dict(status='submitting',array='0-8%9',dependency=DEPENDENCY,
               wall_seconds_per_task=5400,case_budget_seconds=1200,
               absolute_cutoff_utc='2026-10-05T18:00:00Z',
               plan_sha256=sha(control/'plan.json'),stage_manifest_sha256=sha(ROOT/'stage_manifest.json'),
               no_auto_retry=True,no_auto_extension=True,submit_epoch=time.time())
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    command=['sbatch','--parsable','--time=01:30:00','--array=0-8%9',
             f'--dependency={DEPENDENCY}','--begin=now+5minutes',
             '--job-name=estateglobalcompletion','--export=ALL,ESTATE_RUN_MODE=production',str(ROOT/'launch_torch.sh')]
    result=subprocess.run(command,capture_output=True,text=True,check=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit():raise RuntimeError('Ambiguous Slurm response; submission receipt remains submitting')
    state.update(status='submitted',job_id=job,sbatch_stdout=result.stdout.strip(),command=command)
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps(state))

if __name__=='__main__':main()
