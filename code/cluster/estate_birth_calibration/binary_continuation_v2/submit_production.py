"""Exactly-once guarded production submitter; an ambiguous sbatch outcome stays locked."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path('/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2')
def sha(p):
 import hashlib
 return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    control=ROOT/'control';lock=control/'production_submission.lock';receipt=control/'production_submission.json'
    if receipt.exists() or lock.exists(): raise SystemExit('Submission already attempted; inspect receipt/Slurm before any manual recovery')
    subprocess.run(['/share/apps/anaconda3/2025.06/bin/python',str(ROOT/'verify_stage.py'),'--host'],check=True)
    subprocess.run(['/share/apps/anaconda3/2025.06/bin/python',str(ROOT/'verify_starts.py'),str(control/'starts.json')],check=True)
    subprocess.run(['/share/apps/anaconda3/2025.06/bin/python',str(ROOT/'record_smoke_gate.py'),'--verify-existing',str(control/'smoke_gate.json')],check=True)
    stat=os.statvfs('/scratch/td2248')
    free_bytes=stat.f_bavail*stat.f_frsize
    reserve_bytes=10200*1024**3
    if free_bytes<reserve_bytes: raise SystemExit(f'Insufficient scratch reserve: {free_bytes/1024**3:.1f} GiB available; require 10,200 GiB')
    lock.mkdir()
    state=dict(status='submitting',array='0-19%20',stage_inventory_sha256=sha(ROOT/'inventory.json'),starts_sha256=sha(control/'starts.json'),
      smoke_gate_sha256=sha(control/'smoke_gate.json'),task_count=20,cores_per_task=1,memory_GiB_per_task=24,
      task_wall_seconds=43200,objective_max_calls_per_task=500,case_planning_cap_GiB=1,
      scratch_free_bytes_at_submit=free_bytes,scratch_free_GiB_at_submit=round(free_bytes/1024**3,3),
      scratch_reserve_bytes=reserve_bytes,scratch_reserve_GiB=10200,native_reserve_seconds=1800,no_auto_retry=True)
    tmp=control/'production_submission.tmp';tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(receipt)
    result=subprocess.run(['sbatch','--hold','--parsable','--array=0-19%20','--export=ALL,ESTATE_RUN_MODE=production',str(ROOT/'launch_torch.sh')],capture_output=True,text=True,check=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit(): raise RuntimeError('Unparseable sbatch result; submission receipt remains SUBMITTING for manual queue audit')
    state.update(status='submitted',array_job_id=job,submitted_epoch=time.time(),sbatch_stdout=result.stdout.strip(),release_status='held_until_receipt_pinned')
    tmp=control/'production_submission.tmp';tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(receipt)
    # Workers cannot start until the submitted receipt with their array ID is durable.
    released=subprocess.run(['scontrol','release',job],capture_output=True,text=True,check=True)
    state.update(release_status='released',released_epoch=time.time(),scontrol_stdout=released.stdout.strip())
    tmp=control/'production_submission.tmp';tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(receipt)
    print(json.dumps(state))
if __name__=='__main__':main()
