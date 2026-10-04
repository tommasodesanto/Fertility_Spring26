"""Submit exactly one non-retried two-call smoke task and pin its receipt."""
import json,subprocess,sys,time
from pathlib import Path
ROOT=Path('/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2')
def main():
    chain=int(sys.argv[1])
    if chain!=0:raise SystemExit('the smoke gate is pinned to chain 0')
    control=ROOT/'control';lock=control/'smoke_submission.lock';receipt=control/'smoke_submission.json'
    if receipt.exists() or lock.exists():raise SystemExit('Smoke submission already attempted; inspect Slurm before manual recovery')
    subprocess.run(['/share/apps/anaconda3/2025.06/bin/python',str(ROOT/'verify_stage.py'),'--host'],check=True)
    subprocess.run(['/share/apps/anaconda3/2025.06/bin/python',str(ROOT/'verify_starts.py'),str(control/'starts.json')],check=True)
    lock.mkdir()
    result=subprocess.run(['sbatch','--parsable','--time=1:30:00',f'--array={chain}-{chain}%1','--export=ALL,ESTATE_RUN_MODE=smoke',str(ROOT/'launch_torch.sh')],capture_output=True,text=True,check=True)
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit():raise RuntimeError('Unparseable sbatch result; lock remains for manual queue audit')
    receipt.write_text(json.dumps(dict(status='submitted',job_id=job,chain=chain,array=f'{chain}-{chain}%1',
      starts_sha256=__import__('hashlib').sha256((control/'starts.json').read_bytes()).hexdigest(),
      stage_inventory_sha256=__import__('hashlib').sha256((ROOT/'inventory.json').read_bytes()).hexdigest(),
      submitted_epoch=time.time(),no_auto_retry=True),indent=2)+'\n')
    print(json.dumps(dict(status='submitted',job_id=job,chain=chain)))
if __name__=='__main__':main()
