"""Enforce the author's September13 noon cutoff on explicit owned Slurm jobs."""
import argparse
import json
from pathlib import Path
import subprocess
import time

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);args=ap.parse_args()
    plan=json.loads(args.plan.read_text());out=Path(plan['output']);out.mkdir(parents=True,exist_ok=True)
    def save(name,value):
        p=out/name;t=p.with_suffix('.tmp');t.write_text(json.dumps(value,indent=2)+'\n');t.replace(p)
    for stage in plan['stages']:
        ids=stage['job_ids'];assert ids and all(str(x).isdigit() for x in ids)
        while time.time()<stage['deadline_unix']:
            save('heartbeat.json',dict(next_stage=stage['name'],remaining_seconds=stage['deadline_unix']-time.time()))
            time.sleep(min(60,max(0,stage['deadline_unix']-time.time())))
        before=subprocess.run(['squeue','-h','-j',','.join(map(str,ids)),'-o','%i|%u|%j|%T'],capture_output=True,text=True)
        lines=[x.split('|') for x in before.stdout.splitlines() if x.strip()]
        assert all(len(x)==4 and x[1]=='td2248' for x in lines),'Unexpected job owner; no cancellation'
        result=subprocess.run(['scancel',*map(str,ids)],capture_output=True,text=True)
        save(stage['name']+'.json',dict(stage=stage,observed_epoch=time.time(),queue_before=lines,
            queue_stderr=before.stderr,cancel_returncode=result.returncode,cancel_stdout=result.stdout,cancel_stderr=result.stderr,
            meaning='Author deadline; unfinished numerical results are not promoted. Completed checkpoints remain.'))
    save('summary.json',dict(status='cutoff_actions_finished',observed_epoch=time.time()))

if __name__=='__main__':main()
