"""One approved subprocess launch; external 8GiB RSS/600s guard on macOS.

RSS is resident memory, distinct from unsupported RLIMIT_AS virtual memory.
Independent ps polling remains responsive during Numba execution holding GIL.
No numerical retry, changed economic input, or relaxed native gate.
"""
from pathlib import Path
import argparse, json, os, signal, subprocess, sys, time
HERE=Path(__file__).resolve().parent
CAP_KIB=8*1024**2; SECONDS=600; POLL=.25

def write(v): (HERE/'phi095_supervision.json').write_text(json.dumps(v,indent=2)+'\n')

def terminate_group(proc):
    try: os.killpg(proc.pid,signal.SIGTERM)
    except ProcessLookupError: return
    try: proc.wait(timeout=1)
    except subprocess.TimeoutExpired:
        try: os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError: pass
        proc.wait(timeout=2)

def main():
    global HERE
    ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['credit095','price110'],default='credit095');args=ap.parse_args()
    prefix='phi095' if args.kind=='credit095' else 'price110'; case='phi_095_run1' if args.kind=='credit095' else 'price_110_run1'
    if (HERE/case).exists(): raise RuntimeError('No retry: case directory exists')
    if (HERE/(prefix+'_supervision.json')).exists(): raise RuntimeError('No retry: supervisor receipt exists')
    env=dict(os.environ)
    for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        if env.get(key)!='1': raise RuntimeError('All thread caps must be 1')
    env['CREDIT_RSS_SUPERVISED']='1'
    command=[sys.executable,str(HERE/'paired_phi095.py'),'--run','--supervised-rss','--kind',args.kind]
    start=time.monotonic(); maximum=0; samples=0; next_poll=start
    def write(v): (HERE/(prefix+'_supervision.json')).write_text(json.dumps(v,indent=2)+'\n')
    with (HERE/(prefix+'_run1.log')).open('w') as log:
        proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        record=dict(status='running',pid=proc.pid,command=command,cap_kib=CAP_KIB,cap_meaning='8GiB subprocess resident memory RSS, not virtual address space',deadline_seconds=SECONDS,poll_target_seconds=POLL,time_epoch=time.time())
        write(record)
        failure=None
        try:
            while proc.poll() is None:
                elapsed=time.monotonic()-start
                if elapsed>SECONDS: failure='600-second supervisor deadline'; break
                raw=subprocess.run(['ps','-o','rss=','-p',str(proc.pid)],capture_output=True,text=True,timeout=.4)
                if proc.poll() is not None: break
                if raw.returncode!=0 or not raw.stdout.strip(): failure='RSS measurement unavailable'; break
                rss=int(raw.stdout.strip()); maximum=max(maximum,rss); samples+=1
                if rss>CAP_KIB: failure='8GiB RSS cap exceeded'; break
                if samples%20==0: write({**record,'elapsed_seconds':elapsed,'peak_rss_kib':maximum,'samples':samples})
                next_poll+=POLL; time.sleep(max(0,next_poll-time.monotonic()))
            if failure: terminate_group(proc)
            else: proc.wait(timeout=1)
        except Exception as exc:
            failure=repr(exc); terminate_group(proc)
        record.update(status='supervisor_failed' if failure else ('completed' if proc.returncode==0 else 'child_failed'),error=failure,exit_code=proc.returncode,elapsed_seconds=time.monotonic()-start,peak_rss_kib=maximum,samples=samples)
        write(record)
    print(json.dumps(record,indent=2))
    if record['status']!='completed': sys.exit(1)

if __name__=='__main__': main()
