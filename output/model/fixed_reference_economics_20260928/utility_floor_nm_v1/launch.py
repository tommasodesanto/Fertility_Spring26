"""Six detached supervisors; two threads and four GiB RSS per native process."""
import argparse, json, os, shutil, signal, subprocess, time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];PYTHON=ROOT/'code/model/.venv/bin/python'
THREAD_KEYS=('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');tmp.replace(path)
def kill(child):
    try: os.killpg(child.pid,signal.SIGTERM)
    except ProcessLookupError: return
    try: child.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try: os.killpg(child.pid,signal.SIGKILL)
        except ProcessLookupError: pass
def supervise(chain,deadline,fast=False):
    folder=HERE/f'chain_{chain}';folder.mkdir(exist_ok=False);env=os.environ.copy()
    for key in THREAD_KEYS: env[key]='2'
    for name in ('numba_cache','matplotlib'): (folder/name).mkdir()
    env.update(NUMBA_CACHE_DIR=str(folder/'numba_cache'),MPLCONFIGDIR=str(folder/'matplotlib'),PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
    cmd=[str(PYTHON),str(HERE/'bootstrap.py'),'--chain',str(chain),'--out',str(folder/'results'),'--deadline-epoch',str(deadline)]
    if fast: cmd.append('--fast-objective')
    with (folder/'native.log').open('w') as log:
        child=subprocess.Popen(cmd,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    started=time.time();peak=0;reason=None
    write(folder/'launch.json',dict(pid=child.pid,supervisor_pid=os.getpid(),chain=chain,start_epoch=started,deadline_epoch=deadline,threads=2,rss_cap_gib=4,maximum_objective_calls=500,command=cmd,fast_objective=fast))
    stages=[]
    for stage in ('search','selected_postcheck'):
        next_progress=0
        while child.poll() is None:
            raw=subprocess.run(['ps','-o','rss=','-p',str(child.pid)],capture_output=True,text=True).stdout.strip();rss=int(raw or 0)*1024;peak=max(peak,rss);now=time.time();free=shutil.disk_usage(HERE).free
            if now>=next_progress:
                latest=folder/'results/latest.json';best=folder/'results/best_so_far.json'
                write(folder/'watchdog.json',dict(stage=stage,pid=child.pid,supervisor_pid=os.getpid(),rss_bytes=rss,peak_rss_bytes=peak,free_disk_bytes=free,elapsed_seconds=now-started,deadline_epoch=deadline,latest=json.loads(latest.read_text()) if latest.exists() else None,best_summary={'path':str(best)} if best.exists() else None));next_progress=now+60
            if rss>4*1024**3: reason='RSS_above_4GiB'
            elif free<20*1024**3: reason='free_disk_below_20GiB'
            elif now>=deadline: reason='four_hour_hard_stop'
            if reason: kill(child);break
            time.sleep(3)
        stages.append(dict(stage=stage,pid=child.pid,exit_code=child.wait()))
        if stage=='selected_postcheck' or reason or child.returncode!=0: break
        selected_file=folder/'results/search_completed.json'
        if not selected_file.exists() or json.loads(selected_file.read_text())['selected'] is None: break
        if time.time()>=deadline: reason='four_hour_hard_stop';break
        verify=[str(PYTHON),str(HERE/'bootstrap.py'),'--chain',str(chain),'--out',str(folder/'postcheck_results'),'--deadline-epoch',str(deadline),'--verify-only',str(selected_file)]
        with (folder/'postcheck.log').open('w') as log:
            child=subprocess.Popen(verify,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        write(folder/'postcheck_launch.json',dict(pid=child.pid,supervisor_pid=os.getpid(),deadline_epoch=deadline,command=verify,threads=2,rss_cap_gib=4))
    write(folder/'terminal.json',dict(exit_code=child.returncode,stages=stages,termination_reason=reason,elapsed_seconds=time.time()-started,peak_rss_bytes=peak,no_auto_retry=True))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--supervise',type=int,choices=range(6));ap.add_argument('--deadline-epoch',type=float);ap.add_argument('--fast-objective',action='store_true');args=ap.parse_args()
    if args.supervise is not None:
        supervise(args.supervise,args.deadline_epoch,args.fast_objective);return
    if (HERE/'launch.json').exists(): raise RuntimeError('Refusing duplicate launch')
    deadline=time.time()+14400;items=[]
    for chain in range(6):
        cmd=[str(PYTHON),str(HERE/'launch.py'),'--supervise',str(chain),'--deadline-epoch',str(deadline)]
        if args.fast_objective: cmd.append('--fast-objective')
        with (HERE/f'supervisor_{chain}.log').open('w') as log:
            process=subprocess.Popen(cmd,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        items.append(dict(chain=chain,supervisor_pid=process.pid))
    write(HERE/'launch.json',dict(start_epoch=time.time(),deadline_epoch=deadline,chains=items,threads_per_chain=2,total_threads=12,total_rss_cap_gib=24,maximum_objective_calls_per_chain=500,maximum_total_full_GE=3000,final_postchecks_additional=6,fast_objective=args.fast_objective,no_baseline_gate=True,no_auto_retry=True))
if __name__=='__main__': main()
