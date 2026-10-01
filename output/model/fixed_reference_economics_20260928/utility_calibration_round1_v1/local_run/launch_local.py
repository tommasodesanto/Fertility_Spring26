"""Detached, bounded local utility smokes with measured progress and RSS caps."""
import argparse,json,os,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
p=argparse.ArgumentParser();p.add_argument('--arm',choices=['floor','no_A','constant_alpha'],required=True);a=p.parse_args()
out=HERE/a.arm;out.mkdir(exist_ok=False);start=time.time();deadline=min(start+1800,1790819100.);env=os.environ.copy()
for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='2'
env.update(PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg',NUMBA_CACHE_DIR=str(out/'numba_cache'),MPLCONFIGDIR=str(out/'matplotlib'))
(out/'numba_cache').mkdir();(out/'matplotlib').mkdir()

def write(name,r):
 f=out/name;t=f.with_suffix('.tmp');t.write_text(json.dumps(r,indent=2)+'\n');t.replace(f)
with (out/'native.log').open('w') as log:
 cmd=[str(ROOT/'code/model/.venv/bin/python'),str(HERE/'bootstrap.py'),'--mode','smoke','--arm',a.arm,'--start','0','--out',str(out/'smoke'),'--deadline-seconds','1800','--deadline-epoch',str(deadline)]
 child=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 write('launch.json',dict(pid=child.pid,supervisor_pid=os.getpid(),command=cmd,start_epoch=start,deadline_epoch=deadline,threads=2,rss_cap_gib=12,source_overlay='overlay_receipt.json',runtime_scope='Two independent full utility GEs; no search'))
 peak=0;reason=None
 while child.poll() is None:
  raw=subprocess.check_output(['ps','-o','rss=','-p',str(child.pid)],text=True).strip();rss=int(raw or 0)*1024;peak=max(peak,rss)
  now=time.time();write('watchdog.json',dict(pid=child.pid,elapsed_seconds=now-start,rss_bytes=rss,peak_rss_bytes=peak,deadline_epoch=deadline))
  if rss>12*1024**3 or now>=deadline:
   reason='RSS cap exceeded' if rss>12*1024**3 else 'Hard wall deadline reached';os.killpg(child.pid,signal.SIGTERM)
   try:child.wait(timeout=10)
   except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL)
   break
  time.sleep(3)
 code=child.wait();write('terminal.json',dict(exit_code=code,reason=reason,elapsed_seconds=time.time()-start,peak_rss_bytes=peak,finished_epoch=time.time(),no_auto_retry=True))
 sys.exit(code if code>=0 else 1)
