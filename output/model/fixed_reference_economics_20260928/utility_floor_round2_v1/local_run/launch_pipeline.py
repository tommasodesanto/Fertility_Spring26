"""Fresh authenticated gate, then two bounded local floor calibration chains."""
import json,os,signal,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4];PYTHON=ROOT/'code/model/.venv/bin/python';START=time.time();DEADLINE=START+14400
ENV=os.environ.copy()
for k in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):ENV[k]='2'
ENV.update(PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
def write(p,r):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.tmp');t.write_text(json.dumps(r,indent=2)+'\n');t.replace(p)
def spawn(mode,folder,start=0,receipt=None):
 out=HERE/folder;out.mkdir(exist_ok=False);env=ENV.copy();env['NUMBA_CACHE_DIR']=str(out/'numba_cache');env['MPLCONFIGDIR']=str(out/'matplotlib');(out/'numba_cache').mkdir();(out/'matplotlib').mkdir()
 cmd=[str(PYTHON),str(HERE/'bootstrap.py'),'--mode',mode,'--arm','floor','--start',str(start),'--out',str(out/'results'),'--deadline-seconds','1800' if mode=='smoke' else '14400','--deadline-epoch',str(DEADLINE)]
 if receipt:cmd+=['--smoke-receipt',str(receipt)]
 log=(out/'native.log').open('w');child=subprocess.Popen(cmd,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);log.close();r=dict(pid=child.pid,supervisor_pid=os.getpid(),mode=mode,start_index=start,start_epoch=time.time(),deadline_epoch=DEADLINE,threads=2,rss_cap_gib=12,command=cmd);write(out/'launch.json',r);return out,child,r

def monitor(items):
 peaks={c.pid:0 for _,c,_ in items}
 while items:
  for out,c,r in list(items):
   code=c.poll()
   if code is not None:write(out/'terminal.json',dict(exit_code=code,elapsed_seconds=time.time()-r['start_epoch'],peak_rss_bytes=peaks[c.pid],no_auto_retry=True));items.remove((out,c,r));continue
   raw=subprocess.run(['ps','-o','rss=','-p',str(c.pid)],capture_output=True,text=True).stdout.strip();rss=int(raw or 0)*1024;peaks[c.pid]=max(peaks[c.pid],rss);now=time.time();write(out/'watchdog.json',dict(pid=c.pid,watchdog_pid=os.getpid(),rss_bytes=rss,peak_rss_bytes=peaks[c.pid],elapsed_seconds=now-r['start_epoch'],deadline_epoch=DEADLINE))
   if rss>12*1024**3 or now>=DEADLINE:
    os.killpg(c.pid,signal.SIGTERM)
    try:c.wait(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(c.pid,signal.SIGKILL)
  time.sleep(3)
write(HERE/'pipeline_launch.json',dict(supervisor_pid=os.getpid(),start_epoch=START,deadline_epoch=DEADLINE,threads_per_chain=2,total_search_threads=4,stages='fresh2GEsmoke then two concurrent80GE8roundsearches',seed0='verified floor_s2 winner',seed1='historical floor plus documented concave trial; current fixed economy'))
gate=spawn('smoke','fresh_gate');monitor([gate]);require=gate[1].returncode==0
if not require:raise RuntimeError('Fresh smoke gate failed; both local searches blocked')
receipt=gate[0]/'results/smoke_receipt.json';chains=[spawn('run','winner_chain',0,receipt),spawn('run','historical_chain',1,receipt)];write(HERE/'search_launches.json',dict(chains=[r for _,_,r in chains],smoke_receipt=str(receipt)));monitor(chains)
write(HERE/'pipeline_terminal.json',dict(finished_epoch=time.time(),elapsed_seconds=time.time()-START,no_auto_retry=True))
