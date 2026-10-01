"""Attach RSS/wall watchdog to already detached floor process; no model run."""
import json,os,signal,subprocess,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;out=HERE/'floor';r=json.loads((out/'launch.json').read_text());pid=r['pid'];peak=0;reason=None
while True:
 raw=subprocess.run(['ps','-o','rss=','-p',str(pid)],capture_output=True,text=True).stdout.strip()
 if not raw:break
 now=time.time();rss=int(raw)*1024;peak=max(peak,rss)
 f=out/'watchdog.json';t=f.with_suffix('.tmp');t.write_text(json.dumps(dict(pid=pid,watchdog_pid=os.getpid(),elapsed_seconds=now-r['start_epoch'],rss_bytes=rss,peak_rss_bytes=peak,deadline_epoch=r['deadline_epoch']),indent=2)+'\n');t.replace(f)
 if rss>12*1024**3 or now>=r['deadline_epoch']:
  reason='RSS cap exceeded' if rss>12*1024**3 else 'Hard wall deadline reached';os.killpg(pid,signal.SIGTERM);time.sleep(10)
  try:os.killpg(pid,signal.SIGKILL)
  except ProcessLookupError:pass
  break
 time.sleep(3)
(out/'watchdog_terminal.json').write_text(json.dumps(dict(reason=reason,elapsed_seconds=time.time()-r['start_epoch'],peak_rss_bytes=peak,model_exit_code='See smoke completed.json or failure.json'),indent=2)+'\n')
