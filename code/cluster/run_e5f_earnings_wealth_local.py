"""Run a frozen local smoke, then its dependent bounded earnings search.

The caller keeps this supervisor awake with caffeinate. No stage is retried.
"""
from __future__ import annotations
import argparse, datetime, json, os, signal, subprocess, sys, time
from pathlib import Path


def write(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');temp.replace(path)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,required=True)
    args=ap.parse_args();root=args.bundle.resolve();plan=json.loads((root/'plan.json').read_text())
    sys.path.insert(0,str(root/'tools'))
    from run_e5f_earnings_wealth_smoke import _kill_group
    from run_e5f_earnings_wealth_candidate import verify_plan
    verify_plan(plan)
    supervision=plan.get('local_supervision',{})
    stages=[('smoke',int(supervision.get('smoke_seconds',7000))),('search',int(supervision.get('search_seconds',21000)))]
    total_limit=int(supervision.get('total_seconds',28000))
    if min([total_limit]+[limit for _,limit in stages])<=0:raise ValueError('positive supervision budgets required')
    out=root/'output';out.mkdir(exist_ok=False)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1',NUMBA_DISABLE_JIT='0',MPLBACKEND='Agg',PYTHONUNBUFFERED='1',NUMBA_CACHE_DIR=str(root/'numba_cache'))
    receipt={'schema':'earnings_local_overnight_v1','pid':os.getpid(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bundle':str(root),'status':'running','stages':[],'maximum_total_seconds':total_limit,'retry_policy':'no automatic retries'}
    started=time.monotonic();write(root/'execution.json',receipt)
    signal.signal(signal.SIGTERM,lambda *_: (_ for _ in ()).throw(KeyboardInterrupt('supervisor terminated')))
    for mode,limit in stages:
        command=[plan['python'],'-B',str(root/'tools/run_e5f_earnings_wealth_search.py'),'--mode',mode,'--plan',str(root/'plan.json'),'--output',str(out/mode)]
        if mode=='search':command+=['--verified-smoke',str(out/'smoke/smoke_receipt.json')]
        phase={'stage':mode,'status':'running','seconds_limit':limit,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};receipt['stages'].append(phase);write(root/'execution.json',receipt)
        t=time.monotonic();proc=None;last=0.
        try:
            with (root/f'{mode}.log').open('w') as log:
                proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                phase['pid']=proc.pid;write(root/'execution.json',receipt)
                while proc.poll() is None:
                    if time.monotonic()-t>limit or time.monotonic()-started>total_limit:raise TimeoutError('declared stage/total allocation exhausted')
                    if time.monotonic()-last>=30:
                        write(root/'heartbeat.json',{'status':'running','stage':mode,'elapsed_seconds':time.monotonic()-t,'updated_unix':time.time(),'pid':proc.pid});last=time.monotonic()
                    time.sleep(1)
                if proc.returncode:raise RuntimeError(f'{mode} exited {proc.returncode}; see {mode}.log and native receipts')
            expected=out/mode/('smoke_receipt.json' if mode=='smoke' else 'summary.json')
            result=json.loads(expected.read_text());required='verified_smoke' if mode=='smoke' else 'verified_selection'
            if result['status']!=required:raise RuntimeError(f'{mode} did not produce {required}')
            phase.update(status='completed',elapsed_seconds=time.monotonic()-t,receipt=str(expected))
        except BaseException as exc:
            if proc is not None and proc.poll() is None:_kill_group(proc)
            phase.update(status='failed_or_stopped',elapsed_seconds=time.monotonic()-t,error=str(exc));receipt['status']='stopped';write(root/'execution.json',receipt)
            write(root/'heartbeat.json',{'status':'stopped','stage':mode,'updated_unix':time.time()});return 2
        write(root/'execution.json',receipt)
    receipt.update(status='verified_selection_pending_lead_review',elapsed_seconds=time.monotonic()-started);write(root/'execution.json',receipt)
    write(root/'heartbeat.json',{'status':'completed','updated_unix':time.time()});return 0

if __name__=='__main__':raise SystemExit(main())
