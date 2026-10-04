#!/usr/bin/env python3
"""Detached, bounded manager for two authorized one-core transition workers."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[3]
DRIVER = ROOT/'code/model/experiments/birth_count_choice/transition_panel.py'
DEFAULT_PYTHON = ROOT/'output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python'
WALL_SECONDS = 21600
GIB = 1024**3
THREAD_KEYS = ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')


def write(path, value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def absolute_executable(path):
    # Resolving a venv Python symlink loses its pyvenv.cfg/site-package context.
    return str(Path(path).absolute())


def validate(plan_path, config_path, indices):
    plan=json.loads(Path(plan_path).read_text());config=json.loads(Path(config_path).read_text())
    if tuple(indices)!=(5,7):raise ValueError('Exactly the authorized local indices 5,7 are required')
    if config['plan']['sha256']!=sha(plan_path) or config['identity']!=plan['identity']:
        raise ValueError('Panel/plan identity or SHA differs')
    if config['panel_source']['sha256']!=sha(DRIVER):raise ValueError('Panel source SHA differs')
    if not 0<plan['budget']['total_seconds']<=21480:raise ValueError('Plan exceeds external six-hour allowance')
    guesses=[]
    for index,ratio in zip(indices,(.94,.98)):
        row=config['guesses'][index]
        if row['index']!=index or not math.isclose(row['psi'],plan['initial_psi']*ratio,rel_tol=0,abs_tol=1e-15):
            raise ValueError('Authorized index/preference ratio differs')
        guesses.append(dict(index=index,psi=float(row['psi'])))
    return guesses


def process_snapshot():
    # macOS ps reports RSS in KiB. Include children and the owned session group.
    output=subprocess.check_output(['ps','-axo','pid=,ppid=,rss=,pgid='],text=True)
    return {int(pid):(int(parent),int(rss)*1024,int(group))
            for pid,parent,rss,group in (line.split() for line in output.splitlines() if line.strip())}


def owned_processes(snapshot, root):
    owned={pid for pid,row in snapshot.items() if pid==root or row[2]==root}
    while True:
        descendants={pid for pid,row in snapshot.items() if row[0] in owned}
        if descendants.issubset(owned):return owned
        owned|=descendants


def memory_action(first_running, second_running, second_paused, combined):
    if not second_running:return None
    if combined>30*GIB:return 'terminate_second'
    if second_paused and (not first_running or combined<20*GIB):return 'resume_second'
    if first_running and not second_paused and combined>24*GIB:return 'pause_second'
    return None


def signal_worker(worker, sig):
    if worker['process'].poll() is not None:return
    try:os.killpg(worker['process'].pid,sig)
    except ProcessLookupError:pass


def terminate_worker(worker, reason):
    if worker['process'].poll() is not None:return
    signal_worker(worker,signal.SIGTERM)
    signal_worker(worker,signal.SIGCONT)
    worker['paused']=False;worker['termination_reason']=reason
    if worker.get('terminate_at') is None:worker['terminate_at']=time.monotonic()


def artifacts(folder):
    result={}
    for label,name in (('newest_completed','latest_completed.json'),('best_so_far','best_so_far.json')):
        paths=list(Path(folder).rglob(name))
        result[label]=str(max(paths,key=lambda p:p.stat().st_mtime)) if paths else None
    return result


def manage(out, nonce):
    spec=json.loads((out/'invocation.json').read_text())
    if spec['nonce']!=nonce:raise ValueError('Detached manager nonce differs')
    lock=os.open(out/'manager.lock',os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    os.write(lock,str(os.getpid()).encode());os.close(lock)
    guesses=validate(spec['plan'],spec['panel_config'],spec['indices'])
    if sha(spec['plan'])!=spec['plan_sha256'] or sha(spec['panel_config'])!=spec['panel_config_sha256']:
        raise ValueError('Input files changed after launch preparation')
    manager_start=time.monotonic();workers=[];stop_requested=False
    def stop(_signum,_frame):
        nonlocal stop_requested
        stop_requested=True
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    try:
        for guess in guesses:
            folder=out/f"worker_{guess['index']}";folder.mkdir()
            cache=folder/'numba_cache';cache.mkdir();mpl=folder/'matplotlib';mpl.mkdir()
            env=dict(os.environ,**{key:'1' for key in THREAD_KEYS},NUMBA_CACHE_DIR=str(cache),
                MPLCONFIGDIR=str(mpl),PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
            args=[spec['python'],str(DRIVER),'--plan',spec['plan'],'--panel-config',spec['panel_config'],
                '--panel-config-sha256',spec['panel_config_sha256'],'--index',str(guess['index']),
                '--psi',repr(guess['psi']),'--output',str(folder/'run')]
            with (folder/'worker.log').open('ab') as log:
                process=subprocess.Popen(args,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,
                    stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            worker=dict(process=process,folder=folder,args=args,index=guess['index'],psi=guess['psi'],
                started=time.monotonic(),start_epoch=time.time(),paused=False,termination_reason=None,terminate_at=None)
            workers.append(worker)
            write(folder/'launcher_start.json',dict(pid=process.pid,args=args,start_epoch=worker['start_epoch'],
                deadline_epoch=worker['start_epoch']+WALL_SECONDS,threads=1,wall_seconds=WALL_SECONDS,
                plan_sha256=spec['plan_sha256'],panel_config_sha256=spec['panel_config_sha256']))
        while True:
            now=time.monotonic();snapshot=process_snapshot()
            running=[worker['process'].poll() is None for worker in workers]
            rss=[sum(snapshot[pid][1] for pid in owned_processes(snapshot,w['process'].pid)) if active else 0
                 for w,active in zip(workers,running)]
            action=memory_action(*running,workers[1]['paused'],sum(rss))
            if action=='pause_second':signal_worker(workers[1],signal.SIGSTOP);workers[1]['paused']=True
            elif action=='resume_second':signal_worker(workers[1],signal.SIGCONT);workers[1]['paused']=False
            elif action=='terminate_second':terminate_worker(workers[1],'combined_owned_RSS_above_30_GiB')
            for worker,active in zip(workers,running):
                if not active:continue
                if stop_requested or now-worker['started']>=WALL_SECONDS or now-manager_start>=WALL_SECONDS+45:
                    terminate_worker(worker,'manager_stop_requested' if stop_requested else 'six_hour_wall_deadline')
                if worker['terminate_at'] is not None and now-worker['terminate_at']>=10:
                    signal_worker(worker,signal.SIGKILL);worker['forced_kill_after_grace']=True
            rows=[]
            for worker,amount in zip(workers,rss):
                rc=worker['process'].poll()
                rows.append(dict(pid=worker['process'].pid,index=worker['index'],psi=worker['psi'],args=worker['args'],
                    start_epoch=worker['start_epoch'],elapsed_seconds=now-worker['started'],rss_bytes=amount,
                    status='exited' if rc is not None else ('paused' if worker['paused'] else 'running'),
                    returncode=rc,termination_reason=worker['termination_reason'],
                    forced_kill_after_grace=worker.get('forced_kill_after_grace',False),**artifacts(worker['folder']/'run')))
            write(out/'status.json',dict(manager_pid=os.getpid(),updated_epoch=time.time(),workers=rows,
                combined_owned_RSS_bytes=sum(rss),memory_action=action,
                RSS_guards_GiB=dict(pause_second=24,resume_second_below=20,terminate_second=30),
                unrelated_processes_signaled=False,per_worker_wall_seconds=WALL_SECONDS,
                manager_maximum_seconds=WALL_SECONDS+60,paused_wall_time_counts=True))
            if all(w['process'].poll() is not None for w in workers):break
            if now-manager_start>=WALL_SECONDS+50:
                for worker in workers:signal_worker(worker,signal.SIGKILL)
                break
            time.sleep(1 if any(w['terminate_at'] is not None for w in workers) else 15)
    finally:
        for worker in workers:terminate_worker(worker,'manager_finalization')
        # Explicit bounded cleanup, including workers stopped by the RSS guard.
        deadline=min(time.monotonic()+10,manager_start+WALL_SECONDS+60)
        while time.monotonic()<deadline and any(w['process'].poll() is None for w in workers):time.sleep(.2)
        for worker in workers:
            if worker['process'].poll() is None:signal_worker(worker,signal.SIGKILL)
        write(out/'manager_terminal.json',dict(pid=os.getpid(),finished_epoch=time.time(),
            workers=[dict(index=w['index'],pid=w['process'].pid,returncode=w['process'].poll(),
                termination_reason=w['termination_reason']) for w in workers]))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path);parser.add_argument('--panel-config',type=Path)
    parser.add_argument('--output',type=Path);parser.add_argument('--indices',default='5,7')
    parser.add_argument('--python',type=Path,default=DEFAULT_PYTHON)
    parser.add_argument('--manage',action='store_true',help=argparse.SUPPRESS)
    parser.add_argument('--nonce',help=argparse.SUPPRESS)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args(argv)
    if args.self_test:
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp)/'base_python';base.write_text('mock executable')
            link=Path(temp)/'venv/bin/python';link.parent.mkdir(parents=True);link.symlink_to(base)
            assert absolute_executable(link)==str(link.absolute())
            assert absolute_executable(link)!=str(link.resolve())
        assert memory_action(True,True,False,25*GIB)=='pause_second'
        assert memory_action(True,True,True,19*GIB)=='resume_second'
        assert memory_action(False,True,True,23*GIB)=='resume_second'
        assert memory_action(True,True,True,31*GIB)=='terminate_second'
        assert memory_action(True,True,False,23*GIB) is None
        assert owned_processes({10:(1,2,10),11:(10,3,10),12:(11,4,12),99:(1,8,99)},10)=={10,11,12}
        print('PASS: mock memory lifecycle and owned descendant selection; no workers launched');return
    if args.output is None:parser.error('--output required')
    out=args.output.resolve()
    if args.manage:return manage(out,args.nonce)
    if args.plan is None or args.panel_config is None:parser.error('--plan and --panel-config required')
    plan,config=args.plan.resolve(),args.panel_config.resolve();indices=[int(i) for i in args.indices.split(',')]
    validate(plan,config,indices)
    if not args.python.is_file() or not os.access(args.python,os.X_OK):raise ValueError('Worker Python is not executable')
    out.mkdir(parents=True,exist_ok=False)  # Refuses prior output and duplicate launch, even with a dead PID.
    nonce=uuid.uuid4().hex
    write(out/'invocation.json',dict(nonce=nonce,plan=str(plan),panel_config=str(config),indices=indices,
        python=absolute_executable(args.python),plan_sha256=sha(plan),panel_config_sha256=sha(config),prepared_epoch=time.time()))
    with (out/'manager.log').open('ab') as log:
        process=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--manage','--output',str(out),'--nonce',nonce],
            cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    write(out/'manager_start.json',dict(pid=process.pid,started_epoch=time.time(),wall_seconds=WALL_SECONDS+60))
    print(json.dumps(dict(manager_pid=process.pid,output=str(out),worker_indices=indices)))


if __name__=='__main__':main()
