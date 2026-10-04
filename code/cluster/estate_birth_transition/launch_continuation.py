#!/usr/bin/env python3
"""Detached, bounded launcher for one authorized estate-transition continuation."""
from __future__ import annotations
import argparse, hashlib, json, os, signal, subprocess, sys, time
from pathlib import Path

GIB = 1024**3
THREADS = ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, value):
    path=Path(path); tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n'); tmp.replace(path)

def lexical_abs(path):
    # Preserve the venv/bin/python spelling; resolving a symlink loses venv context.
    p=Path(path)
    return str(p if p.is_absolute() else Path.cwd()/p)

def snapshot():
    # macOS ps RSS is KiB. A worker owns its new process group and descendants.
    raw=subprocess.check_output(['ps','-axo','pid=,ppid=,rss=,pgid='],text=True)
    return {int(p):(int(pp),int(rss)*1024,int(pg)) for p,pp,rss,pg in
            (line.split() for line in raw.splitlines() if line.strip())}

def owned(snap, root):
    result={p for p,row in snap.items() if p==root or row[2]==root}
    while True:
        more={p for p,row in snap.items() if row[0] in result}
        if more <= result:return result
        result |= more

def limit_reason(now, deadline, rss, memory_bytes):
    if now >= deadline:return 'wall_timeout'
    if rss > memory_bytes:return 'owned_process_RSS_limit'
    return None

def validate_pin(value):
    pin=json.loads(value) if isinstance(value,str) else value
    if not isinstance(pin,dict) or set(pin)!= {'path','sha256'}:
        raise ValueError('smoke receipt pin must contain exactly path and sha256')
    if not isinstance(pin['path'],str) or not Path(pin['path']).is_file():
        raise ValueError('pinned smoke receipt path must exist')
    if not isinstance(pin['sha256'],str) or sha(pin['path'])!=pin['sha256']:
        raise ValueError('pinned smoke receipt SHA-256 differs')
    return {'path':pin['path'],'sha256':pin['sha256']}

def self_test():
    fake={10:(1,2,10),11:(10,3,10),12:(11,4,12),99:(1,8,99)}
    assert owned(fake,10)=={10,11,12}
    assert lexical_abs('/tmp/venv/bin/python')=='/tmp/venv/bin/python'
    start=100.0; assert 106.0-start>=5.0 and 104.0-start<5.0
    assert limit_reason(106,105,1,10)=='wall_timeout'
    assert limit_reason(104,105,11,10)=='owned_process_RSS_limit'
    assert limit_reason(104,105,10,10) is None
    assert len(THREADS)==6 and GIB==1024**3
    print('PASS: owned process selection, lexical venv path, monotonic deadline arithmetic; zero children')

def manager(spec_path, nonce):
    out=Path(spec_path).parent; spec=json.loads(Path(spec_path).read_text())
    if spec['nonce']!=nonce: raise ValueError('manager nonce mismatch')
    # Claim this invocation before any setup; direct or duplicate manager calls fail closed.
    fd=os.open(out/'manager.lock',os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    os.write(fd,str(os.getpid()).encode()); os.close(fd)
    driver,manifest=spec['driver'],spec['manifest']; start=time.monotonic(); epoch=time.time()
    stop=False; child=None; reason='completed'; forced=False
    def on_stop(_sig,_frame):
        nonlocal stop
        stop=True
    signal.signal(signal.SIGTERM,on_stop); signal.signal(signal.SIGINT,on_stop)
    child=None; rc=None
    try:
        if sha(driver)!=spec['driver_sha256'] or sha(manifest)!=spec['manifest_sha256']:
            raise ValueError('driver or manifest changed before manager launch')
        if spec['smoke_receipt_pin'] is not None:
            validate_pin(spec['smoke_receipt_pin'])
        run=out/'run'; cache=out/'numba_cache'; mpl=out/'matplotlib'
        cache.mkdir(); mpl.mkdir()
        env=dict(os.environ,**{k:'1' for k in THREADS}, NUMBA_CACHE_DIR=str(cache),
                 MPLCONFIGDIR=str(mpl), PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
        args=[spec['python'],driver,'--manifest',manifest,'--'+spec['mode'],'--output',str(run)]
        if spec['mode']=='run':
            args.extend(['--smoke-receipt-pin',json.dumps(spec['smoke_receipt_pin'],separators=(',',':'))])
        with (out/'worker.log').open('ab') as log:
            child=subprocess.Popen(args,cwd=spec['cwd'],env=env,stdin=subprocess.DEVNULL,
                                   stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        child_start=time.monotonic(); deadline=child_start+spec['wall_seconds']; term_at=None
        write(out/'launcher_start.json',dict(manager_pid=os.getpid(),worker_pid=child.pid,
             command=args,started_epoch=time.time(),wall_seconds=spec['wall_seconds'],
             memory_gib=spec['memory_gib'],driver_sha256=spec['driver_sha256'],
             manifest_sha256=spec['manifest_sha256']))
        while child.poll() is None:
            now=time.monotonic(); snap=snapshot()
            pids=owned(snap,child.pid); rss=sum(snap[p][1] for p in pids if p in snap)
            guard=limit_reason(now,deadline,rss,spec['memory_gib']*GIB)
            if (stop or guard) and term_at is None:
                reason='manager_termination' if stop else guard
                try: os.killpg(child.pid,signal.SIGTERM)
                except ProcessLookupError: pass
                term_at=now
            if term_at is not None and now-term_at>=10:
                try: os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError: pass
                forced=True
            write(out/'status.json',dict(manager_pid=os.getpid(),worker_pid=child.pid,
                 updated_epoch=time.time(),elapsed_seconds=now-child_start,deadline_epoch=epoch+spec['wall_seconds'],
                 owned_pids=sorted(pids),owned_rss_bytes=rss,memory_limit_bytes=spec['memory_gib']*GIB,
                 status='terminating' if term_at is not None else 'running',reason=reason))
            time.sleep(1 if term_at is not None else 15)
        rc=child.wait()
        if rc and reason=='completed': reason='worker_exit_nonzero'
    except BaseException as exc:
        reason='manager_error: '+repr(exc); rc=child.poll() if child else None
        if child and rc is None:
            try: os.killpg(child.pid,signal.SIGTERM)
            except ProcessLookupError: pass
            try: child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try: os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError: pass
                forced=True
            rc=child.poll()
    finally:
        write(out/'manager_terminal.json',dict(manager_pid=os.getpid(),finished_epoch=time.time(),
             worker_pid=child.pid if child else None,exit_code=rc if child else None,
             reason=reason,forced_kill_after_grace=forced))

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--python'); p.add_argument('--driver'); p.add_argument('--manifest')
    p.add_argument('--mode',choices=('smoke','run')); p.add_argument('--output')
    p.add_argument('--wall-seconds',type=int); p.add_argument('--memory-gib',type=int)
    p.add_argument('--smoke-receipt-pin',help='JSON object with exactly path and sha256 (required for run)')
    p.add_argument('--self-test',action='store_true'); p.add_argument('--_manage',action='store_true',help=argparse.SUPPRESS)
    p.add_argument('--_nonce',help=argparse.SUPPRESS)
    a=p.parse_args(argv)
    if a.self_test:return self_test()
    if a._manage:return manager(a.output,a._nonce)
    if not all((a.python,a.driver,a.manifest,a.mode,a.output,a.wall_seconds,a.memory_gib)):
        p.error('all launcher arguments are required')
    if not Path(a.python).is_absolute():p.error('--python must be a lexical absolute venv path')
    if not Path(a.python).is_file() or not os.access(a.python,os.X_OK):p.error('--python is not executable')
    if not Path(a.driver).is_file() or not Path(a.manifest).is_file():p.error('driver and manifest must exist')
    if a.mode=='run' and a.smoke_receipt_pin is None:p.error('--smoke-receipt-pin is required for --mode run')
    if a.mode=='smoke' and a.smoke_receipt_pin is not None:p.error('--smoke-receipt-pin is only valid for --mode run')
    pin=None
    if a.mode=='run':
        try: pin=validate_pin(a.smoke_receipt_pin)
        except (ValueError,TypeError,json.JSONDecodeError) as exc:p.error(str(exc))
    maximum=21600
    if not 1<=a.wall_seconds<=maximum:p.error(f'--wall-seconds must be 1..{maximum} for {a.mode}')
    if not 1<=a.memory_gib<=24:p.error('--memory-gib must be 1..24')
    driver,manifest=str(Path(a.driver).absolute()),str(Path(a.manifest).absolute())
    out=Path(a.output).absolute(); out.mkdir(parents=True,exist_ok=False)
    import uuid
    nonce=uuid.uuid4().hex
    spec=dict(nonce=nonce,python=lexical_abs(a.python),driver=driver,manifest=manifest,
              driver_sha256=sha(driver),manifest_sha256=sha(manifest),mode=a.mode,
              wall_seconds=a.wall_seconds,memory_gib=a.memory_gib,cwd=str(Path.cwd()),
              smoke_receipt_pin=pin)
    write(out/'invocation.json',spec)
    with (out/'manager.log').open('ab') as log:
        proc=subprocess.Popen([sys.executable,str(Path(__file__).absolute()),'--_manage','--output',str(out/'invocation.json'),
             '--_nonce',nonce],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    write(out/'manager_start.json',dict(manager_pid=proc.pid,started_epoch=time.time(),
         maximum_seconds=a.wall_seconds+60,mode=a.mode,output=str(out)))
    # Bounded assertion prevents the host sleeping while this manager is active.
    awake_args=['caffeinate','-i','-w',str(proc.pid),'-t',str(a.wall_seconds+60)]
    awake=subprocess.Popen(awake_args,stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL,start_new_session=True)
    write(out/'keep_awake.json',dict(pid=awake.pid,command=awake_args,
         manager_pid=proc.pid,maximum_seconds=a.wall_seconds+60,mode=a.mode))
    print(json.dumps(dict(manager_pid=proc.pid,output=str(out),mode=a.mode,
                          wall_seconds=a.wall_seconds,memory_gib=a.memory_gib)))

if __name__=='__main__':main()
