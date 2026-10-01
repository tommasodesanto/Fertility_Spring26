"""Detached one-child local launcher with 20-minute, 12-GiB RSS supervision."""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, signal, subprocess, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
PYTHON = ROOT / 'code/model/.venv/bin/python'
THREAD_KEYS = ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
               'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
RSS_LIMIT = 12 * 1024**3
WALL_LIMIT = 1200


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    temp.replace(path)


def kill_group(pid: int):
    try:
        os.killpg(pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    time.sleep(5)
    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--go-reviewed', action='store_true', required=True,
                    help='Lead approval is required before starting the one full GE')
    ap.add_argument('--out', type=Path, required=True,
                    help='Fresh isolated result directory')
    ap.add_argument('--deadline-epoch', type=float,
                    help='Reuse a previously approved absolute stop time')
    args = ap.parse_args()
    if not (ROOT / 'code/model/.venv/bin/python').is_file() or not (ROOT / 'CALIBRATION_STATUS.md').is_file():
        raise RuntimeError('Repository root or model Python runtime not found')
    started = time.time()
    deadline = min(started + WALL_LIMIT, args.deadline_epoch or float('inf'))
    out = args.out.resolve()
    if out.exists():
        raise RuntimeError('Refusing existing result directory: ' + str(out))
    out.mkdir(parents=True)
    runtime = out / 'runtime'
    (runtime / 'numba_cache').mkdir(parents=True)
    (runtime / 'matplotlib').mkdir()
    env = os.environ.copy()
    for key in THREAD_KEYS:
        env[key] = '2'
    env.update(NUMBA_CACHE_DIR=str(runtime / 'numba_cache'),
               MPLCONFIGDIR=str(runtime / 'matplotlib'),
               PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    cmd = [str(PYTHON), str(HERE / 'prepare.py'), '--run-once', '--go-reviewed',
           '--deadline-epoch', str(deadline), '--out', str(out / 'results')]
    log_path = out / 'native.log'
    with log_path.open('w') as log:
        child = subprocess.Popen(cmd, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
                                stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    write(out / 'launch.json', dict(pid=child.pid, supervisor_pid=os.getpid(),
                                    start_epoch=started, deadline_epoch=deadline,
                                    maximum_wall_seconds=WALL_LIMIT, threads=2,
                                    rss_cap_gib=12, retries=0, command=cmd))
    peak = 0
    reason = None
    while child.poll() is None:
        raw = subprocess.run(['ps', '-o', 'rss=', '-p', str(child.pid)],
                             capture_output=True, text=True).stdout.strip()
        rss = int(raw or 0) * 1024
        peak = max(peak, rss)
        now = time.time()
        write(out / 'watchdog.json', dict(pid=child.pid, supervisor_pid=os.getpid(),
                                          rss_bytes=rss, peak_rss_bytes=peak,
                                          elapsed_seconds=now-started,
                                          deadline_epoch=deadline))
        if rss > RSS_LIMIT:
            reason = 'RSS_above_12GiB'
        elif now >= deadline:
            reason = '20_minute_hard_stop'
        if reason:
            kill_group(child.pid)
            break
        time.sleep(3)
    code = child.wait()
    write(out / 'terminal.json', dict(pid=child.pid, exit_code=code,
                                      termination_reason=reason,
                                      elapsed_seconds=time.time()-started,
                                      peak_rss_bytes=peak, no_auto_retry=True))


if __name__ == '__main__':
    main()
