#!/usr/bin/env python3
"""One-shot, receipt-backed Fable analysis dispatch. No model code is run here."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import uuid

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
CLAUDE = Path('/Users/tommasodesanto/.local/bin/claude')
MAX_SECONDS = 90 * 60
NEW_WEALTH = ROOT / 'output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/collection'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def atomic_json(path, value):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)

def utc(epoch):
    return dt.datetime.fromtimestamp(epoch, dt.timezone.utc).isoformat()

def supervise():
    record = json.loads((HERE / 'launch_receipt.json').read_text())
    deadline = record['deadline_epoch']
    argv = [str(CLAUDE), '-p', '--model', 'claude-fable-5-1', '--effort', 'high',
            '--session-id', record['session_id'], '--output-format', 'json',
            '--permission-mode', 'dontAsk', '--permission-prompts', 'none',
            '--tools', 'Read,Glob,Grep,Write,Edit,Bash',
            '--allowedTools', 'Read,Glob,Grep,Write,Edit,Bash(python3 *),Bash(rg *),Bash(cat *),Bash(ls *),Bash(find *),Bash(head *),Bash(sed *),Bash(stat *),Bash(shasum *)',
            ]
    result = {'supervisor_pid': os.getpid(), 'session_id': record['session_id'],
              'start_utc': utc(time.time()), 'deadline_utc': utc(deadline)}
    atomic_json(HERE / 'running.json', result)
    with (HERE / 'stdout.json').open('wb') as out, (HERE / 'stderr.log').open('wb') as err:
        try:
            proc = subprocess.Popen(argv, cwd=ROOT, stdin=subprocess.PIPE, stdout=out, stderr=err,
                                    start_new_session=True)
            proc.stdin.write((HERE / 'PROMPT.md').read_bytes())
            proc.stdin.close()
            result['claude_pid'] = proc.pid
            atomic_json(HERE / 'running.json', result)
            while proc.poll() is None and time.time() < deadline:
                time.sleep(10)
            if proc.poll() is None:
                result['timed_out'] = True
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
            result['exit_code'] = proc.returncode
        except Exception as exc:
            result['supervisor_error'] = repr(exc)
        result['end_utc'] = utc(time.time())
        result['stdout_sha256'] = digest(HERE / 'stdout.json')
        result['stderr_sha256'] = digest(HERE / 'stderr.log')
        atomic_json(HERE / 'exit_receipt.json', result)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--supervise', action='store_true', help=argparse.SUPPRESS)
    ap.add_argument('--confirm-ready', action='store_true')
    ap.add_argument('--retry-input-fix', action='store_true', help=argparse.SUPPRESS)
    ap.add_argument('--source-file', action='append', default=[], metavar='LABEL=ABS_PATH')
    ap.add_argument('--cohort-status', choices=['complete', 'partial'])
    args = ap.parse_args()
    if args.supervise:
        supervise(); return
    if args.retry_input_fix:
        old = HERE / 'exit_receipt.json'
        if not old.is_file() or json.loads(old.read_text()).get('exit_code') != 1:
            ap.error('one-shot input fix requires the failed first-attempt exit receipt')
        if 'Input must be provided' not in (HERE / 'stderr.log').read_text():
            ap.error('first-attempt error does not match the known CLI input failure')
        archive = HERE / 'attempt1_cli_input_failure'
        archive.mkdir(exist_ok=False)
        for name in ('exit_receipt.json', 'running.json', 'stdout.json', 'stderr.log', 'supervisor.log'):
            p = HERE / name
            if p.exists():
                p.replace(archive / name)
        receipt = json.loads((HERE / 'launch_receipt.json').read_text())
        receipt['retry_reason'] = 'CLI variadic allowedTools consumed positional prompt; pass prompt on stdin'
        receipt['retry_utc'] = utc(time.time())
        atomic_json(HERE / 'launch_receipt.json', receipt)
        with (HERE / 'supervisor.log').open('wb') as log:
            proc = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--supervise'],
                                    cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                    start_new_session=True)
        receipt['supervisor_pid'] = proc.pid
        atomic_json(HERE / 'launch_receipt.json', receipt)
        print(json.dumps(receipt, indent=2))
        return
    if not args.confirm_ready or not args.cohort_status:
        ap.error('Lead confirmation and explicit complete/partial cohort status required')
    if not CLAUDE.is_file():
        ap.error(f'Claude executable missing: {CLAUDE}')
    if (HERE / 'launch_receipt.json').exists():
        ap.error('Existing launch receipt; duplicate dispatch refused')
    source = {}
    for item in args.source_file:
        label, sep, raw = item.partition('=')
        if not sep or not label or not raw:
            ap.error('source-file must be LABEL=ABS_PATH')
        p = Path(raw).resolve()
        if not p.is_file() or not p.is_absolute() or p.stat().st_size > 50_000_000:
            ap.error(f'missing or oversized source file: {raw}')
        source[label] = p
    for arm in ('original_receipt', 'original_fit', 'original_parameters',
                'alternative_receipt', 'alternative_fit', 'alternative_parameters'):
        if arm not in source:
            ap.error(f'missing required source-file label {arm}')
    snap = HERE / 'input_snapshot'
    snap.mkdir(exist_ok=False)
    manifest = {'cohort_status': args.cohort_status, 'snapshot_utc': utc(time.time()),
                'sources': {}}
    for label, src in source.items():
        dest = snap / f'{label}{src.suffix}'
        shutil.copy2(src, dest)
        os.chmod(dest, 0o444)
        manifest['sources'][label] = {'source': str(src), 'snapshot': str(dest),
                                      'sha256': digest(dest), 'bytes': dest.stat().st_size}
    wealth_dest = snap / 'new_wealth'
    wealth_dest.mkdir()
    for name in ('RESULTS.md', 'verification.json', 'winner_target_fit.csv',
                 'winner_parameters.csv'):
        src = NEW_WEALTH / name
        if src.is_file():
            dest = wealth_dest / name
            shutil.copy2(src, dest)
            os.chmod(dest, 0o444)
            manifest['sources'][f'new_wealth/{name}'] = {'source': str(src),
                'snapshot': str(dest), 'sha256': digest(dest), 'bytes': dest.stat().st_size}
    atomic_json(snap / 'manifest.json', manifest)
    os.chmod(snap / 'manifest.json', 0o444)
    start = time.time()
    receipt = {'session_id': str(uuid.uuid4()), 'launch_epoch': start,
               'deadline_epoch': start + MAX_SECONDS, 'launch_utc': utc(start),
               'deadline_utc': utc(start + MAX_SECONDS),
               'prompt_sha256': digest(HERE / 'PROMPT.md'),
               'input_manifest_sha256': digest(snap / 'manifest.json'),
               'cohort_status': args.cohort_status, 'model': 'claude-fable-5-1',
               'effort': 'high', 'max_seconds': MAX_SECONDS}
    atomic_json(HERE / 'launch_receipt.json', receipt)
    with (HERE / 'supervisor.log').open('wb') as log:
        proc = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--supervise'],
                                cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    receipt['supervisor_pid'] = proc.pid
    atomic_json(HERE / 'launch_receipt.json', receipt)
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
