#!/usr/bin/env python3
"""Single same-session Fable revision with an absolute deadline."""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parents[4]
CLAUDE = '/Users/tommasodesanto/.local/bin/claude'

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()

def write_json(p, x):
    t = p.with_suffix(p.suffix + '.tmp')
    t.write_text(json.dumps(x, indent=2, sort_keys=True) + '\n')
    t.replace(p)

def iso(t):
    return dt.datetime.fromtimestamp(t, dt.timezone.utc).isoformat()

def supervise():
    r = json.loads((HERE / 'launch_receipt.json').read_text())
    argv = [CLAUDE, '-p', '--resume', r['session_id'], '--model', 'claude-fable-5-1',
            '--effort', 'high', '--output-format', 'json', '--permission-mode', 'dontAsk',
            '--permission-prompts', 'none', '--tools', 'Read,Glob,Grep,Write,Edit,Bash',
            '--allowedTools', 'Read,Glob,Grep,Write,Edit,Bash(rg *),Bash(cat *),Bash(sed *),Bash(head *),Bash(python3 -c *)']
    state = {'supervisor_pid': os.getpid(), 'session_id': r['session_id'],
             'start_utc': iso(time.time()), 'deadline_utc': r['deadline_utc']}
    write_json(HERE / 'running.json', state)
    with (HERE / 'stdout.json').open('wb') as out, (HERE / 'stderr.log').open('wb') as err:
        try:
            p = subprocess.Popen(argv, cwd=ROOT, stdin=subprocess.PIPE,
                                 stdout=out, stderr=err, start_new_session=True)
            p.stdin.write((HERE / 'REVISION_PROMPT.md').read_bytes())
            p.stdin.close()
            state['claude_pid'] = p.pid
            write_json(HERE / 'running.json', state)
            while p.poll() is None and time.time() < r['deadline_epoch']:
                time.sleep(10)
            if p.poll() is None:
                state['timed_out'] = True
                os.killpg(p.pid, signal.SIGTERM)
                try:
                    p.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                    p.wait()
            state['exit_code'] = p.returncode
        except Exception as exc:
            state['supervisor_error'] = repr(exc)
        state['end_utc'] = iso(time.time())
        state['stdout_sha256'] = sha(HERE / 'stdout.json')
        state['stderr_sha256'] = sha(HERE / 'stderr.log')
        write_json(HERE / 'exit_receipt.json', state)

if __name__ == '__main__':
    if len(sys.argv) == 2 and sys.argv[1] == '--supervise':
        supervise()
    elif len(sys.argv) == 1:
        if (HERE / 'launch_receipt.json').exists():
            raise SystemExit('duplicate revision refused')
        parent = json.loads((BASE / 'launch_receipt.json').read_text())
        now = time.time()
        deadline = min(now + 20*60, parent['deadline_epoch'])
        if deadline <= now + 60:
            raise SystemExit('insufficient time before original deadline')
        r = {'session_id': parent['session_id'], 'start_epoch': now,
             'start_utc': iso(now), 'deadline_epoch': deadline,
             'deadline_utc': iso(deadline), 'max_seconds': 20*60,
             'original_deadline_utc': parent['deadline_utc'],
             'prompt_sha256': sha(HERE / 'REVISION_PROMPT.md'),
             'original_memo_sha256': json.loads((HERE / 'pre_revision_manifest.json').read_text())['files']['ECONOMIC_MEMO.md']['sha256'],
             'model': 'claude-fable-5-1', 'effort': 'high'}
        write_json(HERE / 'launch_receipt.json', r)
        with (HERE / 'supervisor.log').open('wb') as log:
            p = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--supervise'],
                                 cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        r['supervisor_pid'] = p.pid
        write_json(HERE / 'launch_receipt.json', r)
        print(json.dumps(r, indent=2))
    else:
        raise SystemExit('usage: run_revision.py [--supervise]')
