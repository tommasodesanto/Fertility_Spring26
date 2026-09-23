"""Bounded, first-party Claude Max review; one explicit phase per invocation."""
from pathlib import Path
import argparse, datetime, json, os, signal, subprocess, time

ROOT = Path(__file__).resolve().parents[6]
BASE = Path(__file__).resolve().parent
CLI = '/Users/tommasodesanto/.local/bin/claude'

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--phase', required=True)
    p.add_argument('--seconds', type=int, default=4500)
    p.add_argument('--max-turns', type=int, default=100)
    p.add_argument('--resume')
    a = p.parse_args()
    out = BASE / a.phase
    if out.parent != BASE or not (out / 'prompt.md').is_file():
        raise SystemExit('Use an existing direct-child phase directory with prompt.md')
    if (out / 'receipt.json').exists():
        raise SystemExit('Phase already attempted; preserve evidence and use a new reviewed phase')
    env = os.environ.copy()
    for k in ('ANTHROPIC_API_KEY', 'ANTHROPIC_AUTH_TOKEN', 'ANTHROPIC_BASE_URL'):
        env.pop(k, None)
    auth = subprocess.run([CLI, 'auth', 'status'], env=env, capture_output=True, text=True, check=True)
    info = json.loads(auth.stdout)
    if not all((info.get('loggedIn'), info.get('authMethod') == 'claude.ai',
                info.get('apiProvider') == 'firstParty', info.get('subscriptionType') == 'max')):
        raise SystemExit('Claude Max first-party authentication required; no fallback')
    cmd = [CLI, '--model', 'claude-opus-5-5', '--effort', 'max', '--print',
           '--output-format', 'stream-json', '--verbose', '--max-turns', str(a.max_turns),
           '--tools', 'Read,Glob,Grep,WebSearch,WebFetch', '--allowedTools',
           'Read,Glob,Grep,WebSearch,WebFetch', '--permission-mode', 'dontAsk', '--strict-mcp-config']
    if a.resume:
        cmd += ['--resume', a.resume]
    receipt = dict(status='running', started_utc=now(), model_requested='claude-opus-5-5',
                   auth='claude.ai Max firstParty', seconds=a.seconds, max_turns=a.max_turns,
                   resumed_session=a.resume, supervisor_pid=os.getpid())
    def save():
        tmp = out / 'receipt.tmp'
        tmp.write_text(json.dumps(receipt, indent=2) + '\n')
        tmp.replace(out / 'receipt.json')
    save()
    start = time.monotonic()
    with (out / 'stream.jsonl').open('w') as stream, (out / 'stderr.txt').open('w') as err:
        child = subprocess.Popen(cmd, cwd=ROOT, env=env, stdin=subprocess.PIPE, stdout=stream,
                                 stderr=err, text=True, start_new_session=True)
        receipt['pid'] = child.pid
        save()
        child.stdin.write((out / 'prompt.md').read_text())
        child.stdin.close()
        while child.poll() is None and time.monotonic() - start < a.seconds:
            receipt.update(heartbeat_utc=now(), elapsed_seconds=round(time.monotonic()-start, 1))
            save()
            time.sleep(15)
        if child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            receipt['status'] = 'timed_out'
        else:
            receipt['status'] = 'completed' if child.returncode == 0 else 'failed'
    final = None
    for line in (out / 'stream.jsonl').read_text().splitlines():
        try:
            item = json.loads(line)
        except json.JSONDecodeError:
            continue
        if item.get('type') == 'system' and item.get('subtype') == 'init':
            receipt.update(session_id=item.get('session_id'), model_actual=item.get('model'))
        if item.get('type') == 'result':
            final = item.get('result')
            receipt.update(result_subtype=item.get('subtype'), is_error=item.get('is_error'),
                           usage=item.get('usage'))
    if final:
        (out / 'final.md').write_text(final + '\n')
    receipt.update(exit_code=child.returncode, finished_utc=now(),
                   elapsed_seconds=round(time.monotonic()-start, 1))
    if receipt.get('is_error'):
        receipt['status'] = 'failed'
    save()
    print(json.dumps(receipt, indent=2))

if __name__ == '__main__':
    main()
