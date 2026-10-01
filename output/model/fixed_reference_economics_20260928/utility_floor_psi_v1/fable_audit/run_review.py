"""Bounded read-only Fable review; adapted from the retained target-review supervisor."""
from pathlib import Path
import datetime
import json
import os
import signal
import subprocess
import time

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
OUT = Path(__file__).resolve().parent
CLI = '/Users/tommasodesanto/.local/bin/claude'

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)

def main():
    started = time.time()
    deadline = started + 1800
    record = {'status': 'authenticating', 'started_utc': now(),
              'started_epoch': started, 'deadline_epoch': deadline,
              'time_limit_seconds': 1800, 'supervisor_pid': os.getpid(),
              'model_requested': 'fable', 'effort': 'max', 'max_turns': 60,
              'allowed_tools': ['Read', 'Grep', 'Glob'],
              'prompt': 'prompt.md', 'retries': 0}
    save(OUT / 'launch.json', record)
    env = os.environ.copy()
    for key in ['ANTHROPIC_API_KEY', 'ANTHROPIC_AUTH_TOKEN', 'ANTHROPIC_BASE_URL']:
        env.pop(key, None)
    auth = subprocess.run([CLI, 'auth', 'status'], env=env, capture_output=True,
                          text=True, timeout=30, check=True)
    authentication = json.loads(auth.stdout)
    if not (authentication.get('loggedIn') and
            authentication.get('authMethod') == 'claude.ai' and
            authentication.get('subscriptionType') == 'max' and
            authentication.get('apiProvider') == 'firstParty'):
        record.update(status='blocked_auth', finished_utc=now())
        save(OUT / 'completion.json', record)
        return
    command = [CLI, '--model', 'fable', '--effort', 'max', '--print',
               '--output-format', 'stream-json', '--verbose', '--max-turns', '60',
               '--tools', 'Read,Grep,Glob', '--allowedTools', 'Read,Grep,Glob',
               '--permission-mode', 'dontAsk', '--strict-mcp-config',
               '--safe-mode', '--restricted', '--no-chrome']
    record.update(auth='claude.ai Max firstParty', command=command, status='starting')
    save(OUT / 'launch.json', record)
    assistant_text = []
    final = None
    parsed = 0
    with (OUT / 'stream.jsonl').open('w') as stream, (OUT / 'stderr.txt').open('w') as error:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdin=subprocess.PIPE,
                                   stdout=stream, stderr=error, text=True,
                                   start_new_session=True)
        record.update(pid=process.pid, status='running')
        save(OUT / 'launch.json', record)
        process.stdin.write((OUT / 'prompt.md').read_text())
        process.stdin.close()
        while True:
            lines = (OUT / 'stream.jsonl').read_text().splitlines()
            for line in lines[parsed:]:
                try:
                    item = json.loads(line)
                except json.JSONDecodeError:
                    break
                parsed += 1
                if item.get('type') == 'system' and item.get('subtype') == 'init':
                    record.update(session_id=item.get('session_id'),
                                  model_actual=item.get('model'),
                                  actual_tools=item.get('tools'), initialized_utc=now())
                    save(OUT / 'launch.json', record)
                if item.get('type') == 'assistant':
                    for block in item.get('message', {}).get('content', []):
                        if block.get('type') == 'text':
                            assistant_text.append(block.get('text', ''))
                    if assistant_text:
                        (OUT / 'partial.md').write_text('\n\n'.join(assistant_text) + '\n')
                if item.get('type') == 'result':
                    final = item.get('result')
                    record.update(result_subtype=item.get('subtype'), is_error=item.get('is_error'))
            record.update(last_checked_utc=now(), elapsed_seconds=round(time.time() - started, 1),
                          stream_records=parsed, stream_bytes=(OUT / 'stream.jsonl').stat().st_size)
            save(OUT / 'progress.json', record)
            if process.poll() is not None:
                record['status'] = 'completed' if process.returncode == 0 and not record.get('is_error') else 'failed'
                break
            if time.time() >= deadline:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                record['status'] = 'timed_out'
                break
            time.sleep(min(5, max(.1, deadline - time.time())))
    record.update(exit_code=process.returncode, finished_utc=now(),
                  elapsed_seconds=round(time.time() - started, 1))
    if final:
        (OUT / 'final.md').write_text(final + '\n')
    elif assistant_text:
        (OUT / 'final.md').write_text('Partial assessment; review stopped before a final result.\n\n' + '\n\n'.join(assistant_text) + '\n')
    save(OUT / 'completion.json', record)
    save(OUT / 'progress.json', record)

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        save(OUT / 'failure.json', {'status': 'supervisor_failure', 'error': repr(exc),
                                  'finished_utc': now(), 'supervisor_pid': os.getpid()})
        raise
