#!/usr/bin/env python3
"""Compact read-only monitoring of explicitly registered overnight lanes."""
import argparse
import json
import os
import signal
import subprocess
import time
from pathlib import Path


def read(path):
    return json.loads(Path(path).read_text())


def inspect(root):
    root = Path(root)
    result = {'path': str(root), 'exists': root.exists()}
    for name in ('complete', 'heartbeat', 'best_so_far', 'latest_completed', 'waiting_for_acceptance'):
        path = root / (name + '.json')
        if path.exists():
            data = read(path)
            if name == 'complete':
                result[name] = {k: data[k] for k in ('status', 'error_type', 'error', 'integrity_compromised') if k in data}
            elif name in ('best_so_far', 'latest_completed'):
                result[name] = {k: data[k] for k in ('case', 'status', 'loss', 'case_path', 'point', 'error') if k in data}
            else:
                result[name] = data
            result[name + '_age_seconds'] = time.time() - path.stat().st_mtime
    checkpoint = root / 'checkpoint.json'
    if checkpoint.exists():
        data = read(checkpoint); counts = {}
        for row in data.get('records', []):
            counts[row['status']] = counts.get(row['status'], 0) + 1
        result['counts'] = counts
    result['stale'] = result.get('heartbeat_age_seconds', 0) > 1800 and 'complete' not in result
    # Only requests in this registered controller, never an archive-wide scan.
    active = []
    for request in root.glob('*.request.json'):
        folder = root / request.name.removesuffix('.request.json')
        if (folder / 'success.json').exists() or (folder / 'failure.json').exists():
            continue
        ledger = folder / 'case/stationary_solves.json'
        if ledger.exists():
            rows = read(ledger)
            active.append({'case': folder.name, 'trials': len(rows),
                           'latest': rows[-1] if rows else None,
                           'ledger_age_seconds': time.time() - ledger.stat().st_mtime})
    result['active'] = active
    return result


def collect(registry, remote=True):
    run = read(registry)
    report = {'epoch': time.time(), 'status': run['status'],
              'absolute_end_epoch': run['absolute_end_epoch'], 'local': {}, 'remote': {}}
    for lane, record in run.get('local_lanes', {}).items():
        report['local'][lane] = {stage: inspect(record[stage]) for stage in ('smoke', 'search') if stage in record}
    ps = subprocess.run(['ps', '-Ao', 'pid,ppid,rss,%cpu,lstart,command'], capture_output=True, text=True, check=True)
    processes = []
    for line in ps.stdout.splitlines()[1:]:
        if (run.get('local_root','MISSING_ROOT') not in line or
                'run_e5f_utility_overnight_calibration.py' not in line or '--stage evaluate' not in line):
            continue
        parts = line.split(None, 9)
        if len(parts) == 10:
            processes.append({'pid': int(parts[0]), 'parent': int(parts[1]),
                              'rss_gib': int(parts[2]) / 1048576, 'cpu_percent': float(parts[3]),
                              'start_time': ' '.join(parts[4:9]), 'command': parts[9]})
    report['local_evaluators'] = processes
    report['total_local_evaluator_rss_gib'] = sum(p['rss_gib'] for p in processes)
    report['local_memory_review_required'] = report['total_local_evaluator_rss_gib'] > 32
    # Remote collection is opt-in per registry and never attempts credentials.
    if remote and run.get('remote_status_command'):
        command = ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', 'torch', run['remote_status_command']]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=40)
            report['remote'] = json.loads(result.stdout) if result.returncode == 0 else {'error': result.stderr.strip(), 'returncode': result.returncode}
        except (subprocess.TimeoutExpired, ValueError) as exc:
            report['remote'] = {'error': str(exc)}
    return report


def watch_memory(registry, output):
    """Terminate only an authenticated owned evaluator if local RSS exceeds32GiB.

    The strict controller records a fatal interruption and stops new search;
    this guard never converts a resource stop into a valid model observation.
    """
    registry=Path(registry);output=Path(output);high=0
    cutoff=read(registry)['absolute_end_epoch']
    while time.time()<cutoff:
        state=collect(registry,remote=False)
        state.pop('remote',None)
        temp=output.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(output)
        high=high+1 if state['total_local_evaluator_rss_gib']>32 else 0
        if high>=2 and state['local_evaluators']:
            process=max(state['local_evaluators'],key=lambda p:p['rss_gib'])
            pid=process['pid'];command=process['command']
            # The command contains an absolute owned path with no spaces; the
            # controller's authenticated startup additionally binds PID/parent.
            if ' --output ' not in command:raise RuntimeError('Missing owned output path')
            folder=Path(command.split(' --output ',1)[1].split(' --request ',1)[0])
            run=read(registry);assert folder.resolve().is_relative_to(Path(run['local_root']).resolve())
            startup=read(folder/'startup.json')
            if startup['pid']!=pid or startup['parent_pid']!=process['parent']:
                high=0;continue
            # Reject exited/reused PIDs using fresh identity and process birth
            # time immediately before the signal; artifacts alone are not liveness.
            live=subprocess.run(['ps','-p',str(pid),'-o','ppid=,lstart=,command='],
                                capture_output=True,text=True).stdout.strip().split(None,6)
            if (len(live)!=7 or int(live[0])!=process['parent'] or
                    ' '.join(live[1:6])!=process['start_time'] or live[6]!=command):
                high=0;continue
            try:
                if os.getpgid(pid)!=pid:
                    high=0;continue
                os.killpg(pid,signal.SIGTERM)
                outcome='signal_sent_stop_requested'
            except ProcessLookupError:
                outcome='process_already_exited_no_signal'
            event={'epoch':time.time(),'event':'local_memory_resource_intervention',
                   'outcome':outcome,
                   'total_rss_gib':state['total_local_evaluator_rss_gib'],
                   'threshold_gib':32,'process':process,'case_path':str(folder),
                   'classification':'external_resource_interruption_not_model_inadmissibility'}
            with (output.parent/'resource_interventions.jsonl').open('a') as f:
                f.write(json.dumps(event)+'\n');f.flush()
            high=0
        time.sleep(20)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--registry', type=Path)
    p.add_argument('--lane', type=Path, action='append')
    p.add_argument('--output', type=Path)
    p.add_argument('--watch-local-memory', action='store_true')
    a = p.parse_args()
    if a.watch_local_memory:
        assert a.registry and a.output and not a.lane
        watch_memory(a.registry,a.output)
        raise SystemExit(0)
    result = {str(path): inspect(path) for path in a.lane} if a.lane else collect(a.registry)
    text = json.dumps(result, indent=2, sort_keys=True)
    if a.output:
        temp = a.output.with_suffix('.tmp'); temp.write_text(text+'\n'); temp.replace(a.output)
    print(text)
