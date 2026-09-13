"""Resume one failed original-queue endpoint arm without changing its numerics.

Preparation uses only the standard library on the login node. The compute-node
run imports the original frozen runner and solver, verifies its first mapping
against the saved best point, and retains the original absolute stage deadlines.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import threading
import time


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def prepare(spec_path, output, prior_job):
    spec = read(spec_path)
    dispatch = read(Path(spec['batch']) / 'dispatch.json')
    jobs = [j for j in dispatch['jobs'] if j['job_id'] == prior_job]
    if len(jobs) != 1 or jobs[0]['mode'] != 'terminal' or jobs[0]['count'] != 100:
        raise ValueError('Expected exactly the original terminal-100 job')
    env = dict(os.environ, TZ='UTC')
    result = subprocess.run(['sacct', '-j', str(prior_job), '-X', '-n', '-P',
                             '-o', 'JobID,State,Start'], env=env,
                            capture_output=True, text=True, check=True)
    rows = [r.split('|') for r in result.stdout.splitlines() if r.strip()]
    if len(rows) != 1 or rows[0][0] != str(prior_job) or rows[0][1] != 'FAILED':
        raise ValueError('Original terminal job must have failed and stopped')
    start = datetime.datetime.fromisoformat(rows[0][2]).replace(
        tzinfo=datetime.timezone.utc).timestamp()
    if time.time() >= start + 3600 - 120:
        raise TimeoutError('Original one-hour endpoint window nearly exhausted')
    original = Path(jobs[0]['output'])
    receipt = original / 'endpoint/root_receipt.json'
    previous_contract = original / 'experiment_contract.json'
    contract_path = output / 'restart_contract.json'
    if contract_path.exists():
        raise ValueError('Refusing to duplicate a prepared restart')
    output.mkdir(parents=True, exist_ok=True)
    contract = dict(spec=str(spec_path), output=str(output), prior_job=prior_job,
                    prior_receipt=str(receipt), prior_contract=str(previous_contract),
                    prior_start_unix=start, terminal_deadline_unix=start+3600,
                    run_deadline_unix=min(start+jobs[0]['seconds'],
                                          spec['absolute_deadline_unix']),
                    file_sha256={str(p): sha(p) for p in
                        [Path(__file__).resolve(), spec_path, receipt, previous_contract]},
                    changed='Numerical starting coordinates and saved Jacobian only',
                    production_eligible=False)
    save(contract_path, contract)
    python = '/share/apps/anaconda3/2025.06/bin/python'
    argv = [python, str(Path(__file__).resolve()), '--run', str(contract_path)]
    script = output / 'restart.sbatch'
    script.write_text('\n'.join([
        '#!/bin/bash', '#SBATCH --job-name=e5f_orig_terminal_resume',
        '#SBATCH --account=torch_pr_570_general', '#SBATCH --cpus-per-task=1',
        '#SBATCH --mem=32G', '#SBATCH --time=360',
        f'#SBATCH --output={output}/restart_%j.log', 'set -euo pipefail',
        'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1',
        shlex.join(argv), '']))
    result = subprocess.run(['sbatch', '--parsable', str(script)],
                            capture_output=True, text=True, check=True)
    job = int(result.stdout.strip().split(';')[0])
    save(output / 'submission.json', dict(job_id=job, replaces=prior_job,
         restart_contract_sha256=sha(contract_path), script_sha256=sha(script)))
    print(json.dumps(dict(job_id=job, output=str(output),
                         terminal_deadline_unix=contract['terminal_deadline_unix'])))


def run(contract_path):
    contract = read(contract_path)
    for path, digest in contract['file_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Changed restart input: ' + path)
    spec_path = Path(contract['spec'])
    spec = read(spec_path)
    sys.path.insert(0, str(Path(spec['batch']) / 'source'))
    import run_e5f_original_queue_experiments as runner
    import numpy as np
    from unittest.mock import patch
    c = runner.load_context(spec_path)
    c.spec_path = spec_path
    smoke = read(spec['smoke_summary'])
    if smoke.get('status') != 'passed' or smoke['spec_sha256'] != sha(spec_path):
        raise ValueError('Original exact-loop smoke does not match')
    previous = read(contract['prior_contract'])
    if (previous['spec_sha256'] != sha(spec_path) or previous['mode'] != 'terminal'
            or previous['psi_permanent'] != spec['permanent_psi']):
        raise ValueError('Saved numerical seed has a different experiment contract')
    receipt = read(contract['prior_receipt'])
    best = receipt['best']
    if receipt['verified'] or not best['mapping_valid']:
        raise ValueError('Expected an unverified but admissible saved terminal point')
    controls = dict(c.controls, initial_jacobian=receipt['final_jacobian'],
                    damping=receipt['final_damping'])
    output = Path(contract['output'])
    end = min(contract['run_deadline_unix'], spec['absolute_deadline_unix'])
    deadline = time.monotonic() + end - time.time()
    terminal_deadline = time.monotonic() + contract['terminal_deadline_unix'] - time.time()
    if min(deadline, terminal_deadline) <= time.monotonic() + 120:
        raise TimeoutError('Original deadline expired before restart')
    stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            save(output/'controller_heartbeat.json', dict(remaining_seconds=end-time.time()))
            if time.time() >= end:
                save(output/'controller_failure.json', dict(error='Original hard deadline'))
                os._exit(124)
    threading.Thread(target=heartbeat, daemon=True).start()
    import e5f_original_queue_terminal as terminal
    original_evaluate = terminal._evaluate_trial
    first = True
    def checked_evaluate(**kwargs):
        nonlocal first
        candidate = original_evaluate(**kwargs)
        if first:
            first = False
            error = float(np.max(np.abs(np.asarray(candidate.residual)-best['residual'])))
            # The retained log-price root may round-trip the seed by one ulp.
            if not np.allclose(kwargs['coordinates'], best['prices'], rtol=2e-14, atol=0) or error > 2e-10:
                raise ValueError('Saved best terminal mapping did not reproduce exactly')
            save(output/'seed_reproduction.json', dict(passed=True, residual_max_abs=error))
        return candidate
    try:
        with c.queue.original_queue_adapter(), runner.original_receipts(c), \
                runner.capture_root_paths(c), c.cache.policy_cache(c.joined.pf, max_bytes=12*1024**3):
            with patch.object(terminal, '_evaluate_trial', checked_evaluate):
                endpoint = terminal.solve_terminal(old=c.old, psi=spec['permanent_psi'],
                    audit=c.audit, controls=controls, deadline=min(deadline, terminal_deadline),
                    folder=output/'endpoint', start=best['prices'])
            if not endpoint.verified:
                raise RuntimeError('Restart did not verify the stationary endpoint')
            with runner.capture_paths(c, output/'transition'):
                runner.fixed_terminal_path(c, endpoint, output/'transition', 100,
                                           spec['permanent_psi'], deadline)
        c.driver.verify_pins(spec['file_sha256'])
        save(output/'controller_complete.json', dict(completed=True, production_eligible=False))
    except BaseException as exc:
        save(output/'controller_failure.json', dict(error_type=type(exc).__name__, error=str(exc)))
        raise
    finally:
        stop.set()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', type=Path)
    p.add_argument('--spec', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--prior-job', type=int)
    args = p.parse_args()
    if args.run:
        run(args.run)
    elif args.spec and args.output and args.prior_job:
        prepare(args.spec, args.output, args.prior_job)
    else:
        p.error('Provide --run, or --spec, --output and --prior-job')


if __name__ == '__main__':
    main()
