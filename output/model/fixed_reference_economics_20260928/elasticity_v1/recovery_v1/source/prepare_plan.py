#!/usr/bin/env python3
"""Create the immutable two-case recovery contract from authenticated v2 history."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

STAGE = Path('/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929')
OLD = Path('/scratch/td2248/projects/fixed_reference_elasticity_v2_20260929')
NAMES = ('grid_control', 'grid_control_repeat', 'credit', 'credit_repeat',
         'reference_990', 'credit_990')
PLAN_SHA = '6354bafd069fbfc4eda21271fc91747a22e3aeb8c90a57293f7dfde01badf28f'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    source = STAGE/'source'
    result = OLD/'results_v2'
    plan_path = STAGE/'plan_recovery.json'
    require(not plan_path.exists() and not (STAGE/'solve_v1').exists(), 'Recovery already prepared or started')
    require(sha(result/'plan_v2.json') == PLAN_SHA, 'Original v2 plan changed')
    require(sha(OLD/'source_v2/run_elasticity.py') ==
            '5d461d4999c428c77438ee43dd63db0e69ee0dd10356b0e09d65c62880adc5ca',
            'Original v2 driver changed')
    launch = read(result/'solve_v2/launch.json')
    require(str(launch['slurm_job']) == '18801318' and
            launch['deadline_epoch'] == 1790700222.6872504, 'Original job/deadline differs')
    status = subprocess.run(['sacct', '-j', '18803216', '-n', '-P', '--format=JobID,State,ExitCode'],
                            capture_output=True, text=True, check=True, timeout=20)
    exact = [line.split('|') for line in status.stdout.splitlines() if line.split('|')[0] == '18803216']
    require(len(exact) == 1 and exact[0][1] == 'FAILED', 'Prior continuation not terminal FAILED')
    latest_path = result/'solve_v4/latest_completed.json'
    latest = read(latest_path)
    require(latest['lifecycle_solves'] == 6 and
            [r['case'] for r in latest['completed']] == list(NAMES), 'Prior six-case order differs')
    failure_path = result/'solve_v4/failure.json'
    failure = read(failure_path)
    require(failure['status'] == 'failed' and 'reference_1010' in failure['error'] and
            not (result/'solve_v4/reference_1010/receipt.json').exists(),
            'Prior +1 reference has a certified result or other failure')
    old = read(result/'plan_v2.json')
    old.pop('maximum_lifecycle_solves', None)
    old.pop('stop_criteria', None)
    old['schema'] = 'block0506_three_price_recovery_v1'
    old['maximum_new_lifecycle_solves'] = 2
    old['case_seconds'] = 900
    old['total_seconds'] = 2100
    old['factors'] = [.99, 1., 1.01]
    old['driver_sha256'] = sha(source/'run_recovery.py')
    old['adapter_sha256'] = sha(source/'natural_credit.py')
    old['launcher_sha256'] = sha(source/'launch.sh')
    old['preparer_sha256'] = sha(source/'prepare_plan.py')
    old['original_plan_path'] = '/work/elasticity_results/plan_v2.json'
    old['original_plan_sha256'] = PLAN_SHA
    old['v2_driver_sha256'] = sha(OLD/'source_v2/run_elasticity.py')
    old['original_job'] = '18801318'
    old['original_deadline_epoch'] = launch['deadline_epoch']
    old['failed_job'] = '18803216'
    old['failed_job_state'] = 'FAILED'
    old['failed_case'] = 'reference_1010'
    old['failed_controller_receipt_path'] = '/work/elasticity_results/solve_v4/failure.json'
    old['failed_controller_receipt_sha256'] = sha(failure_path)
    old['failed_case_receipt_path'] = '/work/elasticity_results/solve_v4/reference_1010/receipt.json'
    old['prior_latest_path'] = '/work/elasticity_results/solve_v4/latest_completed.json'
    old['prior_latest_sha256'] = sha(latest_path)
    old['prior_attempts'] = 7
    old['prior_passes'] = 6
    old['automatic_retries'] = 0
    old['author_approval'] = '2026-09-29 explicit two-case recovery authorization'
    old['stop_criteria'] = ['Any model or scientific gate failure stops without retry',
        'Each new fresh-process case completes within 900 seconds including checkpoint and plots',
        'Two new solves and total 2100 seconds from this recovery job start are hard limits',
        'Only .99, 1.00 and 1.01 prices can form a completed comparison; no +/-2% claim']
    prior = []
    for record in latest['completed']:
        name = record['case']
        directory = ('solve_v2' if name in NAMES[:4] else 'solve_v4') + '/' + name
        receipt_path = result/directory/'receipt.json'
        receipt = read(receipt_path)
        require(sha(receipt_path) == record['receipt_sha256'] and
                receipt['status'] == 'passed' and receipt['case'] == name and
                receipt['plan_sha256'] == PLAN_SHA and receipt['standard_plot_count'] == 17,
                'Prior receipt differs: ' + name)
        prior.append(dict(case=name, factor=record['factor'], regime=record['regime'],
            receipt_path='/work/elasticity_results/'+directory+'/receipt.json',
            receipt_sha256=record['receipt_sha256'],
            checkpoint_path='/work/elasticity_results/'+directory+'/conditional_cohort_state.pkl.gz',
            checkpoint_sha256=receipt['checkpoint']['sha256'],
            fit_sha256=sha(receipt_path.parent/'target_fit.csv'),
            parameters_sha256=sha(receipt_path.parent/'parameters.csv'),
            standard_plot_sha256={f.name:sha(f) for f in
                sorted((receipt_path.parent/'standard_diagnostics').glob('*.png'))}))
    old['prior_cases'] = prior
    text = json.dumps(old, indent=2, sort_keys=True, allow_nan=False)+'\n'
    temporary = plan_path.with_suffix('.json.tmp')
    with temporary.open('w') as f:
        f.write(text); f.flush(); os.fsync(f.fileno())
    temporary.replace(plan_path)
    plan_path.chmod(0o444)
    print('Recovery plan', plan_path, sha(plan_path))


if __name__ == '__main__':
    main()
