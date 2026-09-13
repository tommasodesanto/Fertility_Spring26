"""Freeze and smoke-gate one announced four-step, original-queue transition.

Preparation/dispatch use only the standard library on the Torch login node.
All scientific evaluations run on compute nodes. Existing batches are immutable.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shlex
import subprocess
import time

BASE = Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches')
ORIGINAL = BASE / 'afternoon_original_queue_20260913a'
ENDPOINT = BASE / 'afternoon_original_queue_20260913a_terminal_restart_v1/endpoint'
PYTHON = '/share/apps/anaconda3/2025.06/bin/python'


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def submit(manifest_path, mode):
    m = read(manifest_path)
    batch = Path(m['output']).parent
    argv = [PYTHON, str(batch / 'source/run_e5f_announced_original_queue.py'),
            '--manifest', str(manifest_path), '--mode', mode]
    lines = ['#!/bin/bash', f'#SBATCH --job-name=e5f_announced_{mode}',
             '#SBATCH --account=torch_pr_570_general', '#SBATCH --cpus-per-task=1',
             '#SBATCH --mem=32G', f'#SBATCH --time={25 if mode == "smoke" else 370}',
             f'#SBATCH --output={batch}/logs/{mode}_%j.log', 'set -euo pipefail',
             'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1',
             shlex.join(argv)]
    if mode == 'smoke':
        lines.append(shlex.join(['/usr/bin/python3', str(Path(__file__).resolve()),
                                 '--dispatch', str(manifest_path)]))
    script = batch / f'{mode}.sbatch'
    if script.exists():
        raise ValueError('Refusing duplicate submission: ' + str(script))
    script.write_text('\n'.join(lines) + '\n')
    result = subprocess.run(['sbatch', '--parsable', str(script)],
                            text=True, capture_output=True, check=True)
    job = result.stdout.strip().split(';')[0]
    if not job.isdigit():
        raise RuntimeError('Unexpected Slurm response: ' + result.stdout)
    return dict(job_id=int(job), mode=mode, script=str(script), script_sha256=sha(script),
                manifest_sha256=sha(manifest_path), submitted_unix=time.time())


def prepare(batch):
    batch = batch.resolve()
    manifest_path = batch / 'manifest.json'
    if manifest_path.exists():
        raise ValueError('Refusing to overwrite a prepared experiment')
    source = batch / 'source'
    spec_path = ORIGINAL / 'spec.json'
    spec = read(spec_path)
    fits_path = Path(spec['permanent_shock_source'])
    fits = read(fits_path)
    if sha(fits_path) != spec['permanent_shock_source_sha256']:
        raise ValueError('Fitted shock source fingerprint changed')
    if [r['year'] for r in fits] != [2007, 2011, 2015, 2019]:
        raise ValueError('Expected the approved four fitted shocks')
    levels = [float(r['psi']) for r in fits]
    if levels[-1] != spec['permanent_psi']:
        raise ValueError('Final shock differs from the verified endpoint exercise')
    endpoint_pickle = ENDPOINT / 'terminal.pkl.gz'
    endpoint_receipt = ENDPOINT / 'root_receipt.json'
    if not read(endpoint_receipt)['verified']:
        raise ValueError('Endpoint receipt is not verified')
    seed_path = source / 'initial_seed.json'
    seed = read(seed_path)
    if (not seed.get('mapping_valid') or seed.get('evaluation') != 5
            or len(seed.get('prices', [])) != 300):
        raise ValueError('Require the frozen fifth 100-period numerical seed')
    files = [spec_path, fits_path, endpoint_pickle, endpoint_receipt, seed_path,
             source / 'run_e5f_announced_original_queue.py',
             source / 'build_e5f_stationary_shock_figures.py', Path(__file__).resolve()]
    if not all(p.is_file() for p in files):
        raise ValueError('Missing frozen source or input')
    m = dict(spec=str(spec_path), endpoint_pickle=str(endpoint_pickle),
             endpoint_receipt=str(endpoint_receipt), output=str(batch / 'output'),
             psi_levels=levels, shock_years=[2007, 2011, 2015, 2019],
             initial_seed_json=str(seed_path), periods_after_four=100,
             absolute_deadline_unix=time.time() + 7*3600, seconds=6*3600,
             smoke_seconds=1200, max_root_evaluations=8,
             plotter=str(source / 'build_e5f_stationary_shock_figures.py'),
             file_sha256={str(p): sha(p) for p in files},
             information='All four preference levels known at time zero; final level permanent',
             scientific_contract=dict(population_law='original_household_birth_vintage_queue',
                 birth_to_entry_conversion=1/2.1, no_immigration=True,
                 historical_age_conditioning=False, equal_property_tax_rebate=True,
                 annual_property_tax=.01, payroll_tax=.179,
                 structural_parameters_reestimated=False, shocks_reestimated=False,
                 terminal_reconstructed_before_dispatch=True, production_eligible=False),
             budget=dict(total_dates=104, expected_seconds_per_full_mapping=1955,
                 maximum_long_root_evaluations=8, maximum_long_root_rounds=1,
                 expected_full_stage_hours=[3, 5], long_stage_wall_hours=6,
                 smoke_wall_minutes=20, absolute_wall_hours=7,
                 stop='Any failed smoke gate prevents dispatch; full root stops at eight mappings or deadline'))
    (batch / 'logs').mkdir(parents=True, exist_ok=True)
    save(manifest_path, m)
    receipt = submit(manifest_path, 'smoke')
    save(batch / 'submission.json', receipt)
    print(json.dumps(receipt))


def dispatch(manifest_path):
    m = read(manifest_path)
    for path, digest in m['file_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Pinned input changed: ' + path)
    batch = Path(m['output']).parent
    if (batch / 'dispatch.json').exists():
        raise ValueError('Long job already dispatched')
    smoke = read(Path(m['output']) / 'smoke_summary.json')
    if not smoke.get('passed') or smoke.get('manifest_sha256') != sha(manifest_path):
        raise ValueError('Exact-loop and endpoint smoke has not passed')
    if m['absolute_deadline_unix'] - time.time() < 3600:
        raise TimeoutError('Insufficient remaining experiment budget')
    receipt = submit(manifest_path, 'run')
    save(batch / 'dispatch.json', receipt)
    print(json.dumps(receipt))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', type=Path)
    group.add_argument('--dispatch', type=Path)
    args = parser.parse_args()
    prepare(args.prepare) if args.prepare else dispatch(args.dispatch)
