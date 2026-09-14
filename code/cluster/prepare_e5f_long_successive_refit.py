"""Freeze the successive-surprise contract and prepare its bounded Torch chain.

Run on the login node after staging the driver and numerical helper sources.
Preparation performs no model solve; --submit dispatches the exact-loop smoke.
The compute driver dispatches later stages only after their prerequisites pass.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

BASE = Path('/scratch/td2248/projects/Fertility_Spring26_candidate_path_20260911a/batches')
PYTHON = '/share/apps/anaconda3/2025.06/bin/python'
TARGET_SHA = '8945aa2427e26157e2e01b44b327e374daae5078daa458727f025938569a2f74'
YEARS = [2007, 2011, 2015, 2019]
TARGETS = [1.974875, 1.861, 1.755375, 1.64575]
SEEDS = [0.12891531457859182, 0.11696608375682901,
         0.10564290922456478, 0.09221854783921073]


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


def prepare(batch):
    batch = Path(batch).resolve()
    manifest = batch / 'manifest.json'
    if manifest.exists():
        raise ValueError('Refuse to replace an existing experiment manifest')
    source = batch / 'source'
    driver = source / 'run_e5f_long_successive_refit.py'
    numerics = source
    helpers = [source / 'e5f_ssj_scaled_step_root.py',
               source / 'e5f_ssj_toeplitz_jacobian.py']
    spec = BASE / 'afternoon_original_queue_20260913a/spec.json'
    jacobian = BASE / 'afternoon_original_queue_20260913a_ssj_toeplitz_10/derivative/derivative_receipt.json'
    targets = batch / 'inputs/empirical_blocks.csv'
    if sha(targets) != TARGET_SHA:
        raise ValueError('Empirical target fingerprint differs')
    rows = list(csv.DictReader(targets.open()))
    if [(int(r['decision_year']), float(r['period_tfr_arithmetic_mean'])) for r in rows] != list(zip(YEARS, TARGETS)):
        raise ValueError('Four target rows differ')
    for year, row in zip(YEARS, rows):
        if (int(row['birth_year_start']), int(row['birth_year_end'])) != (year + 1, year + 4):
            raise ValueError('Fertility observation clock differs')
    native = read(spec)
    base = read(native['base_spec'])
    runtime = read(base['runtime_manifest'])
    prior = read(runtime['prior_plan'])
    initial_score = Path(prior['initial_score_path'])
    frozen_refs = [driver, *helpers, spec, Path(native['base_spec']),
                   Path(base['runtime_manifest']), Path(runtime['prior_plan']),
                   initial_score, targets, jacobian,
                   source / Path(__file__).name]
    staged_receipt = batch / 'inputs/source_manifest.json'
    if staged_receipt.exists():
        frozen_refs.append(staged_receipt)
    script_paths = {key: str(batch / 'scripts' / (key + '.sbatch'))
                    for key in ['smoke', 'stage_0', 'stage_1', 'stage_2', 'stage_3', 'policy']}
    m = dict(schema='e5f_long_successive_refit_v1', output=str(batch / 'output'),
        spec=str(spec), source_driver=str(driver), numerics_dir=str(numerics),
        jacobian_source=str(jacobian), empirical_blocks=str(targets),
        target_sha256=TARGET_SHA, initial_score=str(initial_score),
        initial_target_fingerprint=prior['target_fingerprint'],
        shock_years=YEARS, year_targets=dict(zip(map(str, YEARS), TARGETS)),
        psi_seeds=SEEDS, bounds_relative=[-.20, .02], fit_tolerance=.005,
        max_trials=6, max_rounds=4, max_root_evaluations=8,
        candidate_seconds=36000, stage_seconds=86400, policy_seconds=43200,
        smoke_seconds=1800, terminal_seconds=3600, terminal_max_evaluations=24,
        total_deadline_unix=time.time() + 7 * 86400, terminal_year=2423,
        scripts=script_paths,
        file_sha256={str(p): sha(p) for p in frozen_refs},
        contract=dict(information='Successive permanent surprises; later surprises not anticipated',
            structural_parameters='estimated, fixed retained 2007 calibration',
            preference_levels='four estimated dated levels',
            initial_distribution='retained calibrated 2007 original-queue checkpoint',
            later_distributions='endogenous inherited states from accepted preceding shock',
            entry_conversion='externally fixed births/2.1',
            migration='externally fixed zero', population_rescaling=False,
            historical_age_reweighting=False,
            supply_scale='empirically normalized fixed asset-price supply curve',
            supply_elasticity='externally fixed 0.63',
            baseline_tax='externally fixed 1% annual; all receipts equally rebated',
            policy_tax='externally fixed 2% annual; unexpected in inherited 2023',
            payroll_tax='externally fixed 0.179; pension balances PAYGO',
            fertility_measurement='retained household-rate analogue of published female TFR',
            terminal_distance='outstanding until checked at each fitted endpoint',
            horizon_robustness='outstanding; no acceptance inferred from finite root',
            production_eligible=False),
        budget=dict(cpus_per_job=1, memory_gb=32, first_stage_dates=104,
            later_stage_dates=[103, 102, 101], policy_dates=100,
            max_candidate_values_total=24, max_native_root_mappings_per_candidate=32,
            root_mapping_seconds_observed=[1666, 2223],
            stage_wall_hours=24, policy_wall_hours=12, total_calendar_days=7,
            interpretation='Caps, not runtime predictions; stop with best evidence if exhausted'))
    for key, name in script_paths.items():
        mode = 'stage' if key.startswith('stage_') else key
        argv = [PYTHON, str(driver), '--manifest', str(manifest), '--mode', mode]
        if mode == 'stage':
            argv.extend(['--stage', key.split('_')[1]])
        minutes = {'smoke': 35, 'stage': 1450, 'policy': 730}[mode]
        p = Path(name)
        p.parent.mkdir(parents=True, exist_ok=True)
        lines = ['#!/bin/bash', '#SBATCH --account=torch_pr_570_general',
                 '#SBATCH --cpus-per-task=1', '#SBATCH --mem=32G',
                 f'#SBATCH --time={minutes}', f'#SBATCH --job-name=e5f_refit_{key}',
                 f'#SBATCH --output={batch}/slurm_{key}_%j.log', 'set -euo pipefail',
                 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1',
                 shlex.join(argv), '']
        p.write_text('\n'.join(lines))
        m['file_sha256'][str(p)] = sha(p)
    save(manifest, m)
    save(batch / 'prepared.json', dict(manifest_sha256=sha(manifest), submitted=False,
        stage_count=4, candidate_cap=24, smoke_script=script_paths['smoke']))
    return manifest


def submit(manifest):
    m = read(manifest)
    receipt = Path(manifest).parent / 'submission.json'
    if receipt.exists():
        raise ValueError('Smoke submission receipt already exists')
    for path, expected in m['file_sha256'].items():
        if sha(path) != expected:
            raise ValueError('Pinned input changed: ' + path)
    reply = subprocess.run(['sbatch', '--parsable', m['scripts']['smoke']],
                           text=True, capture_output=True, check=True)
    result = dict(job_id=int(reply.stdout.strip().split(';')[0]), mode='smoke',
                  manifest_sha256=sha(manifest), submitted_unix=time.time())
    save(receipt, result)
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--submit', action='store_true')
    args = parser.parse_args()
    manifest = prepare(args.batch)
    if args.submit:
        submit(manifest)
    else:
        print(json.dumps(read(args.batch / 'prepared.json')))


if __name__ == '__main__':
    main()
