"""Validate and summarize copied Torch chain folders; reject mixed contracts.

Pass --fetch to copy the eight selected chain folders from Torch with rsync first.
The collector never adopts an estimate or deletes a remote result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

REMOTE = '/scratch/td2248/projects/soft_timing_calibration_20261002_v2/results'
VALID_MODES = ('mock', 'smoke', 'production')


def tasks(mode: str) -> list[tuple[str, int]]:
    return [('original', 0), ('alternative', 0)] if mode == 'smoke' else [
        (arm, chain) for arm in ('original', 'alternative') for chain in range(4)]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def collect_one(root: Path, mode: str, arm: str, chain: int) -> dict:
    folder = root / f'{mode}_{arm}_chain_{chain}'
    start = read(folder / 'launcher_start.json')
    terminal = read(folder / 'launcher_terminal.json')
    contract = read(folder / 'run/start_contract.json')
    completed = read(folder / 'run/completed.json')
    assert start['arm'] == terminal['arm'] == contract['arm'] == arm
    assert int(start['chain']) == int(terminal['chain']) == int(contract['chain']) == chain
    assert start['mode'] == terminal['mode'] == mode
    assert terminal['exit_code'] == 0, folder
    assert contract['objective_calls_max'] == 250 and contract['reserve_seconds'] == 1800
    if mode == 'mock':
        assert completed['status'] == 'mock_loop_passed_zero_solves'
        assert completed['objective_calls'] == 2 and completed['lifecycle_solves'] == 0
    else:
        latest = read(folder / 'run/latest_completed.json')
        best = read(folder / 'run/best_so_far.json')
        search = read(folder / 'run/search_completed.json')
        inputs = read(folder / 'run/input_contract.json')
        assert latest['completed_full_ge'] >= 1 and best['best'] is not None, folder
        if mode == 'smoke':
            assert search['objective_calls'] == 2, folder
        assert inputs['target_fingerprint'] == contract['target_fingerprint'], folder
        assert inputs['weight_fingerprint'] == contract['weight_fingerprint'], folder
        assert completed['status'] == 'selected_numerically_verified', folder
        assert completed['selected_postcheck']['status'] == 'passed', folder
        assert completed['target_fingerprint'] == contract['target_fingerprint']
        assert completed['weight_fingerprint'] == contract['weight_fingerprint']
        report = Path(completed['selected_postcheck']['report'])
        # A copied report may retain the remote absolute path. Resolve by its path
        # below the chain's run directory, without accepting an unrelated folder.
        parts = report.parts
        try:
            report = folder / 'run' / Path(*parts[parts.index('run') + 1:])
        except ValueError:
            raise AssertionError(f'Unexpected report path: {report}')
        with (report / 'target_fit.csv').open(newline='') as stream:
            fits = list(csv.DictReader(stream))
        with (report / 'parameters.csv').open(newline='') as stream:
            params = list(csv.DictReader(stream))
        assert len(fits) == 14 and len(params) == 31, folder
        assert len(list((report / 'standard_diagnostics').glob('*.png'))) == 17, folder
        actual_loss = sum(float(row['loss_contribution'] or 0) for row in fits)
        assert abs(actual_loss - completed['native_loss']) < 1e-8, folder
    return {'arm': arm, 'chain': chain, 'mode': mode, 'status': completed['status'],
            'target_fingerprint': contract['target_fingerprint'],
            'weight_fingerprint': contract['weight_fingerprint'],
            'selected_source_sha256': contract['selected_source_sha256'],
            'native_loss': completed.get('native_loss'),
            'objective_calls': completed.get('objective_calls'),
            'launcher_terminal_sha256': sha(folder / 'launcher_terminal.json'),
            'completed_sha256': sha(folder / 'run/completed.json')}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=VALID_MODES, required=True)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--fetch', action='store_true')
    args = ap.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    if args.fetch:
        for arm, chain in tasks(args.mode):
            name = f'{args.mode}_{arm}_chain_{chain}'
            subprocess.run(['rsync', '-a', '--', f'torch:{REMOTE}/{name}/', str(args.root / name) + '/'], check=True)
    rows = [collect_one(args.root, args.mode, arm, chain) for arm, chain in tasks(args.mode)]
    for key in ('target_fingerprint', 'weight_fingerprint', 'selected_source_sha256'):
        if len({row[key] for row in rows}) != 1:
            raise SystemExit(f'Mixed {key} across chains')
    receipt = {'status': f'all_{len(rows)}_verified', 'mode': args.mode, 'chains': rows,
               'target_fingerprint': rows[0]['target_fingerprint'],
               'weight_fingerprint': rows[0]['weight_fingerprint'],
               'selected_source_sha256': rows[0]['selected_source_sha256']}
    (args.root / f'{args.mode}_collection.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': receipt['status'], 'mode': args.mode,
                      'chains': len(rows), 'collection': str(args.root / f'{args.mode}_collection.json')}))


if __name__ == '__main__':
    main()
