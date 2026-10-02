"""Collect authenticated postchecked hard/quarter candidates without model solves."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
OUT = HERE / 'readout'
REMOTE = '/scratch/td2248/projects/purchase_rules_overnight_v1'
RESTART = '/scratch/td2248/projects/purchase_restart_controller_v2'


def run(*args, input=None, env=None):
    return subprocess.run(args, input=input, env=env, text=True, check=True, capture_output=True).stdout


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select_arm(chains, arm):
    valid = [row for row in chains if row['arm'] == arm and row['status'] == 'postchecked']
    return min(valid, key=lambda row: (row['loss'], row['chain'], row['source_run'])) if valid else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch-reports', action='store_true', help='Copy selected 17-plot native reports to local readout')
    parser.add_argument('--fetch-arrays', action='store_true', help='Also copy large selected native NPZ files')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    source = (HERE / 'scan_remote.py').read_text()
    snapshot = json.loads(run('ssh', '-o', 'BatchMode=yes', 'torch',
                              '/share/apps/anaconda3/2025.06/bin/python -', input=source))
    for row in snapshot['chains']:
        row['origin'] = 'torch'
    restarted = json.loads(run('ssh', '-o', 'BatchMode=yes', 'torch',
        f'env PURCHASE_SOURCE_KIND=restart PURCHASE_RESULTS_ROOT={RESTART}/results '
        f'PURCHASE_PARENT_RESULTS_ROOT={REMOTE}/results '
        '/share/apps/anaconda3/2025.06/bin/python -', input=source))
    if (restarted['target_fingerprint'] != snapshot['target_fingerprint'] or
            restarted['weight_fingerprint'] != snapshot['weight_fingerprint']):
        raise RuntimeError('Restart target or weight fingerprint drift')
    for row in restarted['chains']:
        row['origin'] = 'torch_restart'
    snapshot['chains'].extend(restarted['chains'])
    snapshot['errors'].extend(restarted['errors'])
    local_root = PACKET / 'local_runtime/runs/local10_v1'
    if local_root.is_dir():
        env = dict(os.environ, PURCHASE_RESULTS_ROOT=str(local_root), PURCHASE_PACKET_ROOT=str(PACKET),
                   PURCHASE_CHAIN_FIRST='48', PURCHASE_CHAIN_LAST='58')
        local = json.loads(run(sys.executable, str(HERE / 'scan_remote.py'), env=env))
        if local['target_fingerprint'] != snapshot['target_fingerprint'] or local['weight_fingerprint'] != snapshot['weight_fingerprint']:
            raise RuntimeError('Local/Torch target or weight fingerprint drift')
        for row in local['chains']:
            row['origin'] = 'local'
        snapshot['chains'].extend(local['chains'])
        snapshot['errors'].extend(local['errors'])
        local_restart_root = PACKET / 'local_runtime/restart_v2/runs'
        if local_restart_root.is_dir():
            env.update(PURCHASE_RESULTS_ROOT=str(local_restart_root), PURCHASE_SOURCE_KIND='local_restart',
                       PURCHASE_PARENT_RESULTS_ROOT=str(local_root))
            local_restart = json.loads(run(sys.executable, str(HERE / 'scan_remote.py'), env=env))
            if (local_restart['target_fingerprint'] != snapshot['target_fingerprint'] or
                    local_restart['weight_fingerprint'] != snapshot['weight_fingerprint']):
                raise RuntimeError('Local restart target or weight fingerprint drift')
            for row in local_restart['chains']:
                row['origin'] = 'local_restart'
            snapshot['chains'].extend(local_restart['chains'])
            snapshot['errors'].extend(local_restart['errors'])
    plan = json.loads((PACKET / 'plan.json').read_text())
    if snapshot['target_fingerprint'] != plan['target_fingerprint'] or snapshot['weight_fingerprint'] != plan['weight_fingerprint']:
        raise RuntimeError('Mixed target or weight fingerprint')
    (OUT / 'snapshot.json').write_text(json.dumps(snapshot, indent=2, sort_keys=True) + '\n')
    if snapshot['errors']:
        raise RuntimeError('Invalid postchecked chains; see readout/snapshot.json')
    winners = {}
    for arm in ('hard', 'quarter'):
        best = select_arm(snapshot['chains'], arm)
        if best is None:
            continue
        winners[arm] = best
        (OUT / f'selected_{arm}.json').write_text(json.dumps(best, indent=2, sort_keys=True) + '\n')
        if args.fetch_reports:
            destination = OUT / arm / 'selected_root'
            destination.mkdir(parents=True, exist_ok=True)
            source_report = (f'torch:{best["remote_report"]}/' if best['origin'].startswith('torch')
                             else best['remote_report'] + '/')
            run('rsync', '-az', source_report, str(destination) + '/')
            native_root = Path(best['remote_report']).parent
            other = OUT / arm
            if best['origin'].startswith('torch'):
                chain_root = f'torch:{best["remote_root"]}/postcheck/'
                run('rsync', '-az', chain_root + 'completed.json', str(other / 'completed.json'))
                run('rsync', '-az', chain_root + 'input_contract.json', str(other / 'input_contract.json'))
                for repeat in ('selected_repeat', 'selected_repeat_final'):
                    dest = other / repeat
                    dest.mkdir(parents=True, exist_ok=True)
                    run('rsync', '-az', '--exclude=*.npz', f'torch:{native_root}/{repeat}/', str(dest) + '/')
            else:
                shutil.copy2(Path(best['remote_root']) / 'postcheck/completed.json', other / 'completed.json')
                shutil.copy2(Path(best['remote_root']) / 'postcheck/input_contract.json', other / 'input_contract.json')
                for repeat in ('selected_repeat', 'selected_repeat_final'):
                    dest = other / repeat
                    dest.mkdir(parents=True, exist_ok=True)
                    run('rsync', '-az', '--exclude=*.npz', str(native_root / repeat) + '/', str(dest) + '/')
            for name, expected in best['report_sha256'].items():
                if sha(destination / name) != expected:
                    raise RuntimeError(f'Collected {arm} report hash mismatch: {name}')
        if args.fetch_arrays:
            destination = OUT / arm / 'selected_repeat' / 'stage'
            destination.mkdir(parents=True, exist_ok=True)
            source_arrays = (f'torch:{best["remote_arrays"]}' if best['origin'].startswith('torch')
                             else best['remote_arrays'])
            run('rsync', '-az', source_arrays, str(destination / 'solution_arrays.npz'))
            if (destination / 'solution_arrays.npz').stat().st_size != best['native_arrays_bytes']:
                raise RuntimeError(f'Collected {arm} array size mismatch')
    fit = []
    params = []
    for arm, best in winners.items():
        for row in best['target_fit']:
            fit.append(dict(rule=arm, chain=best['chain'], **row))
        for row in best['parameters']:
            params.append(dict(rule=arm, chain=best['chain'], **row))
    if fit:
        write_csv(OUT / 'target_fit.csv', fit)
        write_csv(OUT / 'parameters.csv', params)
    status = dict(status='postchecked_selection_snapshot', submitted_original_torch=48,
        scanned_restart_torch=48, postchecked_by_arm={
        arm: sum(row['arm'] == arm and row['status'] == 'postchecked' for row in snapshot['chains'])
        for arm in ('hard', 'quarter')}, selected={arm: dict(chain=row['chain'], loss=row['loss'],
        price=row['price'], H0=row['H0'], remote_arrays=row['remote_arrays'],
        source_run=row['source_run'], origin=row['origin'], remote_root=row['remote_root'],
        parent_remote_root=row.get('parent_remote_root')) for arm, row in winners.items()},
        reports_fetched=args.fetch_reports, arrays_fetched=args.fetch_arrays,
        target_fingerprint=plan['target_fingerprint'], weight_fingerprint=plan['weight_fingerprint'])
    (OUT / 'summary.json').write_text(json.dumps(status, indent=2, sort_keys=True) + '\n')
    print(json.dumps(status, sort_keys=True))


if __name__ == '__main__':
    main()
