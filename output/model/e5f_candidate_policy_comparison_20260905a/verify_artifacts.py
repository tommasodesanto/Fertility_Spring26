#!/usr/bin/env python3
"""Verify completed output receipts on Torch, then their collected local copies."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def verify(root, local):
    manifest_path = root / 'artifact_manifest.json'
    if local:
        manifest = read(manifest_path)
        missing_checkpoints = []
        checked = 0
        for relative, digest in manifest['files'].items():
            path = root / relative
            if not path.exists() and relative.endswith('.pkl.gz'):
                missing_checkpoints.append(relative)
                continue
            if sha(path) != digest:
                raise RuntimeError(f'Collected artifact differs: {relative}')
            checked += 1
        result = dict(status='pass', scope='local collected artifacts',
            locally_verified_files=checked, remote_only_checkpoints=missing_checkpoints,
            remote_manifest_sha256=sha(manifest_path),
            remote_validation=manifest['validation'])
        write(root / 'collection_receipt.json', result)
        print(json.dumps({k:v for k,v in result.items() if k!='remote_only_checkpoints'}, indent=2))
        return
    receipts = []
    for stage, expected in [('smoke_diagnostics', 10), ('full_diagnostics', 55)]:
        found = sorted((root / stage).glob('*/*/receipt.json'))
        if len(found) != expected:
            raise RuntimeError(f'{stage}: expected {expected} completed dates, found {len(found)}')
        receipts.extend(found)
        observers = sorted((root / stage).glob('*/observer_receipt.json'))
        if len(observers) != 5 or any(read(p)['status'] != 'complete' for p in observers):
            raise RuntimeError(f'{stage}: observer loop incomplete')
    tax = sorted((root / 'rebate_standard_diagnostics').glob('*/receipt.json'))
    if len(tax) != 11:
        raise RuntimeError('Three tax equilibria and eight component cells required')
    receipts.extend(tax)
    original_hash_count = 0
    for path in receipts:
        receipt = read(path)
        if receipt['status'] != 'complete':
            raise RuntimeError(f'Incomplete receipt: {path}')
        hashes = receipt.get('artifact_hashes', receipt.get('hashes'))
        if len(hashes) != 18:
            raise RuntimeError(f'Expected checkpoint and seventeen graphs: {path}')
        for relative, digest in hashes.items():
            if sha(path.parent / relative) != digest:
                raise RuntimeError(f'Recorded hash differs: {path.parent / relative}')
            original_hash_count += 1
        for bounds in receipt['policy_arrays']['probabilities'].values():
            if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
                raise RuntimeError(f'Invalid probability: {path}')
    if read(root / 'saving_flag_check/summary.json')['status'] != 'complete':
        raise RuntimeError('Fixed-price check incomplete')
    if read(root / 'full/report/summary.json')['status'] != 'complete_post2023_policy_mechanism_comparison':
        raise RuntimeError('Full collector incomplete')
    if read(root / 'shapley/summary.json')['status'] != 'complete':
        raise RuntimeError('Shapley decomposition incomplete')
    files = {}
    folders = ('smoke', 'smoke_diagnostics', 'full', 'full_diagnostics',
               'rebated_impact', 'rebate_standard_diagnostics', 'rebate_diagnostics',
               'shapley', 'saving_flag_check')
    for folder in folders:
        for path in sorted((root / folder).rglob('*')):
            if path.is_file() and path.suffix in ('.json', '.csv', '.png', '.pdf', '.gz'):
                files[str(path.relative_to(root))] = sha(path)
    validation = dict(status='pass', diagnostic_packets=len(receipts),
        original_recorded_hashes_verified=original_hash_count,
        full_policy_dates=55, smoke_policy_dates=10, tax_equilibria_and_cells=11,
        fixed_price_solve_packets=6,
        standard_graphs=17*(len(receipts)+6), manifest_files=len(files))
    write(manifest_path, dict(validation=validation, files=files))
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', action='store_true')
    args = parser.parse_args()
    verify(Path(__file__).resolve().parent, args.local)
