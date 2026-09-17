"""Read-only pension accounting on pinned historical endpoint checkpoints.

No Bellman, market, population or calibration solve is performed. The inherited
snapshots supply serialized distributions; their old equilibrium certificates
do not certify the Social Security budget computed here.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import sys


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--helper', type=Path, required=True)
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ[key] = '1'
    if digest(args.contract) != args.contract_sha256 or digest(args.helper) != args.helper_sha256:
        raise ValueError('Fiscal audit contract/helper hash changed')
    c = json.loads(args.contract.read_text())
    root = args.source_root.resolve()
    for name, expected in c['source_sha256'].items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or digest(path) != expected:
            raise ValueError(f'Inherited source hash mismatch: {name}')
    for key in ('normalized_checkpoint', 'terminal_checkpoint'):
        if digest(c[key]) != c[key + '_sha256']:
            raise ValueError(f'Inherited checkpoint hash mismatch: {key}')
    spec = importlib.util.spec_from_file_location('audited_social_security', args.helper)
    accounting = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(accounting)
    sys.path[:0] = [str(root/'code/model/tools'), str(root/'code/model')]
    with gzip.open(c['normalized_checkpoint'], 'rb') as stream:
        initial = pickle.load(stream)
    old = initial['old']
    results = []
    for name, g in (('pre_announcement_stationary', old.stationary_g_pre),
                    ('announced_2007_initial_state', old.initial_state.g_pre)):
        results.append(dict(state=name, **accounting.fiscal_accounts(g, old.parameters)))
    with gzip.open(c['terminal_checkpoint'], 'rb') as stream:
        terminal = pickle.load(stream)
    results.append(dict(state='terminal_person_household_endpoint',
        **accounting.fiscal_accounts(terminal['fixed_point'].g_pre, terminal['parameters'])))
    packet = dict(status='completed_read_only_budget_audit', model_solves=0,
        contract_sha256=args.contract_sha256, helper_sha256=args.helper_sha256,
        inherited_source_files_verified=len(c['source_sha256']),
        source_root=str(root), rows=results,
        interpretation='Old distributions only; implied pension or tax is an accounting diagnostic, not a recomputed equilibrium',
        budget_scope='Payroll tax and pensions per household head; property tax and rebates excluded')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(packet, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(packet, allow_nan=False))


if __name__ == '__main__':
    main()
