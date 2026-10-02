"""Read only completed, fresh-postchecked Torch calibration chains.

This standard-library script is sent to Torch on stdin by collect.py. It does
not solve a model or alter remote outputs.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

REMOTE = Path(os.environ.get('PURCHASE_RESULTS_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/results'))
PACKET = Path(os.environ.get('PURCHASE_PACKET_ROOT', '/scratch/td2248/projects/purchase_rules_overnight_v1/source/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1'))
PLAN = json.loads((PACKET / 'plan.json').read_text())
CONTRACT = PLAN['base_target_contract']


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def check_chain(index):
    arm = 'hard' if index < 24 or 48 <= index <= 52 else 'quarter'
    root = REMOTE / f'chain_{index}'
    completed = root / 'postcheck/completed.json'
    if not completed.exists():
        return dict(chain=index, arm=arm, status='postcheck_pending', search_completed=(root / 'search/search_completed.json').exists())
    report = json.loads(completed.read_text())
    if report.get('status') != 'selected_numerically_verified':
        return dict(chain=index, arm=arm, status=report.get('status', 'invalid_postcheck'), error='Selected postcheck did not pass')
    selected = report['selected']
    checked = report['selected_postcheck']
    if checked['status'] != 'passed' or checked['weight_fingerprint'] != PLAN['weight_fingerprint']:
        raise ValueError('postcheck status or weight fingerprint drift')
    if canonical(CONTRACT) != PLAN['target_fingerprint']:
        raise ValueError('target contract fingerprint drift')
    base = root / 'postcheck/selected_postcheck/phase_b_ge'
    full = base / 'selected_root'
    target = rows(full / 'target_fit.csv')
    parameter = rows(full / 'parameters.csv')
    if len(target) != 14 or len(parameter) != 31:
        raise ValueError('target or parameter row count drift')
    identity = [{k: x[k] for k in ('moment', 'target', 'weight', 'role')} for x in target]
    if identity != CONTRACT or target != checked['base_target_fit']:
        raise ValueError('target table identity or postcheck mismatch')
    if len([x for x in target if x['role'] == 'scored']) != 10:
        raise ValueError('ten scored targets required')
    loss = sum(float(x['loss_contribution']) for x in target if x['loss_contribution'])
    if not math.isclose(loss, float(checked['weighted_loss']), rel_tol=0, abs_tol=1e-8):
        raise ValueError('loss contribution mismatch')
    if not math.isclose(loss, float(selected['loss']), rel_tol=0, abs_tol=1e-8):
        raise ValueError('search/postcheck loss mismatch')
    bound = PLAN['bounds']
    for row in parameter:
        key = row['parameter']
        if key in bound:
            if not (math.isclose(float(row['lower']), float(bound[key][0]), abs_tol=1e-12) and
                    math.isclose(float(row['upper']), float(bound[key][1]), abs_tol=1e-12)):
                raise ValueError('free parameter bound drift: ' + key)
    if len([x for x in parameter if x['parameter'] in bound]) != 10:
        raise ValueError('ten free coordinates required')
    closure = json.loads((full / 'closure.json').read_text())
    if closure['population_scale'] != 1. or closure['normalized_population'] != 1.:
        raise ValueError('population normalization drift')
    if closure['standard_plot_count'] != 17 or len(list((full / 'standard_diagnostics').glob('*.png'))) != 17:
        raise ValueError('standard diagnostic set incomplete')
    if abs(float(closure['renewal_residual'])) > 1e-6 or abs(float(closure['absolute_housing_residual'])) > 1e-6:
        raise ValueError('renewal or housing closure failed')
    arrays = base / 'selected_repeat/stage/solution_arrays.npz'
    if not arrays.is_file():
        raise ValueError('selected native checkpoint missing')
    # NPZ members are checked by the buyer and policy workers when loaded.
    return dict(chain=index, arm=arm, status='postchecked', loss=loss,
                price=float(checked['price']), H0=float(checked['H0_derived']),
                target_fingerprint=PLAN['target_fingerprint'], weight_fingerprint=PLAN['weight_fingerprint'],
                selected_parameters=selected['parameters'], remote_root=str(root),
                remote_report=str(full), remote_arrays=str(arrays), target_fit=target,
                parameters=parameter, closure=closure,
                report_sha256={name: sha(full / name) for name in ('target_fit.csv', 'parameters.csv', 'closure.json')},
                native_arrays_bytes=arrays.stat().st_size)


if __name__ == '__main__':
    result = {'target_fingerprint': PLAN['target_fingerprint'],
              'weight_fingerprint': PLAN['weight_fingerprint'], 'chains': [], 'errors': []}
    first = int(os.environ.get('PURCHASE_CHAIN_FIRST', '0'))
    last = int(os.environ.get('PURCHASE_CHAIN_LAST', '48'))
    for index in range(first, last):
        try:
            result['chains'].append(check_chain(index))
        except Exception as error:
            result['errors'].append(dict(chain=index, error=f'{type(error).__name__}: {error}'))
    if result['errors']:
        result['status'] = 'rejected_invalid_postcheck'
    else:
        result['status'] = 'valid_snapshot'
    print(json.dumps(result, sort_keys=True, allow_nan=False))
