#!/usr/bin/env python3
"""Inspect submitted utility checkpoints and repeats without solving the model."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_e5f_earnings_entry_checkpoint import _load_checkpoint, audit


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def array_record(value):
    a = np.ascontiguousarray(value)
    return {'shape': list(a.shape), 'dtype': str(a.dtype),
            'finite': bool(np.isfinite(a).all()),
            'sha256': hashlib.sha256(a.view(np.uint8)).hexdigest()}


def signature(row, source):
    score_path = Path(row['score_path'])
    score = read(score_path)
    raw = score_path.parents[1] / 'raw/repetition_01'
    summary = read(raw / 'summary.json')
    checkpoint = raw / 'initial_state.pkl.gz'
    checkpoint_sha = sha(checkpoint)
    assert checkpoint_sha == score['checkpoint_sha256'] == summary['checkpoint_sha256']
    packet = _load_checkpoint(checkpoint, source)
    targets = []
    for target in score['target_fit']:
        assert target['model_checkpoint_sha256'] == checkpoint_sha
        targets.append({k: v for k, v in target.items()
                        if k not in ('model_checkpoint_sha256', 'model_source_path')})
    assert len(targets) == 13
    out = {'price': summary['price'],
           'V': array_record(packet['evaluation'].policy.V),
           'g': array_record(packet['evaluation'].g_current),
           'psi': float(packet['parameters'].psi_child),
           'moments': summary['legacy_stationary_moments'],
           'early_measurement': summary['early_measurement'],
           'fiscal': summary['fiscal'], 'household_budget': summary['household_budget'],
           'normalization': {k: v for k, v in score.get('normalization', {}).items()
                             if k != 'stationary_solve_seconds'},
           'target_fit': targets, 'parameters': score['parameters'], 'loss': score['loss']}
    assert out['V']['finite'] and out['g']['finite']
    assert float(row['objective']) == float(score['loss'])
    return out, checkpoint_sha


def nonfinite(value, path='$'):
    result = []
    if isinstance(value, float) and not math.isfinite(value):
        result.append({'path': path, 'kind': 'nan' if math.isnan(value) else str(value)})
    elif isinstance(value, dict):
        for key, item in value.items():
            result.extend(nonfinite(item, path + '.' + key))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            result.extend(nonfinite(item, path + f'[{i}]'))
    return result


def differences(a, b, path='$'):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return [path + ': key sets differ']
        return [p for k in a for p in differences(a[k], b[k], path + '.' + k)]
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [path + ': lengths differ']
        return [p for i, (x, y) in enumerate(zip(a, b))
                for p in differences(x, y, path + f'[{i}]')]
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return []  # Location is separately reported; never assert finite diagnostics.
    return [] if a == b else [path + ': values differ']


def safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return {'nonfinite_value': 'nan' if math.isnan(value) else str(value)}
    if isinstance(value, dict):
        return {k: safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [safe(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.bundle / 'source'
    results = args.bundle / 'results'
    for cell in ('B_floor', 'B_shares', 'D_floor', 'D_shares'):
        selected = read(results / 'rollup' / cell / 'best_so_far.json')
        checkpoint = Path(selected['score_path']).parents[1] / 'raw/repetition_01/initial_state.pkl.gz'
        receipt = audit(checkpoint, source)
        write(args.output / (cell + '_entry_audit.json'), receipt)
        repeats = [results / 'verification' / cell / f'selected_repeat_{i:02d}/status.json'
                   for i in (1, 2)]
        if not all(p.is_file() for p in repeats):
            write(args.output / (cell + '_comparison.json'), {
                'status': 'two_repetitions_unavailable', 'selected_case_id': selected['case_id'],
                'selected_checkpoint_sha256': receipt['checkpoint_sha256'],
                'repeat_status_files_present': [p.is_file() for p in repeats]})
            print(cell + ': entry audited; repetitions unavailable', flush=True)
            continue
        base, base_sha = signature(selected, source)
        comparison = []
        for p in repeats:
            row = read(p)
            assert row['status'] == 'completed'
            repeated, repeated_sha = signature(row, source)
            comparison.append({'repeat_case_id': row['case_id'],
                               'checkpoint_sha256': repeated_sha,
                               'different_paths': differences(base, repeated),
                               'nonfinite_locations': nonfinite(repeated)})
        all_equal = all(not x['different_paths'] for x in comparison)
        write(args.output / (cell + '_comparison.json'), {
            'status': 'finite_values_and_array_bytes_equal_with_reported_nonfinite_locations'
                      if all_equal else 'scientific_difference',
            'selected_case_id': selected['case_id'], 'selected_checkpoint_sha256': base_sha,
            'original_signature': safe(base), 'nonfinite_locations': nonfinite(base),
            'comparisons_against_original_selected': comparison,
            'exclusions': ['checked target provenance fields', 'normalization.stationary_solve_seconds'],
            'limitation': 'Matching NaN locations are explicitly disclosed; their scientific meaning needs review.'})
        print(cell + ': comparisons written; equal=' + str(all_equal), flush=True)


if __name__ == '__main__':
    main()
