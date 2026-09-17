"""Collect hash-pinned signed price-path finite differences; no model solves."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

SCHEMA = 'e5f_matched_pf_price_jacobian_v1'
STATUS = 'passed_conditional_historical_price_probe'
SHARED_FIELDS = (
    'arm', 'years', 'anchor_prices', 'target_fingerprint', 'path_date_count',
    'normalized_checkpoint_sha256', 'terminal_checkpoint_sha256',
    'checkpoint_sha256', 'selected_summary_sha256', 'demographic_sources',
    'source_sha256', 'initial_price_rule', 'terminal_preference_rule',
    'psi_path', 'transfer_path', 'terminal_price', 'supply_rule', 'probe_log_step',
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_hash(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError(f'Invalid SHA-256 pin: {name}')


def vector(value, n, name, *, positive=False):
    result = np.asarray(value, dtype=float)
    if (result.shape != (n,) or not np.isfinite(result).all()
            or (positive and np.any(result <= 0))):
        raise ValueError(f'Invalid {name} vector')
    return result


def equal_vectors(actual, expected, name):
    if not np.allclose(actual, expected, rtol=0., atol=1e-12):
        raise ValueError(f'Inconsistent {name}')


def gates_pass(gates):
    if not isinstance(gates, dict) or not gates:
        raise ValueError('A nonempty numerical-gate receipt is required')
    for name, gate in gates.items():
        if type(gate) is bool:
            passed = gate
        elif isinstance(gate, dict):
            passed = gate.get('passed') is True
            if 'value' not in gate or 'tolerance' not in gate:
                raise ValueError(f'Incomplete numerical gate: {name}')
            actual, tolerance = float(gate['value']), float(gate['tolerance'])
            passed = passed and np.isfinite([actual, tolerance]).all() and tolerance >= 0 and abs(actual) <= tolerance
        else:
            passed = False
        if not passed:
            raise ValueError(f'Failed numerical gate: {name}')


def load_case(folder):
    folder = Path(folder).resolve()
    contract_path, summary_path, path_file = (folder / name for name in
                                             ('contract.json', 'summary.json', 'transition_path.csv'))
    c, s = json.loads(contract_path.read_text()), json.loads(summary_path.read_text())
    if s.get('status') != STATUS or s.get('mapping_valid') is not True:
        raise ValueError(f'Price probe did not pass its mapping gates: {folder}')
    gates_pass(s.get('gates'))
    if c.get('contract_sha256') != s.get('contract_sha256') or not c.get('contract_sha256'):
        raise ValueError('Summary/originating contract identity differs')
    shared = {name: c[name] for name in SHARED_FIELDS}
    # Any additional input receipt hashes must also agree. The originating
    # per-probe contract hash is the only intentionally different SHA field.
    shared.update({name: value for name, value in c.items()
                   if name.endswith('_sha256') and name != 'contract_sha256'})
    shared['mapping_gate_names'] = sorted(s['gates'])
    for name, value in shared.items():
        if name.endswith('_sha256') and name != 'source_sha256':
            require_hash(value, name)
    require_hash(c['contract_sha256'], 'contract_sha256')
    if not isinstance(c['source_sha256'], dict) or not c['source_sha256']:
        raise ValueError('A nonempty source fingerprint map is required')
    for source, value in c['source_sha256'].items():
        require_hash(value, source)
    if set(c['demographic_sources']) != {'population_mid', 'births_mid', 'survival', 'vintage_2025', 'acs_headship'}:
        raise ValueError('All five demographic source pins are required')
    for name, entry in c['demographic_sources'].items():
        require_hash(entry['sha256'], name)
        if not Path(entry['path']).is_absolute():
            raise ValueError(f'Demographic source path must be absolute: {name}')
    if c['arm'] not in ('sequential', 'nested') or s['arm'] != c['arm']:
        raise ValueError('Invalid or inconsistent arm')
    n, coordinate = c['path_date_count'], c['probe_coordinate']
    if type(n) is not int or n < 1 or type(coordinate) is not int or not -1 <= coordinate < n:
        raise ValueError('Invalid path count or probe coordinate')
    if c['probe_log_step'] != .01:
        raise ValueError('Expected the contracted .01 log-price finite difference')
    years = c['years']
    if years != list(2007 + 4 * np.arange(n)) or s['years'] != years:
        raise ValueError('Invalid or inconsistent dated year vector')
    if s['target_fingerprint'] != c['target_fingerprint']:
        raise ValueError('Summary/contract target fingerprint differs')
    anchor = vector(c['anchor_prices'], n, 'anchor prices', positive=True)
    equal_vectors(vector(s['anchor_prices'], n, 'summary anchor', positive=True), anchor, 'anchor prices')
    expected = anchor.copy()
    if coordinate >= 0:
        expected[coordinate] *= np.exp(.01)
    prices = vector(s['prices'], n, 'actual prices', positive=True)
    equal_vectors(prices, vector(c['prices'], n, 'contract prices', positive=True), 'contract/summary prices')
    equal_vectors(prices, expected, 'single-coordinate prices')
    residual = vector(s['residual'], n, 'signed residual')
    csv_hash = digest(path_file)
    if s['artifact_sha256'].get('transition_path.csv') != csv_hash:
        raise ValueError('Transition path hash does not match the successful receipt')
    with path_file.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != n or [int(row['calendar_year']) for row in rows] != years:
        raise ValueError('Transition CSV omits, duplicates, or changes dates')
    demand = vector([float(row['housing_demand']) for row in rows], n, 'demand')
    supply = vector([float(row['housing_supply']) for row in rows], n, 'supply', positive=True)
    if np.any(demand < 0):
        raise ValueError('Housing demand must be nonnegative')
    signed = (demand - supply) / supply
    equal_vectors(residual, signed, 'signed (demand-supply)/supply residual')
    equal_vectors(vector([float(row['asset_price']) for row in rows], n, 'CSV prices', positive=True), prices, 'CSV/summary prices')
    if not np.isclose(float(s['maximum_market_residual']), np.max(np.abs(signed)), rtol=0, atol=1e-12):
        raise ValueError('Maximum residual receipt is inconsistent')
    return dict(coordinate=coordinate, shared=shared, prices=prices, residual=signed,
        provenance=dict(directory=str(folder),
            contract_sha256=digest(contract_path), summary_sha256=digest(summary_path),
            transition_path_sha256=csv_hash, originating_contract_sha256=c['contract_sha256']))


def collect(anchor_dir, probe_dirs):
    anchor = load_case(anchor_dir)
    if anchor['coordinate'] != -1:
        raise ValueError('Anchor must have probe_coordinate=-1')
    n = len(anchor['prices'])
    columns = {}
    for directory in probe_dirs:
        case = load_case(directory)
        j = case['coordinate']
        if j < 0 or j in columns:
            raise ValueError(f'Duplicate or invalid probe coordinate: {j}')
        if set(case['shared']) != set(anchor['shared']):
            raise ValueError('Mixed probe contract: shared input field set')
        for name in anchor['shared']:
            if case['shared'][name] != anchor['shared'][name]:
                raise ValueError(f'Mixed probe contract: {name}')
        columns[j] = case
    if set(columns) != set(range(n)):
        raise ValueError(f'Missing probe coordinates: {sorted(set(range(n)) - set(columns))}')
    J = np.column_stack([(columns[j]['residual'] - anchor['residual']) / .01 for j in range(n)])
    if J.shape != (n, n) or not np.isfinite(J).all():
        raise ValueError('Finite-difference Jacobian must be finite NxN')
    condition = float(np.linalg.cond(J))
    return dict(schema=SCHEMA, status='complete_validated_finite_difference_jacobian',
        prices=anchor['prices'].tolist(), residual=anchor['residual'].tolist(), jacobian=J.tolist(),
        derivative='d signed((housing_demand-housing_supply)/housing_supply) / d log(asset_price)',
        log_step=.01, condition_number=condition if np.isfinite(condition) else None,
        condition_number_finite=bool(np.isfinite(condition)),
        arm=anchor['shared']['arm'], years=anchor['shared']['years'],
        target_fingerprint=anchor['shared']['target_fingerprint'], shared_contract=anchor['shared'],
        provenance=dict(anchor=anchor['provenance'],
            columns=[dict(coordinate=j, **columns[j]['provenance']) for j in range(n)]),
        equilibrium_verified=False, calibrated_history=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--anchor', required=True, type=Path)
    parser.add_argument('--probes', required=True, nargs='+', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f'Refusing to overwrite {args.output}')
    result = collect(args.anchor, args.probes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + '.tmp')
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(args.output)


if __name__ == '__main__':
    main()
