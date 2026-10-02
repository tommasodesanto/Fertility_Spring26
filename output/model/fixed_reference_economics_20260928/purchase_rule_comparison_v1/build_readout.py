"""Read saved four-case purchase-rule diagnostics; never solve the model."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CASES = ('hard80', 'hard100', 'quarter80', 'quarter100')
PRICE = 0.7152515073815459
H0 = 6.778473404808042


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def approx(a, b, tol=1e-9):
    return math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol)


def load_case(path, label):
    if not path.is_file():
        return None
    obj = read(path)
    if 'cases' in obj:
        assert obj.get('status') in ('fixed_price_diagnostic_passed',
                                     'fixed_price_quarter_saving_solvency_retry_passed',
                                     'fixed_price_quarter_saving_diagnostic_passed'), path
        matches = [c for c in obj['cases'] if c['label'] == label]
        assert len(matches) == 1, (path, label)
        return matches[0]
    if 'target_fit' in obj and 'parameters' in obj:
        assert obj.get('status') in ('passed_fixed_price_diagnostic',
                                      'fixed_coordinate_strict_purchase_experiment_passed'), path
        return obj
    raise ValueError(f'No complete native report in {path}')


def rows_for(case):
    assert len(case['target_fit']) == 14 and len(case['parameters']) == 31
    assert approx(case['price'], PRICE)
    if 'H0' in case:
        assert approx(case['H0'], H0)
    if 'H0_derived' in case:
        assert approx(case['H0_derived'], H0)
    for row in case['target_fit']:
        if row['loss_contribution']:
            assert approx(float(row['model']) - float(row['target']), row['gap'], 1e-7)
            assert approx(float(row['weight']) * float(row['gap'])**2,
                          row['loss_contribution'], 1e-7)
    loss = sum(float(r['loss_contribution'] or 0) for r in case['target_fit'])
    assert approx(loss, case['loss'], 1e-7)
    return loss


def write_csv(path, columns, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hard100', type=Path, default=BASE/'hard_zero_down_solvency_v2/collected/run/completed.json')
    ap.add_argument('--quarter100', type=Path, default=BASE/'quarter_saving_solvency_v2/collected/run/completed.json')
    ap.add_argument('--out', type=Path, default=HERE/'readout')
    args = ap.parse_args()
    paths = {
        'hard80': BASE/'strict_purchase_sandbox_v1/collected/run/completed.json',
        'hard100': args.hard100,
        'quarter80': BASE/'quarter_saving_benchmarks_v1/collected/run/latest_completed.json',
        'quarter100': args.quarter100,
    }
    loaded = {name: load_case(paths[name], name) for name in CASES}
    assert loaded['hard80'] is not None and loaded['quarter80'] is not None
    if 'target_fingerprint' in loaded['hard80']:
        old = loaded['hard80']
        for name in ('hard100', 'quarter100'):
            if paths[name].is_file():
                raw = read(paths[name])
                if 'target_fingerprint' in raw:
                    assert raw['target_fingerprint'] == old['target_fingerprint']
                    assert raw['weight_fingerprint'] == old['weight_fingerprint']
    point = read(BASE/'hard_zero_down_v1/input/strict80/incumbent.json')['parameters']
    assert len(point) == 10
    loss = {name: rows_for(case) for name, case in loaded.items() if case is not None}
    targets = [r['moment'] for r in loaded['hard80']['target_fit']]
    params = [r['parameter'] for r in loaded['hard80']['parameters']]
    for name, case in loaded.items():
        if case is None:
            continue
        assert [r['moment'] for r in case['target_fit']] == targets
        assert [r['parameter'] for r in case['parameters']] == params
        parameter_map = {r['parameter']: r for r in case['parameters']}
        for key, value in point.items():
            assert approx(parameter_map[key]['estimate'], value, 1e-10), (name, key)
        assert approx(parameter_map['delta_alpha_jump']['estimate'], 0)
        assert approx(parameter_map['financed_share']['estimate'],
                      1.0 if name.endswith('100') else 0.8)
        for old, new in zip(loaded['hard80']['target_fit'], case['target_fit']):
            assert old['role'] == new['role']
            for key in ('target', 'weight'):
                assert old[key] == new[key] or (old[key] and new[key] and approx(old[key], new[key]))
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    target_cols = ['moment', 'role', 'target', 'weight']
    target_cols += [f'{name}_{k}' for name in CASES for k in ('model', 'gap', 'loss_contribution')]
    target_rows = []
    for i, base in enumerate(loaded['hard80']['target_fit']):
        row = {key: base[key] for key in ('moment', 'role', 'target', 'weight')}
        for name in CASES:
            value = None if loaded[name] is None else loaded[name]['target_fit'][i]
            for key in ('model', 'gap', 'loss_contribution'):
                row[f'{name}_{key}'] = '' if value is None else value[key]
        target_rows.append(row)
    write_csv(out/'target_fit.csv', target_cols, target_rows)
    param_cols = ['parameter'] + [f'{name}_{key}' for name in CASES
                                 for key in ('estimate', 'lower', 'upper', 'near_bound', 'status')]
    param_rows = []
    for i, parameter in enumerate(params):
        row = {'parameter': parameter}
        for name in CASES:
            value = None if loaded[name] is None else loaded[name]['parameters'][i]
            for key in ('estimate', 'lower', 'upper', 'near_bound', 'status'):
                row[f'{name}_{key}'] = '' if value is None else value[key]
            if value is not None:
                if parameter in point:
                    row[f'{name}_status'] = 'fixed at strict-80 estimate for this diagnostic'
                elif parameter == 'delta_alpha_jump':
                    row[f'{name}_lower'] = '0.0'
                    row[f'{name}_upper'] = '0.0'
                    row[f'{name}_near_bound'] = ''
                    row[f'{name}_status'] = 'externally fixed at zero for this diagnostic'
                elif parameter == 'financed_share':
                    row[f'{name}_status'] = 'externally set policy input for this diagnostic'
        param_rows.append(row)
    write_csv(out/'parameters.csv', param_cols, param_rows)
    statuses = {name: ('passed' if loaded[name] is not None else 'pending_native_report')
                for name in CASES}
    receipt = dict(status='four_case_readout_complete' if all(loaded.values()) else 'partial_readout',
                   cases=statuses, loss=loss, target_rows=len(target_rows),
                   parameter_rows=len(param_rows), fixed_price=PRICE, fixed_H0=H0,
                   sources={name: dict(path=str(path.resolve()), sha256=sha(path)) if path.is_file()
                            else dict(path=str(path.resolve()), available=False)
                            for name, path in paths.items()},
                   prior_rejected_phi1_attempts=dict(hard='hard_zero_down_v1 job 19006866: negative-estate gate',
                                                     quarter='quarter_saving_benchmarks_v1 job 19007303: negative-estate gate'))
    (out/'verification.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    lines = ['# Four fixed-price purchase-rule diagnostics', '',
             'All cases hold the strict-80 price, $H_0$, $N_0=1$, and ten calibrated coordinates fixed. This is a partial-equilibrium comparison, not recalibration.', '',
             '| Case | Purchase eligibility | Saving rule | Gate | Loss |',
             '| --- | --- | --- | --- | ---: |']
    for name, rule, saving in (('hard80','strict wealth-only; $\\phi=0.8$','standard'),
                               ('hard100','strict wealth-only; $\\phi=1$','standard'),
                               ('quarter80','strict wealth-only; $\\phi=0.8$','quarter-saving constraint'),
                               ('quarter100','strict wealth-only; $\\phi=1$','quarter-saving constraint')):
        lines.append(f'| {name} | {rule} | {saving} | {statuses[name]} | {loss.get(name, "unavailable")} |')
    lines += ['', 'The ten estimated coordinates are fixed at their strict-80 values in all four cases. `delta_alpha_jump` is externally fixed at zero; any inherited source-table status calling it free is corrected in the display table. The failed first $\\phi=1$ attempts reached the negative-estate production gate before target-fit and parameter reports were serialized. Their missing cells are blank, not estimated. Successful isolated solvency retries populate those columns when collected.', '',
              'Complete 14-moment target fit: [target_fit.csv](target_fit.csv). Complete 31-parameter comparison: [parameters.csv](parameters.csv). Source hashes and gate states: [verification.json](verification.json).', '']
    (out/'RESULTS.md').write_text('\n'.join(lines))
    print(receipt['status'], statuses)


if __name__ == '__main__':
    main()
