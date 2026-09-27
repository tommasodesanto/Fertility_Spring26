#!/usr/bin/env python3
"""Describe warm/cold numerical acceptance evidence without promoting a run."""
import argparse
import csv
import json
import math
from pathlib import Path


def read(path):
    return json.loads(Path(path).read_text())


def rows(path, key):
    with Path(path).open(newline='') as f:
        result = list(csv.DictReader(f))
    mapped = {r[key]: r for r in result}
    assert len(mapped) == len(result)
    return mapped


def cases(root):
    done = read(root / 'complete.json')
    assert done['status'] == 'exact_loop_smoke_passed'
    records = done['records']
    assert len(records) == 2 and all(r['status'] == 'success' for r in records)
    return [Path(r['case_path']) for r in records]


def compare(warm_root, cold_root, output):
    warm, cold = cases(warm_root), cases(cold_root)
    out = Path(output); out.mkdir(parents=True, exist_ok=False)
    a, b = read(warm[0] / 'receipt.json'), read(cold[0] / 'receipt.json')
    assert a['point'] == b['point'], 'Different structural parameter points'
    assert a['normalization_inputs']['warm_price'] is True
    assert b['normalization_inputs']['warm_price'] is False
    for field in ('initial_psi', 'initial_step', 'maximum_stationary_solves'):
        assert a['normalization_inputs'][field] == b['normalization_inputs'][field]
    comparisons = {}
    for filename, key in [('target_fit.csv', 'moment'), ('parameters.csv', 'parameter')]:
        x, y = rows(warm[0] / filename, key), rows(cold[0] / filename, key)
        assert x.keys() == y.keys(), filename
        result = []
        for name in x:
            if key == 'moment':
                for field in ('target', 'weight'):
                    assert x[name][field] == y[name][field], (name, field)
                field = 'model'
            else:
                for field in ('lower', 'upper'):
                    assert x[name][field] == y[name][field], (name, field)
                field = 'estimate'
            u, v = float(x[name][field]), float(y[name][field])
            assert math.isfinite(u) and math.isfinite(v)
            result.append({key: name, 'warm': u, 'cold': v, 'difference': u-v,
                           'relative_difference': (u-v)/max(abs(v), 1e-12),
                           'target': x[name].get('target', ''),
                           'weight': x[name].get('weight', ''),
                           'warm_loss_contribution': x[name].get('loss_contribution', ''),
                           'cold_loss_contribution': y[name].get('loss_contribution', '')})
        with (out / filename).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(result[0]))
            writer.writeheader(); writer.writerows(result)
        comparisons[filename] = result
    totals = {}
    for label, paths in [('warm', warm), ('cold', cold)]:
        totals[label] = []
        for path in paths:
            r = read(path / 'receipt.json')
            ledger = read(path / 'stationary_solves.json')
            assert len(ledger) == r['objective_stationary_solves']
            assert all(row['status'] == 'completed' for row in ledger)
            assert r['market_residual'] <= 2e-4
            assert all(row['market_residual'] < 2.5e-5 for row in ledger)
            assert abs(r['normalization']['completed_fertility']-2.1) <= 5e-4
            assert len(list((path / 'standard_diagnostics').glob('*.png'))) == 17
            totals[label].append(dict(case=str(path), loss=r['loss'], price=r['price'],
                psi_child=r['normalization']['psi_child'], wall_seconds=r.get('wall_seconds'),
                solve_seconds=sum(row['seconds'] for row in ledger),
                stationary_solves=len(ledger), price_evaluations=sum(row['price_evaluations'] for row in ledger)))
    result = dict(status='comparison_prepared_for_lead_review', cases=totals,
                  comparisons=comparisons, promoted=False,
                  note='All model-moment and parameter differences are reported. '
                       'Equilibrium/fertility gates do not themselves bound every moment difference.')
    (out / 'comparison.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--warm-smoke', type=Path, required=True)
    p.add_argument('--cold-smoke', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    result = compare(a.warm_smoke, a.cold_smoke, a.output)
    print(json.dumps({'status': result['status'], 'cases': result['cases']}, indent=2))
