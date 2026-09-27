#!/usr/bin/env python3
"""Prepare a new, explicit overnight lane from the reviewed integration contract.

No model solve or submission. Model source and empirical targets stay pinned;
identity weighting is a separately named diagnostic objective in native units.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def pin(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': sha(path)}


def canonical(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def write(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def prepare(a):
    c = read(a.reference_contract)
    if sha(c['objective']['path']) != c['objective']['sha256']:
        raise ValueError('Reviewed objective changed')
    objective = read(c['objective']['path'])
    provenance = read(c['files']['target_provenance']['path'])
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    old_tools = Path(c['runtime_tools'])
    new_tools = a.runtime_tools.resolve(strict=True)
    for key, record in c['files'].items():
        p = Path(record['path'])
        if p.is_relative_to(old_tools):
            c['files'][key] = pin(new_tools / p.relative_to(old_tools))
        elif sha(p) != record['sha256']:
            raise ValueError('Reference changed: ' + str(p))
    inventory = {'root': str(new_tools), 'files': {
        str(p.relative_to(new_tools)): sha(p) for p in sorted(new_tools.rglob('*'))
        if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc', '.pyo')}}
    inventory['file_count'] = len(inventory['files'])
    write(out / 'runtime_inventory.json', inventory)
    c['runtime_tools'] = str(new_tools)
    c['files']['runtime_inventory'] = pin(out / 'runtime_inventory.json')
    for rel in inventory['files']:
        c['files']['runtime_file:' + rel] = pin(new_tools / rel)
    c['files']['overnight_generator'] = pin(__file__)
    c['files']['reviewed_reference_contract'] = pin(a.reference_contract)
    lane = a.weighting
    rationale = ('Retained working minimum-distance weights, including early-fertility weight 100.'
                 if lane == 'primary' else
                 'Diagnostic identity weights on all 13 scored raw moment gaps in native units; '
                 'not unit-invariant or an adopted optimal weighting matrix.')
    if lane == 'identity':
        for row in objective['target_rows']:
            if row['actual_weight'] is not None:
                row['actual_weight'] = 1.0
                row['weight_status'] = rationale
        for row in provenance['target_rows']:
            if row['weight'] is not None:
                row['weight'] = 1.0
                row['weight_status'] = rationale
    objective.update(contract_id='supervised_calibration_20260927_' + lane,
                     objective_name='supervised_calibration_20260927_' + lane,
                     status='author-authorized provisional overnight exploration',
                     weighting_experiment=lane, weighting_rationale=rationale)
    provenance.update(status='author-authorized provisional overnight exploration',
                      weighting_experiment=lane,
                      fitted_parameter_count=10, outer_search_parameter_count=9)
    write(out / 'target_provenance.json', provenance)
    c['files']['target_provenance'] = pin(out / 'target_provenance.json')
    objective['current_target_provenance'] = pin(out / 'target_provenance.json')
    write(out / 'objective.json', objective)
    c.update(objective=pin(out / 'objective.json'),
             objective_canonical_sha256=canonical(objective),
             target_weight_fingerprint=canonical(objective['target_rows']),
             seed=a.seed, execution={'kind': 'slurm'},
             overnight_lane=lane, weighting_rationale=rationale)
    c['normalization']['warm_price'] = not a.cold
    c['budget'].update(workers=a.workers, points_per_round=a.workers,
                       absolute_end_epoch=a.absolute_end_epoch)
    c['approval'] = {'authority': 'Tommaso September27: eight-hour supervised run; '
                     '24 cluster and10 local workers; numerical repairs and labeled weight experiments',
                     'production_authorized': True,
                     'condition': 'Exact-loop smoke and lead numerical acceptance before search'}
    c['status'] = 'reviewed_smoke'
    c['production_blockers'] = ['Exact-loop smoke and numerical acceptance pending']
    c['identification'].update(fitted_parameters=10, outer_search_parameters=9)
    for change in c['economic_changes']:
        if change['object'] in ('child_benefit_curvature', 'early_fertility'):
            change['status'] = 'lead-selected provisional setting under author overnight delegation'
    c['economic_changes'].append({'object': 'weighting',
        'status': 'primary retained' if lane == 'primary' else 'separate diagnostic only',
        'change': rationale})
    write(out / 'contract.json', c)
    return {'contract': pin(out / 'contract.json'), 'weighting': lane,
            'workers': a.workers, 'launch_performed': False}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-contract', type=Path, required=True)
    p.add_argument('--runtime-tools', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, required=True, choices=range(1, 25))
    p.add_argument('--weighting', choices=['primary', 'identity'], required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--absolute-end-epoch', type=float, required=True)
    p.add_argument('--cold', action='store_true')
    print(json.dumps(prepare(p.parse_args()), indent=2))
