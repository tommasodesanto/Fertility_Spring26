#!/usr/bin/env python3
"""Prepare a new, explicit overnight lane from the reviewed integration contract.

No model solve or submission. Model source and empirical targets stay pinned;
identity weighting is a separately named diagnostic objective in native units.
"""
import argparse
import copy
import hashlib
import json
import math
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


def weighting_update(objective, provenance, lane):
    descriptions = {
        'primary': 'Retained working minimum-distance weights, including early-fertility weight 100.',
        'identity': 'Diagnostic identity weights on all 13 scored raw moment gaps in native units; not unit-invariant or an adopted optimal weighting matrix.',
        'early_fertility_3000': 'Diagnostic early-fertility weight 3000; all other primary weights retained. No target dropped; not an estimated optimal weight.'}
    rationale = descriptions[lane]
    for rows, key, name in [(objective['target_rows'], 'actual_weight', 'restriction_id'),
                            (provenance['target_rows'], 'weight', 'id')]:
        early = 0
        for row in rows:
            if row[name] == 'early_fertility':
                early += 1
            if row[key] is not None and (lane == 'identity' or
                    lane == 'early_fertility_3000' and row[name] == 'early_fertility'):
                row[key] = 1.0 if lane == 'identity' else 3000.0
                row['weight_status'] = rationale
        if lane == 'early_fertility_3000' and early != 1:
            raise ValueError('Expected exactly one early-fertility row')
    return rationale


def initial_case(c, objective, folder):
    folder = Path(folder).resolve(strict=True)
    receipt = read(folder / 'receipt.json')
    if (receipt.get('status') != 'verified_provisional_calibration_point'
            or receipt['target_system_sha256'] != c['objective']['sha256']
            or receipt['source_manifest_sha256'] != c['source_manifest']['sha256']
            or sha(folder / 'initial_state.pkl.gz') != receipt['case_checkpoint_sha256']):
        raise ValueError('Initial candidate provenance/checkpoint differs from reference')
    restrictions = {r['parameter']: r for r in objective['parameter_restrictions']}
    point = receipt['point']
    if set(point) != set(restrictions):
        raise ValueError('Initial candidate parameter set differs')
    for name, value in point.items():
        r = restrictions[name]
        if not math.isfinite(value) or not r['lower'] <= value <= r['upper']:
            raise ValueError('Initial candidate outside bounds: ' + name)
    return point, {'initial_case_receipt': pin(folder / 'receipt.json'),
                   'initial_case_checkpoint': pin(folder / 'initial_state.pkl.gz')}


def prepare(a):
    if a.reference_sha256 and sha(a.reference_contract) != a.reference_sha256:
        raise ValueError('Reference contract changed')
    c = read(a.reference_contract)
    if sha(c['objective']['path']) != c['objective']['sha256']:
        raise ValueError('Reviewed objective changed')
    objective = read(c['objective']['path'])
    provenance = read(c['files']['target_provenance']['path'])
    if bool(a.local_hostname) != bool(a.local_authorization_id):
        raise ValueError('Local hostname and authorization must be explicit together')
    if a.initial_case:
        c['initial_point'], extra_pins = initial_case(c, objective, a.initial_case)
        c['files'].update(extra_pins)
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
    rationale = weighting_update(objective, provenance, lane)
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
    if a.local_hostname:
        c['execution'] = {'kind': 'local', 'hostname': a.local_hostname,
                          'authorization_id': a.local_authorization_id}
    c['normalization']['warm_price'] = not a.cold
    c['budget'].update(workers=a.workers, points_per_round=a.workers,
                       absolute_end_epoch=a.absolute_end_epoch)
    if a.rounds is not None:
        c['budget']['rounds'] = a.rounds
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
    p.add_argument('--reference-sha256')
    p.add_argument('--runtime-tools', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, required=True, choices=range(1, 25))
    p.add_argument('--weighting', choices=['primary', 'identity', 'early_fertility_3000'], required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--absolute-end-epoch', type=float, required=True)
    p.add_argument('--cold', action='store_true')
    p.add_argument('--initial-case', type=Path)
    p.add_argument('--local-hostname')
    p.add_argument('--local-authorization-id')
    p.add_argument('--rounds', type=int, choices=range(1, 31))
    print(json.dumps(prepare(p.parse_args()), indent=2))
