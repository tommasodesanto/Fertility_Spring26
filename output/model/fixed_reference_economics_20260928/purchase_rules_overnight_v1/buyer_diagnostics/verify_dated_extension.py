"""Read-only authentication of accepted quarter-rule T48 date-zero packets."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path('/scratch/td2248/projects/purchase_mechanism_horizon_extension_v1')
STORE = Path('/scratch/td2248/projects/purchase_mechanism_v1')
MANIFEST_SHA = '19db8dfcc928c4cdea5810f70c4e62a9de65372bd038646d0bc2d49532913879'
CASES = {'case_06_quarter_control_h48': ('quarter', 'control', 48),
         'case_07_quarter_temporary_h48': ('quarter', 'temporary', 48)}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def authenticate(root: Path, store: Path, case: str, date: str) -> dict:
    require(case in CASES and date == 'date_000', 'Only accepted quarter T48 date zero is routed')
    arm, kind, horizon = CASES[case]
    manifest_path = store / 'selection/manifest.json'
    require(sha(manifest_path) == MANIFEST_SHA, 'Immutable selection manifest differs')
    selected = read(manifest_path)['arms'][arm]
    chosen = store / 'selection' / f'selected_{arm}.json'
    require(sha(chosen) == selected['selected_json_sha256'], 'Selected arm JSON differs')
    selected_json = read(chosen)
    require(selected_json['status'] == 'postchecked' and selected_json['arm'] == arm
            and int(selected_json['chain']) == int(selected['chain']), 'Selected arm identity differs')
    snapshot = store / 'selected_postchecks' / f"chain_{selected['chain']}"
    require(selected['snapshot_remote_root'] == str(snapshot)
            and selected_json['snapshot_remote_root'] == str(snapshot)
            and selected_json['origin'] == selected['origin']
            and selected_json['remote_root'] == selected['physical_remote_root']
            and selected_json.get('parent_remote_root') == selected['parent_remote_root'],
            'Selected source provenance differs')
    selected_completed = snapshot / 'postcheck/completed.json'
    require(sha(selected_completed) == selected['completed_sha256']
            and read(selected_completed)['status'] == 'selected_numerically_verified',
            'Selected postcheck differs')
    run = root / 'results' / case / 'run'
    require(read(run.parent / 'launcher_terminal.json')['exit_code'] == 0,
            'Mechanism launcher did not finish successfully')
    contract = read(run / 'run_contract.json')
    require(contract['arm'] == arm and contract['kind'] == kind
            and int(contract['horizon']) == horizon
            and contract['selected_json_sha256'] == selected['selected_json_sha256']
            and contract['selected_sha256'] == selected['completed_sha256']
            and contract['fixed_H0'] is True,
            'Dated contract differs from selected fit')
    completed = run / ('completed.json' if kind == 'control' else 'dated_path/completed.json')
    accepted = read(completed)
    require(accepted['status'] == 'passed' and accepted['arm'] == arm
            and accepted['kind'] == kind and int(accepted['horizon']) == horizon
            and accepted['terminal']['all_checks_pass'] is True
            and len(accepted['rows']) == len(accepted['phi_path']) == horizon,
            'Dated completion or terminal gate failed')
    if kind == 'temporary':
        dated_root = read(run / 'dated_path/root.json')
        require(dated_root['converged'] is True and all(dated_root['gates'].values()),
                'Dated market/fiscal root failed')
    accepted_path = Path(accepted['accepted_mapping'])
    require(accepted_path.is_relative_to('/work/results/run'), 'Accepted mapping escaped run')
    mapping_rel = accepted_path.relative_to('/work/results/run')
    require(mapping_rel == Path('control') if kind == 'control'
            else mapping_rel.parent == Path('dated_path') and mapping_rel.name.startswith('mapping_'),
            'Accepted mapping directory differs')
    mapping = read(run / mapping_rel / 'mapping.json')
    require(all(mapping['gates'].values()) and len(mapping['rows']) == horizon,
            'Accepted native mapping gates failed')
    packets = [p for p in mapping['diagnostic_packets'] if int(p['period']) == 0]
    require(len(packets) == 1, 'Unique accepted date-zero packet missing')
    pin = packets[0]
    expected = accepted_path / date / 'diagnostic_packet.pkl.gz'
    require(Path(pin['path']) == expected, 'Saved mapping does not pin exact date-zero path')
    packet_rel = expected.relative_to('/work/results/run')
    packet = run / packet_rel
    require(packet.is_file() and sha(packet) == pin['sha256'], 'Accepted date-zero packet SHA differs')
    phi = float(accepted['phi_path'][0])
    require(phi == (0.8 if kind == 'control' else 1.0), 'Date-zero financed share differs')
    return {'case': case, 'relative_packet': str(packet_rel), 'observed_phi': phi,
            'packet_sha256': pin['sha256'], 'completed_sha256': sha(completed),
            'mapping_sha256': sha(run / mapping_rel / 'mapping.json'),
            'selection_manifest_sha256': MANIFEST_SHA, 'no_model_solve': True}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, default=ROOT)
    ap.add_argument('--store', type=Path, default=STORE)
    ap.add_argument('--case', required=True)
    ap.add_argument('--date', default='date_000')
    a = ap.parse_args()
    print(json.dumps(authenticate(a.root, a.store, a.case, a.date), sort_keys=True))
