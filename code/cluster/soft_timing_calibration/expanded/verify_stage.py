"""Verify the staged source and the mounted execution contract without model solves."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

REMOTE = Path('/scratch/td2248/projects/soft_timing_calibration_20261002_v3')
CONTAINER_STAGE = Path('/work/deployment')
REPO = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
NORMAL = Path('output/model/fixed_reference_economics_20260928/normalized_calibration_v2')
TIMING = Path('output/model/fixed_reference_economics_20260928/purchase_timing_sandbox_v1')
SOFT = Path('output/model/fixed_reference_economics_20260928/soft_timing_review_v1/soft_selected.json')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    stage = CONTAINER_STAGE if sys.argv[1:] == ['--container'] else REMOTE
    inventory = json.loads((stage / 'inventory.json').read_text())
    mounted = sys.argv[1:] == ['--container']
    for rel, expected in inventory['files'].items():
        path = (REPO if mounted else stage / 'source') / rel
        if not path.is_file() or sha(path) != expected:
            raise SystemExit(f'Stage source mismatch: {rel}')
    for name, expected in inventory['entrypoints'].items():
        path = stage / name
        if not path.is_file() or sha(path) != expected:
            raise SystemExit(f'Stage entrypoint mismatch: {name}')
    if mounted:
        plan = json.loads((REPO / NORMAL / 'plan.json').read_text())
        timing = json.loads((REPO / TIMING / 'manifest.json').read_text())
        selected = json.loads((REPO / SOFT).read_text())
        for key in ('target_fingerprint', 'weight_fingerprint'):
            if len({inventory[key], plan[key], timing[key]}) != 1:
                raise SystemExit(f'Mixed {key}')
        if selected['selected']['weight_fingerprint'] != inventory['weight_fingerprint']:
            raise SystemExit('Selected soft weight fingerprint mismatch')
        selected_contract = [{k: row[k] for k in ('moment', 'target', 'weight', 'role')}
                             for row in selected['selected']['target_fit']]
        if selected_contract != plan['base_target_contract']:
            raise SystemExit('Selected soft target contract mismatch')
        source = REPO / selected['source']
        if sha(source) != inventory['selected_source_sha256'] or sha(source) != selected['source_sha256']:
            raise SystemExit('Selected checkpoint source mismatch')
        for rel, expected in json.loads((REPO / NORMAL / 'source_pins.json').read_text()).items():
            if sha(REPO / rel) != expected:
                raise SystemExit(f'Executed source pin mismatch: {rel}')
        for pair, pins in timing['source_pairs'].items():
            a, b = pair.split('|')
            if sha(REPO / a) != pins['original'] or sha(REPO / b) != pins['sandbox']:
                raise SystemExit(f'Timing pair mismatch: {pair}')
    print(json.dumps({'status': 'passed', 'mode': 'container' if mounted else 'host',
                      'source_files': len(inventory['files']),
                      'target_fingerprint': inventory['target_fingerprint'],
                      'weight_fingerprint': inventory['weight_fingerprint']}))


if __name__ == '__main__':
    main()
