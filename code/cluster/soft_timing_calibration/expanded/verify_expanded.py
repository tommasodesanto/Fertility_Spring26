"""Zero-solve gate for the 40 additional matched starts on the v3 stage."""
from __future__ import annotations

import ast
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

PARENT = Path('/scratch/td2248/projects/soft_timing_calibration_20261002_v2')
PACKET = Path('output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1')
DRIVER = Path('code/model/experiments/purchase_timing_sandbox/calibrate.py')
CHECKED_FUNCTIONS = ('checked_inputs', 'install_timing_observer', 'starts',
                     'completion_receipt', 'selected_repeat_path', 'objective')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def function_ast(path: Path, name: str) -> str:
    nodes = [n for n in ast.walk(ast.parse(path.read_text()))
             if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name]
    if len(nodes) != 1:
        raise AssertionError(f'Expected one {name} in {path}')
    return ast.dump(nodes[0], include_attributes=False)


def verify(stage: Path, parent: Path = PARENT) -> dict:
    inv = json.loads((stage / 'inventory.json').read_text())
    parent_inv = json.loads((parent / 'inventory.json').read_text())
    if sha(parent / 'stage.tar.gz') != inv['parent_stage_archive_sha256']:
        raise AssertionError('Passed v2 stage archive identity changed')
    parent_archive = stage / 'parent_calibrate.py'
    if sha(parent_archive) != inv['parent_driver_sha256']:
        raise AssertionError('Pinned v2 driver snapshot changed')
    if sha(parent / 'source' / DRIVER) != inv['parent_driver_sha256']:
        raise AssertionError('Passed v2 smoke driver does not match pinned snapshot')
    parent_files = parent_inv['files']
    for rel, expected in parent_files.items():
        if rel in (str(DRIVER), str(PACKET / 'driver_plan.json')):
            continue
        if inv['files'].get(rel) != expected:
            raise AssertionError(f'Unreviewed change to parent source: {rel}')
    if set(inv['files']) - set(parent_files) != {str(PACKET / 'expanded_start_plan.json')}:
        raise AssertionError('Unexpected new source input')
    for name in CHECKED_FUNCTIONS:
        if function_ast(parent_archive, name) != function_ast(stage / 'source' / DRIVER, name):
            raise AssertionError(f'Core driver function changed: {name}')
    if inv['target_fingerprint'] != parent_inv['target_fingerprint'] or \
            inv['weight_fingerprint'] != parent_inv['weight_fingerprint'] or \
            inv['selected_source_sha256'] != parent_inv['selected_source_sha256']:
        raise AssertionError('Economic source or objective fingerprint drift')

    table_path = stage / 'source' / PACKET / 'expanded_start_plan.json'
    if sha(table_path) != inv['expanded_start_sha256'] or \
            sha(table_path) != inv['files'][str(PACKET / 'expanded_start_plan.json')]:
        raise AssertionError('Expanded table SHA mismatch')
    table = json.loads(table_path.read_text())
    plan = json.loads((stage / 'source' / PACKET / 'driver_plan.json').read_text())
    for key in ('target_fingerprint', 'weight_fingerprint'):
        if table[key] != inv[key] or plan[key] != inv[key]:
            raise AssertionError(f'{key} mismatch')
    if table['source_checkpoint_sha256'] != inv['selected_source_sha256'] or \
            plan['selected_source_sha256'] != inv['selected_source_sha256']:
        raise AssertionError('Checkpoint mismatch')
    bounds, starts = table['bounds'], table['starts']
    if bounds != plan['bounds'] or len(bounds) != 10 or len(starts) != 24:
        raise AssertionError('Expanded starts/bounds mismatch')
    if starts[:4] != plan['starts'] or len(table['start_provenance']) != 24:
        raise AssertionError('Legacy starts or provenance mismatch')
    vectors = []
    for idx, row in enumerate(starts):
        if set(row) != set(bounds):
            raise AssertionError(f'Start {idx} coordinates differ from bounds')
        vals = tuple(float(row[key]) for key in sorted(bounds))
        for key, value in row.items():
            lo, hi = bounds[key]
            if not math.isfinite(value) or not lo <= value <= hi:
                raise AssertionError(f'Start {idx} violates {key} bounds')
        vectors.append(vals)
    if len(set(vectors)) != 24:
        raise AssertionError('Duplicate start vectors')
    groups = [item['group'] for item in table['start_provenance']]
    if {g: groups.count(g) for g in set(groups)} != \
            {'legacy': 4, 'historical': 7, 'near_selected': 8, 'broad': 5}:
        raise AssertionError('Start group counts changed')
    tasks = [('original', task + 4) if task < 20 else ('alternative', task - 16)
             for task in range(40)]
    if len(set(tasks)) != 40 or \
            {chain for arm, chain in tasks if arm == 'original'} != set(range(4, 24)) or \
            {chain for arm, chain in tasks if arm == 'alternative'} != set(range(4, 24)):
        raise AssertionError('Expanded task map incomplete')

    gate = subprocess.run([sys.executable, str(parent / 'verify_smoke_gate.py'), str(parent)],
                          check=True, capture_output=True, text=True)
    parent_smoke = json.loads(gate.stdout.strip().splitlines()[-1])
    if parent_smoke['status'] != 'both_smokes_passed_on_current_stage':
        raise AssertionError('Parent native smokes did not pass')
    return {'status': 'expanded_zero_solve_gate_passed', 'tasks': 40,
            'new_chains_per_arm': 20, 'total_chains_including_v2': 48,
            'expanded_start_sha256': sha(table_path),
            'parent_driver_sha256': sha(parent_archive),
            'parent_smoke_status': parent_smoke['status'],
            'target_fingerprint': inv['target_fingerprint'],
            'weight_fingerprint': inv['weight_fingerprint']}


if __name__ == '__main__':
    stage = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    parent = Path(sys.argv[2]) if len(sys.argv) > 2 else PARENT
    print(json.dumps(verify(stage, parent), sort_keys=True))
