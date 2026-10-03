"""Read-only production gate: require the best-start native smoke receipt on staged source."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from collect_torch import collect_one, read, tasks

REMOTE = Path('/scratch/td2248/projects/soft_timing_continuation_20261003_v1')
NORMAL = 'output/model/fixed_reference_economics_20260928/normalized_calibration_v2/source_pins.json'
TIMING = 'output/model/fixed_reference_economics_20260928/purchase_timing_sandbox_v1/manifest.json'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(stage: Path) -> dict:
    inventory = read(stage / 'inventory.json')
    rows = [collect_one(stage / 'results', 'smoke', arm, chain) for arm, chain in tasks('smoke')]
    for key in ('target_fingerprint', 'weight_fingerprint', 'selected_source_sha256'):
        if {row[key] for row in rows} != {inventory[key]}:
            raise AssertionError(f'Smoke/stage {key} mismatch')
    normal_digest = sha(stage / 'source' / NORMAL)
    timing_digest = sha(stage / 'source' / TIMING)
    stage_digest = sha(stage / 'inventory.json')
    for arm, chain in tasks('smoke'):
        launcher = stage / 'results' / f'smoke_{arm}_chain_{chain}'
        folder = launcher / 'run'
        start = read(launcher / 'launcher_start.json')
        contract = read(folder / 'start_contract.json')
        search = read(folder / 'search_contract.json')
        completed = read(folder / 'completed.json')
        child = read(folder / 'native_postcheck/completed.json')
        heartbeat = read(folder / 'heartbeat.json')
        cases = read(folder / 'cases.json')
        if start['wall_seconds'] != 5400 or start['deadline_epoch'] - start['start_epoch'] != 5400:
            raise AssertionError(f'{arm} smoke deadline differs from Slurm allocation')
        if start['stage_inventory_sha256'] != stage_digest:
            raise AssertionError(f'{arm} smoke ran under a different stage inventory')
        if start['maximum_objective_calls'] != 500 or start['final_native_reserve_seconds'] != 1800:
            raise AssertionError('Launcher search budget drift')
        if contract['objective_calls_max'] != 500 or contract['starts_count'] != 10:
            raise AssertionError('Continuation search/start contract drift')
        if contract['continuation_checkpoint_sha256'] != inventory['selected_continuation_checkpoint_sha256']:
            raise AssertionError('Best-start checkpoint drift')
        if contract['starts_file_sha256'] != inventory['start_plan_sha256']:
            raise AssertionError('Continuation start plan SHA drift')
        if contract['bounds']['beta_annual'] != [0.94, 0.99]:
            raise AssertionError('Original beta bound drift')
        if heartbeat['status'] != 'completed' or len(cases) != 2:
            raise AssertionError(f'{arm} smoke checkpoints incomplete')
        if contract['normalized_source_pins_sha256'] != normal_digest:
            raise AssertionError(f'{arm} original source pin identity mismatch')
        if arm == 'alternative' and contract['timing_manifest_sha256'] != timing_digest:
            raise AssertionError('Alternative timing manifest identity mismatch')
        if search['max_objective_calls'] != 2:
            raise AssertionError(f'{arm} smoke was not two-call loop')
        if completed['objective_calls'] != 2:
            raise AssertionError(f'{arm} smoke objective count mismatch')
        if completed['selected_postcheck']['status'] != 'passed':
            raise AssertionError(f'{arm} native selected postcheck failed')
        if child['status'] != 'full_native_postcheck_passed':
            raise AssertionError(f'{arm} fresh-child postcheck did not pass')
        if child['search_receipt_sha256'] != sha(folder / 'search_completed.json'):
            raise AssertionError(f'{arm} fresh child used a different search receipt')
        if child['target_fingerprint'] != inventory['target_fingerprint'] or child['weight_fingerprint'] != inventory['weight_fingerprint']:
            raise AssertionError(f'{arm} fresh-child target or weight drift')
        comparison = completed['smoke_fast_full_comparison']
        if comparison['status'] != 'matched_saved_full_baseline' or comparison['absolute_tolerance'] != 1e-10:
            raise AssertionError(f'{arm} fast/full comparison missing or changed')
        if comparison['checks']['target_fit.csv']['rows'] != 14 or comparison['checks']['parameters.csv']['rows'] != 31:
            raise AssertionError(f'{arm} fast/full row count drift')
    return {'status': 'continuation_smoke_passed_on_current_stage',
            'arms': [row['arm'] for row in rows],
            'target_fingerprint': inventory['target_fingerprint'],
            'weight_fingerprint': inventory['weight_fingerprint'],
            'selected_source_sha256': inventory['selected_source_sha256']}


def main() -> None:
    stage = Path(sys.argv[1]) if len(sys.argv) == 2 else REMOTE
    print(json.dumps(verify(stage), sort_keys=True))


if __name__ == '__main__':
    main()
