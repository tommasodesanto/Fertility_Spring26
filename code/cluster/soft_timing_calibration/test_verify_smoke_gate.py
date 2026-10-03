"""Synthetic fixture for the read-only production gate; contains no model data."""
import csv
import hashlib
import json
import tempfile
from pathlib import Path

from verify_smoke_gate import verify


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + '\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(root):
    normal = root / 'source/output/model/fixed_reference_economics_20260928/normalized_calibration_v2/source_pins.json'
    timing = root / 'source/output/model/fixed_reference_economics_20260928/purchase_timing_sandbox_v1/manifest.json'
    put(normal, {'synthetic': True})
    put(timing, {'synthetic': True})
    put(root / 'inventory.json', dict(target_fingerprint='synthetic_target',
        weight_fingerprint='synthetic_weight', selected_source_sha256='synthetic_selection'))
    for arm in ('original', 'alternative'):
        top = root / 'results' / f'smoke_{arm}_chain_0'
        run = top / 'run'
        report = run / 'native_postcheck/selected_postcheck/phase_b_ge/selected_root'
        put(top / 'launcher_start.json', dict(arm=arm, chain=0, mode='smoke',
            wall_seconds=5400, start_epoch=100, deadline_epoch=5500,
            stage_inventory_sha256=sha(root / 'inventory.json')))
        put(top / 'launcher_terminal.json', dict(arm=arm, chain=0, mode='smoke', exit_code=0))
        put(run / 'start_contract.json', dict(arm=arm, chain=0, objective_calls_max=250,
            reserve_seconds=1800, target_fingerprint='synthetic_target',
            weight_fingerprint='synthetic_weight', selected_source_sha256='synthetic_selection',
            normalized_source_pins_sha256=sha(normal),
            timing_manifest_sha256=sha(timing) if arm == 'alternative' else None))
        put(run / 'input_contract.json', dict(target_fingerprint='synthetic_target',
            weight_fingerprint='synthetic_weight'))
        put(run / 'search_contract.json', {'max_objective_calls': 2})
        put(run / 'search_completed.json', {'objective_calls': 2})
        put(run / 'latest_completed.json', {'completed_full_ge': 2})
        put(run / 'best_so_far.json', {'best': {'loss': 0}})
        put(run / 'heartbeat.json', {'status': 'completed'})
        put(run / 'cases.json', [{}, {}])
        put(run / 'completed.json', dict(status='selected_numerically_verified',
            objective_calls=2, selected_postcheck=dict(status='passed',
            report='/work/results/run/native_postcheck/selected_postcheck/phase_b_ge/selected_root'),
            target_fingerprint='synthetic_target', weight_fingerprint='synthetic_weight',
            native_loss=0))
        report.mkdir(parents=True)
        with (report / 'target_fit.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=['moment', 'loss_contribution'])
            writer.writeheader()
            writer.writerows({'moment': str(i), 'loss_contribution': '0'} for i in range(14))
        with (report / 'parameters.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=['parameter', 'estimate'])
            writer.writeheader()
            writer.writerows({'parameter': str(i), 'estimate': '0'} for i in range(31))
        plots = report / 'standard_diagnostics'
        plots.mkdir()
        for i in range(17):
            (plots / f'{i}.png').write_bytes(b'synthetic')


def expect_failure(root, label):
    try:
        verify(root)
    except (AssertionError, KeyError, FileNotFoundError):
        print('synthetic rejection passed:', label)
    else:
        raise AssertionError('Gate accepted ' + label)


with tempfile.TemporaryDirectory(prefix='synthetic_soft_gate_') as directory:
    root = Path(directory)
    fixture(root)
    assert verify(root)['status'] == 'both_smokes_passed_on_current_stage'
    print('synthetic two-smoke fixture passed')
    bad = root / 'results/smoke_alternative_chain_0/run/completed.json'
    saved = bad.read_text()
    data = json.loads(saved)
    data['weight_fingerprint'] = 'wrong'
    put(bad, data)
    expect_failure(root, 'mixed fingerprint')
    bad.write_text(saved)
    plot = root / 'results/smoke_alternative_chain_0/run/native_postcheck/selected_postcheck/phase_b_ge/selected_root/standard_diagnostics/16.png'
    plot.unlink()
    expect_failure(root, 'missing standard plot')
    plot.write_bytes(b'synthetic')
    launcher = root / 'results/smoke_original_chain_0/launcher_start.json'
    changed = json.loads(launcher.read_text())
    changed['stage_inventory_sha256'] = 'wrong'
    put(launcher, changed)
    expect_failure(root, 'stage changed after smoke')
