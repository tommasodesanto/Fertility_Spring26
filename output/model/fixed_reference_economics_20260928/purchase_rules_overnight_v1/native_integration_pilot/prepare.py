"""Freeze provisional verified checkpoints for the two one-date native pilots."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

PACKET = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
SCAN = PACKET / 'collection/scan_remote.py'
TORCH_PYTHON = '/share/apps/anaconda3/2025.06/bin/python -'


def run(command, *, input=None, env=None):
    return subprocess.run(command, input=input, env=env, text=True, capture_output=True, check=True).stdout


def read_remote(extra=''):
    return json.loads(run(['ssh', '-o', 'BatchMode=yes', 'torch',
                           (extra + ' ' + TORCH_PYTHON).strip()], input=SCAN.read_text()))


def read_local(root, kind):
    env = dict(os.environ, PURCHASE_PACKET_ROOT=str(PACKET), PURCHASE_RESULTS_ROOT=str(root),
               PURCHASE_SOURCE_KIND=kind, PURCHASE_CHAIN_FIRST='48', PURCHASE_CHAIN_LAST='58')
    if kind == 'local_restart':
        env['PURCHASE_PARENT_RESULTS_ROOT'] = str(PACKET / 'local_runtime/runs/local10_v1')
    return json.loads(run([sys.executable, str(SCAN)], env=env))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if HERE.joinpath('manifest.json').exists():
        raise RuntimeError('Pilot snapshot already exists; never replace a selected checkpoint')
    sources = [
        ('torch', read_remote()),
        ('torch_restart', read_remote(
            'env PURCHASE_SOURCE_KIND=restart '
            'PURCHASE_RESULTS_ROOT=/scratch/td2248/projects/purchase_restart_controller_v2/results '
            'PURCHASE_PARENT_RESULTS_ROOT=/scratch/td2248/projects/purchase_rules_overnight_v1/results')),
        ('local', read_local(PACKET / 'local_runtime/runs/local10_v1', 'original')),
        ('local_restart', read_local(PACKET / 'local_runtime/restart_v2/runs', 'local_restart')),
    ]
    plan = json.loads((PACKET / 'plan.json').read_text())
    rows = []
    for origin, snapshot in sources:
        if (snapshot['errors'] or snapshot['target_fingerprint'] != plan['target_fingerprint']
                or snapshot['weight_fingerprint'] != plan['weight_fingerprint']):
            raise RuntimeError('Invalid candidate snapshot: ' + origin)
        rows.extend(dict(row, origin=origin) for row in snapshot['chains']
                    if row['status'] == 'postchecked')
    manifest = dict(schema='provisional_native_integration_pilot_v1', created_epoch=time.time(),
                    final_calibration=False, production_policy=False,
                    target_fingerprint=plan['target_fingerprint'],
                    weight_fingerprint=plan['weight_fingerprint'], arms={})
    for arm in ('hard', 'quarter'):
        candidate = min((row for row in rows if row['arm'] == arm),
                        key=lambda row: (row['loss'], row['chain'], row['source_run']))
        selected = HERE / ('selected_' + arm + '.json')
        selected.write_text(json.dumps(candidate, indent=2, sort_keys=True) + '\n')
        completed = Path(candidate['remote_root']) / 'postcheck/completed.json'
        native_hash = (digest(completed) if candidate['origin'].startswith('local') else
                       run(['ssh', '-o', 'BatchMode=yes', 'torch',
                            'sha256sum ' + str(completed)]).split()[0])
        manifest['arms'][arm] = dict(chain=candidate['chain'], origin=candidate['origin'],
            source_run=candidate['source_run'], provisional_loss=candidate['loss'],
            physical_root=candidate['remote_root'], selected_sha256=digest(selected),
            completed_sha256=native_hash,
            local_upload_required=candidate['origin'].startswith('local'))
    (HERE / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps({arm: dict(chain=row['chain'], origin=row['origin'], loss=row['provisional_loss'])
                      for arm, row in manifest['arms'].items()}, sort_keys=True))


if __name__ == '__main__':
    main()
