"""Stage the collector's two postchecked winners for the dated-policy jobs."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
READOUT = PACKET / 'collection/readout'
CALIB_REMOTE = '/scratch/td2248/projects/purchase_rules_overnight_v1'
MECH_REMOTE = '/scratch/td2248/projects/purchase_mechanism_v1'


def run(*args):
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def remote_sha(path):
    return run('ssh', '-o', 'BatchMode=yes', 'torch', f'sha256sum {path}').split()[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', help='Upload selected files after validation')
    args = parser.parse_args()
    plan = json.loads((PACKET / 'plan.json').read_text())
    manifest = dict(schema='purchase_policy_selected_snapshot_v1',
                    target_fingerprint=plan['target_fingerprint'],
                    weight_fingerprint=plan['weight_fingerprint'], arms={})
    staged = HERE / 'selection_snapshot'
    staged.mkdir(parents=True, exist_ok=True)
    for arm in ('hard', 'quarter'):
        source = READOUT / f'selected_{arm}.json'
        selected = json.loads(source.read_text())
        if selected['status'] != 'postchecked' or selected['arm'] != arm:
            raise RuntimeError('No authenticated selected ' + arm + ' candidate')
        if selected['target_fingerprint'] != plan['target_fingerprint'] or selected['weight_fingerprint'] != plan['weight_fingerprint']:
            raise RuntimeError('Mixed target/weight fingerprints')
        chain = int(selected['chain'])
        if selected['origin'] == 'local':
            local = PACKET / f'local_runtime/runs/local10_v1/chain{chain}/postcheck'
            if not local.is_dir():
                local = PACKET / f'local_runtime/runs/local10_v1/chain_{chain}/postcheck'
            completed = local / 'completed.json'
            if not completed.is_file() or json.loads(completed.read_text())['status'] != 'selected_numerically_verified':
                raise RuntimeError('Local postcheck not complete')
            if args.apply:
                destination = f'{CALIB_REMOTE}/results/chain_{chain}'
                run('ssh', '-o', 'BatchMode=yes', 'torch', f'mkdir -p {destination}/postcheck')
                run('rsync', '-az', str(local) + '/', f'torch:{destination}/postcheck/')
            digest = sha(completed)
        elif selected['origin'] == 'torch':
            digest = remote_sha(f'{CALIB_REMOTE}/results/chain_{chain}/postcheck/completed.json')
        else:
            raise RuntimeError('Unknown selected origin')
        if args.apply:
            remote_completed = f'{CALIB_REMOTE}/results/chain_{chain}/postcheck/completed.json'
            if remote_sha(remote_completed) != digest:
                raise RuntimeError('Uploaded postcheck hash differs')
        copy = dict(selected)
        if selected['origin'] == 'local':
            copy['origin'] = 'local_uploaded_to_torch'
            copy['remote_root'] = f'{CALIB_REMOTE}/results/chain_{chain}'
            copy['remote_report'] = copy['remote_root'] + '/postcheck/selected_postcheck/phase_b_ge/selected_root'
            copy['remote_arrays'] = copy['remote_root'] + '/postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz'
        target = staged / f'selected_{arm}.json'
        target.write_text(json.dumps(copy, indent=2, sort_keys=True) + '\n')
        manifest['arms'][arm] = dict(chain=chain, selected_json_sha256=sha(target),
                                    completed_sha256=digest, original_origin=selected['origin'])
    target = staged / 'manifest.json'
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    if args.apply:
        run('scp', '-q', *(str(staged / name) for name in ('selected_hard.json', 'selected_quarter.json', 'manifest.json')),
            f'torch:{MECH_REMOTE}/selection/')
        if remote_sha(f'{MECH_REMOTE}/selection/manifest.json') != sha(target):
            raise RuntimeError('Remote selection manifest hash differs')
    print(json.dumps(dict(status='uploaded' if args.apply else 'validated_local_selection',
                          manifest=str(target), arms=manifest['arms']), sort_keys=True))


if __name__ == '__main__':
    main()
