"""Stage the collector's two postchecked winners for the dated-policy jobs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shlex
import subprocess
from pathlib import Path, PurePosixPath

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
READOUT = PACKET / 'collection/readout'
MECH_REMOTE = '/scratch/td2248/projects/purchase_mechanism_v1'
ORIGINS = {
    'torch': '/scratch/td2248/projects/purchase_rules_overnight_v1/results',
    'torch_restart': '/scratch/td2248/projects/purchase_restart_controller_v2/results',
    'local': str(PACKET / 'local_runtime/runs/local10_v1'),
    'local_restart': str(PACKET / 'local_runtime/restart_v2/runs'),
}
REPORT = 'selected_postcheck/phase_b_ge/selected_root'
ARRAY = 'selected_postcheck/phase_b_ge/selected_repeat/stage/solution_arrays.npz'


def run(*args, input=None):
    return subprocess.run(args, check=True, text=True, input=input, capture_output=True).stdout


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def files_local(folder):
    root = Path(folder)
    return {str(p.relative_to(root)): sha(p) for p in sorted(root.rglob('*')) if p.is_file()}


REMOTE_FILES = '''import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]); out={}
for p in sorted(root.rglob('*')):
 if p.is_file():
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1048576),b''): h.update(b)
  out[str(p.relative_to(root))]=h.hexdigest()
print(json.dumps(out,sort_keys=True))
'''


def files_remote(folder):
    cmd = '/share/apps/anaconda3/2025.06/bin/python - ' + shlex.quote(str(folder))
    return json.loads(run('ssh', '-o', 'BatchMode=yes', 'torch', cmd, input=REMOTE_FILES))


def read_source(origin, path):
    if origin.startswith('torch'):
        return run('ssh', '-o', 'BatchMode=yes', 'torch', 'cat ' + shlex.quote(path))
    return Path(path).read_text()


def physical_source(selected):
    origin = selected.get('origin')
    if origin not in ORIGINS:
        raise RuntimeError('Unknown selected origin')
    chain = selected.get('chain')
    if type(chain) is not int or chain < 0:
        raise RuntimeError('Invalid original chain ID')
    raw = selected.get('remote_root')
    if not isinstance(raw, str) or not raw.startswith('/') or '..' in PurePosixPath(raw).parts:
        raise RuntimeError('Missing or unsafe physical remote_root')
    root = PurePosixPath(raw)
    if root.parent != PurePosixPath(ORIGINS[origin]) or re.fullmatch(r'chain_?' + str(chain), root.name) is None:
        raise RuntimeError('Selected physical root differs from origin and chain')
    if selected.get('remote_report') != str(root / 'postcheck' / REPORT) or selected.get('remote_arrays') != str(root / 'postcheck' / ARRAY):
        raise RuntimeError('Collector report or arrays path differs from physical source')
    if origin.endswith('restart') and not selected.get('parent_remote_root'):
        raise RuntimeError('Restart parent provenance missing')
    return origin, str(root)


def publish_snapshot(origin, root, chain, files, identity):
    base = MECH_REMOTE + '/selected_postchecks'
    final = f'{base}/chain_{chain}'
    exists = run('ssh', '-o', 'BatchMode=yes', 'torch',
                 f'test -e {shlex.quote(final)} && echo yes || echo no').strip() == 'yes'
    if exists:
        receipt = json.loads(run('ssh', '-o', 'BatchMode=yes', 'torch',
                                 'cat ' + shlex.quote(final + '/source_identity.json')))
        if receipt != identity or files_remote(final + '/postcheck') != files:
            raise RuntimeError('Immutable snapshot differs; use a new named version')
        return
    temp = f'{base}/.chain_{chain}_{os.getpid()}'
    run('ssh', '-o', 'BatchMode=yes', 'torch', f'mkdir -p {shlex.quote(base)} && mkdir {shlex.quote(temp)}')
    try:
        if origin.startswith('torch'):
            run('ssh', '-o', 'BatchMode=yes', 'torch',
                f'cp -a {shlex.quote(root + "/postcheck")} {shlex.quote(temp + "/postcheck")}')
        else:
            run('rsync', '-a', root + '/postcheck/', f'torch:{temp}/postcheck/')
        if files_remote(temp + '/postcheck') != files:
            raise RuntimeError('Copied snapshot file hashes differ')
        receipt = json.dumps(identity, indent=2, sort_keys=True) + '\n'
        writer = 'import sys\nopen(sys.argv[1],"w").write(' + repr(receipt) + ')\n'
        run('ssh', '-o', 'BatchMode=yes', 'torch',
            '/share/apps/anaconda3/2025.06/bin/python - ' + shlex.quote(temp + '/source_identity.json'), input=writer)
        run('ssh', '-o', 'BatchMode=yes', 'torch',
            f'test ! -e {shlex.quote(final)} && mv -T {shlex.quote(temp)} {shlex.quote(final)}')
    finally:
        run('ssh', '-o', 'BatchMode=yes', 'torch', 'rm -rf ' + shlex.quote(temp))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', help='Upload selected files after validation')
    args = parser.parse_args()
    plan = json.loads((PACKET / 'plan.json').read_text())
    manifest = dict(schema='purchase_policy_selected_snapshot_v2',
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
        chain = selected['chain']
        origin, root = physical_source(selected)
        postcheck = root + '/postcheck'
        files = files_remote(postcheck) if origin.startswith('torch') else files_local(postcheck)
        for required in ('completed.json', 'input_contract.json', ARRAY):
            if required not in files:
                raise RuntimeError('Missing selected postcheck file: ' + required)
        if set(selected['report_sha256']) != {'target_fit.csv', 'parameters.csv', 'closure.json'}:
            raise RuntimeError('Incomplete collector report hash set')
        for name, expected in selected['report_sha256'].items():
            if files.get(REPORT + '/' + name) != expected:
                raise RuntimeError('Collector report hash differs: ' + name)
        report = postcheck + '/' + REPORT
        target_rows = list(csv.DictReader(io.StringIO(read_source(origin, report + '/target_fit.csv'))))
        parameter_rows = list(csv.DictReader(io.StringIO(read_source(origin, report + '/parameters.csv'))))
        closure = json.loads(read_source(origin, report + '/closure.json'))
        if (len(target_rows) != 14 or target_rows != selected['target_fit'] or
                len(parameter_rows) != 31 or parameter_rows != selected['parameters'] or
                closure != selected['closure'] or
                len([name for name in files if name.startswith(REPORT + '/standard_diagnostics/') and name.endswith('.png')]) != 17):
            raise RuntimeError('Full native target, parameter, closure or 17-plot identity differs')
        completed_text = read_source(origin, postcheck + '/completed.json')
        if json.loads(completed_text).get('status') != 'selected_numerically_verified':
            raise RuntimeError('Selected postcheck is not numerically verified')
        identity = dict(original_selected_json_sha256=sha(source), origin=origin, remote_root=root,
                        parent_remote_root=selected.get('parent_remote_root'), chain=chain)
        snapshot = f'{MECH_REMOTE}/selected_postchecks/chain_{chain}'
        if args.apply:
            publish_snapshot(origin, root, chain, files, identity)
        copy = dict(selected)
        copy['snapshot_remote_root'] = snapshot
        target = staged / f'selected_{arm}.json'
        target.write_text(json.dumps(copy, indent=2, sort_keys=True) + '\n')
        manifest['arms'][arm] = dict(chain=chain, origin=origin, physical_remote_root=root,
            parent_remote_root=selected.get('parent_remote_root'), snapshot_remote_root=snapshot,
            selected_json_sha256=sha(target), original_selected_json_sha256=sha(source),
            completed_sha256=files['completed.json'],
            report_sha256={name: files[REPORT + '/' + name] for name in ('target_fit.csv','parameters.csv','closure.json')},
            native_arrays_sha256=files[ARRAY], postcheck_files=files)
    target = staged / 'manifest.json'
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    if args.apply:
        existing = files_remote(f'{MECH_REMOTE}/selection').get('manifest.json')
        if existing is not None and existing != sha(target):
            raise RuntimeError('Published selection differs; use a new named version')
        run('scp', '-q', *(str(staged / name) for name in ('selected_hard.json', 'selected_quarter.json', 'manifest.json')),
            f'torch:{MECH_REMOTE}/selection/')
        if files_remote(f'{MECH_REMOTE}/selection').get('manifest.json') != sha(target):
            raise RuntimeError('Remote selection manifest hash differs')
    print(json.dumps(dict(status='uploaded' if args.apply else 'validated_local_selection',
                          manifest=str(target), arms={arm: {k:v for k,v in row.items() if k != 'postcheck_files'}
                                                     for arm,row in manifest['arms'].items()}), sort_keys=True))


if __name__ == '__main__':
    main()
