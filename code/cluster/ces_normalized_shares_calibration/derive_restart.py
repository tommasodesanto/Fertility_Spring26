"""Derive an immutable local/remote restart by hardlinking authenticated v3 bytes.

No network, scheduler, or model actions. The caller must run verify_stage --host
on the completed stage before use; inherited source bytes are not rescanned here.
"""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

PLAN = 'output/model/experiments/ces_normalized_shares/overnight_v1/start_plan.json'
ADAPTER = 'code/model/experiments/ces_normalized_shares/adapter.py'
CALIBRATE = 'code/model/experiments/ces_normalized_shares/calibrate.py'
ENTRYPOINTS = ('launch_torch.sh', 'submit_torch.sh', 'verify_stage.py',
               'verify_smoke_gate.py', 'collect_torch.py')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(value):
    path = Path(value)
    if path.is_absolute() or not path.parts or '..' in path.parts:
        raise ValueError('Unsafe relative path: ' + str(value))
    return path


def fresh_write(path, data, mode=0o644):
    """Replace the directory entry, never write through an inherited hardlink."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(data)
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def derive(parent_deployment, out, patches, patch_inventory, psi_child=None,
           dependency_overlay=None, seed_overrides=None):
    parent = Path(parent_deployment).resolve()
    out = Path(out).resolve()
    patches = Path(patches).resolve()
    if out.exists() or out == parent or out.is_relative_to(parent):
        raise ValueError('Output must be a fresh deployment outside parent')
    if seed_overrides is not None and psi_child is not None:
        raise ValueError('Specify seed_overrides or psi_child, not both')
    overrides = seed_overrides if seed_overrides is not None else ({'psi_child': psi_child} if psi_child is not None else None)
    if not isinstance(overrides, dict) or not overrides:
        raise ValueError('Explicit nonempty seed override mapping is required')
    prior = json.loads((parent / 'inventory.json').read_text())
    parent_sha = sha(parent / 'inventory.json')
    overlay = (Path(dependency_overlay).resolve() if dependency_overlay else
               parent / 'dependency_overlay')
    dependency = json.loads((overlay / 'inventory.json').read_text())
    if dependency['source_v3_inventory_sha256'] != parent_sha:
        raise ValueError('Dependency overlay does not authenticate this parent inventory')
    artifacts = dependency['files']
    if len(artifacts) != 3 or {Path(name).name for name in artifacts} != {
            'initial_state.pkl.gz', 'receipt.json', 'parameters.csv'} or len({
            str(Path(name).parent) for name in artifacts}) != 1:
        raise ValueError('Expected precisely the three missing reference-case artifacts')
    plan = json.loads((parent / 'source' / PLAN).read_text())
    if sha(parent / 'source' / PLAN) != prior['files'][PLAN]:
        raise ValueError('Parent plan digest mismatch')
    if prior['start_plan_sha256'] != prior['files'][PLAN]:
        raise ValueError('Parent plan metadata mismatch')
    allowed = set(plan['coordinates']) - {'delta_alpha_jump', 'delta_alpha'}
    if not set(overrides) <= allowed:
        raise ValueError('Seed overrides may name only inherited free coordinates')
    for name, value in overrides.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError('Seed override must be a finite number: ' + name)
        lower, upper = plan['bounds'][name]
        if not lower <= value <= upper:
            raise ValueError('Seed override outside inherited bounds: ' + name)
    if len(plan['starts']) != 4 or len(plan['coordinates']) != 11:
        raise ValueError('Parent start/coordinate contract mismatch')
    if plan['target_fingerprint'] != prior['target_fingerprint'] or plan['weight_fingerprint'] != prior['weight_fingerprint']:
        raise ValueError('Parent target metadata mismatch')
    expected = json.loads(Path(patch_inventory).read_text())
    if 'files' in expected:
        expected = expected['files']
    patch_paths = {ADAPTER: patches / 'adapter.py', CALIBRATE: patches / 'calibrate.py',
                   **{n: patches / n for n in ENTRYPOINTS}}
    if set(expected) != set(patch_paths):
        raise ValueError('Patch inventory must name adapter and calibrate relative paths and exactly five entrypoints')
    for name, path in patch_paths.items():
        if sha(path) != expected[name]:
            raise ValueError('Lead-reviewed patch digest mismatch: ' + name)
    for name, record in artifacts.items():
        safe_relative(name)
        if name in prior['files']:
            raise ValueError('Dependency artifact must be absent from parent: ' + name)
        if (overlay / name).is_symlink():
            raise ValueError('Dependency must be a regular immutable file: ' + name)
        if sha(overlay / name) != record['sha256'] or (overlay / name).stat().st_size != record['bytes']:
            raise ValueError('Dependency artifact digest/size mismatch: ' + name)
    checkpoint = next(name for name in artifacts if name.endswith('/initial_state.pkl.gz'))
    receipt_name = next(name for name in artifacts if name.endswith('/receipt.json'))
    if json.loads((overlay / receipt_name).read_text())['case_checkpoint_sha256'] != artifacts[checkpoint]['sha256']:
        raise ValueError('Dependency receipt/checkpoint mismatch')
    updated = copy.deepcopy(plan)
    for start in updated['starts']:
        if set(start) != set(plan['coordinates']):
            raise ValueError('Unexpected start coordinate set')
        start.update(overrides)
    updated['status'] = 'prepared_not_launched'
    updated['seed_adjustment'] = dict(overrides=overrides,
        previous_values=[{name: start[name] for name in overrides} for start in plan['starts']],
        parent_inventory_sha256=parent_sha,
        reason='The failed v3 initial candidate found no fertility root within the approved price caps; explicit caller-selected diagnostic coordinate overrides, not a feasibility certification or economic adoption',
        no_adoption=True)
    # Every inherited source entry is linked independently; do not copy caches or results.
    out.mkdir(parents=True)
    source = out / 'source'
    source.mkdir()
    files = dict(prior['files'])
    for name in files:
        relative = safe_relative(name)
        origin = parent / 'source' / relative
        if origin.is_symlink() or not origin.is_file():
            raise ValueError('Parent source must be a regular immutable file: ' + name)
        destination = source / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(origin, destination)
    fresh_write(source / PLAN, encoded(updated))
    for name in (ADAPTER, CALIBRATE):
        if name not in prior['files']:
            raise ValueError('Reviewed replacement missing from parent inventory: ' + name)
        fresh_write(source / name, patch_paths[name].read_bytes())
    for name, record in artifacts.items():
        destination = source / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(overlay / name, destination)
        files[name] = record['sha256']
    files[PLAN] = sha(source / PLAN)
    for name in (ADAPTER, CALIBRATE):
        files[name] = sha(source / name)
    entries = {}
    old_remote = prior['remote_root']
    for name in ENTRYPOINTS:
        data = patch_paths[name].read_bytes().replace(old_remote.encode(), str(out).encode())
        fresh_write(out / name, data, 0o755 if name.endswith('.sh') else 0o644)
        entries[name] = sha(out / name)
    inventory = copy.deepcopy(prior)
    inventory.update(files=files, entrypoints=entries, remote_root=str(out),
        source_prefix='source/', parent_inventory_sha256=parent_sha,
        start_plan_sha256=files[PLAN])
    # Preserve archive history but make no claim that this derived stage has an archive.
    inventory['derivation_method'] = 'authenticated_parent_hardlinks_plus_reviewed_replacements'
    deps = inventory.setdefault('authenticated_dependencies', {})
    deps['declared_artifact_count'] = deps.get('declared_artifact_count', 0) + 3
    deps['unresolved_pointers'] = []
    deps['dependency_overlay_inventory_sha256'] = sha(overlay / 'inventory.json')
    whitelist = sorted([PLAN, ADAPTER, CALIBRATE, *artifacts])
    changed = sorted(name for name, digest in files.items() if prior['files'].get(name) != digest)
    if changed != whitelist:
        raise ValueError('Source diff differs from exact permitted whitelist')
    fresh_write(out / 'inventory.json', encoded(inventory))
    receipt = dict(status='prepared_no_submission', parent_deployment=str(parent),
        remote_root=str(out), parent_inventory_sha256=parent_sha,
        inventory_sha256=sha(out / 'inventory.json'), patch_inventory_sha256=sha(patch_inventory),
        expected_source_diff=whitelist, actual_source_diff=changed,
        inherited_source_validation='Parent was authenticated by preflight; hardlinks preserve inode bytes; full derived verify_stage required before execution',
        source_files=len(files), target_fingerprint=updated['target_fingerprint'],
        weight_fingerprint=updated['weight_fingerprint'], seed_overrides=overrides,
        archive_created=False, required_next_check='verify_stage.py --host')
    fresh_write(out / 'derivation_receipt.json', encoded(receipt))
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-deployment', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--patches', required=True)
    parser.add_argument('--patch-inventory', required=True)
    seeds = parser.add_mutually_exclusive_group(required=True)
    seeds.add_argument('--seed-overrides', type=json.loads,
        help='JSON object mapping inherited free coordinates to explicit diagnostic values')
    seeds.add_argument('--psi-child', type=float, help='Legacy single-coordinate override')
    parser.add_argument('--dependency-overlay')
    args = parser.parse_args()
    print(json.dumps(derive(args.parent_deployment, args.out, args.patches,
        args.patch_inventory, args.psi_child, args.dependency_overlay, args.seed_overrides), sort_keys=True))


if __name__ == '__main__':
    main()
