#!/usr/bin/env python3
"""Build and verify the isolated v9 diagnostic-output-routing package.

This utility deliberately has no submit action.  It authenticates the v7
inventory, hardlinks its immutable source tree, detaches exactly one source
file, applies the reviewed routing delta, and refreshes the handoff, plans,
and panel pins that explicitly contain that source identity.
"""
from __future__ import annotations

import argparse
import copy
import difflib
import gzip
import hashlib
import io
import json
import os
import shutil
import tarfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
TASK = ROOT / 'output/model/transition_readiness_v1/current_baseline_20261003'
V7 = TASK / 'deployment_v7'
V9 = TASK / 'deployment_v9'
PLANS7 = TASK / 'plans_v7'
PLANS5 = TASK / 'plans_v5'
PLANS9 = TASK / 'plans_v9'
REMOTE = '/scratch/td2248/projects/current_estate_transition_20261003_v9'
REMOTE_PARENT = '/scratch/td2248/projects/current_estate_transition_20261003_v7'
V7_SHA = '5f1ee713bd1c5ab62d2a1776667f6191d61e69794cfe9a1fd22e6a6841397e92'
TARGET = 'code/model/experiments/transition_readiness/floor_runtime.py'
OLD = "        folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);P=copy.deepcopy(self.P);P.psi_child=float(psi)\n"
ADD = """        saved_output=getattr(P,'native_inherited_distribution_evidence_dir',None)
        P.native_inherited_distribution_evidence_dir=str(folder/'inherited_state_evidence')
        write(folder/'output_override.json',dict(field='native_inherited_distribution_evidence_dir',
            saved=saved_output,effective=P.native_inherited_distribution_evidence_dir,economic_change=False))
"""


def sha_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def load(path: Path):
    return json.loads(path.read_text())


def dump(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def v7_inventory():
    inv = load(V7 / 'inventory.json')
    require(sha(V7 / 'inventory.json') == V7_SHA, 'Authenticated v7 inventory SHA differs')
    for rel, digest in inv['files'].items():
        require(sha(V7 / 'source' / rel) == digest, 'v7 source drift: ' + rel)
    return inv


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def patch_source(before: bytes) -> bytes:
    text = before.decode()
    require(text.count(OLD) == 1, 'Reviewed stationary prefix is not unique')
    require(ADD not in text, 'Routing patch unexpectedly already present')
    return text.replace(OLD, OLD + ADD).encode()


def plan_pin(path: Path) -> dict:
    return {'path': str(path), 'sha256': sha(path)}


def verify_reference_handoff(identity: dict, handoff_path: Path) -> None:
    require(identity.get('reference_sha256') == sha(handoff_path),
            'Reference identity does not match selected handoff SHA-256')


def replace_identity(value, new_identity):
    """Replace only identity dictionaries that carry the source-pin vector."""
    if isinstance(value, dict):
        return {key: (copy.deepcopy(new_identity) if key == 'identity' else replace_identity(item, new_identity))
                for key, item in value.items()}
    if isinstance(value, list):
        return [replace_identity(item, new_identity) for item in value]
    return value


def build(_: argparse.Namespace) -> dict:
    inv7 = v7_inventory()
    require(not V9.exists() or not any(V9.iterdir()), 'Refusing nonempty deployment_v9')
    require(not PLANS9.exists() or not any(PLANS9.iterdir()), 'Refusing nonempty plans_v9')
    V9.mkdir(parents=True, exist_ok=True)
    PLANS9.mkdir(parents=True, exist_ok=True)
    shutil.copytree(V7 / 'source', V9 / 'source', copy_function=os.link)
    source = V9 / 'source' / TARGET
    before = source.read_bytes()
    # copytree used hardlinks; unlink before the sole source write.
    source.unlink()
    source.write_bytes(patch_source(before))
    new_hash = sha(source)
    require(sha(V7 / 'source' / TARGET) == inv7['files'][TARGET], 'Detached write changed v7 source')

    handoff = load(PLANS5 / 'handoff.json')
    require(handoff['source_pins'][TARGET] == inv7['files'][TARGET], 'v5 handoff source pin differs from v7')
    handoff['source_pins'][TARGET] = new_hash
    handoff_path = PLANS9 / 'handoff.json'
    dump(handoff_path, handoff)
    identity = None
    fit_old = load(PLANS7 / 'fit_plan.json')
    identity = copy.deepcopy(fit_old['identity'])
    require(identity['source_pins'][TARGET] == inv7['files'][TARGET], 'v7 fit source pin differs')
    identity['source_pins'][TARGET] = new_hash
    identity['reference_sha256'] = sha(handoff_path)
    verify_reference_handoff(identity, handoff_path)

    smoke = replace_identity(load(PLANS5 / 'smoke_plan.json'), identity)
    smoke['handoff'] = plan_pin(handoff_path)
    smoke_path = PLANS9 / 'smoke_plan.json'
    dump(smoke_path, smoke)
    fit = replace_identity(fit_old, identity)
    fit['handoff'] = plan_pin(handoff_path)
    fit_path = PLANS9 / 'fit_plan.json'
    dump(fit_path, fit)
    panel = replace_identity(load(PLANS7 / 'panel_config.json'), identity)
    panel['plan'] = plan_pin(fit_path)
    panel_path = PLANS9 / 'panel_config.json'
    dump(panel_path, panel)
    # Add the refreshed plan objects to the staged source set; all are genuine
    # runtime inputs, including the copied handoff required by validate_handoff.
    for path in (handoff_path, smoke_path, fit_path, panel_path):
        target = V9 / 'source' / rel(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())

    inv = copy.deepcopy(inv7)
    inv['remote_root'] = REMOTE
    inv['files'][TARGET] = new_hash
    for path in (handoff_path, smoke_path, fit_path, panel_path):
        inv['files'][rel(path)] = sha(V9 / 'source' / rel(path))
    inv['identity'] = identity
    inv['plans'] = {'smoke': {'path': rel(smoke_path), 'sha256': sha(smoke_path)},
                    'fit': {'path': rel(fit_path), 'sha256': sha(fit_path)}}
    inv['panel_config'] = {'path': rel(panel_path), 'sha256': sha(panel_path)}
    inv['output_routing_delta'] = {'source': TARGET, 'old_sha256': inv7['files'][TARGET],
        'new_sha256': new_hash, 'economic_change': False, 'requires_new_native_smoke': True,
        'parent_inventory_sha256': V7_SHA}
    # Stage-local entrypoints are only remote-root/output-name retargets.
    entrypoints = {}
    for name in ('deploy.py', 'launch_torch.sh', 'panel_launch_torch_v7.sh'):
        data = (V7 / name).read_bytes()
        prior = b'current_estate_transition_20261003_v5' if name == 'launch_torch.sh' else b'current_estate_transition_20261003_v7'
        require(data.count(prior) >= 1, 'Expected inherited remote root absent: ' + name)
        data = data.replace(prior, b'current_estate_transition_20261003_v9')
        out = V9 / name
        out.write_bytes(data)
        if name.endswith('.sh'):
            out.chmod(0o755)
        entrypoints[name] = sha(out)
    inv['entrypoints'] = entrypoints
    dump(V9 / 'inventory.json', inv)
    entries = {'inventory.json': (V9 / 'inventory.json').read_bytes(),
               'source/' + TARGET: source.read_bytes()}
    for path in (handoff_path, smoke_path, fit_path, panel_path):
        entries['source/' + rel(path)] = (V9 / 'source' / rel(path)).read_bytes()
    for name in entrypoints:
        entries[name] = (V9 / name).read_bytes()
    with (V9 / 'stage.tar.gz').open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode='w') as archive:
            for name, data in sorted(entries.items()):
                info = tarfile.TarInfo(name); info.size = len(data); info.mtime = 0
                info.mode = 0o755 if name.endswith('.sh') else 0o644
                archive.addfile(info, io.BytesIO(data))
    receipt = verify(argparse.Namespace())
    receipt.update(status='prepared_local_no_submission', archive_sha256=sha(V9 / 'stage.tar.gz'),
                   inventory_sha256=sha(V9 / 'inventory.json'), remote_parent=REMOTE_PARENT,
                   remote_parent_inventory_sha256=V7_SHA)
    dump(V9 / 'stage_receipt.json', receipt)
    verify(argparse.Namespace())
    return receipt


def verify(_: argparse.Namespace) -> dict:
    inv7 = v7_inventory(); inv = load(V9 / 'inventory.json')
    require(inv['remote_root'] == REMOTE, 'v9 remote root differs')
    require(inv['output_routing_delta']['parent_inventory_sha256'] == V7_SHA, 'Parent inventory pin differs')
    for key, path in (('smoke', PLANS9 / 'smoke_plan.json'), ('fit', PLANS9 / 'fit_plan.json')):
        require(inv['plans'][key] == {'path': rel(path), 'sha256': sha(path)}, 'Plan pin differs: ' + key)
    require(inv['panel_config'] == {'path': rel(PLANS9 / 'panel_config.json'), 'sha256': sha(PLANS9 / 'panel_config.json')}, 'Panel pin differs')
    for file_rel, digest in inv['files'].items():
        require(sha(V9 / 'source' / file_rel) == digest, 'v9 source drift: ' + file_rel)
    changed = {key for key in set(inv['files']) | set(inv7['files']) if inv['files'].get(key) != inv7['files'].get(key)}
    expected = {TARGET, *(rel(p) for p in (PLANS9 / 'handoff.json', PLANS9 / 'smoke_plan.json', PLANS9 / 'fit_plan.json', PLANS9 / 'panel_config.json'))}
    require(changed == expected, 'Unexpected source/data delta: ' + repr(sorted(changed ^ expected)))
    require(inv['identity']['source_pins'][TARGET] == sha(V9 / 'source' / TARGET), 'Inventory identity missing patched hash')
    for path in (PLANS9 / 'handoff.json', PLANS9 / 'smoke_plan.json', PLANS9 / 'fit_plan.json', PLANS9 / 'panel_config.json'):
        data = load(path)
        require(data['identity']['source_pins'][TARGET] == sha(V9 / 'source' / TARGET) if 'identity' in data else data['source_pins'][TARGET] == sha(V9 / 'source' / TARGET), 'Transitive pin differs: ' + str(path))
    handoff_path = PLANS9 / 'handoff.json'
    require(inv['identity']['reference_sha256'] == sha(handoff_path), 'Inventory identity handoff SHA differs')
    verify_reference_handoff(inv['identity'], handoff_path)
    fit = load(PLANS9 / 'fit_plan.json'); smoke = load(PLANS9 / 'smoke_plan.json'); panel = load(PLANS9 / 'panel_config.json')
    for data, name in ((fit, 'fit'), (smoke, 'smoke'), (panel, 'panel')):
        verify_reference_handoff(data['identity'], handoff_path)
        require(data['identity'] == inv['identity'], 'Plan identity differs from inventory: ' + name)
    require(fit['handoff'] == plan_pin(handoff_path) and smoke['handoff'] == plan_pin(handoff_path), 'Handoff pin differs')
    require(panel['plan'] == plan_pin(PLANS9 / 'fit_plan.json'), 'Panel plan pin differs')
    # Assert numerical controls did not move; only source identity/handoff path changes are allowed.
    old_fit = load(PLANS7 / 'fit_plan.json')
    for key in ('budget', 'fit', 'path', 'endpoint', 'initial_psi', 'psi_bound_ratios', 'target_contract', 'gates', 'horizons'):
        require(fit[key] == old_fit[key], 'Economic/numerical plan drift: ' + key)
    receipt_path = V9 / 'stage_receipt.json'
    if receipt_path.exists():
        receipt = load(receipt_path)
        require(receipt['remote_parent'] == REMOTE_PARENT and receipt['remote_parent_inventory_sha256'] == V7_SHA,
                'Stage receipt parent pin differs')
    require(patch_source((V7 / 'source' / TARGET).read_bytes()) == (V9 / 'source' / TARGET).read_bytes(),
            'Patched source differs from the exact reviewed four-line delta')
    return {'verified_parent_inventory_sha256': V7_SHA, 'source_files': len(inv['files']),
            'changed_source_files': sorted(changed), 'unchanged_v7_source_files': len(inv7['files']) - 1,
            'reference_handoff_sha256': sha(handoff_path),
            'new_native_smoke_required': True, 'no_submission_performed': True}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'verify'))
    args = parser.parse_args()
    print(json.dumps(build(args) if args.action == 'build' else verify(args), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
