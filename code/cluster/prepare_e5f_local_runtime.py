#!/usr/bin/env python3
"""Relocate authenticated E5F Torch inputs without changing model-source bytes.

The original downloaded tree is read-only input. Only the two ancestor runtime
BASE literals, their manifest/lock pins, and the recovery runner's corresponding
PARENT_LOCK literal change. Every change is recorded; original evidence remains.
This prepares files and optionally loads a seed, never solves or launches jobs.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import sys

TORCH_ROOT = '/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a'
SNAPSHOT = 'calibration_code_integration_20260927_v2'
PAIR = 'nightpair_20260925_v1'
PARENT_LOCK = '6443195fa3f7de0dce5cc8a4c05e2709b99d586421c96a35dba3c07e93e061a1'
CONTRACT_HASH = '399abb6e9eab0d447dca627f94e3de4e6a8006d8920d02241fadac21f6c6ebae'
SEED_BASE = 'utility_overnight_20260923_v1/results/production/B_floor/worker09_proposal16'
SEED = SEED_BASE + '/result/evaluation/raw/repetition_01/initial_state.pkl.gz'
SEED_SHA = '83a28e46b36e2fbe30338d366611f3ec209f0c5a68309ee4ee9fa8523b66adee'
PLAN_SHA = '9cad55d0b8ec186b95d3cb4e9f866add69ba78d5e894b0cdf780ce8191f1889c'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + '.relocation.tmp')
    tmp.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


class CompatibleUnpickler(pickle.Unpickler):
    """Only namespace aliases used by Python 3.13 / NumPy 2 seed pickles."""
    def find_class(self, module, name):
        if module == 'pathlib._local':
            module = 'pathlib'
        if module.startswith('numpy._core'):
            module = module.replace('numpy._core', 'numpy.core', 1)
        return super().find_class(module, name)


def load_checkpoint(path):
    """Read the exact authenticated seed; never convert arrays or parameters."""
    if sha(path) != SEED_SHA:
        raise RuntimeError('Portable checkpoint differs from the pinned Torch seed')
    with gzip.open(path, 'rb') as stream:
        return CompatibleUnpickler(stream).load()


def verify_original(root):
    contract_path = root / SNAPSHOT / 'launch_v1/contract.json'
    assert sha(contract_path) == CONTRACT_HASH, 'Original v2 contract changed'
    c = read(contract_path)
    def original_path(path):
        assert str(path).startswith(TORCH_ROOT + '/')
        return root / str(path)[len(TORCH_ROOT) + 1:]
    for item in list(c['files'].values()) + [c['base_contract'], c['objective']]:
        assert sha(original_path(item['path'])) == item['sha256'], item['path']
    for name, manifest_path in ((SNAPSHOT, 'launch_v1/source_inventory.json'),
                                (PAIR, 'inputs/source_manifest.json')):
        inventory = read(root / name / manifest_path)
        source = root / name / 'source'
        for relative, expected in inventory['files'].items():
            path = (source / relative).resolve()
            assert path.is_relative_to(source.resolve()) and sha(path) == expected, relative
    lock_path = root / PAIR / 'inputs/launch_lock.json'
    assert sha(lock_path) == PARENT_LOCK
    lock = read(lock_path)
    for relative, expected in lock['runtime_file_sha256'].items():
        assert sha(root / PAIR / relative) == expected, relative
    for relative, key in [('ancestor_commute.py', 'ancestor_sha256'),
                          ('inputs/source_manifest.json', 'source_manifest_sha256'),
                          ('inputs/objective.json', 'objective_sha256'),
                          ('inputs/proposal_bank.json', 'proposal_bank_sha256')]:
        assert sha(root / PAIR / relative) == lock[key], relative
    assert sha(root / 'paygo_tax_comparison_20260924/run_paygo_two_rate.py') == lock['tax_driver_sha256']
    assert sha(root / SEED) == SEED_SHA
    assert sha(root / SEED_BASE / 'plan.json') == PLAN_SHA
    return c


def prepare(original, destination):
    original, destination = original.resolve(strict=True), destination.resolve()
    if destination.exists():
        raise RuntimeError('Refusing existing destination; original inputs are never overwritten')
    c = verify_original(original)
    # Hard links are safe only with replacement writes below; never open copied
    # source files in writable mode. Models, plan and checkpoint stay identical.
    shutil.copytree(original, destination, copy_function=os.link)
    changes = []
    def text_change(relative, before, after):
        path = destination / relative
        text = path.read_text()
        if text.count(before) != 1:
            raise RuntimeError('Expected one explicit substitution: ' + relative)
        old = sha(path)
        tmp = path.with_name(path.name + '.relocation.tmp')
        tmp.write_text(text.replace(before, after))
        tmp.replace(path)
        changes.append(dict(path=relative, kind='literal_only', before=before,
                            after=after, original_sha256=old, relocated_sha256=sha(path)))
    def json_change(relative, value, reason):
        path = destination / relative
        old = sha(path)
        write(path, value)
        changes.append(dict(path=relative, kind='json_path_or_pin_only', reason=reason,
                            original_sha256=old, relocated_sha256=sha(path)))
    for relative in (PAIR + '/run_pair.py', PAIR + '/ancestor_commute.py'):
        text_change(relative, TORCH_ROOT, str(destination))
    manifest_path = PAIR + '/inputs/source_manifest.json'
    manifest = read(destination / manifest_path)
    manifest['source_root'] = str(destination / PAIR / 'source')
    json_change(manifest_path, manifest, 'Relocated source root; files hash map unchanged')
    lock_path = PAIR + '/inputs/launch_lock.json'
    lock = read(destination / lock_path)
    lock['ancestor_sha256'] = sha(destination / PAIR / 'ancestor_commute.py')
    lock['source_manifest_sha256'] = sha(destination / manifest_path)
    lock['runtime_file_sha256']['run_pair.py'] = sha(destination / PAIR / 'run_pair.py')
    json_change(lock_path, lock, 'Repin two path-only runtime/manifest edits')
    text_change(SNAPSHOT + '/tools/run_e5f_utility_comparison.py', PARENT_LOCK,
                sha(destination / lock_path))
    def relocate(value):
        if isinstance(value, dict):
            return {key: relocate(item) for key, item in value.items()}
        if isinstance(value, list):
            return [relocate(item) for item in value]
        if isinstance(value, str) and value.startswith(TORCH_ROOT + '/'):
            return str(destination) + value[len(TORCH_ROOT):]
        return value
    base_relative = 'utility_four_arm_preparation_20260925_v2/launch_v1/contract.json'
    base = read(destination / base_relative)
    # setup only needs reference_root; preserve every historical input record.
    base['reference_root'] = str(destination / PAIR)
    json_change(base_relative, base, 'Only runtime ancestry root relocated')
    inventory_path = SNAPSHOT + '/launch_v1/runtime_inventory.json'
    inventory = read(destination / inventory_path)
    inventory['files']['run_e5f_utility_comparison.py'] = sha(destination / SNAPSHOT / 'tools/run_e5f_utility_comparison.py')
    json_change(inventory_path, inventory, 'Repin path-only runner lock literal')
    relocated = relocate(copy.deepcopy(c))
    for item in list(relocated['files'].values()) + [relocated['base_contract'], relocated['objective']]:
        item['sha256'] = sha(item['path'])
    relocated['source_manifest']['sha256'] = sha(relocated['source_manifest']['path'])
    relocated['portable_ancestry'] = dict(original_contract_sha256=CONTRACT_HASH,
        original_parent_lock_sha256=PARENT_LOCK, relocated_parent_lock_sha256=sha(destination / lock_path),
        model_source_bytes_unchanged=True, checkpoint_bytes_unchanged=True,
        relocation_receipt=str(destination / 'relocation_receipt.json'))
    json_change(SNAPSHOT + '/launch_v1/contract.json', relocated,
                'Execution path and relocated runtime pins only; target/objective bytes unchanged')
    receipt = dict(status='prepared_not_executed', original_root=str(original),
        local_root=str(destination), original_contract_sha256=CONTRACT_HASH,
        local_contract_sha256=sha(destination / SNAPSHOT / 'launch_v1/contract.json'),
        seed_sha256=SEED_SHA, changes=changes,
        model_sources_byte_identical=True, native_solve_count=0,
        note='Controller still requires explicit local execution authorization and compatible seed loader; no fake SLURM.')
    write(destination / 'relocation_receipt.json', receipt)
    # Recheck originals after all writes: hard-linked source/input bytes survived.
    verify_original(original)
    return receipt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original', type=Path, required=True)
    p.add_argument('--destination', type=Path, required=True)
    a = p.parse_args()
    r = prepare(a.original, a.destination)
    print(json.dumps({k: r[k] for k in ('status', 'local_root', 'local_contract_sha256', 'native_solve_count')}))


if __name__ == '__main__':
    main()
