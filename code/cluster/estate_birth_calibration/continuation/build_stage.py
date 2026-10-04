"""Build a separate, pinned derivative of the reviewed Estate-A stage; no solves."""
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT_ARCHIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/stage.tar.gz'
OUT = ROOT / 'output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1/deployment'
REMOTE = '/scratch/td2248/projects/estate_birth_continuation_20261004_v1'
DRIVER = 'code/model/experiments/birth_count_choice/cluster_calibrate.py'
PARENT_INVENTORY_SHA = 'd14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def main():
    with tarfile.open(PARENT_ARCHIVE) as archive:
        parent_inventory_bytes = archive.extractfile('inventory.json').read()
        assert sha(parent_inventory_bytes) == PARENT_INVENTORY_SHA, 'Parent stage inventory drift'
        parent = json.loads(parent_inventory_bytes)
        source = {name.removeprefix('source/'): archive.extractfile(name).read()
                  for name in archive.getnames() if name.startswith('source/')}
    assert {name: sha(blob) for name, blob in source.items()} == parent['files'], 'Parent stage source drift'
    original = source[DRIVER].decode()
    begin = original.index('def checked_plan(')
    end = original.index('def experimental_parameter_rows(', begin)
    replacement = (HERE / 'checked_plan.py.txt').read_text()
    derivative = original[:begin] + replacement + original[end:]
    assert derivative.replace(replacement, original[begin:end], 1) == original
    assert derivative.count('def checked_plan(') == 1
    derivative = derivative.replace('start+21600', 'start+43200', 1)
    derivative = derivative.replace("wall_seconds=21600,deterministic_seed", "wall_seconds=43200,deterministic_seed", 1)
    derivative = derivative.replace('plan,anchor=checked_plan(args.starts_file,args.starts_file_sha256)',
        "plan,anchor=checked_plan(args.starts_file,args.starts_file_sha256)\n    require(plan['continuation_arm']==args.arm, 'Continuation arm/plan mismatch')", 1)
    assert derivative.count('start+43200') == 1
    assert derivative.count('wall_seconds=43200,deterministic_seed') == 1
    assert derivative.count('Continuation arm/plan mismatch') == 1
    source[DRIVER] = derivative.encode()
    scripts = {p.name: p.read_bytes() for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py', '.sh') and p.name != 'build_stage.py'}
    inventory = dict(parent_inventory_sha256=PARENT_INVENTORY_SHA,
                     parent_driver_sha256=parent['files'][DRIVER],
                     derivative_driver_sha256=sha(source[DRIVER]),
                     derivative_scope='checked_plan and 12-hour ceiling only',
                     files={name: sha(blob) for name, blob in sorted(source.items())},
                     entrypoints={name: sha(blob) for name, blob in scripts.items()},
                     target_fingerprint=parent['target_fingerprint'],
                     weight_fingerprint=parent['weight_fingerprint'],
                     selected_source_sha256=parent['selected_source_sha256'],
                     remote_root=REMOTE, parent_job_id='19127370', no_auto_retry=True)
    OUT.mkdir(parents=True, exist_ok=True)
    inventory_bytes = (json.dumps(inventory, sort_keys=True, indent=2) + '\n').encode()
    (OUT / 'inventory.json').write_bytes(inventory_bytes)
    entries = {'source/' + name: blob for name, blob in source.items()}
    entries.update(scripts)
    entries['inventory.json'] = inventory_bytes
    with (OUT / 'stage.tar.gz').open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode='w') as archive:
            for name, blob in sorted(entries.items()):
                info = tarfile.TarInfo(name)
                info.size = len(blob)
                info.mtime = 0
                info.mode = 0o755 if name.endswith('.sh') else 0o644
                archive.addfile(info, io.BytesIO(blob))
    receipt = dict(status='prepared_no_submission', archive=str(OUT / 'stage.tar.gz'),
                   archive_sha256=sha((OUT / 'stage.tar.gz').read_bytes()),
                   inventory_sha256=sha(inventory_bytes), source_files=len(source),
                   parent_inventory_sha256=PARENT_INVENTORY_SHA, remote_root=REMOTE)
    (OUT / 'stage_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
