"""Create a byte-pinned count-three search overlay from the reviewed continuation."""
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE_ARCHIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1/deployment/stage.tar.gz'
BASE_INVENTORY_SHA = 'a3f82c0b40d64911ad1354e77c777553eb9e27260abfb4c8e88a55a323f9d6cc'
OUT = ROOT / 'output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment'
REMOTE = '/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1'
DRIVER = 'code/model/experiments/birth_count_choice/cluster_calibrate.py'


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def main():
    with tarfile.open(BASE_ARCHIVE) as archive:
        base_bytes = archive.extractfile('inventory.json').read()
        assert sha(base_bytes) == BASE_INVENTORY_SHA, 'Reviewed continuation inventory drift'
        base = json.loads(base_bytes)
        source = {name.removeprefix('source/'): archive.extractfile(name).read()
                  for name in archive.getnames() if name.startswith('source/')}
    assert {name: sha(blob) for name, blob in source.items()} == base['files'], 'Base source drift'
    original = source[DRIVER].decode()
    begin, end = original.index('def checked_plan('), original.index('def experimental_parameter_rows(')
    replacement = (HERE / 'checked_plan.py.txt').read_text()
    derivative = original[:begin] + replacement + original[end:]
    assert derivative.replace(replacement, original[begin:end], 1) == original
    assert derivative.count('def checked_plan(') == 1
    assert 'start+43200' in derivative and 'Continuation arm/plan mismatch' in derivative
    source[DRIVER] = derivative.encode()
    scripts = {p.name: p.read_bytes() for p in HERE.iterdir()
               if p.is_file() and p.suffix in ('.py', '.sh') and p.name != 'build_stage.py'}
    inventory = dict(base_inventory_sha256=BASE_INVENTORY_SHA,
                     original_parent_inventory_sha256=base['parent_inventory_sha256'],
                     base_driver_sha256=base['files'][DRIVER],
                     derivative_driver_sha256=sha(source[DRIVER]),
                     derivative_scope='count3 expansion checked_plan only',
                     files={name: sha(blob) for name, blob in sorted(source.items())},
                     entrypoints={name: sha(blob) for name, blob in scripts.items()},
                     target_fingerprint=base['target_fingerprint'],
                     weight_fingerprint=base['weight_fingerprint'],
                     selected_source_sha256=base['selected_source_sha256'],
                     remote_root=REMOTE, parent_array_job_id='19127370',
                     parent_controller_job_id='19136605', no_auto_retry=True)
    OUT.mkdir(parents=True, exist_ok=True)
    inv_blob = (json.dumps(inventory, sort_keys=True, indent=2) + '\n').encode()
    (OUT / 'inventory.json').write_bytes(inv_blob)
    entries = {'source/' + name: blob for name, blob in source.items()}
    entries.update(scripts)
    entries['inventory.json'] = inv_blob
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
                   inventory_sha256=sha(inv_blob), source_files=len(source),
                   base_inventory_sha256=BASE_INVENTORY_SHA, remote_root=REMOTE)
    (OUT / 'stage_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
