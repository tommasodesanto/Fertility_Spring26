"""Deterministic compact source archive; no SSH, solves or giant checkpoint copy."""
from pathlib import Path
import gzip,hashlib,io,json,tarfile
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[6]
PACKET=HERE.parent.parent
runner=json.loads((HERE.parent/'source_hashes.json').read_text())
INPUT_PACKET=PACKET.parent
prep=json.loads((INPUT_PACKET/'preflight.json').read_text())
pins=dict(prep['source_pins'],**runner)
extra=[INPUT_PACKET/'reference_parameters.csv',INPUT_PACKET/'reference_target_fit.csv',INPUT_PACKET/'proposed_120x9/arrays.npz',HERE.parent/'source_hashes.json',HERE.parent/'README.md']
extra+=list((ROOT/'code/model/refactor_lab').rglob('*.py'))
for path in extra:pins[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
for rel,digest in pins.items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
assert sum((ROOT/p).stat().st_size for p in pins)<5_000_000,'Unexpected large staging payload'
(HERE/'inventory.json').write_text(json.dumps(dict(files=pins,bytes=sum((ROOT/p).stat().st_size for p in pins),includes_checkpoints=False,includes_original_bundle=False),indent=2,sort_keys=True)+'\n')
# Package-owned folders overlay as directories; other exact pinned files individually.
directories=['code/model/refactor_lab','output/model/publication_refactor_20260929/grid_resolution_v1','output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed']
mounts=directories+[p for p in sorted(pins) if not any(p.startswith(d+'/') for d in directories)]
(HERE/'mounts.txt').write_text('\n'.join(mounts)+'\n')
archive=HERE/'grid_resolution_credit053_stage.tar.gz'
with archive.open('wb') as raw:
    with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped:
        with tarfile.open(fileobj=zipped,mode='w') as tar:
            for rel in sorted(pins):
                blob=(ROOT/rel).read_bytes();info=tarfile.TarInfo('source/'+rel);info.size=len(blob);info.mode=0o644;info.mtime=0;tar.addfile(info,io.BytesIO(blob))
            for name in ['inventory.json','mounts.txt','launch_torch.sh']:
                blob=(HERE/name).read_bytes();info=tarfile.TarInfo(name);info.size=len(blob);info.mode=0o755 if name.endswith('.sh') else 0o644;info.mtime=0;tar.addfile(info,io.BytesIO(blob))
print(json.dumps(dict(archive=str(archive),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_bytes=archive.stat().st_size,source_files=len(pins),source_bytes=sum((ROOT/p).stat().st_size for p in pins)),indent=2))
