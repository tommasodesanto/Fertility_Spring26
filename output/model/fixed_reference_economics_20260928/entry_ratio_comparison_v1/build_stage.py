"""Compact overlay only; reuse authenticated previous remote stage."""
from pathlib import Path
import gzip,hashlib,io,json,tarfile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
names=['run_comparison.py','source_pins.json','README.md','launch_torch.sh','build_stage.py','plan.json']
pins={str((HERE/n).relative_to(ROOT)):hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
(HERE/'inventory.json').write_text(json.dumps(dict(files=pins,includes_checkpoints=False,reused_remote_stage='/scratch/td2248/projects/grid_resolution_credit053_v2'),indent=2,sort_keys=True)+'\n')
archive=HERE/'entry_ratio_comparison_v1_stage.tar.gz'
with archive.open('wb') as raw:
    with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped:
        with tarfile.open(fileobj=zipped,mode='w') as tar:
            for rel in sorted(pins):
                blob=(ROOT/rel).read_bytes();info=tarfile.TarInfo('source/'+rel);info.size=len(blob);info.mode=0o644;info.mtime=0;tar.addfile(info,io.BytesIO(blob))
            for n in ('inventory.json','launch_torch.sh'):
                blob=(HERE/n).read_bytes();info=tarfile.TarInfo(n);info.size=len(blob);info.mode=0o755 if n.endswith('.sh') else 0o644;info.mtime=0;tar.addfile(info,io.BytesIO(blob))
receipt=dict(archive=str(archive),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),bytes=archive.stat().st_size,source_files=len(pins))
(HERE/'stage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
