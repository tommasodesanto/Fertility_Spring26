"""Compact read-only Torch collection; selected diagnostic PNGs, no arrays."""
import io,json,shlex,subprocess,tarfile,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
script=r'''
import io,json,tarfile,hashlib,sys
from pathlib import Path
remote=Path('/scratch/td2248/projects/grid_resolution_credit053_v2');root=remote/'results'
paths=set(root.glob('full/*.json'))|set(root.glob('full/*.csv'))|set(root.glob('preflight/*.json'))
for arm in ['control_160x15','proposal_120x9']:
    prefix=root/'full'/arm
    paths |= set(prefix.glob('*.json'))
    paths |= set(prefix.glob('seed/*.json'))
    paths |= set(prefix.glob('phase_b_ge/*/closure.json'))
    paths |= set(prefix.glob('phase_b_ge/*/stage/summary.json'))
    paths |= set(prefix.glob('phase_b_ge/*.json'))
    for final in ['selected_root','selected_repeat_final']:
        folder=prefix/'phase_b_ge'/final
        paths |= set(folder.glob('*.json'))|set(folder.glob('*.csv'))
    paths |= set((prefix/'phase_b_ge/selected_root/standard_diagnostics').glob('*.png'))
paths={p for p in paths if p.is_file()}
assert sum(p.stat().st_size for p in paths)<40_000_000
files={str(p.relative_to(root)):dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in sorted(paths)}
repeatplots={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for arm in ['control_160x15','proposal_120x9'] for p in (root/'full'/arm/'phase_b_ge/selected_repeat_final/standard_diagnostics').glob('*.png')}
inv=json.loads((remote/'inventory.json').read_text());source_verified={rel:hashlib.sha256((remote/'source'/rel).read_bytes()).hexdigest()==digest for rel,digest in inv['files'].items()};assert all(source_verified.values())
receipt=dict(remote_root=str(remote),files=files,repeat_plot_hashes=repeatplots,source_inventory_verified=source_verified,archive_sha256=hashlib.sha256((remote/'grid_resolution_credit053_v2_stage.tar.gz').read_bytes()).hexdigest(),remote_writes=False,arrays_downloaded=False)
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as archive:
    for p in sorted(paths):archive.add(p,arcname=str(p.relative_to(root)))
    blob=json.dumps(receipt,indent=2).encode();info=tarfile.TarInfo('remote_hash_receipt.json');info.size=len(blob);archive.addfile(info,io.BytesIO(blob))
'''
result=subprocess.run(['ssh','torch','/share/apps/anaconda3/2025.06/bin/python -c '+shlex.quote(script)],capture_output=True,check=True)
with tarfile.open(fileobj=io.BytesIO(result.stdout),mode='r:gz') as archive:archive.extractall(HERE,filter='data')
receipt=json.loads((HERE/'remote_hash_receipt.json').read_text())
for rel,row in receipt['files'].items():assert hashlib.sha256((HERE/rel).read_bytes()).hexdigest()==row['sha256'],rel
print(json.dumps(dict(files=len(receipt['files']),bytes=sum(v['bytes'] for v in receipt['files'].values()),selected_pngs=sum(p.endswith('.png') for p in receipt['files']),repeat_png_hashes=len(receipt['repeat_plot_hashes']),verified_source_files=len(receipt['source_inventory_verified']),archive_sha256=receipt['archive_sha256'])))
