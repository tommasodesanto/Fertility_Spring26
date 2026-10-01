"""Build a source-only, hash-pinned overlay for normalized calibration v2."""
from __future__ import annotations
import argparse, gzip, hashlib, io, json, tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
ROOT=HERE.parents[4]
REMOTE_BASE='/scratch/td2248/projects/grid_resolution_credit053_v2'

def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()

def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive',type=Path,default=HERE/'normalized_calibration_v2_stage.tar.gz')
    a=ap.parse_args()
    for name in ('run_psi.py','normalized_objective.py','plan.json','source_pins.json'):
        if not (PACKET/name).is_file(): raise SystemExit(f'Missing finalized source: {name}')
    pins=json.loads((PACKET/'source_pins.json').read_text())
    for rel,digest in pins.items():
        p=ROOT/rel
        if not p.is_file() or sha(p)!=digest: raise SystemExit(f'Source pin drift: {rel}')
    sources=set(PACKET.glob('*.py'))|set(HERE.glob('*.py'))|set(HERE.glob('*.sh'))
    sources|={PACKET/'plan.json',PACKET/'source_pins.json',PACKET/'README.md'}
    sources|={ROOT/rel for rel in pins if (ROOT/rel).suffix in {'.py','.sh','.json','.csv','.md','.npz'} and (ROOT/rel).stat().st_size<=2_000_000}
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(sources) if p.is_file()}
    missing=set(pins)-set(files)
    if missing: raise SystemExit('Pinned files omitted from package: '+', '.join(sorted(missing)))
    inv={'files':files,'includes_checkpoints':False,'includes_caches':False,'reused_remote_stage':REMOTE_BASE,'source_prefix':'source/'}
    inv_path=HERE/'inventory.json';inv_path.write_text(json.dumps(inv,indent=2,sort_keys=True)+'\n')
    entries=[('source/'+rel,ROOT/rel) for rel in sorted(files)]
    entries += [('inventory.json',inv_path)]
    entries += [(name,HERE/name) for name in ('launch_torch.sh','launch_smoke_torch.sh','submit_torch.sh')]
    for _,p in entries:
        if not p.is_file():raise SystemExit(f'Missing deployment artifact: {p}')
    with a.archive.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped:
            with tarfile.open(fileobj=zipped,mode='w') as tar:
                for name,path in entries:
                    blob=path.read_bytes();info=tarfile.TarInfo(name)
                    info.size=len(blob);info.mtime=0;info.mode=0o755 if path.suffix=='.sh' else 0o644
                    tar.addfile(info,io.BytesIO(blob))
    receipt={'archive':str(a.archive),'sha256':sha(a.archive),'bytes':a.archive.stat().st_size,'source_files':len(files),'source_pins':len(pins),'reused_remote_stage':REMOTE_BASE}
    (HERE/'stage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':main()
