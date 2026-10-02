"""Build a source-only, hash-pinned overlay for purchase-rule overnight calibration."""
from __future__ import annotations
import argparse, gzip, hashlib, io, json, tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
ROOT=HERE.parents[4]
REMOTE_BASE='/scratch/td2248/projects/normalized_floor_calibration_v1'

def sha(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()

def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive',type=Path,default=HERE/'purchase_rules_overnight_v1_stage.tar.gz')
    a=ap.parse_args()
    for name in ('run_psi.py','normalized_objective.py','plan.json','source_pins.json','engine_pins.json','manifest.json'):
        if not (PACKET/name).is_file(): raise SystemExit(f'Missing finalized source: {name}')
    pins=json.loads((PACKET/'engine_pins.json').read_text())
    for rel,digest in pins.items():
        p=PACKET/rel
        if not p.is_file() or sha(p)!=digest: raise SystemExit(f'Engine pin drift: {rel}')
    sources=set(PACKET.glob('*.py'))|set(PACKET.glob('*.json'))|set(HERE.glob('*.py'))|set(HERE.glob('*.sh'))
    sources|={PACKET/'README.md'}|set((PACKET/'center').glob('*'))
    sources|={PACKET/rel for rel in pins}
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(sources) if p.is_file()}
    inv={'files':files,'includes_checkpoints':False,'includes_caches':False,'reused_remote_stage':REMOTE_BASE,'source_prefix':'source/'}
    inv_path=HERE/'inventory.json';inv_path.write_text(json.dumps(inv,indent=2,sort_keys=True)+'\n')
    entries=[('source/'+rel,ROOT/rel) for rel in sorted(files)]
    entries += [('inventory.json',inv_path)]
    entries += [(name,HERE/name) for name in ('launch_torch.sh','submit_torch.sh')]
    for _,p in entries:
        if not p.is_file():raise SystemExit(f'Missing deployment artifact: {p}')
    with a.archive.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped:
            with tarfile.open(fileobj=zipped,mode='w') as tar:
                for name,path in entries:
                    blob=path.read_bytes();info=tarfile.TarInfo(name)
                    info.size=len(blob);info.mtime=0;info.mode=0o755 if path.suffix=='.sh' else 0o644
                    tar.addfile(info,io.BytesIO(blob))
    receipt={'archive':str(a.archive),'sha256':sha(a.archive),'bytes':a.archive.stat().st_size,'source_files':len(files),'engine_source_pins':len(pins),'reused_remote_stage':REMOTE_BASE}
    (HERE/'stage_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':main()
