#!/usr/bin/env python3
"""Package exact-path native transition inputs; never solve or submit.

The package is mounted at original_root with Apptainer, preserving scientific
JSON/source bytes and absolute-path authentication. This is not a relocation
rewrite. Run again with a fresh output after the lead pins a new root plan.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[2]
PORTABLE = ROOT / 'tmp/e5f_overnight_local_20260927/portable'

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def build(plan, output, archive=False):
    plan = plan.resolve(strict=True)
    if output.exists():
        raise ValueError('Fresh output directory required')
    files = set()
    def add(p):
        p = Path(p).absolute()
        if not p.is_relative_to(ROOT):
            raise ValueError('Dependency outside project requires explicit review: ' + str(p))
        if p.is_file(): files.add(p)
    add(plan)
    # Small ancestry bundles needed by frozen driver.verify/pair_runtime.
    for name in ('nightpair_20260925_v1', 'calibration_code_integration_20260927_v2',
                 'utility_overnight_20260923_v1', 'utility_four_arm_preparation_20260925_v2',
                 'paygo_tax_comparison_20260924', 'tools_v4'):
        directory = PORTABLE / name
        if not directory.is_dir(): raise FileNotFoundError(directory)
        for p in directory.rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts: add(p)
    base = PORTABLE / 'night_launch_v4/primary_continuation'
    for p in base.iterdir():
        if p.is_file(): add(p)
    for p in (base/'search/de_0093/case').iterdir():
        if p.is_file(): add(p)
    # Current executable sources; no compiled Mac environment or generated results.
    for p in (ROOT/'code/model').rglob('*.py'):
        if not any(x in p.parts for x in ('.venv','__pycache__')): add(p)
    # Follow file-valued JSON references, including source-pin map keys. Directory
    # strings are not traversed: an unseen closure dependency must fail startup.
    seen = set()
    while True:
        pending = [p for p in files if p.suffix == '.json' and p not in seen]
        if not pending: break
        for p in pending:
            seen.add(p)
            value = json.loads(p.read_text())
            def walk(v):
                if isinstance(v,dict):
                    for k,item in v.items():
                        if k == 'path' and 'sha256' in v: walk(item)
                        elif isinstance(k,str) and k.startswith(str(ROOT)+'/') and isinstance(item,str) and len(item)==64: walk(k)
                        elif k.endswith('_receipt_path'): walk(item)
                        elif isinstance(item,(dict,list)): walk(item)
                elif isinstance(v,list):
                    for item in v: walk(item)
                elif isinstance(v,str) and v.startswith(str(ROOT)+'/'):
                    add(Path(v))
            walk(value)
    records = {str(p.relative_to(ROOT)): {'sha256':sha(p),'bytes':p.stat().st_size}
               for p in sorted(files)}
    total = sum(r['bytes'] for r in records.values())
    if total > 2 * 1024**3:
        raise ValueError('Dependency package exceeds explicit 2GiB preparation cap')
    output.mkdir(parents=True)
    receipt = dict(status='prepared_not_runtime_verified', original_root=str(ROOT),
                   plan_path=str(plan), plan_sha256=sha(plan), files=records,
                   total_bytes=total, relocation='exact original path via Apptainer bind; no byte rewriting',
                   limitations=['Requires Linux dependency import/startup check',
                                'Requires authenticated cluster exact-loop numerical smoke',
                                'Does not authorize submission or imply equilibrium certification'])
    (output/'manifest.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    if archive:
        with tarfile.open(output/'project_inputs.tar','w',dereference=True) as tar:
            for relative,record in records.items():
                p=ROOT/relative
                if sha(p)!=record['sha256']: raise RuntimeError('Source changed during packaging: '+relative)
                tar.add(p,arcname=relative,recursive=False)
    return receipt

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--archive',action='store_true')
    args=parser.parse_args()
    receipt=build(args.plan,args.output,args.archive)
    print(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2))
