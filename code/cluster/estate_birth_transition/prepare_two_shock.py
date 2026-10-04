#!/usr/bin/env python3
"""Authenticate and stage an isolated, root-preserving Estate-A source package."""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, sys, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from two_shock_source_overlay import prepare_overlay

REPO = Path(__file__).resolve().parents[3]
BASE = REPO / 'output/model/transition_readiness_v1/current_baseline_20261003'
V9 = BASE / 'deployment_v9'
DEST = BASE / 'two_shock_v1'
EXPECTED_INVENTORY = '15cc5f6036d57b985ea2ea9098f015dd527d2ce11ec94bd458a49913f017d977'
EXPECTED_INPUTS = {
    'output/model/publication_refactor_20260929/local_export_v1/inputs/bundle.json': '427e67a3d9dd663cd23c3f8533c55a1a64b4f9350d396c97b5c5bd4700bc90b7',
    'output/model/publication_refactor_20260929/local_export_v1/inputs/arrays.npz': 'a48ecb71055f979284e7d63ccfa19bca8570ff8a93635fd352918356c69ee0b0',
    'output/model/fertility_identification_20260928/fixed_reference_manifest.json': '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4',
    'output/model/fixed_reference_economics_20260928/sources/fixed_price_v1/run_fixed_price.py': '96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44',
}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def dump(path, data):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

def main():
    ap=argparse.ArgumentParser()
    args=ap.parse_args()
    invpath=V9/'inventory.json'; invhash=sha(invpath)
    if invhash != EXPECTED_INVENTORY: raise SystemExit(f'v9 inventory identity mismatch: {invhash}')
    inv=json.loads(invpath.read_text()); src=V9/'source'; out=DEST/'frozen/source'
    if len(inv['files']) != 709: raise SystemExit('unexpected authenticated source count')
    for rel,digest in inv['files'].items():
        p=src/rel
        if not p.is_file() or sha(p)!=digest: raise SystemExit(f'v9 source mismatch: {rel}')
    for rel,digest in inv['entrypoints'].items():
        p=V9/rel
        if not p.is_file() or sha(p)!=digest: raise SystemExit(f'v9 entrypoint mismatch: {rel}')
    if out.exists():
        for rel,digest in inv['files'].items():
            if not (out/rel).is_file() or sha(out/rel)!=digest:
                raise SystemExit(f'existing frozen source does not match v9: {rel}')
    else:
        out.mkdir(parents=True)
        for rel in inv['files']:
            dest=out/rel; dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(src/rel,dest)
            if sha(dest)!=inv['files'][rel]: raise SystemExit(f'copy verification failed: {rel}')
    aux_receipts=[]
    for rel,digest in EXPECTED_INPUTS.items():
        original=REPO/rel
        if not original.is_file() or sha(original)!=digest:
            raise SystemExit(f'auxiliary input is not authenticated by the frozen parent pins: {rel}')
        dest=out/rel; dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():
            if sha(dest)!=digest: raise SystemExit(f'existing auxiliary input pin mismatch: {rel}')
        else: shutil.copyfile(original,dest)
        if sha(dest)!=digest: raise SystemExit(f'auxiliary input copy mismatch: {rel}')
        aux_receipts.append({'relative_path':rel,'sha256':digest,'authenticated_by':
            'v9 launch_torch.sh fixed bundle/arrays SHA pins' if rel.endswith(('bundle.json','arrays.npz')) else
            'v9 frozen source code/model/refactor_lab/inputs.py MANIFEST_SHA256' if rel.endswith('fixed_reference_manifest.json') else
            'v9 frozen small_credit_replication driver.py AUTH_SHA'})
    overlay=prepare_overlay(REPO,DEST/'frozen',V9/'inventory.json')
    plans_src=src/'output/model/transition_readiness_v1/current_baseline_20261003/plans_v9'
    planroot=DEST/'inputs'
    oldroot=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
    newroot=out
    for planname,outputname in [('fit_plan.json','base_fit_plan.json'),('smoke_plan.json','base_smoke_plan.json')]:
        d=json.loads((plans_src/planname).read_text())
        # Relocate path pins and prove each underlying file retains its v9 hash.
        for item in d['source_files'].values():
            rel=Path(item['path']).relative_to(oldroot)
            item['path']=str(newroot/rel)
        hp=Path(d['handoff']['path']).relative_to(oldroot)
        d['handoff']['path']=str(newroot/hp)
        # Contract target paths relocate; target values and all model controls stay byte-for-byte semantic.
        for key in ('annual','blocks'):
            rel=Path(d['target_contract'][key]['path']).relative_to(oldroot)
            d['target_contract'][key]['path']=str(newroot/rel)
        for key,item in d['source_files'].items():
            if sha(item['path']) != item['sha256']: raise SystemExit(f'source content hash changed during relocation: {key}')
        if sha(d['handoff']['path']) != d['handoff']['sha256']: raise SystemExit('handoff content hash changed during relocation')
        for key in ('annual','blocks'):
            if sha(d['target_contract'][key]['path']) != d['target_contract'][key]['sha256']:
                raise SystemExit(f'target input content hash changed during relocation: {key}')
        dump(planroot/outputname,d)
    receipt={'schema':'two_shock_frozen_source_package_v1','source_inventory_sha256':invhash,
        'stage_receipt_sha256':sha(V9/'stage_receipt.json'),'source_file_count':len(inv['files']),
        'source_bytes':sum((out/r).stat().st_size for r in inv['files']),
        'entrypoints_sha256':inv['entrypoints'],'auxiliary_inputs':aux_receipts,
        'legacy_source_overlay':{'manifest_sha256':overlay['manifest_sha256'],
            'manifest_file_count':overlay['manifest_file_count'],'recovered_from_v9':overlay['recovered_from_v9'],
            'overlay_metadata':'frozen/source_overlay/overlay.json'},
        'plan_sources':{'base_fit_plan.json':sha(planroot/'base_fit_plan.json'),'base_smoke_plan.json':sha(planroot/'base_smoke_plan.json')},
        'source_authentication':'every file and entrypoint checked against immutable deployment_v9 inventory before physical copy; destination rehashed after copy',
        'native_calls':0,'model_solves':0}
    dump(DEST/'frozen_source_completion_receipt.json',receipt)
    print(json.dumps(receipt,indent=2))
if __name__=='__main__': main()
