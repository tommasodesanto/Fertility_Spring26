"""Verify continuation source pins without running the model."""
import hashlib,json,sys
from pathlib import Path
REMOTE=Path('/scratch/td2248/projects/estate_birth_binary_continuation_20261004_v2')
REPO=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(stage,container=False):
    inv=json.loads((stage/'inventory.json').read_text());root=REPO if container else stage/'source'
    for rel,digest in inv['files'].items():
        if not (root/rel).is_file() or sha(root/rel)!=digest: raise RuntimeError('Source pin drift: '+rel)
    for rel,digest in inv['entrypoints'].items():
        if not (stage/rel).is_file() or sha(stage/rel)!=digest: raise RuntimeError('Entrypoint pin drift: '+rel)
    assert inv['parent_inventory_sha256']=='974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22'
    assert inv['files']['code/model/experiments/birth_count_choice/cluster_calibrate.py']==inv['derivative_driver_sha256']
    return dict(status='passed_zero_solves',source_files=len(inv['files']),target_fingerprint=inv['target_fingerprint'],weight_fingerprint=inv['weight_fingerprint'])
if __name__=='__main__':
    mounted='--container' in sys.argv
    print(json.dumps(verify(Path('/work/deployment') if mounted else REMOTE,mounted)))
