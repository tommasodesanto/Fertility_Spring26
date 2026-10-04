"""Verify the pinned expansion source without model solves."""
import hashlib,json,sys
from pathlib import Path
REMOTE=Path('/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1')
REPO=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    mounted=sys.argv[1:]==['--container'];stage=Path('/work/deployment') if mounted else REMOTE
    root=REPO if mounted else stage/'source';inventory=json.loads((stage/'inventory.json').read_text())
    for rel,digest in inventory['files'].items():
        if not (root/rel).is_file() or sha(root/rel)!=digest:raise RuntimeError('Source pin drift: '+rel)
    for rel,digest in inventory['entrypoints'].items():
        if sha(stage/rel)!=digest:raise RuntimeError('Entrypoint drift: '+rel)
    assert inventory['base_inventory_sha256']=='a3f82c0b40d64911ad1354e77c777553eb9e27260abfb4c8e88a55a323f9d6cc'
    assert inventory['original_parent_inventory_sha256']=='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
    assert inventory['files']['code/model/experiments/birth_count_choice/cluster_calibrate.py']==inventory['derivative_driver_sha256']
    print(json.dumps(dict(status='passed_zero_solves',source_files=len(inventory['files']))))
if __name__=='__main__':main()
