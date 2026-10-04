"""Verify exact staged source and both matched-arm inputs without model solves."""
import hashlib,json,sys
from pathlib import Path
REMOTE=Path('/scratch/td2248/projects/estate_birth_continuation_20261004_v1')
REPO=Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
STARTS='output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    mounted=sys.argv[1:]==['--container'];stage=Path('/work/deployment') if mounted else REMOTE
    root=REPO if mounted else stage/'source';inventory=json.loads((stage/'inventory.json').read_text())
    for rel,digest in inventory['files'].items():
        if not (root/rel).is_file() or sha(root/rel)!=digest:raise RuntimeError('Source pin drift: '+rel)
    for rel,digest in inventory['entrypoints'].items():
        if sha(stage/rel)!=digest:raise RuntimeError('Entrypoint drift: '+rel)
    assert inventory['parent_inventory_sha256']=='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
    assert inventory['files']['code/model/experiments/birth_count_choice/cluster_calibrate.py']==inventory['derivative_driver_sha256']
    print(json.dumps(dict(status='passed_zero_solves',source_files=len(inventory['files']))))
if __name__=='__main__':main()
