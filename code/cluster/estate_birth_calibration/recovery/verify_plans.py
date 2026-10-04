"""Verify dynamic plan bytes and pinned parent checkpoint source, without model solves."""
import hashlib,json
from pathlib import Path
ROOT=Path('/scratch/td2248/projects/estate_birth_recovery_20261004_v1')
PARENT=Path('/scratch/td2248/projects/estate_birth_calibration_20261003_v3')
REL='output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    receipt=json.loads((ROOT/'control/starts_receipt.json').read_text())
    assert receipt['status']=='provisional_parent_saved_best_only'
    assert sha(PARENT/'inventory.json')==receipt['parent_inventory_sha256']=='d14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
    for arm in ('binary','count3'):
        path=ROOT/'source'/REL/f'plan_{arm}.json'
        assert sha(path)==receipt[f'plan_{arm}_sha256']
        plan=json.loads(path.read_text())
        assert plan['recovery_arm']==arm and len(plan['starts'])==(5 if arm=='binary' else 10)
    print(json.dumps(dict(status='plans_verified_zero_solves',binary_sha=receipt['plan_binary_sha256'],count3_sha=receipt['plan_count3_sha256'])))
if __name__=='__main__':main()
