"""Verify arm-specific dynamic plans generated only after all parent gates pass."""
import hashlib
import json
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(stage):
    gate = json.loads((stage / 'control/gate/parent_gate.json').read_text())
    inventory = json.loads((stage / 'inventory.json').read_text())
    assert gate['status'] == 'all_ten_parent_endpoints_verified'
    assert gate['parent_inventory_sha256'] == inventory['parent_inventory_sha256']
    assert {r['task'] for r in gate['parent_receipts']} == set(range(10))
    for arm in ('binary', 'count3'):
        path = stage / 'source/output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1' / f'plan_{arm}.json'
        plan = json.loads(path.read_text())
        assert sha(path) == gate[f'plan_{arm}_sha256']
        assert plan['continuation_arm'] == arm
        assert plan['parent_receipts'] == gate['parent_receipts']
        assert plan['starts'] == [r['parameters'] for r in gate['parent_receipts'] if r['arm'] == arm]
        assert plan['target_fingerprint'] == inventory['target_fingerprint']
        assert plan['weight_fingerprint'] == inventory['weight_fingerprint']
    return dict(status='both_continuation_plans_verified', parent_job_id=gate['parent_job_id'])


if __name__ == '__main__':
    print(json.dumps(verify(Path(sys.argv[1]))))
