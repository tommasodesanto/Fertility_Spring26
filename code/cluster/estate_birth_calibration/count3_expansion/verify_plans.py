"""Verify the single dynamic count-three expansion plan."""
import hashlib
import json
import sys
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(stage):
    receipt = json.loads((stage / 'control/start_plan_receipt.json').read_text())
    inventory = json.loads((stage / 'inventory.json').read_text())
    path = stage / 'source/output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/start_plan.json'
    plan = json.loads(path.read_text())
    assert receipt['status'] == 'five_count3_starts_prepared_zero_solves'
    assert sha(path) == receipt['plan_sha256']
    assert receipt['base_inventory_sha256'] == inventory['base_inventory_sha256']
    assert receipt['starts'] == plan['starts'] and len(plan['starts']) == 5
    assert plan['continuation_arm'] == 'count3' and plan['source_tasks'] == receipt['source_tasks']
    assert plan['target_fingerprint'] == inventory['target_fingerprint']
    assert plan['weight_fingerprint'] == inventory['weight_fingerprint']
    return dict(status='count3_expansion_plan_verified', plan_sha256=sha(path))


if __name__ == '__main__':
    print(json.dumps(verify(Path(sys.argv[1]))))
