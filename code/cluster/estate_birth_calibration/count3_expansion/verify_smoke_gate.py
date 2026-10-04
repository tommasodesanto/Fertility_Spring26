"""Gate five production tasks on a fresh two-call native count-three smoke."""
import hashlib
import json
import sys
from pathlib import Path

REMOTE = Path('/scratch/td2248/projects/estate_birth_count3_expansion_20261004_v1')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def verify(stage):
    inv = read(stage / 'inventory.json')
    fingerprint = sha(stage / 'inventory.json')
    plan_path = stage / 'source/output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/start_plan.json'
    plan_sha = sha(plan_path)
    assert plan_sha == read(stage / 'control/start_plan_receipt.json')['plan_sha256']
    launcher = stage / 'results/smoke_count3_chain_0'
    run = launcher / 'run'
    started, contract, done = (read(launcher / 'launcher_start.json'), read(run / 'start_contract.json'), read(run / 'completed.json'))
    child, search = read(run / 'native_postcheck/completed.json'), read(run / 'search_completed.json')
    assert started['stage_inventory_sha256'] == fingerprint and started['wall_seconds'] == 5400
    assert started['deadline_epoch'] - started['start_epoch'] == 5400
    assert contract['arm'] == done['arm'] == 'count3' and contract['chain'] == done['chain'] == 0
    assert contract['birth_cap'] == 3 and contract['starts_count'] == 5 and contract['bounds']['beta_annual'] == [.93, .99]
    assert contract['starts_file_sha256'] == done['starts_file_sha256'] == child['starts_file_sha256'] == plan_sha
    assert all(contract[k] == done[k] == child[k] == inv[k] for k in ('target_fingerprint', 'weight_fingerprint'))
    assert done['status'] == 'selected_numerically_verified' and done['objective_calls'] == 2
    assert len(read(run / 'cases.json')) == 2 and read(run / 'heartbeat.json')['status'] == 'completed'
    assert child['status'] == 'full_native_postcheck_passed' and child['search_receipt_sha256'] == sha(run / 'search_completed.json')
    assert len(done['target_fit']) == 14 and len(done['parameters']) == 31
    assert done['repeat']['status'] == 'exact_full_ge_repeat_passed'
    assert len(done['repeat']['standard_plot_hashes']) == 17 and done['repeat']['experimental_target_fit_exact']
    assert done['selected_postcheck']['status'] == 'passed'
    assert done['smoke_fast_full_comparison']['status'] == 'search_full_new_target_exact'
    assert search['search_evaluator'] == 'full_native' and search['selected_evaluator'] == 'full_native'
    return dict(status='count3_expansion_native_smoke_passed', stage_inventory_sha256=fingerprint,
                plan_sha256=plan_sha, target_fingerprint=inv['target_fingerprint'],
                weight_fingerprint=inv['weight_fingerprint'])


if __name__ == '__main__':
    print(json.dumps(verify(Path(sys.argv[1]) if len(sys.argv) == 2 else REMOTE)))
