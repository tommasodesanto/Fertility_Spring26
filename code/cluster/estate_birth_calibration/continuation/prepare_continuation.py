"""Fail-closed parent gate and verified endpoint plan for the two birth menus."""
import argparse
import hashlib
import json
from pathlib import Path

PARENT_SHA = 'd14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
JOB = '19127370'
OLD_PLAN = 'output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json'
NEW_DIR = 'output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prepare(parent, stage, receipt_root, sacct_path):
    inventory = read(parent / 'inventory.json')
    require(sha(parent / 'inventory.json') == PARENT_SHA, 'Parent inventory bytes drift')
    require(read(stage / 'inventory.json')['parent_inventory_sha256'] == PARENT_SHA, 'Continuation parent pin drift')
    old_plan_path = parent / 'source' / OLD_PLAN
    require(sha(old_plan_path) == inventory['start_plan_sha256'], 'Parent plan bytes drift')
    old = read(old_plan_path)
    require(old['target_fingerprint'] == inventory['target_fingerprint'] and
            old['weight_fingerprint'] == inventory['weight_fingerprint'], 'Parent target/weight drift')
    require(old['arms'] == {'binary': 1, 'count3': 3} and old['bounds']['beta_annual'] == [.93, .99], 'Parent economics drift')
    sacct = {}
    for line in sacct_path.read_text().splitlines():
        fields = line.split('|')
        if len(fields) != 4 or not fields[0].startswith(JOB + '_') or '.' in fields[0]:
            continue
        task = int(fields[0].split('_', 1)[1])
        require(task not in sacct, f'Duplicate accounting row for task {task}')
        require(fields[2] == 'COMPLETED' and fields[3] == '0:0', f'Parent task {task} accounting failure')
        sacct[task] = fields[1]
    require(set(sacct) == set(range(10)), 'Parent accounting task set incomplete')
    rows = []
    for task in range(10):
        arm = 'binary' if task < 5 else 'count3'
        chain = task % 5
        folder = parent / 'results' / f'production_{arm}_chain_{chain}'
        launched = read(folder / 'launcher_start.json')
        terminal = read(folder / 'launcher_terminal.json')
        contract = read(folder / 'run/start_contract.json')
        completed_path = folder / 'run/completed.json'
        completed = read(completed_path)
        child = read(folder / 'run/native_postcheck/completed.json')
        require(launched['slurm_job_id'] == sacct[task] and launched['stage_inventory_sha256'] == PARENT_SHA,
                f'Parent task {task} launch identity failed')
        require(terminal['slurm_job_id'] == sacct[task] and terminal['exit_code'] == 0,
                f'Parent task {task} did not exit successfully')
        require(contract['arm'] == arm and contract['chain'] == chain and contract['birth_cap'] == old['arms'][arm],
                f'Parent task {task} arm/chain drift')
        require(contract['starts_file_sha256'] == inventory['start_plan_sha256'],
                f'Parent task {task} starts drift')
        require(completed['status'] == 'selected_numerically_verified' and completed['arm'] == arm and completed['chain'] == chain,
                f'Parent task {task} lacks verified endpoint')
        require(completed['repeat']['status'] == 'exact_full_ge_repeat_passed' and
                len(completed['repeat']['standard_plot_hashes']) == 17,
                f'Parent task {task} lacks exact repeat')
        require(len(completed['target_fit']) == 14 and len(completed['parameters']) == 31 and
                completed['selected_postcheck']['status'] == 'passed' and child['status'] == 'full_native_postcheck_passed',
                f'Parent task {task} native report failed')
        require(child['search_receipt_sha256'] == sha(folder / 'run/search_completed.json'),
                f'Parent task {task} selected search drift')
        require(all(completed[key] == old[key] == contract[key] for key in ('target_fingerprint', 'weight_fingerprint')),
                f'Parent task {task} target/weight drift')
        point = completed['selected']['parameters']
        reported = {r['parameter']: float(r['estimate']) for r in completed['parameters'] if r['parameter'] in point}
        require(point == reported, f'Parent task {task} selected point/report drift')
        require(set(point) == set(old['bounds']) and
                all(lo <= point[key] <= hi for key, (lo, hi) in old['bounds'].items()),
                f'Parent task {task} endpoint outside bounds')
        rows.append(dict(task=task, arm=arm, chain=chain, slurm_task_job_id=sacct[task], parameters=point,
                         completed_sha256=sha(completed_path), native_loss=completed['native_loss']))
    require(len(rows) == 10, 'Ten parent endpoints required')
    result = dict(status='all_ten_parent_endpoints_verified', parent_job_id=JOB,
                  parent_inventory_sha256=PARENT_SHA, parent_receipts=rows,
                  no_failure_restart=True, no_scientific_adoption=True)
    receipt_root.mkdir(parents=True, exist_ok=False)
    (receipt_root / 'parent_gate.json').write_text(json.dumps(result, indent=2) + '\n')
    out_dir = stage / 'source' / NEW_DIR
    out_dir.mkdir(parents=True, exist_ok=False)
    for arm in ('binary', 'count3'):
        plan = dict(old)
        plan['starts'] = [r['parameters'] for r in rows if r['arm'] == arm]
        plan['continuation_arm'] = arm
        plan['parent_job_id'] = JOB
        plan['parent_inventory_sha256'] = PARENT_SHA
        plan['parent_receipt_root'] = '/work/parent'
        plan['parent_receipts'] = rows
        plan['provisional_seed_source'] = '/work/parent/inventory.json'
        plan['provisional_seed_source_sha256'] = PARENT_SHA
        plan['provisional_seed'] = dict(source='own_arm_parent_verified_endpoints', verified=True)
        plan['deterministic_seed'] = 20261004
        plan['start_generation'] = 'Each chain starts at its own numerically verified parent endpoint; new Nelder-Mead simplex, no optimizer-state resume.'
        plan['per_chain'] = dict(cpus=1, memory_GiB=24, max_objective_calls=500,
                                 max_lifecycle_per_GE=32, final_native_reserve_seconds=1800,
                                 absolute_stop='2026-10-04T10:00:00-04:00')
        plan['no_auto_retry'] = True
        plan['no_parameter_promotion'] = True
        plan['production_release_required'] = False
        plan['economic_changes']['continuation'] = 'No further economic change; verified endpoint seeds and fresh optimizer starts.'
        path = out_dir / f'plan_{arm}.json'
        path.write_text(json.dumps(plan, sort_keys=True, indent=2) + '\n')
        result[f'plan_{arm}_sha256'] = sha(path)
    (receipt_root / 'parent_gate.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--receipt-root', type=Path, required=True)
    parser.add_argument('--sacct', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.parent, args.stage, args.receipt_root, args.sacct)))
