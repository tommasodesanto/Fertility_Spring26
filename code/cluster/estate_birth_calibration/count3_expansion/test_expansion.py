"""Synthetic authenticated receipts test seed propagation and exact two-call mock loop; no model solves."""
import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT_ARCHIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/stage.tar.gz'
BASE_ARCHIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1/deployment/stage.tar.gz'
ARCHIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/deployment/stage.tar.gz'
SMOKE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/smoke_collection_reviewed.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + '\n')


def extract(archive, name, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        path.write_bytes(tar.extractfile(name).read())


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    prepare = load('prepare_expansion', HERE / 'prepare_expansion.py')
    verify = load('verify_plans', HERE / 'verify_plans.py')
    collect = load('collect_torch', HERE / 'collect_torch.py')
    with tempfile.TemporaryDirectory() as temporary:
        home = Path(temporary)
        parent, base, stage, repo = (home / p for p in ('parent', 'base', 'stage', 'repo'))
        extract(PARENT_ARCHIVE, 'inventory.json', parent / 'inventory.json')
        extract(PARENT_ARCHIVE, 'source/output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json',
                parent / 'source/output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json')
        extract(BASE_ARCHIVE, 'inventory.json', base / 'inventory.json')
        extract(ARCHIVE, 'inventory.json', stage / 'inventory.json')
        driver = 'code/model/experiments/birth_count_choice/cluster_calibrate.py'
        anchor = 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
        extract(ARCHIVE, 'source/' + driver, repo / driver)
        extract(ARCHIVE, 'source/' + anchor, repo / anchor)
        old = json.loads((parent / 'source/output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json').read_text())
        samples = json.loads(SMOKE.read_text())['chains']
        rows = []
        for task in range(10):
            arm = 'binary' if task < 5 else 'count3'
            chain = task % 5
            sample = copy.deepcopy(samples[0 if task < 5 else 1])
            sample['arm'], sample['chain'] = arm, chain
            if arm == 'count3':
                point = sample['selected']['parameters']
                point['beta_annual'] = .955 + chain * .001
                point['first_birth_fixed_cost'] = .5
                point['theta0'] = .2
                point['h_P'] = 2.4
                for parameter_row in sample['parameters']:
                    if parameter_row['parameter'] in point:
                        parameter_row['estimate'] = str(point[parameter_row['parameter']])
            sample['native_loss'] += chain * .1
            sample['status'] = 'selected_numerically_verified'
            folder = parent / 'results' / f'production_{arm}_chain_{chain}'
            write(folder / 'run/completed.json', sample)
            job_id = str(900000 + task)
            write(folder / 'launcher_terminal.json', dict(exit_code=0, slurm_job_id=job_id))
            rows.append(dict(task=task, arm=arm, chain=chain, slurm_task_job_id=job_id,
                             completed_sha256=sha(folder / 'run/completed.json'),
                             parameters=sample['selected']['parameters'], native_loss=sample['native_loss']))
        base_plan = dict(old, continuation_arm='count3', parent_receipts=rows,
                         parent_job_id='19127370', parent_inventory_sha256=prepare.PARENT_SHA)
        base_plan_path = base / 'source' / prepare.SOURCE_PLAN
        write(base_plan_path, base_plan)
        gate_path = base / 'control/gate/parent_gate.json'
        write(gate_path, dict(status='all_ten_parent_endpoints_verified',
                              parent_inventory_sha256=prepare.PARENT_SHA,
                              parent_receipts=rows, plan_count3_sha256=sha(base_plan_path)))
        write(base / 'control/controller_terminal.json', dict(exit_code=0, controller_job_id='19136605'))
        write(base / 'control/smoke_gate.json', dict(status='matched_both_arm_smoke_passed',
                                                     stage_inventory_sha256=prepare.BASE_SHA))
        write(base / 'control/production_submission.json', dict(controller_job_id='19136605', array='0-9%10', cores_max=10))
        receipt = prepare.prepare(parent, base, stage)
        assert receipt['status'] == 'five_count3_starts_prepared_zero_solves'
        assert verify.verify(stage)['status'] == 'count3_expansion_plan_verified'
        plan_path = stage / 'source' / prepare.NEW_PLAN
        plan = json.loads(plan_path.read_text())
        assert plan['source_tasks'] == {'binary_best': 0, 'count3_best': 5}
        assert len(plan['starts']) == 5 and len({json.dumps(s, sort_keys=True) for s in plan['starts']}) == 5
        assert plan['starts'][0] == rows[0]['parameters']
        assert plan['starts'][1] == {k: (rows[0]['parameters'][k] + rows[5]['parameters'][k]) / 2 for k in old['bounds']}
        local_plan = copy.deepcopy(plan)
        local_plan['parent_receipt_root'], local_plan['base_controller_root'] = str(parent), str(base)
        local_plan_path = repo / prepare.NEW_PLAN
        write(local_plan_path, local_plan)
        command = [sys.executable, str(repo / driver), '--arm', 'count3', '--chain', '0',
                   '--out', str(home / 'mock_run'), '--deadline-epoch', str(time.time() + 4000),
                   '--starts-file', str(local_plan_path), '--starts-file-sha256', sha(local_plan_path), '--mock-smoke']
        completed = subprocess.run(command, capture_output=True, text=True)
        assert completed.returncode == 0, completed.stderr
        result = json.loads((home / 'mock_run/completed.json').read_text())
        assert result['status'] == 'mock_loop_passed_zero_solves' and result['objective_calls'] == 2
        bad = copy.deepcopy(local_plan)
        bad['starts'][0]['beta_annual'] += .001
        bad_path = repo / prepare.NEW_PLAN.replace('start_plan.json', 'bad_plan.json')
        write(bad_path, bad)
        bad_command = command[:]
        bad_command[bad_command.index('--out') + 1] = str(home / 'bad_run')
        bad_command[bad_command.index('--starts-file') + 1] = str(bad_path)
        bad_command[bad_command.index('--starts-file-sha256') + 1] = sha(bad_path)
        failure = subprocess.run(bad_command, capture_output=True, text=True)
        assert failure.returncode != 0 and 'Expansion seed construction drift' in failure.stderr
        sample = copy.deepcopy(samples[1])
        sample['starts_file_sha256'] = sha(plan_path)
        sample['arm'], sample['chain'] = 'count3', 0
        folder = stage / 'results/smoke_count3_chain_0'
        inv = json.loads((stage / 'inventory.json').read_text())
        write(folder / 'launcher_start.json', dict(stage_inventory_sha256=sha(stage / 'inventory.json'), mode='smoke', arm='count3', chain=0))
        write(folder / 'run/start_contract.json', dict(arm='count3', chain=0, birth_cap=3,
             experiment_flags=dict(birth_count_choice_enabled=True, birth_count_choice_cap=3,
                                   bequest_net_of_selling_cost=True, estate_flow_net_of_selling_cost=True),
             selected_source_sha256=inv['selected_source_sha256'], starts_file_sha256=sha(plan_path),
             seed=plan['starts'][0], all_starts=plan['starts'], bounds=plan['bounds'],
             provisional_seed_source_sha256=receipt['parent_gate_sha256'],
             target_fingerprint=inv['target_fingerprint'], weight_fingerprint=inv['weight_fingerprint']))
        write(folder / 'run/search_completed.json', dict(starts_file_sha256=sha(plan_path)))
        write(folder / 'run/native_postcheck/completed.json', dict(status='full_native_postcheck_passed',
             search_receipt_sha256=sha(folder / 'run/search_completed.json'), starts_file_sha256=sha(plan_path)))
        write(folder / 'run/completed.json', sample)
        data = collect.collect(stage, 'smoke')
        assert len(data['chains']) == 1 and data['chains'][0]['status'] == 'selected_numerically_verified'
        plan_path.write_text(plan_path.read_text() + ' ')
        try:
            collect.collect(stage, 'smoke')
            raise RuntimeError('Tampered plan escaped collector')
        except AssertionError:
            pass
    print(json.dumps(dict(status='passed_zero_model_solves', count3_mock_calls=2,
                          synthetic_parent_receipts=10, negatives=['seed_propagation', 'collector_plan_sha'])))


if __name__ == '__main__':
    main()
