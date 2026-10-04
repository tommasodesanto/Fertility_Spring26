"""Synthetic receipts exercise the exact two-call loop and fail-closed gates; zero model solves."""
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
ARCHIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/stage.tar.gz'
DERIVATIVE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1/deployment/stage.tar.gz'
SMOKE = ROOT / 'output/model/experiments/birth_count_choice/estate_a_calibration_v1/deployment/attempt3/smoke_collection_reviewed.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def extract(archive, member, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        path.write_bytes(tar.extractfile(member).read())


def import_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    prepare = import_file('prepare_continuation', HERE / 'prepare_continuation.py')
    collector = import_file('continuation_collector', HERE / 'collect_torch.py')
    with tempfile.TemporaryDirectory() as temporary:
        home = Path(temporary)
        parent, stage, repo = home / 'parent', home / 'stage', home / 'repo'
        extract(ARCHIVE, 'inventory.json', parent / 'inventory.json')
        extract(ARCHIVE, 'source/' + prepare.OLD_PLAN, parent / 'source' / prepare.OLD_PLAN)
        extract(DERIVATIVE, 'inventory.json', stage / 'inventory.json')
        driver_rel = 'code/model/experiments/birth_count_choice/cluster_calibrate.py'
        extract(DERIVATIVE, 'source/' + driver_rel, repo / driver_rel)
        anchor_rel = 'output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json'
        extract(DERIVATIVE, 'source/' + anchor_rel, repo / anchor_rel)
        smoke = json.loads(SMOKE.read_text())['chains']
        assert len(smoke) == 2
        old_inv = json.loads((parent / 'inventory.json').read_text())
        accounting = []
        for task in range(10):
            arm = 'binary' if task < 5 else 'count3'
            chain = task % 5
            sample = copy.deepcopy(smoke[0 if arm == 'binary' else 1])
            sample['chain'] = chain
            sample['status'] = 'selected_numerically_verified'
            folder = parent / 'results' / f'production_{arm}_chain_{chain}'
            job_id = str(900000 + task)
            write(folder / 'launcher_start.json', dict(slurm_job_id=job_id, stage_inventory_sha256=prepare.PARENT_SHA))
            write(folder / 'launcher_terminal.json', dict(slurm_job_id=job_id, exit_code=0))
            write(folder / 'run/start_contract.json', dict(arm=arm, chain=chain, birth_cap=1 if arm == 'binary' else 3,
                 starts_file_sha256=old_inv['start_plan_sha256'], target_fingerprint=old_inv['target_fingerprint'],
                 weight_fingerprint=old_inv['weight_fingerprint']))
            write(folder / 'run/search_completed.json', dict(status='provisional_search_finished'))
            write(folder / 'run/native_postcheck/completed.json', dict(status='full_native_postcheck_passed',
                 search_receipt_sha256=sha(folder / 'run/search_completed.json')))
            write(folder / 'run/completed.json', sample)
            accounting.append(f'{prepare.JOB}_{task}|{job_id}|COMPLETED|0:0')
        sacct = home / 'sacct.txt'
        sacct.write_text('\n'.join(accounting) + '\n')
        gate = prepare.prepare(parent, stage, stage / 'control/gate', sacct)
        assert gate['status'] == 'all_ten_parent_endpoints_verified'
        for arm in ('binary', 'count3'):
            path = stage / 'source' / prepare.NEW_DIR / f'plan_{arm}.json'
            plan = json.loads(path.read_text())
            plan['parent_receipt_root'] = str(parent)
            local_plan = repo / prepare.NEW_DIR / f'plan_{arm}.json'
            write(local_plan, plan)
            command = [sys.executable, str(repo / driver_rel), '--arm', arm, '--chain', '0',
                       '--out', str(home / f'mock_run_{arm}'), '--deadline-epoch', str(time.time() + 4000),
                       '--starts-file', str(local_plan), '--starts-file-sha256', sha(local_plan), '--mock-smoke']
            subprocess.run(command, check=True, capture_output=True, text=True)
            result = json.loads((home / f'mock_run_{arm}/completed.json').read_text())
            assert result['status'] == 'mock_loop_passed_zero_solves' and result['objective_calls'] == 2
            if arm == 'binary':
                bad = copy.deepcopy(plan)
                bad['parent_receipts'][0]['completed_sha256'] = '0' * 64
                bad_path = repo / prepare.NEW_DIR / 'bad_plan.json'
                write(bad_path, bad)
                bad_command = command[:]
                bad_command[bad_command.index('--out') + 1] = str(home / 'bad_run')
                failure = subprocess.run(bad_command[:bad_command.index('--starts-file') + 1] +
                    [str(bad_path), '--starts-file-sha256', sha(bad_path), '--mock-smoke'],
                    capture_output=True, text=True)
                assert failure.returncode != 0 and 'Parent completion bytes drift' in failure.stderr
        # Native-shaped smoke receipts exercise collector with separate arm plan hashes.
        inventory_sha = sha(stage / 'inventory.json')
        inv = json.loads((stage / 'inventory.json').read_text())
        for index, arm in enumerate(('binary', 'count3')):
            plan_path = stage / 'source' / prepare.NEW_DIR / f'plan_{arm}.json'
            plan = json.loads(plan_path.read_text())
            plan_sha = sha(plan_path)
            folder = stage / 'results' / f'smoke_{arm}_chain_0'
            sample = copy.deepcopy(smoke[index])
            sample['starts_file_sha256'] = plan_sha
            write(folder / 'launcher_start.json', dict(stage_inventory_sha256=inventory_sha,
                 mode='smoke', arm=arm, chain=0))
            write(folder / 'run/start_contract.json', dict(arm=arm, chain=0, birth_cap=1 if arm == 'binary' else 3,
                 experiment_flags=dict(birth_count_choice_enabled=True,birth_count_choice_cap=1 if arm == 'binary' else 3,
                    bequest_net_of_selling_cost=True, estate_flow_net_of_selling_cost=True),
                 selected_source_sha256=inv['selected_source_sha256'], starts_file_sha256=plan_sha,
                 seed=plan['starts'][0], all_starts=plan['starts'], bounds=plan['bounds'],
                 provisional_seed_source_sha256=inv['parent_inventory_sha256'],
                 target_fingerprint=inv['target_fingerprint'], weight_fingerprint=inv['weight_fingerprint']))
            write(folder / 'run/search_completed.json', dict(starts_file_sha256=plan_sha))
            write(folder / 'run/native_postcheck/completed.json', dict(status='full_native_postcheck_passed',
                 search_receipt_sha256=sha(folder / 'run/search_completed.json'), starts_file_sha256=plan_sha))
            write(folder / 'run/completed.json', sample)
        result = collector.collect(stage, 'smoke')
        assert len(result['chains']) == 2 and all(r['status'] == 'selected_numerically_verified' for r in result['chains'])
        wrong = stage / 'source' / prepare.NEW_DIR / 'plan_binary.json'
        wrong.write_text(wrong.read_text() + ' ')
        try:
            collector.collect(stage, 'smoke')
            raise AssertionError('Tampered plan escaped collection gate')
        except AssertionError as exc:
            assert 'Tampered' not in str(exc)
    print(json.dumps(dict(status='passed_zero_model_solves', mock_objective_calls_per_arm=2,
                          parent_gate='ten_verified_synthetic_receipts',
                          negatives=['parent_completed_sha', 'collector_plan_sha'])))


if __name__ == '__main__':
    main()
