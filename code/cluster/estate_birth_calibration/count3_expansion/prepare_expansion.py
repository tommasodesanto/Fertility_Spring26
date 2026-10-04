"""Prepare five bounded count-three starts from verified parent endpoints; zero solves."""
import argparse
import hashlib
import json
import math
from pathlib import Path

BASE_SHA = 'a3f82c0b40d64911ad1354e77c777553eb9e27260abfb4c8e88a55a323f9d6cc'
PARENT_SHA = 'd14a39bcbb55067060ca492943f1e0067000c10733a48c3a3aff3b06e61e2afe'
SOURCE_PLAN = 'output/model/experiments/birth_count_choice/estate_a_continuation_20261004_v1/plan_count3.json'
NEW_PLAN = 'output/model/experiments/birth_count_choice/estate_a_count3_expansion_20261004_v1/start_plan.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def make_starts(rows, bounds):
    require(len(rows) == 10 and {r['task'] for r in rows} == set(range(10)), 'Ten parent endpoints required')
    binary = min((r for r in rows if r['arm'] == 'binary'), key=lambda r: (r['native_loss'], r['task']))
    count3 = min((r for r in rows if r['arm'] == 'count3'), key=lambda r: (r['native_loss'], r['task']))
    b, c = binary['parameters'], count3['parameters']
    require(set(b) == set(c) == set(bounds) and len(bounds) == 10, 'Ten free coordinates required')
    one = dict(b)
    two = {k: (b[k] + c[k]) / 2 for k in c}
    three = dict(c)
    three['beta_annual'] = .94
    three['first_birth_fixed_cost'] = c['first_birth_fixed_cost'] * .5
    four = dict(c)
    four['kappa_fert'] = c['kappa_fert'] * 3
    four['kappa_fert_continuation'] = c['kappa_fert_continuation'] * 3
    five = dict(c)
    five['theta0'] = c['theta0'] * 2
    five['psi_child'] = c['psi_child'] * 1.5
    five['h_P'] = c['h_P'] - .25
    starts = [{k: min(bounds[k][1], max(bounds[k][0], float(v))) for k, v in row.items()}
              for row in (one, two, three, four, five)]
    hashes = [hashlib.sha256(json.dumps(s, sort_keys=True, separators=(',', ':')).encode()).hexdigest() for s in starts]
    originals = {hashlib.sha256(json.dumps(r['parameters'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
                 for r in rows if r['arm'] == 'count3'}
    require(len(set(hashes)) == 5 and not set(hashes) & originals, 'Expansion starts duplicate each other or original count-three endpoints')
    require(all(all(math.isfinite(v) and bounds[k][0] <= v <= bounds[k][1] for k, v in s.items()) for s in starts), 'Bounded seed failure')
    return starts, dict(binary_best=binary['task'], count3_best=count3['task'])


def prepare(parent, base, stage):
    require(sha(parent / 'inventory.json') == PARENT_SHA, 'Parent inventory drift')
    require(sha(base / 'inventory.json') == BASE_SHA, 'Base continuation inventory drift')
    inventory = read(stage / 'inventory.json')
    require(inventory['base_inventory_sha256'] == BASE_SHA and inventory['original_parent_inventory_sha256'] == PARENT_SHA, 'Expansion source pin drift')
    terminal = read(base / 'control/controller_terminal.json')
    require(terminal['exit_code'] == 0 and terminal['controller_job_id'] == '19136605', 'Base controller did not pass')
    smoke = read(base / 'control/smoke_gate.json')
    require(smoke['status'] == 'matched_both_arm_smoke_passed' and smoke['stage_inventory_sha256'] == BASE_SHA, 'Base two-arm native smoke drift')
    production = read(base / 'control/production_submission.json')
    require(production['controller_job_id'] == '19136605' and production['array'] == '0-9%10' and production['cores_max'] == 10, 'Base production release drift')
    gate_path = base / 'control/gate/parent_gate.json'
    gate = read(gate_path)
    require(gate['status'] == 'all_ten_parent_endpoints_verified' and gate['parent_inventory_sha256'] == PARENT_SHA, 'Parent gate drift')
    rows = gate['parent_receipts']
    require([r['task'] for r in rows] == list(range(10)), 'Parent task order/map incomplete')
    for r in rows:
        folder = parent / 'results' / f"production_{r['arm']}_chain_{r['chain']}"
        require(sha(folder / 'run/completed.json') == r['completed_sha256'], 'Parent selected receipt drift')
        done = read(folder / 'run/completed.json')
        require(done['status'] == 'selected_numerically_verified' and done['selected']['parameters'] == r['parameters'] and done['native_loss'] == r['native_loss'], 'Parent endpoint gate drift')
        require(done['repeat']['status'] == 'exact_full_ge_repeat_passed' and len(done['repeat']['standard_plot_hashes']) == 17, 'Parent exact repeat drift')
        require(len(done['target_fit']) == 14 and len(done['parameters']) == 31, 'Parent full report drift')
        require(done['target_fingerprint'] == inventory['target_fingerprint'] and done['weight_fingerprint'] == inventory['weight_fingerprint'], 'Parent objective drift')
    base_plan_path = base / 'source' / SOURCE_PLAN
    require(sha(base_plan_path) == gate['plan_count3_sha256'], 'Base count-three plan bytes drift')
    base_plan = read(base_plan_path)
    require(base_plan['continuation_arm'] == 'count3' and base_plan['parent_receipts'] == rows, 'Base count-three plan drift')
    original_plan = read(parent / 'source/output/model/experiments/birth_count_choice/estate_a_calibration_v1/start_plan.json')
    require(base_plan['bounds'] == original_plan['bounds'] and base_plan['target_contract'] == original_plan['target_contract'], 'Parent bounds or target contract drift')
    starts, tasks = make_starts(rows, base_plan['bounds'])
    plan = dict(base_plan)
    plan.update(starts=starts, continuation_arm='count3', source_tasks=tasks,
                parent_receipt_root='/work/parent', base_controller_root='/work/basecontroller',
                base_controller_inventory_sha256=BASE_SHA, parent_gate_sha256=sha(gate_path),
                base_count3_plan_sha256=sha(base_plan_path),
                base_smoke_gate_sha256=sha(base / 'control/smoke_gate.json'),
                provisional_seed_source='/work/basecontroller/control/gate/parent_gate.json',
                provisional_seed_source_sha256=sha(gate_path),
                provisional_seed=dict(source='verified_parent_binary_and_count3_endpoints', verified=True),
                deterministic_seed=20261004,
                start_generation='Best verified binary endpoint; binary/count3 midpoint; count3 beta=.94 and first cost/2; count3 both fertility choice scales*3; count3 theta0*2, psi_child*1.5, h_P-.25; all clipped to existing bounds.',
                no_auto_retry=True, no_parameter_promotion=True, production_release_required=False)
    plan['per_chain'] = dict(cpus=1, memory_GiB=24, max_objective_calls=500,
                             max_lifecycle_per_GE=32, final_native_reserve_seconds=1800,
                             absolute_stop='2026-10-04T10:00:00-04:00')
    plan['economic_changes']['count3_expansion'] = 'No new economics, targets, weights, bounds or gates; deterministic search seeds only.'
    path = stage / 'source' / NEW_PLAN
    path.parent.mkdir(parents=True, exist_ok=False)
    path.write_text(json.dumps(plan, indent=2, sort_keys=True) + '\n')
    receipt = dict(status='five_count3_starts_prepared_zero_solves', parent_job_id='19127370',
                   base_controller_job_id='19136605', base_inventory_sha256=BASE_SHA,
                   parent_gate_sha256=sha(gate_path), base_smoke_gate_sha256=sha(base / 'control/smoke_gate.json'),
                   base_count3_plan_sha256=sha(base_plan_path),
                   plan_sha256=sha(path), source_tasks=tasks, starts=starts,
                   target_fingerprint=inventory['target_fingerprint'], weight_fingerprint=inventory['weight_fingerprint'],
                   no_auto_retry=True, no_scientific_adoption=True)
    (stage / 'control').mkdir(exist_ok=True)
    out = stage / 'control/start_plan_receipt.json'
    require(not out.exists(), 'Refusing duplicate start plan receipt')
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--stage', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.parent, args.base, args.stage)))
