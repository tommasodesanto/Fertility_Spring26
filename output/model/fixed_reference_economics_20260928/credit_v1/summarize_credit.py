#!/usr/bin/env python3
"""Zero-solve Torch extraction of the matched-grid borrowing mechanism.

Uses saved policies and one supplied-policy period evaluation per regime.
Negative next-saving positions are net financial debt, not identified mortgages.
"""
import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import signal
import sys
import time

LABEL = '2007 stationary reference — block0506, September 28 verified export'


def require(test, message):
    if not test:
        raise RuntimeError(message)


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def birth_cells(ev, P, grid, model):
    import numpy as np
    require(P.sequential_births and not P.joint_nested_choice and P.n_parity == 4 and
            P.child_state_mode == 'independent_count' and model.readiness_settled_state(P) == 0,
            'Birth extraction requires the authenticated sequential architecture')
    gp, policy = ev.g_pre, ev.policy
    fec = np.asarray(model.get_fecundity_by_age(P))
    rows = {}
    order_totals = np.zeros(3)
    for j in range(P.J):
        if not P.A_f_start <= j + 1 <= P.A_f_end:
            continue
        for n in range(3):
            risk = np.zeros((len(grid), 1 + P.n_house, P.I, P.Nz))
            attempts = np.zeros_like(risk)
            for m in range(n + 1):
                pool = gp[:, :, :, j, :, n, m]
                pr = (policy.fert_probs[:, :, :, j, :, 1] if n == 0 else
                      policy.fert2_probs[:, :, :, j, :, 1, n-1, m])
                risk += pool
                attempts += pool * pr
            born = attempts * fec[j]
            order_totals[n] += born.sum()
            for tenure, ts in (('renter', slice(0, 1)), ('owner', slice(1, None))):
                for wealth, mask in (('nonpositive', grid <= 0), ('positive', grid > 0)):
                    key = (float(P.age_start + P.da*j), tenure, wealth, 'first' if n == 0 else 'subsequent')
                    row = rows.setdefault(key, dict(age_left=key[0], age_right=key[0]+float(P.da),
                        inherited_tenure=tenure, inherited_net_financial_wealth=wealth,
                        birth_order=key[3], at_risk_mass=0., birth_mass=0.))
                    row['at_risk_mass'] += float(risk[mask, ts].sum())
                    row['birth_mass'] += float(born[mask, ts].sum())
    require(abs(float(order_totals.sum()) - float(ev.births)) <= 2e-10, 'Grouped birth flows do not add up')
    return rows, order_totals.tolist()


def debt_rows(ev, P, grid, model, case, scope):
    import numpy as np
    from e5f_overnight_estate_audit import policy_mass_branches
    branches = policy_mass_branches(ev, P)
    if bool(P.native_due_stayer_credit):
        other, stay = branches
    else:
        require(len(branches) == 1, 'Unexpected credit branch structure')
        mass, saving, consumption = branches[0]
        stayer = model.realize_stayer_cross_section(ev.g_post_fertility, ev.policy.loc_probs,
            ev.policy.tenure_choice, ev.policy.tenure_probs)
        other_mass = mass - stayer
        require(float(other_mass.min()) >= -1e-12, 'Stayer mass exceeds current distribution')
        other = (np.maximum(other_mass, 0.), saving, consumption)
        stay = (stayer, saving, consumption)
    groups = [('renter', other[0][:, :1], other[1][:, :1]),
              ('buyer', other[0][:, 1:], other[1][:, 1:]),
              ('owner_stayer', stay[0][:, 1:], stay[1][:, 1:])]
    total = float(ev.g_current.sum())
    require(abs(sum(float(g.sum()) for _, g, _ in groups)-total) <= 2e-10, 'Debt groups do not partition households')
    rows = []
    for name, mass, saving in groups:
        denominator = float(mass.sum())
        indebted = float(mass[saving < -1e-9].sum())
        debt = float(np.sum(mass * np.maximum(-saving, 0.)))
        rows.append(dict(case=case, scope=scope, current_branch=name, branch_mass=denominator,
            next_saving_debt_mass=indebted, debt_share_within_branch=indebted/denominator if denominator else None,
            next_net_financial_debt=debt, debt_per_all_households=debt/total,
            debt_per_branch_household=debt/denominator if denominator else None))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-script-sha', required=True)
    parser.add_argument('--seconds', type=int, default=300)
    args = parser.parse_args()
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm required')
    require(0 < args.seconds <= 300 and int(os.environ.get('SLURM_CPUS_PER_TASK', '1')) == 1, 'Budget differs')
    require(sha(__file__) == args.expected_script_sha, 'Postprocessor source pin differs')
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('Postprocessor time budget reached')))
    signal.alarm(args.seconds)
    started = time.monotonic()
    root, out = args.run_root.resolve(), args.output.resolve()
    launch, completion = read(root / 'launch.json'), read(root / 'completed.json')
    require(completion['status'] == 'passed' and completion['lifecycle_solves'] == 3, 'Three-case run not successful')
    plan = launch['plan']
    require(sha(args.plan) == launch['plan_sha256'] and read(args.plan) == plan,
            'Immutable run plan file does not match successful launch')
    require(plan['reference_label'] == LABEL and completion['plan_sha256'] == launch['plan_sha256'], 'Run plan identity differs')
    for item in plan['files']:
        require(sha(item['path']) == item['sha256'], 'Pinned run source/input differs: ' + item['path'])
    require(not out.exists(), 'Output must be new')
    out.mkdir(parents=True)
    base = load(plan['base_driver_path'], '_credit_summary_reference')
    manifest, _, _, _, prepared, reference = base.authenticate(out)
    import numpy as np
    cal, model = prepared.rt['primitive'].pf.calendar, prepared.rt['model']
    grid_spec = read(plan['credit_grid_path'])
    grid = np.asarray(grid_spec['grid'])
    idx = np.asarray(grid_spec['old_indices'], dtype=np.intp)
    require(np.array_equal(grid[idx], reference['b_grid']), 'Grid embedding differs from reference')
    inherited = np.zeros((len(grid),) + reference['stationary_g_pre'].shape[1:], dtype=reference['stationary_g_pre'].dtype)
    inherited[idx] = reference['stationary_g_pre']
    case_rows, debts, summary, pins = {}, [], {}, {}
    for name in ('grid_control', 'credit'):
        receipt_path = root / name / 'receipt.json'
        record = next(row for row in completion['completed'] if row['case'] == name)
        require(sha(receipt_path) == record['receipt_sha256'], 'Successful case receipt changed')
        receipt = read(receipt_path)
        require(receipt['status'] == 'passed' and receipt['plan_sha256'] == launch['plan_sha256'] and
                receipt['reference_manifest_sha256'] == base.MANIFEST_SHA, 'Case identity differs')
        checkpoint = receipt['checkpoint']
        require(sha(checkpoint['path']) == checkpoint['sha256'], 'Case checkpoint changed')
        observer_path = root / name / 'observers.json'
        require(sha(observer_path) == receipt['artifact_hashes']['observers.json'], 'Saved observers changed')
        with gzip.open(checkpoint['path'], 'rb') as stream:
            packet = pickle.load(stream)
        P = copy.deepcopy(packet['parameters'])
        require(np.array_equal(packet['b_grid'], grid), 'Cases do not use identical refined grids')
        P.native_inherited_distribution_evidence_dir = str(out / (name+'_inherited_diagnostics'))
        shared = model.precompute_shared(P, grid)
        policy = packet['evaluation'].policy
        counter = cal.SolveCounter()
        impact = cal.evaluate_period(policy.price, inherited, P, grid, shared, counter,
            supply_rule=packet['supply_rule'], supplied_policy=policy)
        require(counter.total == 0, 'Postprocessor unexpectedly solved a household problem')
        require(np.array_equal(impact.g_pre, inherited) and impact.feasibility_projection_mass == 0,
                'Impact altered inherited states')
        rows, flows = birth_cells(impact, P, grid, model)
        for actual, key in zip(flows, ('first_births', 'second_births', 'third_bin_entries')):
            require(abs(actual-receipt['baseline_state_impact_summary'][key]) <= 2e-10,
                    'Birth decomposition differs from certified receipt: ' + key)
        case_rows[name] = rows
        debts.extend(debt_rows(impact, P, grid, model, name, 'impact'))
        debts.extend(debt_rows(packet['evaluation'], packet['parameters'], grid, model, name, 'cohort'))
        mean_age = read(observer_path)['fertility']['uniform_birth_time']['moments']['period_mean_age_first_birth']
        summary[name] = dict(completed_fertility=receipt['completed_fertility'],
            cohort_period_mean_age_first_birth=mean_age, impact_births=float(impact.births), impact_births_by_order=flows)
        pins[name] = dict(receipt_sha256=sha(receipt_path), checkpoint=checkpoint, observers_sha256=sha(observer_path))
        del impact, packet, P, shared, policy
    require(set(case_rows['grid_control']) == set(case_rows['credit']), 'Decomposition cells differ')
    differences = []
    for key, baseline in case_rows['grid_control'].items():
        credit = case_rows['credit'][key]
        require(abs(baseline['at_risk_mass']-credit['at_risk_mass']) <= 2e-10, 'Inherited risk pool differs')
        row = {k: v for k, v in baseline.items() if k != 'birth_mass'}
        row.update(reference_births=baseline['birth_mass'], credit_births=credit['birth_mass'],
            birth_difference=credit['birth_mass']-baseline['birth_mass'])
        differences.append(row)
    require(abs(sum(r['birth_difference'] for r in differences)-
        (summary['credit']['impact_births']-summary['grid_control']['impact_births'])) <= 2e-10, 'Differences do not add up')
    base.table(out / 'impact_birth_decomposition.csv', differences)
    base.table(out / 'next_saving_debt.csv', debts)
    base.write(out / 'summary.json', dict(status='passed', reference_label=LABEL, lifecycle_solves=0,
        supplied_policy_period_evaluations=2, elapsed_seconds=time.monotonic()-started,
        script_sha256=sha(__file__), run_plan_sha256=launch['plan_sha256'], case_pins=pins, cases=summary,
        interpretation={'births':'First and subsequent birth flows grouped by inherited age, tenure and net financial wealth; identical inherited risk pools.',
            'debt':'Negative next-saving net financial positions; levels, not new borrowing flows or identified mortgage balances.',
            'mean_age':'Saved flow-weighted first-birth age in the recomputed cohort; not completed-cohort timing or transition.'}))
    signal.alarm(0)


if __name__ == '__main__':
    main()
