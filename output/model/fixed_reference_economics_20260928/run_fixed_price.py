#!/usr/bin/env python3
"""Three-case, fixed-preference price experiment; execute only on Torch.

The plan must pin this file's SHA256 and the existing reference manifest.
Two fresh-process exact controls precede one +10% price/rent case. No
calibration evaluator or normalization method is called. All times include
authentication, solving and reporting; timeout/failure stops the whole loop.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import os
import pickle
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
LABEL = '2007 stationary reference — block0506, September 28 verified export'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
MANIFEST = ROOT / 'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
CASES = (('control_1', 1.0), ('control_2', 1.0), ('price_110', 1.1))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def table(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def csv_rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def serialized(value):
    """Same compact array identity convention as the reference authenticator."""
    if hasattr(value, 'shape') and getattr(value, 'size', 0) > 1024:
        h = hashlib.sha256()
        for start in range(0, value.size, 65536):
            h.update(value.flat[start:start + 65536].tobytes())
        return dict(serialized_array=True, shape=list(value.shape), dtype=str(value.dtype),
                    size=int(value.size), sha256_c_order_bytes=h.hexdigest())
    if hasattr(value, 'tolist'):
        return serialized(value.tolist())
    if isinstance(value, dict):
        return {str(k): serialized(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialized(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {'nonfinite_float': repr(value)}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError('Unsupported serialized value: ' + str(type(value)))


def plan_contract(path):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(),
            'Torch Slurm execution required, including preflight')
    plan = read(path)
    require(plan['schema'] == 'block0506_fixed_price_v1', 'Wrong plan schema')
    require(plan['reference_label'] == LABEL, 'Wrong reference label')
    require(plan['driver_sha256'] == sha(__file__), 'Driver pin mismatch')
    require(plan['reference_manifest_sha256'] == MANIFEST_SHA == sha(MANIFEST), 'Manifest pin mismatch')
    require(plan['cases'] == [dict(name=k, price_factor=v) for k, v in CASES], 'Case loop differs')
    require(plan['maximum_lifecycle_solves'] == 3, 'Solve budget differs')
    require(0 < plan['case_seconds'] <= 360 and 0 < plan['total_seconds'] <= 1200, 'Time budget differs')
    require(plan['threads'] == 1 and plan['memory_gib'] == 16, 'Resource contract differs')
    return plan


def progress(output, phase, **more):
    write(output / 'progress.json', dict(reference_label=LABEL, phase=phase,
                                       time_epoch=time.time(), **more))


def authenticate(output):
    """Register immutable runtime classes, then load the actual frozen export."""
    manifest = read(MANIFEST)
    require(manifest['label'] == LABEL and not manifest['economic_contract']['counterfactual_renormalization'],
            'Reference economic contract differs')
    for key in ('contract', 'objective', 'source_manifest', 'source_contract', 'native_ancestry_contract'):
        pin = manifest[key]
        require(sha(pin['path']) == pin['sha256'], 'Reference pin differs: ' + key)
    contract = read(manifest['contract']['path'])
    for pin in contract['files'].values():
        require(sha(pin['path']) == pin['sha256'], 'Controller ancestry pin differs: ' + pin['path'])
    objective = read(manifest['objective']['path'])
    fingerprint = hashlib.sha256(json.dumps(objective['target_rows'], sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    require(fingerprint == manifest['target_weight_fingerprint'], 'Target/weight fingerprint differs')
    export = Path(manifest['local_export'])
    for name, digest in manifest['artifact_hashes'].items():
        require(sha(export / name) == digest, 'Export artifact differs: ' + name)
    checkpoint = export / 'initial_state.pkl.gz'
    require(sha(checkpoint) == manifest['checkpoint']['sha256'], 'Block0506 checkpoint differs')
    require(csv_rows(export / 'target_fit.csv') == manifest['full_target_table'], 'Reference target table differs')
    require(csv_rows(export / 'parameters.csv') == manifest['full_parameter_table'], 'Reference parameter table differs')
    sys.path[:0] = [str(ROOT / 'code/model/tools'),
                   str(ROOT / 'tmp/e5f_overnight_local_20260927/portable/tools_v4')]
    import e5f_evening_calibration_runtime as runtime
    prepared = runtime.setup(dict(contract, objective=manifest['objective']), objective,
                             output / 'runtime_preparation')
    with gzip.open(checkpoint, 'rb') as stream:
        reference = pickle.load(stream)
    require(serialized(vars(reference['parameters'])) == manifest['actual_serialized_parameters'],
            'Serialized block0506 parameter identity differs')
    return manifest, contract, objective, runtime, prepared, reference


def actual_parameters(prepared, P, grid):
    result = dict(prepared.tax.actual_parameters(P), delta_alpha_jump=P.delta_alpha_jump,
                  child_benefit_curvature=P.child_benefit_curvature, tenure_choice_kappa=P.tenure_choice_kappa)
    result.update(psi_child=P.psi_child, child_benefit_CRRA_coefficient=(1-P.child_benefit_curvature)*P.psi_child,
        theta1=P.theta1, sigma=P.sigma, alpha_cons=P.alpha_cons, delta_alpha=P.delta_alpha, h_P=0.,
        utility_reference_rent=P.utility_reference_rent, q_annual=(1+P.q)**(1/P.period_years)-1,
        financed_share=P.phi[0], housing_supply_elasticity=P.xi_supply[0], payroll_tax=P.tau_pay,
        pension_period=P.pension, annual_depreciation=prepared.ancestor.ANNUAL_DEP,
        period_depreciation=P.delta, annual_property_tax=prepared.ancestor.ANNUAL_PROPERTY_TAX,
        period_property_tax=P.tau_H, selling_cost=P.psi, rental_cap=P.hR_max,
        wealth_grid_nodes=len(grid), income_states=len(P.z_grid))
    return result


def gates(packet, prepared, output, *, stationary):
    import numpy as np
    from e5f_social_security import fiscal_accounts
    P, grid, ev, sd = (packet[k] for k in ('parameters', 'b_grid', 'evaluation', 'shared'))
    rt = prepared.rt
    gate = sys.modules['e5f_evening_calibration_runtime'].require_abs_gate
    require(np.array_equal(ev.g_pre, packet['stationary_g_pre']), 'Inherited pre-choice mass altered')
    gate(ev.feasibility_projection_mass, 0., 'Feasibility projection')
    for name in ('g_pre', 'g_post_fertility', 'g_current'):
        mass = getattr(ev, name)
        require(np.isfinite(mass).all() and np.min(mass) >= 0., 'Invalid mass: ' + name)
        gate(float(mass.sum() - ev.g_pre.sum()), 2e-10, 'Within-period mass conservation')
    budget = rt['primitive'].dated_budget(ev, P, sd, grid, float(P.user_cost_rate*ev.policy.price[0]))
    gate(budget['budget_excess_mass'], 2e-10, 'Household budget')
    purchase = rt['accounting'].audit_purchase_accounting(ev, P, sd, grid, rt['model'])
    gate(purchase['maximum_occupied_transaction_wealth_error'], 1e-9, 'Transaction wealth')
    for key, value in purchase.items():
        if key.endswith('violation_mass') or key in ('transaction_outside_grid_mass', 'negative_estate_exposure_mass', 'saving_outside_grid_mass'):
            gate(value, 2e-10, key)
    estate = prepared.estate.audit(ev, P, grid)
    sys.modules['e5f_evening_calibration_runtime'].require_no_negative_estates(estate)
    require(estate['status'] == 'funded', 'Estate entry funding gate failed')
    arrays = rt['audit'].policy_array_audit(packet, output)
    require(arrays['occupied_negative_steps'] == 0, 'Occupied value monotonicity gate failed')
    require(all(not v['nonfinite'] and 0 <= v['minimum'] <= v['maximum'] <= 1
                for v in arrays['probabilities'].values()), 'Probability gate failed')
    result = dict(household_budget=budget, purchase=purchase, estate=estate, policy_arrays=arrays,
        fiscal=fiscal_accounts(ev.g_current, P), feasibility_projection_mass=float(ev.feasibility_projection_mass),
        housing_market_clearing_required=False, relative_market_residual=float(ev.relative_market_residual),
        housing_demand=np.asarray(ev.demand_by_loc).tolist(), housing_supply=np.asarray(ev.supply_by_loc).tolist(),
        housing_excess=(ev.demand_by_loc-ev.supply_by_loc).tolist())
    # Fixed fiscal inputs are still certified; a failed accounting gate is not hidden.
    result['fiscal_certificate'] = rt['certify_initial_pension'](ev.g_current, P,
        marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    if stationary:
        operator = rt['primitive'].pf.transition.operator_gates(packet['solution'], ev.policy,
            packet['stationary_g_pre'], P, grid, sd)
        for key in ('stationary_post_fertility_nesting_l1', 'one_step_constant_path_nesting_l1',
                    'mature_flow_abs_error', 'birth_flow_abs_error', 'topcode_adjusted_birth_flow_abs_error'):
            gate(operator[key], 5e-9, key)
        gate(operator['zero_entry_mass_accounting_residual'], 2e-8, 'Mass accounting')
        result['stationary_operator'] = operator
    write(output / 'gates.json', rt['primitive'].pf.calendar.jsonable(result))
    return result


def aggregates(ev, P, grid):
    import numpy as np
    from e5f_overnight_estate_audit import policy_mass_branches
    cal = sys.modules['run_dynamic_population_transition']
    g = ev.g_current
    mass = float(g.sum())
    branches = policy_mass_branches(ev, P)
    by_order = [float(ev.g_post_fertility[..., n:, :].sum()-ev.g_pre[..., n:, :].sum())
                for n in range(1, P.n_parity)]
    require(abs(sum(by_order)-float(ev.births)) < 2e-10, 'Birth-order flows do not add up')
    return dict(household_mass=mass, births=float(ev.births), births_per_household=float(ev.births)/mass,
        first_births=by_order[0], second_births=by_order[1], third_bin_entries=by_order[2],
        ownership_rate=float(g[:, 1:].sum())/mass,
        rooms_per_household=float(cal.housing_demand_by_location(g, ev.policy.hR_pol, P).sum())/mass,
        nonhousing_consumption_per_household=sum(float(np.sum(m*c)) for m, _, c in branches)/mass,
        next_liquid_assets_per_household=sum(float(np.sum(m*b)) for m, b, _ in branches)/mass,
        current_liquid_assets_per_household=float(np.sum(g*grid[:, None, None, None, None, None, None]))/mass,
        owner_stayer_mass=float(ev.g_stay_distribution.sum()),
        reference_label=LABEL)


def exact_control(reference, candidate, fits, manifest, prepared, output):
    import e5f_current_transition_runtime as native
    check = native.compare_arrays(reference, candidate)
    write(output / 'control_arrays.json', check)
    bad = {k: v for k, v in check['arrays'].items() if v.get('status') != 'compared' or not v.get('exact') or not v.get('finite')}
    require(not bad, 'Exact reference array replay failed: ' + ', '.join(list(bad)[:15]))
    expected = manifest['full_target_table']
    require(len(fits) == 14 and [r['moment'] for r in fits] == [r['moment'] for r in expected], 'Control target rows differ')
    for actual, old in zip(fits, expected):
        for key in ('target', 'model', 'gap', 'weight', 'loss_contribution'):
            require((actual[key] == '' and old[key] == '') or
                    (actual[key] != '' and old[key] != '' and float(actual[key]) == float(old[key])),
                    'Exact control target differs: ' + actual['moment'] + '/' + key)
    require(float(candidate['evaluation'].relative_market_residual) <= 2e-4, 'Control housing gate failed')
    return dict(status='passed', all_numeric_arrays_exact=True, arrays=check['array_count'],
                all_14_fit_rows_exact=True, all_31_parameters_exact=True)


def child(args, plan):
    import numpy as np
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    # Internal child mode cannot bypass the two-control order.
    launch = read(output.parent / 'launch.json')
    require(launch['plan_sha256'] == sha(args.plan) and launch['plan'] == plan,
            'Child lacks matching controller launch contract')
    index = [name for name, _ in CASES].index(args.case)
    if index:
        previous = read(output.parent / 'latest_completed.json')['completed']
        require([r['case'] for r in previous] == [name for name, _ in CASES[:index]],
                'Preceding controls incomplete')
        for record in previous:
            receipt_path = output.parent / record['case'] / 'receipt.json'
            receipt = read(receipt_path)
            require(sha(receipt_path) == record['receipt_sha256'] and receipt['control']['status'] == 'passed'
                    and receipt['plan_sha256'] == sha(args.plan), 'Preceding control pin/gate differs')
    progress(output, 'authenticate')
    manifest, contract, objective, runtime, prepared, reference = authenticate(output)
    rt = prepared.rt
    model, cal = rt['model'], rt['primitive'].pf.calendar
    P = copy.deepcopy(reference['parameters'])
    grid = np.asarray(reference['b_grid']).copy()
    # This logging destination is not an economic primitive. Authenticate the
    # original serialized P first, then keep new logs out of the reference case.
    inherited_log_destination = P.native_inherited_distribution_evidence_dir
    P.native_inherited_distribution_evidence_dir = str(output / 'inherited_state_diagnostics')
    public_before = serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    require(P.native_due_stayer_credit and P.native_exact_inherited_distribution, 'Reference DUE/inherited-state contract absent')
    actual = actual_parameters(prepared, P, grid)
    require(len(actual) == 31, 'Complete parameter table unavailable')
    params = copy.deepcopy(manifest['full_parameter_table'])
    for row in params:
        require(float(actual[row['parameter']]) == float(row['estimate']), 'Reference parameter mismatch: ' + row['parameter'])
        row['status'] = 'Fixed inherited reference value; no estimation or fertility renormalization'
    factor = dict(CASES)[args.case]
    price = factor*np.asarray(reference['solution'].p_eq)
    sd = model.precompute_shared(P, grid)
    progress(output, 'one_lifecycle_solve', price_factor=factor, deadline_epoch=args.deadline)
    require(time.time() < args.deadline, 'Case deadline reached before solve')
    started = time.monotonic()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=sd, verbose=False, fast_stats=False)
    elapsed = time.monotonic()-started
    require(time.time() < args.deadline, 'Case deadline exceeded in solve')
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(reconstruction['stationary_post_fertility_nesting_l1'], 5e-9, 'Stationary reconstruction')
    runtime.require_abs_gate(reconstruction['stationary_feasibility_projection_mass'], 0., 'Stationary projection')
    supply = cal.HousingSupplyRule('static-elastic', float(price[0]),
        float(P.H0[0]*(P.user_cost_rate*price[0]/P.r_bar[0])**P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol, evaluation=ev,
        stationary_g_pre=pre, supply_rule=supply, demographic_seed=reference.get('demographic_seed'))
    progress(output, 'cohort_audit', lifecycle_solve_seconds=elapsed)
    cohort_gates = gates(packet, prepared, output, stationary=True)
    fertility = {p: rt['observe_initial_fertility'](ev, P, age_projection=p)
                 for p in ('uniform_birth_time', 'constant_post_cell')}
    housing = rt['observe_initial_housing_wealth'](ev, P, grid, sd, diagnostic_enabled=True,
        age_projection='uniform_within_age_cell', diagnostic_allow_family_proxies=True,
        include_wealth=True, include_birth_response=True)
    recent = rt['observe_recent_parent_flow'](ev, P, diagnostic_enabled=True,
        snapshot=rt['SNAPSHOT'], age_projection=rt['AGE_PROJECTION'], diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=args.case, reference_checkpoint_sha256=manifest['checkpoint']['sha256']))
    completed = float(rt['chain'].extract_moments(sol, P)['tfr'])
    fits = runtime.score_targets(objective, fertility, housing, recent['model_value'], completed)
    require(len(fits) == 14, 'All 14 fit rows required')
    control = exact_control(reference, packet, fits, manifest, prepared, output) if factor == 1. else None
    # Change descriptive roles only after exact control comparison.
    for row in fits:
        if row['moment'] == 'initial_normalization':
            row['role'] = 'reference replacement benchmark; not imposed after shock'
    table(output / 'target_fit.csv', fits)
    table(output / 'parameters.csv', params)
    write(output / 'observers.json', cal.jsonable(dict(fertility=fertility, housing_wealth=housing, recent_parent=recent)))
    cohort_summary = aggregates(ev, P, grid)
    # Retain each cohort packet's branch-owned caches while impact is evaluated.
    impact_P = copy.deepcopy(P)
    impact_sd = model.precompute_shared(impact_P, grid)
    impact = cal.evaluate_period(price, reference['stationary_g_pre'], impact_P, grid, impact_sd,
        cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    impact_output = output / 'baseline_state_impact'
    impact_output.mkdir()
    impact_packet = dict(parameters=impact_P, b_grid=grid, shared=impact_sd, evaluation=impact,
                         stationary_g_pre=reference['stationary_g_pre'])
    impact_gates = gates(impact_packet, prepared, impact_output, stationary=False)
    impact_summary = aggregates(impact, impact_P, grid)
    if factor == 1.:
        for name in ('g_pre', 'g_post_fertility', 'g_current', 'g_stay_distribution'):
            require(np.array_equal(getattr(impact, name), getattr(reference['evaluation'], name)),
                    'Control impact array differs: ' + name)
    write(impact_output / 'summary.json', impact_summary)
    progress(output, 'standard_17_plot_rendering')
    rt['audit'].standard_diagnostics(packet, output, validate_production_young=False)
    plots = sorted(p.name for p in (output / 'standard_diagnostics').glob('*.png'))
    require(plots == sorted(manifest['standard_diagnostic_names']), 'Standard 17-plot set differs')
    if factor == 1.:
        for name in plots:
            require(sha(output / 'standard_diagnostics' / name) == manifest['artifact_hashes']['standard_diagnostics/' + name],
                    'Control diagnostic hash differs: ' + name)
    require(serialized({k: v for k, v in vars(P).items() if not k.startswith('_')}) == public_before,
            'Public economic/numerical parameter mutated during case')
    runtime.verify_sources(dict(contract, objective=manifest['objective']))
    require(sha(MANIFEST) == MANIFEST_SHA, 'Reference manifest changed during case')
    checkpoint = None
    if factor != 1.:
        progress(output, 'save_new_shock_checkpoint')
        # Reference binaries are never rewritten or copied.
        checkpoint_path = output / 'conditional_cohort_state.pkl.gz'
        with gzip.open(checkpoint_path, 'wb', compresslevel=1) as stream:
            pickle.dump(packet, stream, protocol=5)
        checkpoint = dict(path=str(checkpoint_path), sha256=sha(checkpoint_path))
    require(time.time() < args.deadline, 'Case deadline exceeded in reporting')
    require(sha(__file__) == plan['driver_sha256'] and sha(args.plan) == launch['plan_sha256'],
            'Driver or plan changed during case')
    receipt = dict(status='passed', reference_label=LABEL, case=args.case, price_factor=factor,
        price=price.tolist(), rent=(P.user_cost_rate*price).tolist(), reference_manifest_sha256=MANIFEST_SHA,
        reference_checkpoint=manifest['checkpoint'], source_manifest=manifest['source_manifest'],
        target_weight_fingerprint=manifest['target_weight_fingerprint'], driver_sha256=sha(__file__),
        plan_sha256=sha(args.plan), lifecycle_solves=1, lifecycle_solve_seconds=elapsed,
        diagnostic_destination_override=dict(reference=inherited_log_destination,
            experiment=P.native_inherited_distribution_evidence_dir, economic_change=False),
        fixed_psi=float(P.psi_child), completed_fertility=completed,
        replacement_gap=float(sol.adult_entry_stationary_relative_gap),
        replacement_signed_residual=float(sol.adult_entry_stationary_residual),
        target_comparison_loss=sum(float(r['loss_contribution']) for r in fits if r['loss_contribution'] != ''),
        normalization_performed=False, control=control, checkpoint=checkpoint,
        cohort_summary=cohort_summary, baseline_state_impact_summary=impact_summary,
        cohort_gates=cohort_gates, impact_gates=impact_gates, standard_plot_count=len(plots),
        artifact_hashes={name: sha(output / name) for name in
            ['target_fit.csv', 'parameters.csv', 'observers.json', 'gates.json'] +
            ['standard_diagnostics/' + name for name in plots]},
        interpretation='Prescribed permanent prices and unchanged fiscal inputs. Impact uses exact baseline pre-choice states. Cohort distribution uses normalized entry; neither is a cleared equilibrium, renewed demographic steady state or transition.',
        economic_changes=['Asset price and implied rent +10%'] if factor != 1. else [],
        limitations=['Standard policy plots include buyer-conditional full-grid states.',
            'Inherited lifecycle_2023.csv contains conditional stationary-cohort diagnostics, not a 2023 transition.',
            'Estate recipient and physical/financial counterparty closure remains provisional.'])
    write(output / 'receipt.json', cal.jsonable(receipt))
    progress(output, 'complete')


def controller(args, plan):
    started = time.time()
    deadline = started + plan['total_seconds']
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / 'launch.json', dict(reference_label=LABEL, plan=plan, plan_sha256=sha(args.plan),
        started_epoch=started, deadline_epoch=deadline, slurm_job=os.environ['SLURM_JOB_ID']))
    completed = []
    write(args.output / 'latest_completed.json', dict(reference_label=LABEL, completed=[],
          lifecycle_solves=0, remaining_cases=3))
    for name, factor in CASES:
        require(read(args.plan) == plan and sha(__file__) == plan['driver_sha256'],
                'Plan or driver changed during loop')
        require(time.time() < deadline, 'Global deadline reached; no further case')
        case_output = args.output / name
        case_deadline = min(deadline, time.time()+plan['case_seconds'])
        command = [sys.executable, str(Path(__file__).resolve()), '--plan', str(args.plan),
                   '--output', str(case_output), '--case', name, '--deadline', str(case_deadline)]
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   NUMEXPR_NUM_THREADS='1', NUMBA_NUM_THREADS='1', MPLBACKEND='Agg')
        with (args.output / (name+'.log')).open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env,
                                       start_new_session=True)
            try:
                while process.poll() is None:
                    progress(args.output, 'case_running', case=name, pid=process.pid,
                             elapsed_seconds=time.time()-started, case_deadline_epoch=case_deadline,
                             completed_cases=len(completed))
                    if time.time() >= case_deadline:
                        raise TimeoutError('Case/global deadline reached: ' + name)
                    time.sleep(2)
                require(process.returncode == 0, 'Case failed; stop without retry: ' + name)
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
        receipt = read(case_output / 'receipt.json')
        require(receipt['status'] == 'passed' and receipt['lifecycle_solves'] == 1, 'Missing successful case receipt')
        if factor == 1.:
            require(receipt['control']['status'] == 'passed', 'Shock blocked: control replay unavailable')
        completed.append(dict(case=name, receipt_sha256=sha(case_output / 'receipt.json'),
                              cohort=receipt['cohort_summary'], impact=receipt['baseline_state_impact_summary']))
        write(args.output / 'latest_completed.json', dict(reference_label=LABEL, completed=completed,
              lifecycle_solves=len(completed), remaining_cases=3-len(completed)))
    comparison = {}
    for scope in ('cohort', 'impact'):
        baseline, shocked = completed[0][scope], completed[2][scope]
        require(completed[1][scope] == baseline, 'The two control summaries differ')
        comparison[scope] = {key: dict(reference=value, price_110=shocked[key],
            difference=shocked[key]-value) for key, value in baseline.items()
            if isinstance(value, (int, float))}
    write(args.output / 'comparison.json', dict(reference_label=LABEL, **comparison,
          interpretation='Impact holds baseline pre-choice mass fixed; cohort allows normalized-entry composition to change. Neither clears housing or supplies an equilibrium transition.'))
    require(time.time() < deadline, 'Global deadline exceeded during finalization')
    write(args.output / 'completed.json', dict(status='passed', reference_label=LABEL,
          completed=completed, lifecycle_solves=3, elapsed_seconds=time.time()-started,
          plan_sha256=sha(args.plan), driver_sha256=sha(__file__), normalization_performed=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', choices=[name for name, _ in CASES])
    parser.add_argument('--deadline', type=float)
    args = parser.parse_args()
    args.plan = args.plan.resolve()
    args.output = args.output.resolve()
    plan = plan_contract(args.plan)
    require(not args.output.exists(), 'Versioned output must not already exist')
    try:
        if args.case:
            require(args.deadline is not None and time.time() < args.deadline <= time.time()+plan['case_seconds'],
                    'Child requires bounded absolute deadline')
            child(args, plan)
        else:
            require(args.deadline is None, 'Controller sets its own immutable deadline')
            controller(args, plan)
    except BaseException as error:
        if args.output.exists():
            write(args.output / 'failure.json', dict(reference_label=LABEL, status='failed',
                  error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc(),
                  failure_ledger=serialized(getattr(error, 'ledger', getattr(error, 'audit', None))),
                  time_epoch=time.time(), retries=0))
        raise


if __name__ == '__main__':
    main()
