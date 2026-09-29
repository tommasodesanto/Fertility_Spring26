#!/usr/bin/env python3
"""Bounded, fixed-price credit comparison; Torch only, with audited overlay.

This controller deliberately has no default credit implementation. The lead's
immutable adapter must supply install(), audit_purchase_accounting(), and
audit_solvency(); its source and a passing review receipt must be pinned in the
plan. One fresh exact control precedes matched-grid baseline and credit solves. No recalibration,
equilibrium iteration, fertility normalization, or automatic retry occurs.

Adapter API_VERSION is ``block0506_credit_adapter_v1``. Its keyword-only
install(prepared, reference, P, grid, output, overlay_root, credit_mode) may
modify only the explicitly declared P fields and runtime solver functions.
audit_purchase_accounting(ev, P, sd, grid, model) preserves transaction/grid/
estate checks while replacing artificial credit limits. audit_solvency(packet,
prepared, output, *, stationary) certifies reachable continuation and repayment
and reports finite-grid exclusions separately. All three hooks are required.
"""
from __future__ import annotations

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
import subprocess
import sys
import time
import traceback
from types import SimpleNamespace

LABEL = '2007 stationary reference — block0506, September 28 verified export'
ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
MANIFEST = ROOT / 'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
CASES = ('control', 'grid_control', 'credit')


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


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def verify_plan(path):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(),
            'All execution and validation require Torch Slurm')
    plan = read(path)
    require(plan['schema'] == 'block0506_fixed_price_credit_v1', 'Unexpected plan schema')
    require(plan['reference_label'] == LABEL, 'Reference label differs')
    require(plan['driver_sha256'] == sha(__file__), 'Driver pin differs')
    require(plan['reference_manifest_sha256'] == MANIFEST_SHA == sha(MANIFEST), 'Manifest differs')
    require(plan['cases'] == list(CASES) and plan['maximum_lifecycle_solves'] == 3,
            'Exactly one exact control, matched-grid baseline, and credit solve are authorized')
    require(0 < plan['case_seconds'] <= 600 and 0 < plan['total_seconds'] <= 1800,
            'Explicit case/total time budget required, maximum 600/1800 seconds')
    require(plan['threads'] == 1 and plan['memory_gib'] == 16, 'Resource contract differs')
    require(plan['price_factor'] == 1 and plan['renormalize_fertility'] is False,
            'This experiment fixes prices and preferences')
    require(plan['credit_mode'] and isinstance(plan['credit_mode'], str), 'Named credit mode required')
    root = Path(plan['overlay_root']).resolve()
    adapter = Path(plan['adapter_path']).resolve()
    require(adapter.is_relative_to(root), 'Adapter must belong to the immutable overlay')
    pins = plan['files']
    for item in pins:
        require(sha(item['path']) == item['sha256'], 'Executed source/input pin differs: ' + item['path'])
    pinned = {str(Path(item['path']).resolve()) for item in pins}
    require(len(pinned) == len(pins), 'Duplicate source pins')
    for item in (adapter, Path(plan['base_driver_path']).resolve(), Path(plan['review_receipt_path']).resolve()):
        require(str(item) in pinned, 'Missing required source/review pin: ' + str(item))
    require(plan.get('credit_grid_path'), 'Matched-grid comparison requires a refined-grid input')
    require(str(Path(plan['credit_grid_path']).resolve()) in pinned, 'Credit grid input is not pinned')
    # No unpinned overlay code may silently enter through a sibling import.
    for item in root.rglob('*.py'):
        require(str(item.resolve()) in pinned, 'Unpinned overlay Python source: ' + str(item))
    review = read(plan['review_receipt_path'])
    require(review['status'] == 'passed' and review['credit_mode'] == plan['credit_mode'],
            'Audited credit implementation is not ready')
    require(review['adapter_sha256'] == sha(adapter), 'Review belongs to another adapter')
    require(review['reference_manifest_sha256'] == MANIFEST_SHA, 'Review used another calibration')
    require(review['lifetime_solvency_and_repayment_retained'] is True,
            'Review must certify retained solvency and debt repayment')
    require(review['artificial_limits_removed'] is True,
            'Review must certify the advertised removal of artificial limits')
    require(review['unresolved_blockers'] == [], 'Unresolved credit implementation blocker')
    return plan


def progress(output, phase, **more):
    write(output / 'progress.json', dict(reference_label=LABEL, phase=phase,
                                       time_epoch=time.time(), **more))


def credit_audits(base, packet, prepared, output, adapter, stationary):
    """Keep universal native gates; replace only the artificial credit ledger."""
    original = prepared.rt['accounting']

    def checked_purchase(ev, P, sd, grid, model):
        result = adapter.audit_purchase_accounting(ev, P, sd, grid, model)
        required = ('maximum_occupied_transaction_wealth_error', 'transaction_outside_grid_mass',
                    'negative_estate_exposure_mass', 'saving_outside_grid_mass')
        require(all(key in result for key in required), 'Credit purchase audit omitted a universal gate')
        require(result.get('artificial_credit_limits_enforced') is False,
                'Credit audit must explicitly replace artificial credit tests')
        return result

    prepared.rt['accounting'] = SimpleNamespace(audit_purchase_accounting=checked_purchase)
    try:
        result = base.gates(packet, prepared, output, stationary=stationary)
    finally:
        prepared.rt['accounting'] = original
    solvency = adapter.audit_solvency(packet, prepared, output, stationary=stationary)
    require(solvency['status'] == 'passed', 'Independent solvency audit failed')
    require(solvency['lifetime_solvency_and_repayment_retained'] is True,
            'Solvency audit did not establish repayment')
    require(0. <= solvency['occupied_unreachable_continuation_mass'] <= 2e-10,
            'Occupied policies enter infeasible continuation support')
    require(0. <= solvency['negative_estate_exposure_mass'] <= 2e-10, 'Death estate solvency failed')
    require(isinstance(solvency['finite_grid'], dict) and solvency['finite_grid'],
            'Finite-grid restrictions must be reported explicitly')
    result['credit_solvency'] = solvency
    base.write(output / 'gates.json', prepared.rt['primitive'].pf.calendar.jsonable(result))
    return result


def aggregates(base, ev, P, grid, model):
    # The reference helper reports owner-stayer mass. With the DUE-specific
    # credit rule disabled, evaluate_period does not populate that cache;
    # reconstruct the same economically meaningful flow solely for reporting.
    if getattr(ev, 'g_stay_distribution', None) is None:
        ev = copy.copy(ev)
        ev.g_stay_distribution = model.realize_stayer_cross_section(
            ev.g_post_fertility, ev.policy.loc_probs, ev.policy.tenure_choice, ev.policy.tenure_probs)
    return base.aggregates(ev, P, grid)


def embed_credit_grid(base, P, reference, path):
    """Add numerical knots without moving or redistributing any inherited mass."""
    import numpy as np
    spec = read(path)
    old = np.asarray(reference['b_grid'])
    grid = np.asarray(spec['grid'], dtype=old.dtype)
    raw_indices = np.asarray(spec['old_indices'])
    require(raw_indices.ndim == 1 and raw_indices.shape == old.shape and
            np.issubdtype(raw_indices.dtype, np.integer), 'Old grid indices must be aligned integers')
    indices = raw_indices.astype(np.intp)
    require(grid.ndim == 1 and len(grid) > len(old) and np.isfinite(grid).all() and
            np.all(np.diff(grid) > 0), 'Credit grid must add finite, strictly increasing knots')
    require(grid[0] == old[0] and grid[-1] == old[-1], 'Credit grid endpoints changed')
    require(np.all(np.diff(indices) > 0) and indices[0] == 0 and indices[-1] == len(grid)-1,
            'Invalid original-node embedding')
    require(np.array_equal(grid[indices], old), 'Refined grid does not exactly contain original nodes')
    require(np.array_equal(P.fixed_reference_entry_grid, old), 'Original entry grid differs')
    weights = np.asarray(P.fixed_reference_entry_conditional)
    require(weights.ndim == 2 and weights.shape[0] == len(old) and np.isfinite(weights).all() and
            np.all(weights >= 0), 'Invalid reference conditional entry law')
    expanded = np.zeros((len(grid),) + weights.shape[1:], dtype=weights.dtype)
    expanded[indices] = weights
    original_pre = np.asarray(reference['stationary_g_pre'])
    require(original_pre.shape[0] == len(old), 'Inherited mass is not on original grid')
    impact_pre = np.zeros((len(grid),) + original_pre.shape[1:], dtype=original_pre.dtype)
    impact_pre[indices] = original_pre
    added = np.ones(len(grid), dtype=bool)
    added[indices] = False
    require(np.array_equal(expanded[indices], weights) and not np.any(expanded[added]),
            'Entry embedding changed point masses')
    require(np.array_equal(impact_pre[indices], original_pre) and not np.any(impact_pre[added]),
            'Impact embedding changed inherited point masses')
    P.fixed_reference_entry_grid = grid.copy()
    P.fixed_reference_entry_conditional = expanded
    P.Nb = len(grid)
    # make_grid() uses this companion field when the authenticated earnings
    # adapter is installed; keep its numerical metadata internally consistent.
    if hasattr(P, 'earnings_transaction_grid'):
        require(np.array_equal(P.earnings_transaction_grid, old), 'Original transaction grid differs')
        P.earnings_transaction_grid = grid.copy()
    fields = ('fixed_reference_entry_grid', 'fixed_reference_entry_conditional', 'Nb',
              'earnings_transaction_grid')
    overrides = {name: base.serialized(getattr(P, name)) for name in fields if hasattr(P, name)}
    receipt = dict(classification='numerical_only_exact_point_mass_embedding', path=str(path), sha256=sha(path),
        reference_nodes=len(old), experiment_nodes=len(grid), added_nodes=int(added.sum()),
        endpoints_unchanged=True, old_indices=indices.tolist(),
        original_entry_and_inherited_point_masses_exact=True, added_node_entry_mass=0., added_node_impact_mass=0.,
        limitation='Borrowing comparison uses this grid for both regimes. The original-to-refined baseline change is separate; a matched grid is not a grid-convergence certificate.')
    return grid, impact_pre, overrides, receipt


def child(args, plan):
    import numpy as np
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    launch = read(out.parent / 'launch.json')
    require(launch['plan'] == plan and launch['plan_sha256'] == sha(args.plan), 'Controller contract absent')
    index = CASES.index(args.case)
    if index:
        previous = read(out.parent / 'latest_completed.json')['completed']
        require([case['case'] for case in previous] == list(CASES[:index]),
                'Required preceding baseline cases did not finish in order')
        for record in previous:
            control_file = out.parent / record['case'] / 'receipt.json'
            require(sha(control_file) == record['receipt_sha256'], 'Baseline receipt changed')
            prior = read(control_file)
            require(prior['status'] == 'passed' and prior['control']['status'] == 'passed' and
                    prior['plan_sha256'] == sha(args.plan), 'Required preceding baseline failed')
    base = load_module(plan['base_driver_path'], '_fixed_credit_reference_helpers')
    progress(out, 'authenticate')
    manifest, contract, objective, runtime, prepared, reference = base.authenticate(out)
    P = copy.deepcopy(reference['parameters'])
    grid = np.asarray(reference['b_grid']).copy()
    baseline_parameters = base.actual_parameters(prepared, P, grid)
    require(len(baseline_parameters) == 31, 'Full reference parameter table missing')
    for row in manifest['full_parameter_table']:
        require(float(baseline_parameters[row['parameter']]) == float(row['estimate']),
                'Reference parameter differs: ' + row['parameter'])
    P.native_inherited_distribution_evidence_dir = str(out / 'inherited_state_diagnostics')
    original_public = base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    inherited_mass_pin = base.serialized(reference['stationary_g_pre'])
    inherited_parameters_pin = base.serialized(vars(reference['parameters']))
    impact_pre = reference['stationary_g_pre']
    numerical_overrides = {}
    grid_receipt = None
    if args.case != 'control':
        grid, impact_pre, numerical_overrides, grid_receipt = embed_credit_grid(
            base, P, reference, plan['credit_grid_path'])
        base.write(out / 'grid_embedding.json', grid_receipt)
    planned_grid_pin = base.serialized(grid)
    planned_impact_pin = base.serialized(impact_pre)
    adapter = None
    install_receipt = None
    if args.case == 'credit':
        adapter = load_module(plan['adapter_path'], '_fixed_credit_audited_overlay')
        require(adapter.API_VERSION == 'block0506_credit_adapter_v1', 'Unknown credit adapter API')
        for name in ('install', 'audit_purchase_accounting', 'audit_solvency'):
            require(callable(getattr(adapter, name, None)), 'Missing credit adapter hook: ' + name)
        install_receipt = adapter.install(prepared=prepared, reference=reference, P=P, grid=grid,
            output=out, overlay_root=Path(plan['overlay_root']), credit_mode=plan['credit_mode'])
        require(install_receipt['status'] == 'installed' and
                install_receipt['credit_mode'] == plan['credit_mode'], 'Credit installation incomplete')
        require(install_receipt['economic_changes'] == plan['economic_changes'], 'Undisclosed economic change')
    actual_public = base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    changed = {k: dict(reference=original_public.get(k), experiment=actual_public.get(k))
               for k in set(original_public) | set(actual_public)
               if original_public.get(k) != actual_public.get(k)}
    expected_changes = dict(plan['parameter_overrides']) if args.case == 'credit' else {}
    require(set(expected_changes).issubset({'native_due_stayer_credit'}),
            'Only the declared credit-regime flag may change outside pinned grid metadata')
    require(not (set(expected_changes) & set(numerical_overrides)), 'Ambiguous numerical override')
    expected_changes.update(numerical_overrides)
    require(set(changed) == set(expected_changes), 'Unplanned parameter additions/removals/changes')
    for key, value in expected_changes.items():
        require(key in actual_public and actual_public[key] == value, 'Parameter override differs: ' + key)
    require(float(P.psi_child) == float(reference['parameters'].psi_child), 'Child preference changed')
    require(base.serialized(grid) == planned_grid_pin, 'Adapter changed planned asset grid')
    require(base.serialized(impact_pre) == planned_impact_pin, 'Adapter changed embedded inherited mass')
    require(base.serialized(reference['stationary_g_pre']) == inherited_mass_pin,
            'Credit installation mutated inherited household states')
    actual_parameters = base.actual_parameters(prepared, P, grid)
    expected_parameters = dict(baseline_parameters, wealth_grid_nodes=len(grid))
    require(actual_parameters == expected_parameters,
            'A calibrated or externally fixed parameter changed')
    require(base.serialized(vars(reference['parameters'])) == inherited_parameters_pin,
            'Original reference parameter object changed')
    rt = prepared.rt
    model, cal = rt['model'], rt['primitive'].pf.calendar
    price = np.asarray(reference['solution'].p_eq).copy()
    sd = model.precompute_shared(P, grid)
    progress(out, 'one_lifecycle_solve', deadline_epoch=args.deadline)
    require(time.time() < args.deadline, 'Deadline reached before solve')
    started = time.monotonic()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=sd, verbose=False, fast_stats=False)
    elapsed = time.monotonic() - started
    require(time.time() < args.deadline, 'Case deadline exceeded during solve')
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(reconstruction['stationary_post_fertility_nesting_l1'], 5e-9, 'Cohort reconstruction')
    runtime.require_abs_gate(reconstruction['stationary_feasibility_projection_mass'], 0., 'Cohort projection')
    supply = cal.HousingSupplyRule('static-elastic', float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(),
        supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol, evaluation=ev,
        stationary_g_pre=pre, supply_rule=supply, demographic_seed=reference.get('demographic_seed'))
    progress(out, 'cohort_audit', lifecycle_solve_seconds=elapsed)
    cohort_gates = (credit_audits(base, packet, prepared, out, adapter, True) if adapter else
                    base.gates(packet, prepared, out, stationary=True))
    fertility = {kind: rt['observe_initial_fertility'](ev, P, age_projection=kind)
                 for kind in ('uniform_birth_time', 'constant_post_cell')}
    housing = rt['observe_initial_housing_wealth'](ev, P, grid, sd, diagnostic_enabled=True,
        age_projection='uniform_within_age_cell', diagnostic_allow_family_proxies=True,
        include_wealth=True, include_birth_response=True)
    recent = rt['observe_recent_parent_flow'](ev, P, diagnostic_enabled=True,
        snapshot=rt['SNAPSHOT'], age_projection=rt['AGE_PROJECTION'], diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=args.case, reference_checkpoint_sha256=manifest['checkpoint']['sha256']))
    completed_fertility = float(rt['chain'].extract_moments(sol, P)['tfr'])
    fits = runtime.score_targets(objective, fertility, housing, recent['model_value'], completed_fertility)
    require(len(fits) == 14, 'Full fit table missing')
    control = None
    if args.case == 'control':
        control = base.exact_control(reference, packet, fits, manifest, prepared, out)
    elif args.case == 'grid_control':
        control = dict(status='passed', kind='matched_grid_baseline', exact_reference_replay=False,
            interpretation='Unchanged reference credit rules and economics on the refined numerical grid; native gates passed.')
    for row in fits:
        if row['moment'] == 'initial_normalization':
            row['role'] = 'reference replacement benchmark; not imposed in credit experiment'
    parameters = copy.deepcopy(manifest['full_parameter_table'])
    for row in parameters:
        row['status'] = 'Fixed reference value; credit regime separately disclosed; no renormalization'
        if row['parameter'] == 'wealth_grid_nodes' and grid_receipt:
            row['estimate'] = len(grid)
            row['status'] = ('Numerical refinement only: ' + str(len(reference['b_grid'])) + ' -> ' +
                             str(len(grid)) + '; exact inherited and entrant point masses preserved')
        if adapter and row['parameter'] == 'financed_share':
            row['status'] = 'Reference value retained for provenance; artificial LTV constraint inactive in this experiment'
    base.table(out / 'target_fit.csv', fits)
    base.table(out / 'parameters.csv', parameters)
    base.write(out / 'observers.json', cal.jsonable(dict(fertility=fertility, housing_wealth=housing, recent_parent=recent)))
    cohort = aggregates(base, ev, P, grid, model)
    impact_P = copy.deepcopy(P)
    impact_sd = model.precompute_shared(impact_P, grid)
    impact = cal.evaluate_period(price, impact_pre, impact_P, grid, impact_sd,
        cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    impact_out = out / 'baseline_state_impact'
    impact_out.mkdir()
    impact_packet = dict(parameters=impact_P, b_grid=grid, shared=impact_sd, evaluation=impact,
        stationary_g_pre=impact_pre, solution=sol)
    impact_gates = (credit_audits(base, impact_packet, prepared, impact_out, adapter, False) if adapter else
                    base.gates(impact_packet, prepared, impact_out, stationary=False))
    impact_summary = aggregates(base, impact, impact_P, grid, model)
    if args.case == 'control':
        for name in ('g_pre', 'g_post_fertility', 'g_current', 'g_stay_distribution'):
            require(np.array_equal(getattr(impact, name), getattr(reference['evaluation'], name)),
                    'Exact control impact differs: ' + name)
    base.write(impact_out / 'summary.json', impact_summary)
    progress(out, 'standard_17_plot_rendering')
    rt['audit'].standard_diagnostics(packet, out, validate_production_young=False)
    plots = sorted(p.name for p in (out / 'standard_diagnostics').glob('*.png'))
    require(plots == sorted(manifest['standard_diagnostic_names']), 'Standard plot set differs')
    if args.case == 'control':
        for name in plots:
            require(sha(out / 'standard_diagnostics' / name) == manifest['artifact_hashes']['standard_diagnostics/' + name],
                    'Exact control plot differs: ' + name)
    require(base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')}) == actual_public,
            'Parameters changed during solution or reporting')
    require(base.serialized(reference['stationary_g_pre']) == inherited_mass_pin,
            'Inherited household states changed during evaluation')
    require(base.serialized(impact_pre) == planned_impact_pin, 'Embedded inherited states changed during evaluation')
    require(base.serialized(grid) == planned_grid_pin, 'Planned grid changed during evaluation')
    require(base.serialized(vars(reference['parameters'])) == inherited_parameters_pin,
            'Original reference parameters changed during evaluation')
    runtime.verify_sources(dict(contract, objective=manifest['objective']))
    require(verify_plan(args.plan) == plan, 'Executed source or plan changed during case')
    checkpoint = None
    if args.case != 'control':
        path = out / 'conditional_cohort_state.pkl.gz'
        with gzip.open(path, 'wb', compresslevel=1) as stream:
            pickle.dump(packet, stream, protocol=5)
        checkpoint = dict(path=str(path), sha256=sha(path))
    require(time.time() < args.deadline, 'Case deadline exceeded in reporting')
    receipt = dict(status='passed', reference_label=LABEL, case=args.case,
        credit_mode=plan['credit_mode'] if adapter else 'reference', price=price.tolist(),
        baseline_role=('exact_original_grid_reference' if args.case == 'control' else
                       'matched_refined_grid_reference' if args.case == 'grid_control' else 'matched_refined_grid_credit'),
        rent=(P.user_cost_rate * price).tolist(), fixed_psi=float(P.psi_child),
        completed_fertility=completed_fertility, replacement_gap=float(sol.adult_entry_stationary_relative_gap),
        reference_manifest_sha256=MANIFEST_SHA, reference_checkpoint=manifest['checkpoint'],
        source_manifest=manifest['source_manifest'], plan_sha256=sha(args.plan), driver_sha256=sha(__file__),
        files=plan['files'], lifecycle_solves=1, lifecycle_solve_seconds=elapsed,
        normalization_performed=False, control=control, checkpoint=checkpoint,
        parameter_changes=changed,
        numerical_parameter_changes={key: value for key, value in changed.items() if key in numerical_overrides},
        credit_parameter_changes={key: value for key, value in changed.items() if key not in numerical_overrides},
        installation=install_receipt, numerical_grid_change=grid_receipt,
        economic_changes=plan['economic_changes'] if adapter else [],
        cohort_summary=cohort, baseline_state_impact_summary=impact_summary,
        cohort_gates=cohort_gates, impact_gates=impact_gates, standard_plot_count=len(plots),
        artifact_hashes={name: sha(out / name) for name in
            ['target_fit.csv', 'parameters.csv', 'observers.json', 'gates.json'] +
            (['grid_embedding.json'] if grid_receipt else []) +
            ['standard_diagnostics/' + name for name in plots]},
        interpretation='Permanent credit change at fixed reference prices, rents, preferences, and fiscal inputs. Impact uses identical inherited states. Cohort allows normalized-entry composition to change. Neither is a cleared equilibrium, transition, or demographic steady state.')
    base.write(out / 'receipt.json', cal.jsonable(receipt))
    progress(out, 'complete')


def controller(args, plan):
    started = time.time()
    deadline = started + plan['total_seconds']
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / 'launch.json', dict(reference_label=LABEL, plan=plan, plan_sha256=sha(args.plan),
        started_epoch=started, deadline_epoch=deadline, slurm_job=os.environ['SLURM_JOB_ID']))
    completed = []
    write(args.output / 'latest_completed.json', dict(reference_label=LABEL, completed=[], lifecycle_solves=0))
    for name in CASES:
        require(verify_plan(args.plan) == plan, 'Pinned plan changed during loop')
        require(time.time() < deadline, 'Global deadline reached; stop without another case')
        case_deadline = min(deadline, time.time() + plan['case_seconds'])
        command = [sys.executable, str(Path(__file__).resolve()), '--plan', str(args.plan),
                   '--output', str(args.output / name), '--case', name, '--deadline', str(case_deadline)]
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   NUMEXPR_NUM_THREADS='1', NUMBA_NUM_THREADS='1', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
        with (args.output / (name + '.log')).open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
            try:
                while process.poll() is None:
                    progress(args.output, 'case_running', case=name, pid=process.pid,
                             elapsed_seconds=time.time() - started, case_deadline_epoch=case_deadline,
                             completed_cases=len(completed))
                    if time.time() >= case_deadline:
                        raise TimeoutError('Case/global deadline reached: ' + name)
                    time.sleep(2)
                require(process.returncode == 0, 'Case failed; no retry: ' + name)
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
        receipt_path = args.output / name / 'receipt.json'
        receipt = read(receipt_path)
        require(receipt['status'] == 'passed' and receipt['lifecycle_solves'] == 1, 'Case receipt unavailable')
        completed.append(dict(case=name, receipt_sha256=sha(receipt_path), cohort=receipt['cohort_summary'],
            impact=receipt['baseline_state_impact_summary'], completed_fertility=receipt['completed_fertility']))
        write(args.output / 'latest_completed.json', dict(reference_label=LABEL, completed=completed,
            lifecycle_solves=len(completed), remaining_cases=3-len(completed)))
    comparison = {}
    grid_change = {}
    for scope in ('cohort', 'impact'):
        original, baseline, credit = (completed[index][scope] for index in (0, 1, 2))
        comparison[scope] = {key: dict(reference=value, credit=credit[key], difference=credit[key]-value)
            for key, value in baseline.items() if isinstance(value, (int, float))}
        grid_change[scope] = {key: dict(original_grid_reference=value, refined_grid_reference=baseline[key],
            difference=baseline[key]-value) for key, value in original.items() if isinstance(value, (int, float))}
    write(args.output / 'comparison.json', dict(reference_label=LABEL, credit_mode=plan['credit_mode'],
        baseline_reference='grid_control: unchanged reference economics on the identical refined grid used for credit',
        original_reference='control: exact original-grid replay of ' + LABEL,
        economic_changes=plan['economic_changes'], **comparison,
        numerical_grid_change=read(args.output / 'credit/receipt.json')['numerical_grid_change'],
        grid_change=dict(**grid_change, completed_fertility=dict(original_grid_reference=completed[0]['completed_fertility'],
            refined_grid_reference=completed[1]['completed_fertility'],
            difference=completed[1]['completed_fertility']-completed[0]['completed_fertility'])),
        completed_fertility=dict(reference=completed[1]['completed_fertility'], credit=completed[2]['completed_fertility'],
            difference=completed[2]['completed_fertility']-completed[1]['completed_fertility']),
        interpretation='Fixed-price borrowing mechanism compares identical refined grids; numerical refinement effect is shown separately. No market-clearing or steady-state claim.'))
    require(time.time() < deadline, 'Global deadline exceeded during finalization')
    write(args.output / 'completed.json', dict(status='passed', reference_label=LABEL, lifecycle_solves=3,
        elapsed_seconds=time.time()-started, completed=completed, plan_sha256=sha(args.plan), normalization_performed=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--deadline', type=float)
    args = parser.parse_args()
    args.plan, args.output = args.plan.resolve(), args.output.resolve()
    plan = verify_plan(args.plan)
    require(not args.output.exists(), 'Versioned output must not already exist')
    try:
        if args.case:
            require(args.deadline is not None and time.time() < args.deadline <= time.time()+plan['case_seconds'],
                    'Child requires a bounded absolute deadline')
            child(args, plan)
        else:
            require(args.deadline is None, 'Controller sets immutable deadline')
            controller(args, plan)
    except BaseException as error:
        if args.output.exists():
            helper = sys.modules.get('_fixed_credit_reference_helpers')
            ledger = getattr(error, 'ledger', getattr(error, 'audit', None))
            ledger = helper.serialized(ledger) if helper is not None else repr(ledger)
            write(args.output / 'failure.json', dict(reference_label=LABEL, status='failed', retries=0,
                error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc(),
                failure_ledger=ledger, time_epoch=time.time()))
        raise


if __name__ == '__main__':
    main()
