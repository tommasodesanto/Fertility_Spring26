#!/usr/bin/env python3
"""Bounded closed stationary GE controller for the block0506 natural-credit endpoint.

Torch Slurm only. Every distinct price is solved in a fresh interpreter. A
candidate need not renew; only a selected root and its exact repeat can certify
GE. No calibration, preference normalization, or inherited-state transition.
"""
from __future__ import annotations
import argparse
import copy
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

LABEL = '2007 stationary reference — block0506, September 28 verified export'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
MANIFEST = ROOT / 'output/model/fertility_identification_20260928/fixed_reference_manifest.json'
MODE = 'natural_solvency_boolean_grid_stationary_ge_v1'
RENEWAL_TOL = 1e-6
PAYGO_TOL = 1e-6


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def progress(out, phase, **details):
    write(Path(out) / 'progress.json', dict(reference_label=LABEL, phase=phase,
                                           time_epoch=time.time(), **details))


def verify_plan(path):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(),
            'Torch Slurm required')
    p = read(path)
    require(p['schema'] == 'block0506_stationary_credit_ge_v1' and
            p['reference_label'] == LABEL, 'Wrong GE contract/reference')
    require(p['driver_sha256'] == sha(__file__) and
            p['reference_manifest_sha256'] == MANIFEST_SHA == sha(MANIFEST),
            'Driver or manifest identity differs')
    require(p['reference_checkpoint_sha256'] == 'b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d',
            'Reference checkpoint pin differs')
    require(p['credit_seed_checkpoint_sha256'] == '7b0b92f0e2ed33cfc80a3c63c2c15d6417b14ab88db3930284e440e94e2e515c' and
            p['credit_seed_receipt_sha256'] == '56fd38695487ae1262b5d772254d77943c046f106f91fcd42a92dc8d452ed18c',
            'Credit seed pin differs')
    require(p['maximum_lifecycle_solves'] == 12 and p['reserve_exact_repeat_solves'] == 1 and
            p['case_seconds'] == 300 and p['total_seconds'] == 2400 and
            p['threads'] == 1 and p['memory_gib'] == 16, 'Resource contract differs')
    require(p['credit_mode'] == MODE and p['renormalize_fertility'] is False and
            p['psi_child_fixed'] == 0.1355551166583114, 'Economic contract differs')
    require(p['upper_price_ratios'] == [1.05, 1.1, 1.2, 1.35] and
            p['lower_price_search_authorized'] is False, 'Price domain differs')
    pinned = {}
    for item in p['files']:
        path = str(Path(item['path']).resolve())
        require(path not in pinned and sha(path) == item['sha256'], 'Missing/changed pin: ' + path)
        pinned[path] = item['sha256']
    for key in ('adapter_path', 'base_driver_path', 'credit_driver_path',
                'preflight_receipt_path', 'credit_seed_receipt_path', 'credit_seed_checkpoint_path',
                'credit_seed_fit_path', 'credit_seed_parameters_path'):
        require(str(Path(p[key]).resolve()) in pinned, 'Unpinned required input: ' + key)
    require(Path(p['adapter_path']).resolve().is_relative_to(Path(p['overlay_root']).resolve()),
            'GE adapter outside isolated overlay')
    for source in Path(p['overlay_root']).glob('*.py'):
        require(str(source.resolve()) in pinned, 'Unpinned source: ' + str(source))
    seed = read(p['credit_seed_receipt_path'])
    require(seed['status'] == 'passed' and seed['case'] == 'credit' and
            seed['checkpoint']['sha256'] == p['credit_seed_checkpoint_sha256'],
            'Prior credit seed receipt differs')
    preflight = read(p['preflight_receipt_path'])
    require(preflight['status'] == 'ready' and preflight['model_solves'] == 0 and
            preflight['q0_grid_reproduces_prior_262_nodes_exactly'] and
            preflight['source_sha256']['natural_credit.py'] == sha(p['adapter_path']),
            'GE preflight or adapter identity differs')
    return p


def closed_accounting(entry, adjusted_births, demand, supply, price):
    values = (entry, adjusted_births, demand, supply, price)
    require(all(math.isfinite(float(x)) for x in values) and
            min(entry, demand, supply, price) > 0 and adjusted_births >= 0,
            'Invalid stationary closure input')
    N = supply / demand
    return dict(price=float(price), entry_per_normalized_household=float(entry),
        adjusted_births_per_normalized_household=float(adjusted_births),
        renewal_residual=float(adjusted_births / (2.1 * entry) - 1),
        population_scale=float(N), normalized_housing_demand=float(demand),
        absolute_housing_supply=float(supply), absolute_housing_demand=float(N * demand),
        absolute_housing_residual=float(N * demand - supply), absolute_entry=float(N * entry),
        absolute_adjusted_births=float(N * adjusted_births), outside_entry=0.,
        retention=1., replacement_conversion=1 / 2.1)


def scaled_native_step(packet, prepared, closure):
    """Replay the native 16/20 entry queue at actual population without resetting it."""
    import numpy as np
    native = prepared.rt['primitive'].pf.transition
    cal = prepared.rt['primitive'].pf.calendar
    P, grid, shared = (packet[key] for key in ('parameters', 'b_grid', 'shared'))
    require(abs(closure['renewal_residual']) <= RENEWAL_TOL, 'Root renewal not certified')
    N = closure['population_scale']
    ev = copy.copy(packet['evaluation'])
    for name in ('g_pre', 'g_post_fertility', 'g_current', 'g_stay_distribution'):
        value = getattr(ev, name, None)
        if value is not None:
            setattr(ev, name, N * value)
    ev.births *= N
    ev.births_by_loc = N * ev.births_by_loc
    ev.demand_by_loc = N * ev.demand_by_loc
    actual = float(native.calendar_topcode_birth_accounting(
        ev.g_pre, ev.g_post_fertility, float(ev.births), P)['topcode_adjusted_birth_children'])
    require(abs(actual - closure['absolute_adjusted_births']) <= max(1., N) * 5e-9,
            'Scaled adjusted births disagree')
    queue = native.SplitBirthEntryQueue.constant_prehistory(actual)
    due, next_queue = queue.step(actual)
    require(len(queue.due_in_16) == 3 and len(queue.due_in_20) == 4 and
            abs(due - actual / 2.1) <= max(1., N) * 2e-10,
            'Native split queue timing/conversion differs')
    half = actual / (2 * 2.1)
    require(all(abs(x-half) <= max(1., N)*2e-10 for x in
                queue.due_in_16 + queue.due_in_20 +
                next_queue.due_in_16 + next_queue.due_in_20) and
            next_queue.due_in_16 == queue.due_in_16 and
            next_queue.due_in_20 == queue.due_in_20,
            'Native 16/20 queue renewal tuple did not replay')
    following, _, deaths, mass_residual = native.advance_sequential_calendar_distribution(
        ev, np.asarray([due]), P, grid, shared)
    target = N * packet['stationary_g_pre']
    raw_l1 = float(np.abs(following - target).sum())
    entry_gap = float(due - closure['absolute_entry'])
    accounted = following - target
    accounted[:, :, :, 0] -= cal.entrant_cohort(np.asarray([entry_gap]), P, grid)
    adjusted_l1 = float(np.abs(accounted).sum())
    require(adjusted_l1 <= max(1., N) * 5e-9 and
            abs(float(mass_residual)) <= max(1., N) * 2e-8,
            'Native scaled stationary operator failed')
    return dict(status='passed', population_scale=N, actual_birth_derived_entry=float(due),
        stationary_required_entry=closure['absolute_entry'], entry_gap=entry_gap,
        raw_distribution_l1=raw_l1, distribution_l1_after_accounted_renewal_error=adjusted_l1,
        mass_residual=float(mass_residual), deaths=float(deaths),
        split_16_20_queue_reproduced=True, next_queue_type=type(next_queue).__name__)


def compare_seed(packet, fits, parameters, p, out, base):
    import e5f_current_transition_runtime as native
    with gzip.open(p['credit_seed_checkpoint_path'], 'rb') as stream:
        seed = pickle.load(stream)
    census = native.compare_arrays(seed, packet)
    base.write(out / 'seed_array_replay.json', census)
    bad = [name for name, row in census['arrays'].items()
           if row.get('status') != 'compared' or not row.get('exact') or not row.get('finite')]
    require(not bad, 'q0 credit seed numeric array mismatch: ' + ', '.join(bad[:12]))
    prior_fit = base.csv_rows(p['credit_seed_fit_path'])
    require([r['moment'] for r in fits] == [r['moment'] for r in prior_fit], 'q0 fit row names differ')
    for actual, old in zip(fits, prior_fit):
        for key in ('target', 'model', 'gap', 'weight', 'loss_contribution'):
            require((actual[key] == '' and old[key] == '') or
                    (actual[key] != '' and old[key] != '' and float(actual[key]) == float(old[key])),
                    'q0 fit differs: ' + actual['moment'] + '/' + key)
    old_parameters = base.csv_rows(p['credit_seed_parameters_path'])
    require(len(parameters) == len(old_parameters) == 31, 'q0 parameter row count differs')
    for row, old in zip(parameters, old_parameters):
        require(row['parameter'] == old['parameter'] and
                float(row['estimate']) == float(old['estimate']),
                'q0 parameter differs: ' + row['parameter'])
    return dict(status='passed', numeric_arrays_exact=census['array_count'],
                fit_rows_exact=14, parameter_estimates_exact=31,
                checkpoint_sha256=p['credit_seed_checkpoint_sha256'])


def compare_repeat(packet, fits, parameters, out, base):
    import e5f_current_transition_runtime as native
    selected = read(out.parent / 'selected.json')
    first = out.parent / selected['case']
    require(sha(first / 'receipt.json') == selected['receipt_sha256'] and
            sha(selected['checkpoint']['path']) == selected['checkpoint']['sha256'],
            'Selected receipt/checkpoint changed before repeat')
    with gzip.open(selected['checkpoint']['path'], 'rb') as stream:
        original = pickle.load(stream)
    census = native.compare_arrays(original, packet)
    base.write(out / 'selected_repeat_arrays.json', census)
    bad = [name for name, row in census['arrays'].items()
           if row.get('status') != 'compared' or not row.get('exact') or not row.get('finite')]
    require(not bad, 'Selected repeat numeric arrays differ: ' + ', '.join(bad[:12]))
    old_fit = base.csv_rows(first / 'target_fit.csv')
    require([r['moment'] for r in fits] == [r['moment'] for r in old_fit],
            'Selected repeat fit rows differ')
    for actual, old in zip(fits, old_fit):
        for key in ('target', 'model', 'gap', 'weight', 'loss_contribution'):
            require((actual[key] == '' and old[key] == '') or
                    (actual[key] != '' and old[key] != '' and float(actual[key]) == float(old[key])),
                    'Selected repeat fit differs: ' + actual['moment'] + '/' + key)
    old_parameters = base.csv_rows(first / 'parameters.csv')
    require(len(parameters) == len(old_parameters) == 31 and
            all(a['parameter'] == b['parameter'] and
                float(a['estimate']) == float(b['estimate'])
                for a, b in zip(parameters, old_parameters)),
            'Selected repeat parameter estimates differ')
    return dict(status='passed', numeric_arrays_exact=census['array_count'],
                fit_rows_exact=14, parameter_estimates_exact=31)


def solve_case(args, p):
    import numpy as np
    import run_credit as credit
    import run_fixed_price as base
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    require(read(out.parent / 'launch.json')['plan_sha256'] == sha(args.plan), 'Controller launch absent')
    progress(out, 'authenticate', price_factor=args.factor)
    manifest, contract, objective, runtime, prepared, reference = base.authenticate(out)
    require(manifest['checkpoint']['sha256'] == p['reference_checkpoint_sha256'], 'Wrong checkpoint')
    P = copy.deepcopy(reference['parameters'])
    P.native_inherited_distribution_evidence_dir = str(out / 'inherited_state_diagnostics')
    model, cal = prepared.rt['model'], prepared.rt['primitive'].pf.calendar
    original_grid = np.asarray(reference['b_grid'])
    q0 = float(np.asarray(reference['solution'].p_eq)[0])
    price = np.asarray([q0 * args.factor])
    baseline_parameters = base.actual_parameters(prepared, P, original_grid)
    require(len(baseline_parameters) == 31, 'Reference parameter table missing')
    for row in manifest['full_parameter_table']:
        require(float(baseline_parameters[row['parameter']]) == float(row['estimate']),
                'Reference parameter differs: ' + row['parameter'])
    original_public = base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    reference_pin = base.serialized(vars(reference['parameters']))
    adapter = credit.load_module(p['adapter_path'], '_block0506_credit_ge_adapter')
    require(adapter.API_VERSION == 'block0506_credit_adapter_v1' and adapter.MODE == MODE,
            'GE adapter API differs')
    grid = adapter.build_grid(model, P, original_grid, float(price[0]))
    indices = np.searchsorted(grid, original_grid)
    require(np.array_equal(grid[indices], original_grid), 'Grid lost original atoms')
    if args.factor == 1.:
        require(len(grid) == 262, 'q0 grid does not reproduce 262 nodes')
    grid_json = out / 'credit_grid.json'
    base.write(grid_json, dict(grid=grid.tolist(), old_indices=indices.tolist(), price=float(price[0])))
    grid, embedded_pre, numerical, embedding = credit.embed_credit_grid(base, P, reference, grid_json)
    installation = adapter.install(prepared=prepared, reference=reference, P=P, grid=grid,
        price=float(price[0]), output=out, overlay_root=Path(p['overlay_root']), credit_mode=MODE)
    require(installation['status'] == 'installed' and installation['credit_mode'] == MODE and
            installation['economic_changes'] == p['economic_changes'], 'Credit installation differs')
    actual_public = base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    changed = {k: dict(reference=original_public.get(k), experiment=actual_public.get(k))
               for k in set(original_public) | set(actual_public)
               if original_public.get(k) != actual_public.get(k)}
    expected = dict(numerical, native_due_stayer_credit=False)
    require(set(changed) == set(expected) and
            all(actual_public[k] == v for k, v in expected.items()),
            'Undeclared parameter mutation')
    require(float(P.psi_child) == p['psi_child_fixed'] and
            base.actual_parameters(prepared, P, grid) == dict(baseline_parameters, wealth_grid_nodes=len(grid)) and
            base.serialized(vars(reference['parameters'])) == reference_pin,
            'Fixed economic parameters changed')
    sd = model.precompute_shared(P, grid)
    progress(out, 'one_lifecycle_solve', price=float(price[0]), grid_nodes=len(grid),
             deadline_epoch=args.deadline)
    require(time.time() < args.deadline, 'Deadline reached before solve')
    started = time.monotonic()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=sd,
                                              verbose=False, fast_stats=False)
    elapsed = time.monotonic() - started
    require(time.time() < args.deadline, 'Case deadline during solve')
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, sd)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, sd)
    runtime.require_abs_gate(reconstruction['stationary_post_fertility_nesting_l1'], 5e-9,
                             'Cohort reconstruction')
    runtime.require_abs_gate(reconstruction['stationary_feasibility_projection_mass'], 0.,
                             'Cohort projection')
    supply = cal.HousingSupplyRule('static-elastic', float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]),
        float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, sd, cal.SolveCounter(),
                             supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=sd, solution=sol,
                  evaluation=ev, stationary_g_pre=pre, supply_rule=supply,
                  demographic_seed=reference.get('demographic_seed'))
    gates = credit.credit_audits(base, packet, prepared, out, adapter, True)
    fiscal = gates['fiscal']['scaled_pension_budget_residual']
    require(abs(float(fiscal)) <= PAYGO_TOL and
            gates['fiscal_certificate']['fiscal_gate'] and
            gates['fiscal_certificate']['marginal_gate'], 'Actual PAYGO or stationary margins fail')
    adjusted = float(sol.adult_entry_adjusted_birth_children)
    closure = closed_accounting(float(sol.entry_rate), adjusted,
        float(np.asarray(ev.demand_by_loc).sum()), float(np.asarray(ev.supply_by_loc).sum()),
        float(price[0]))
    closure['supply_exponent'] = float(P.xi_supply[0])
    require(abs(closure['supply_exponent'] - 0.63) <= 1e-12 and
            abs(closure['absolute_housing_residual']) <= 1e-12,
            'Absolute housing supply contract failed')
    housing_certificate = dict(housing_market_clearing_required=True,
        absolute_demand=closure['absolute_housing_demand'],
        absolute_supply=closure['absolute_housing_supply'],
        absolute_residual=closure['absolute_housing_residual'],
        relative_residual=abs(closure['absolute_housing_residual']) / closure['absolute_housing_supply'],
        gate=abs(closure['absolute_housing_residual']) <= 1e-12,
        interpretation='N times normalized demand equals the unchanged absolute housing supply at q')
    rt = prepared.rt
    fertility = {kind: rt['observe_initial_fertility'](ev, P, age_projection=kind)
                 for kind in ('uniform_birth_time', 'constant_post_cell')}
    housing = rt['observe_initial_housing_wealth'](ev, P, grid, sd,
        diagnostic_enabled=True, age_projection='uniform_within_age_cell',
        diagnostic_allow_family_proxies=True, include_wealth=True, include_birth_response=True)
    recent = rt['observe_recent_parent_flow'](ev, P, diagnostic_enabled=True,
        snapshot=rt['SNAPSHOT'], age_projection=rt['AGE_PROJECTION'],
        diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=args.case, reference_checkpoint_sha256=manifest['checkpoint']['sha256']))
    completed_fertility = float(rt['chain'].extract_moments(sol, P)['tfr'])
    fits = runtime.score_targets(objective, fertility, housing, recent['model_value'],
                                 completed_fertility)
    require(len(fits) == 14, '14 fit rows required')
    for row in fits:
        if row['moment'] == 'initial_normalization':
            row['role'] = 'reference replacement target; no normalization in GE experiment'
    parameters = copy.deepcopy(manifest['full_parameter_table'])
    require(len(parameters) == 31, '31 parameter rows required')
    for row in parameters:
        row['status'] = 'Fixed reference value; no normalization'
        if row['parameter'] == 'wealth_grid_nodes':
            row['estimate'] = len(grid)
            row['status'] = 'Numerical embedding only; original entrant and household atoms preserved'
        if row['parameter'] == 'financed_share':
            row['status'] = 'Reference value retained; artificial LTV limit inactive'
    base.table(out / 'target_fit.csv', fits)
    base.table(out / 'parameters.csv', parameters)
    base.write(out / 'observers.json', cal.jsonable(dict(fertility=fertility,
        housing_wealth=housing, recent_parent=recent)))
    base.write(out / 'grid_embedding.json', embedding)
    cohort = credit.aggregates(base, ev, P, grid, model)
    seed_replay = None
    if args.case == 'q0_smoke':
        seed_replay = compare_seed(packet, fits, parameters, p, out, base)
    repeat_check = compare_repeat(packet, fits, parameters, out, base) if args.role == 'repeat' else None
    selected = abs(closure['renewal_residual']) <= RENEWAL_TOL
    native = scaled_native_step(packet, prepared, closure) if selected else None
    if selected or args.case == 'q0_smoke':
        # The q0 smoke shows the raw normalized-household supply imbalance.
        # Only a renewing endpoint uses Hs/N in the reporting copy.
        report_ev = copy.copy(ev)
        if selected:
            report_ev.supply_by_loc = np.asarray(ev.supply_by_loc) / closure['population_scale']
            report_ev.relative_market_residual = float(np.max(np.abs(
                (ev.demand_by_loc - report_ev.supply_by_loc) / report_ev.supply_by_loc)))
        report_packet = dict(packet, evaluation=report_ev)
        rt['audit'].standard_diagnostics(report_packet, out, validate_production_young=False)
        plots = sorted(path.name for path in (out / 'standard_diagnostics').glob('*.png'))
        require(plots == sorted(manifest['standard_diagnostic_names']), 'Standard 17 plots differ')
        base.write(out / 'reporting_units.json', dict(
            economic_absolute_supply=closure['absolute_housing_supply'],
            reporting_supply_per_normalized_household=float(report_ev.supply_by_loc.sum()),
            normalized_demand=float(ev.demand_by_loc.sum()),
            population_scale=closure['population_scale'], economic_H0_unchanged=True,
            supply_divided_by_population_for_plot=bool(selected),
            q0_smoke_not_a_cleared_economy=bool(args.case == 'q0_smoke' and not selected),
            plotted_quantity_unit='rooms per normalized adult household'))
    else:
        plots = []
    with gzip.open(out / 'conditional_cohort_state.pkl.gz', 'wb', compresslevel=1) as stream:
        pickle.dump(packet, stream, protocol=5)
    checkpoint = dict(path=str(out / 'conditional_cohort_state.pkl.gz'),
                      sha256=sha(out / 'conditional_cohort_state.pkl.gz'))
    require(base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')}) == actual_public and
            base.serialized(vars(reference['parameters'])) == reference_pin and
            base.serialized(grid) == base.serialized(np.asarray(read(grid_json)['grid'])),
            'Model/entry grid mutated during solution or reporting')
    runtime.verify_sources(dict(contract, objective=manifest['objective']))
    require(verify_plan(args.plan) == p and time.time() < args.deadline,
            'Source identity or deadline changed during case')
    receipt = dict(status='passed', reference_label=LABEL, case=args.case, role=args.role,
        lifecycle_solves=1, lifecycle_solve_seconds=elapsed, price=float(price[0]),
        price_factor=args.factor, rent=float(P.user_cost_rate * price[0]),
        completed_fertility=completed_fertility, fixed_psi=float(P.psi_child),
        renewal_residual=closure['renewal_residual'], paygo_residual=float(fiscal),
        closure=closure, ge_housing_certificate=housing_certificate,
        native_scaled_step=native, cohort_summary=cohort,
        seed_replay=seed_replay, repeat_check=repeat_check, gates=gates, checkpoint=checkpoint,
        grid_embedding=embedding, installation=installation,
        parameter_changes=changed, normalization_performed=False,
        reference_checkpoint=manifest['checkpoint'], reference_manifest_sha256=MANIFEST_SHA,
        plan_sha256=sha(args.plan), standard_plot_count=len(plots),
        estate_closure='Provisional reference net-estate funding/residual sink retained; counterparties unresolved')
    base.write(out / 'receipt.json', cal.jsonable(receipt))
    progress(out, 'complete', receipt_sha256=sha(out / 'receipt.json'))


def choose_bracket(records):
    ordered = sorted(records, key=lambda x: x['factor'])
    for left, right in zip(ordered, ordered[1:]):
        if left['residual'] * right['residual'] <= 0:
            return left, right
    return None


def next_log_factor(left, right):
    lo, hi = math.log(left['factor']), math.log(right['factor'])
    fl, fh = left['residual'], right['residual']
    require(lo < hi and fl * fh <= 0 and fl != fh, 'Invalid renewal bracket')
    raw = (lo * fh - hi * fl) / (fh - fl)
    safe = min(max(raw, lo + .1 * (hi - lo)), hi - .1 * (hi - lo))
    return math.exp(safe)


def controller(args, p):
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    deadline = started + p['total_seconds']
    write(out / 'launch.json', dict(reference_label=LABEL, plan=p,
        plan_sha256=sha(args.plan), started_epoch=started,
        deadline_epoch=deadline, slurm_job=os.environ['SLURM_JOB_ID']))
    records = []
    write(out / 'latest_completed.json', dict(completed=[], lifecycle_solves=0))
    write(out / 'best_so_far.json', dict(status='no_case_completed'))

    def run_case(name, factor, role):
        require(verify_plan(args.plan) == p, 'Plan/source changed during loop')
        require(len(records) < p['maximum_lifecycle_solves'] - (1 if role != 'repeat' else 0),
                'Repeat reserve or total solve budget exhausted')
        require(time.time() < deadline, 'Total deadline reached before case')
        case_deadline = min(deadline, time.time() + p['case_seconds'])
        command = [sys.executable, str(Path(__file__).resolve()), '--plan', str(args.plan),
                   '--output', str(out / name), '--case', name, '--factor', repr(factor),
                   '--role', role, '--deadline', repr(case_deadline)]
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', NUMBA_NUM_THREADS='1',
                   MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
        with (out / (name + '.log')).open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                       env=env, start_new_session=True)
            try:
                while process.poll() is None:
                    progress(out, 'case_running', case=name, price_factor=factor,
                             elapsed_seconds=time.time() - started,
                             completed_cases=len(records), deadline_epoch=case_deadline,
                             best_so_far=read(out / 'best_so_far.json'))
                    if time.time() >= case_deadline:
                        raise TimeoutError('Case/global deadline: ' + name)
                    time.sleep(2)
                require(process.returncode == 0, 'Case failed without retry: ' + name)
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
        receipt_path = out / name / 'receipt.json'
        receipt = read(receipt_path)
        require(receipt['status'] == 'passed' and receipt['lifecycle_solves'] == 1 and
                receipt['plan_sha256'] == sha(args.plan), 'Case receipt incomplete')
        row = dict(case=name, role=role, factor=factor, price=receipt['price'],
                   residual=receipt['renewal_residual'], population=receipt['closure']['population_scale'],
                   receipt_sha256=sha(receipt_path), checkpoint=receipt['checkpoint'])
        records.append(row)
        best = min(records, key=lambda x: abs(x['residual']))
        write(out / 'latest_completed.json', dict(completed=records,
            lifecycle_solves=len(records), remaining_total=12-len(records),
            repeat_reserved=(role != 'repeat')))
        write(out / 'best_so_far.json', dict(status='case_completed', **best))
        return row

    q0 = run_case('q0_smoke', 1., 'seed')
    q0_receipt = read(out / 'q0_smoke' / 'receipt.json')
    require(q0_receipt['seed_replay']['status'] == 'passed' and
            q0_receipt['standard_plot_count'] == 17 and
            q0_receipt['checkpoint']['sha256'] == q0['checkpoint']['sha256'],
            'q0 exact-loop smoke incomplete; no price search')
    selected = q0 if abs(q0['residual']) <= RENEWAL_TOL else None
    bracket = None
    if selected is None:
        for index, factor in enumerate(p['upper_price_ratios']):
            row = run_case('upper_' + str(int(round(100 * factor))), factor, 'search')
            bracket = choose_bracket(records)
            if abs(row['residual']) <= RENEWAL_TOL:
                selected = row
                break
            if bracket is not None:
                break
    require(selected is not None or bracket is not None,
            'Upper prices through 1.35q0 did not bracket renewal; no lower search authorized')
    if selected is None:
        left, right = bracket
        for index in range(p['maximum_lifecycle_solves'] - 1 - len(records)):
            factor = next_log_factor(left, right)
            row = run_case('root_' + str(index + 1).zfill(2), factor, 'search')
            if abs(row['residual']) <= RENEWAL_TOL:
                selected = row
                break
            if row['residual'] * left['residual'] > 0:
                left = row
            else:
                right = row
            write(out / 'bracket.json', dict(left=left, right=right,
                width_log_price=math.log(right['factor'] / left['factor'])))
    require(selected is not None, 'Solve cap reached without 1e-6 renewal root')
    selected_receipt = read(out / selected['case'] / 'receipt.json')
    require(selected_receipt['native_scaled_step']['status'] == 'passed' and
            selected_receipt['standard_plot_count'] == 17,
            'Selected root has incomplete GE diagnostics')
    write(out / 'selected.json', selected)
    # The fresh repeat verifies every saved numeric array, fit and parameter
    # inside its authenticated runtime; the controller does no model imports.
    repeat = run_case('selected_repeat', selected['factor'], 'repeat')
    repeat_receipt = read(out / repeat['case'] / 'receipt.json')
    require(repeat_receipt['repeat_check']['status'] == 'passed' and
            selected['residual'] == repeat['residual'] and
            abs(repeat['residual']) <= RENEWAL_TOL and
            repeat_receipt['native_scaled_step']['status'] == 'passed' and
            repeat_receipt['standard_plot_count'] == 17,
            'Repeat has incomplete GE certification')
    # Search selection is reported from the verified repeat, with exact
    # source/fit/parameter provenance and the complete standard graph set.
    write(out / 'completed.json', dict(status='certified_stationary_ge',
        reference_label=LABEL, selected=selected, repeat=repeat,
        selected_receipt_sha256=selected['receipt_sha256'],
        repeat_receipt_sha256=repeat['receipt_sha256'],
        exact_repeat_numeric_arrays=repeat_receipt['repeat_check']['numeric_arrays_exact'],
        total_lifecycle_solves=len(records), elapsed_seconds=time.time()-started,
        plan_sha256=sha(args.plan), normalization_performed=False,
        interpretation='Closed stationary credit GE endpoint; no dated transition or temporary fertility impact.'))
    progress(out, 'complete', selected_case=selected['case'], repeat_case=repeat['case'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case')
    parser.add_argument('--factor', type=float)
    parser.add_argument('--role', choices=('seed', 'search', 'selected', 'repeat'))
    parser.add_argument('--deadline', type=float)
    args = parser.parse_args()
    args.plan, args.output = args.plan.resolve(), args.output.resolve()
    plan = verify_plan(args.plan)
    require(not args.output.exists(), 'Versioned output must not already exist')
    try:
        if args.case:
            require(args.factor is not None and args.factor >= 1 and args.role and
                    args.deadline is not None and
                    time.time() < args.deadline <= time.time() + plan['case_seconds'],
                    'Invalid bounded child case')
            solve_case(args, plan)
        else:
            require(args.deadline is None and args.factor is None and args.role is None,
                    'Controller owns prices and deadline')
            controller(args, plan)
    except BaseException as error:
        if args.output.exists():
            write(args.output / 'failure.json', dict(status='failed', reference_label=LABEL,
                error_type=type(error).__name__, error=str(error),
                traceback=traceback.format_exc(), retries=0, time_epoch=time.time()))
        raise


if __name__ == '__main__':
    main()
