#!/usr/bin/env python3
"""Zero-solve, saved-policy housing-supply closure checks on Torch only."""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import inspect
import json
import os
import pickle
import sys
import time
import traceback
from pathlib import Path

LABEL = '2007 stationary reference — block0506, September 28 verified export'
REFERENCE_SHA = 'b15ba92dc60e3d5590d2beb6e05d36f71d17b20b1a432edc2c2db926a217309d'
GE_SHA = '2a3e85bca910f2057727637075da174f551de8a88dfca4eb6c4dbe6527db7ecd'
MANIFEST_SHA = '147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4'
BASELINE_STOCK = 5.847942339355678
RENEWAL_TOL = 1e-6
PAYGO_TOL = 1e-6


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def table(path, records):
    require(bool(records), 'Cannot write an empty table')
    with Path(path).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def forbid_solves(model, cal):
    def blocked(*_args, **_kwargs):
        raise RuntimeError('ZERO-SOLVE CONTRACT: household solve attempted')
    for obj, name in ((model, 'solve_markov_income_at_prices'),
                      (model, 'solve_bellman_full_markov_income'),
                      (cal, 'solve_policy')):
        require(hasattr(obj, name), 'Missing guarded solve API: ' + name)
        setattr(obj, name, blocked)


def verify_model_independence(model, cal):
    root = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
    model_path = root/'code/model/intergen_eqscale_seq_optimized/solver.py'
    calendar_path = root/'code/model/tools/run_dynamic_population_transition.py'
    require(Path(model.__file__).resolve() == model_path and
            sha(model_path) == 'b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1',
            'Authenticated active optimized solver differs')
    require(Path(inspect.getsourcefile(cal.evaluate_period)).resolve() == calendar_path and
            sha(calendar_path) == '849709afb281193a48ccc44b352ed548c615c28da00b87301f9b7ad212bb1d8c',
            'Authenticated period evaluator differs')
    functions = ('solve_markov_income_at_prices', 'precompute_shared',
                 'solve_bellman_full_markov_income', 'forward_distribution_markov_income')
    for name in functions:
        fn = getattr(model, name)
        require(Path(inspect.getsourcefile(fn)).resolve() == model_path and
                'P.H0' not in inspect.getsource(fn) and
                'P.xi_supply' not in inspect.getsource(fn),
                'Household optimization/distribution calls supply: ' + name)
    require('supply_rule.quantity(policy.price)' in inspect.getsource(cal.evaluate_period),
            'Supply does not enter evaluator as expected')
    return dict(model_path=str(model_path), model_sha256=sha(model_path),
                calendar_path=str(calendar_path), calendar_sha256=sha(calendar_path),
                supply_independent_functions=list(functions),
                reporting_market_stats_rebased_separately=True)


def verify_plan(p, driver):
    require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(),
            'Torch Slurm execution required')
    require(p['reference_label'] == LABEL and p['schema'] == 'block0506_supply_replay_v1',
            'Wrong supply plan')
    require(p['driver_sha256'] == sha(driver) and
            p['reference_checkpoint_sha256'] == REFERENCE_SHA and
            p['credit_ge_checkpoint_sha256'] == GE_SHA and
            p['manifest_sha256'] == MANIFEST_SHA,
            'Source/reference pin differs')
    require(p['maximum_lifecycle_solves'] == 0 and p['threads'] == 1 and
            p['memory_gib'] == 16 and p['wall_seconds'] == 300,
            'Resource or zero-solve contract differs')
    for item in p['files']:
        require(sha(item['path']) == item['sha256'], 'Missing/changed pin: ' + item['path'])
    require(sha(p['reference_checkpoint_path']) == REFERENCE_SHA and
            sha(p['credit_ge_checkpoint_path']) == GE_SHA,
            'Checkpoint bytes differ')


def replay_case(name, source, P, supply_rule, supply_description, prepared, out, fit_rows,
                parameter_rows, checkpoint_sha, original_h0, original_eta, credit_adapter=None):
    import numpy as np
    import run_credit as credit
    import run_fixed_price as base

    cal = prepared.rt['primitive'].pf.calendar
    native = prepared.rt['primitive'].pf.transition
    model = prepared.rt['model']
    grid = np.asarray(source['b_grid'])
    sol = source['solution']
    saved_ev = source['evaluation']
    price = np.asarray(sol.p_eq)
    pre = np.asarray(source['stationary_g_pre'])
    case_dir = out/name
    case_dir.mkdir()
    P.native_inherited_distribution_evidence_dir = str(case_dir/'inherited_state_diagnostics')
    shared = source.get('shared')
    if shared is None:
        shared = model.precompute_shared(P, grid)
    ev = cal.evaluate_period(price, pre, P, grid, shared, cal.SolveCounter(),
                             supply_rule=supply_rule, supplied_policy=saved_ev.policy)
    unchanged = {}
    for key in ('g_pre', 'g_post_fertility', 'g_current', 'births_by_loc', 'demand_by_loc'):
        lhs, rhs = np.asarray(getattr(ev, key)), np.asarray(getattr(saved_ev, key))
        error = float(np.max(np.abs(lhs-rhs)))
        require(error <= 1e-12, 'Supply changes saved household object ' + key)
        unchanged[key] = error
    require(abs(float(ev.births)-float(saved_ev.births)) <= 1e-12,
            'Supply changed saved births')
    demand = float(np.asarray(ev.demand_by_loc).sum())
    supply = float(np.asarray(ev.supply_by_loc).sum())
    N = supply/demand
    entry = float(sol.entry_rate)
    adjusted = float(sol.adult_entry_adjusted_birth_children)
    renewal = adjusted/(2.1*entry)-1.
    require(abs(renewal) <= RENEWAL_TOL and N > 0 and np.isfinite(N),
            'Renewal or population closure failed')
    if name == 'baseline_h0_plus10':
        require(abs(N-1.1) <= 5e-7, 'Baseline +10% scale did not imply N=1.1')
    else:
        require(abs(float(price[0])-0.821517532009075) <= 1e-13,
                'Credit price is not the selected renewal root')
    absolute_residual = N*demand-supply
    require(abs(absolute_residual) <= 1e-12, 'Absolute housing market did not clear')

    # Apply actual household mass before native renewal, estate, fiscal and
    # solvency audits. The household policies themselves are never re-solved.
    actual = copy.copy(ev)
    for key in ('g_pre', 'g_post_fertility', 'g_current', 'g_stay_distribution',
                'inherited_g_pre'):
        value = getattr(actual, key, None)
        if value is not None:
            setattr(actual, key, N*np.asarray(value))
    actual.births = N*ev.births
    actual.births_by_loc = N*np.asarray(ev.births_by_loc)
    actual.demand_by_loc = N*np.asarray(ev.demand_by_loc)
    actual.relative_market_residual = abs(absolute_residual)/supply
    packet = dict(source, parameters=P, b_grid=grid, shared=shared, solution=sol,
                  evaluation=actual, stationary_g_pre=N*pre, supply_rule=supply_rule)
    # The frozen solution packs market statistics using its old H0. Replace
    # only these derived reporting fields; household arrays stay saved.
    reported_solution = copy.copy(sol)
    reported_solution.housing_supply = np.asarray(ev.supply_by_loc)/N
    reported_solution.aggregate_housing_supply = float(reported_solution.housing_supply.sum())
    reported_solution.aggregate_housing_excess = float(demand-reported_solution.aggregate_housing_supply)
    reported_solution.best_max_abs_rel_excess = abs(reported_solution.aggregate_housing_excess)/demand
    reported_solution.best_market_metric = reported_solution.best_max_abs_rel_excess
    reported_solution.converged = True
    packet['solution'] = reported_solution
    if credit_adapter is None:
        gates = base.gates(packet, prepared, case_dir, stationary=False)
    else:
        gates = credit.credit_audits(base, packet, prepared, case_dir, credit_adapter,
                                     stationary=False)
        require(gates['credit_solvency']['status'] == 'passed' and
                gates['credit_solvency']['negative_estate_exposure_mass'] <= 2e-10,
                'Credit solvency failed at actual population')
    fiscal = gates['fiscal']['scaled_pension_budget_residual']
    require(abs(float(fiscal)) <= PAYGO_TOL and
            gates['fiscal_certificate']['fiscal_gate'] and
            gates['fiscal_certificate']['marginal_gate'],
            'Actual-population PAYGO or marginal fiscal gate failed')
    require(gates['estate']['status'] == 'funded', 'Estate funding failed')

    # Native four-year entry clock: half the newborns enter after 16 years,
    # half after 20 years. Reapply its transition operator at actual mass.
    derived = float(native.calendar_topcode_birth_accounting(
        actual.g_pre, actual.g_post_fertility, float(actual.births), P)
        ['topcode_adjusted_birth_children'])
    require(abs(derived-N*adjusted) <= max(1., N)*5e-9,
            'Topcode-adjusted actual births differ')
    queue = native.SplitBirthEntryQueue.constant_prehistory(derived)
    due, next_queue = queue.step(derived)
    half = derived/(2*2.1)
    require(len(queue.due_in_16) == 3 and len(queue.due_in_20) == 4 and
            all(abs(x-half) <= max(1., N)*2e-10 for x in
                queue.due_in_16+queue.due_in_20+next_queue.due_in_16+next_queue.due_in_20) and
            next_queue.due_in_16 == queue.due_in_16 and
            next_queue.due_in_20 == queue.due_in_20 and
            abs(due-derived/2.1) <= max(1., N)*2e-10,
            'Native 16/20 queue failed')
    following, _, deaths, mass_residual = native.advance_sequential_calendar_distribution(
        actual, np.asarray([due]), P, grid, shared)
    target = N*pre
    raw_l1 = float(np.abs(following-target).sum())
    entry_gap = float(due-N*entry)
    accounted = following-target
    accounted[:, :, :, 0] -= cal.entrant_cohort(np.asarray([entry_gap]), P, grid)
    accounted_l1 = float(np.abs(accounted).sum())
    require(accounted_l1 <= max(1., N)*5e-9 and
            abs(float(mass_residual)) <= max(1., N)*2e-8,
            'Native actual-population stationary operator failed')
    np.savez_compressed(case_dir/'actual_arrays.npz',
        g_pre=actual.g_pre, g_post_fertility=actual.g_post_fertility,
        g_current=actual.g_current, births_by_loc=actual.births_by_loc,
        demand_by_loc=actual.demand_by_loc, supply_by_loc=actual.supply_by_loc,
        following_g_pre=following, target_g_pre=target,
        g_stay_distribution=np.asarray(actual.g_stay_distribution)
        if actual.g_stay_distribution is not None else np.empty(0))
    table(case_dir/'target_fit.csv', fit_rows)
    table(case_dir/'parameters.csv', parameter_rows)
    outcome = dict(reference_label=LABEL, case=name, model_solves=0, status='passed',
        source_checkpoint_sha256=checkpoint_sha, price=float(price[0]),
        normalized_demand=demand, absolute_supply=supply, population=N,
        absolute_demand=N*demand, housing_residual=absolute_residual,
        renewal_residual=renewal, adjusted_births_per_normalized_household=adjusted,
        actual_adjusted_births=derived, actual_entrant_flow=due,
        native_queue_16_20_pass=True, native_raw_l1=raw_l1,
        native_accounted_l1=accounted_l1, native_mass_residual=float(mass_residual),
        deaths=float(deaths), paygo_residual=float(fiscal),
        negative_estate_gate_passed=True,
        negative_estate_exposure_mass=float(gates['credit_solvency']['negative_estate_exposure_mass'])
            if credit_adapter else None,
        saved_household_array_max_errors=unchanged, supply_override=supply_description,
        saved_solution_market_stats_rebased_to_supply_per_household=True,
        serialized_H0=float(P.H0[0]), serialized_xi_supply=float(P.xi_supply[0]),
        original_H0=original_h0, original_xi_supply=original_eta,
        arrays_sha256=sha(case_dir/'actual_arrays.npz'),
        fit_rows=14, parameter_rows=31,
        estate_counterparty_caveat='Provisional net-estate entrant funding and residual sink; counterparties unresolved')
    write(case_dir/'verification.json', outcome)
    return outcome


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    p = read(args.plan)
    verify_plan(p, __file__)
    require(not args.output.exists(), 'Output already exists; immutable run only')
    args.output.mkdir(parents=True)
    started = time.monotonic()
    write(args.output/'progress.json', dict(phase='authenticated_inputs', model_solves=0))
    import numpy as np
    import run_fixed_price as base
    import run_credit as credit
    manifest, contract, _, runtime, prepared, reference = base.authenticate(args.output)
    require(sha(p['reference_checkpoint_path']) == REFERENCE_SHA and
            manifest['checkpoint']['sha256'] == REFERENCE_SHA,
            'Reference identity changed during authentication')
    adapter = credit.load_module(p['credit_adapter_path'], '_block0506_credit_ge_adapter')
    with gzip.open(p['credit_ge_checkpoint_path'], 'rb') as f:
        ge = pickle.load(f)
    cal = prepared.rt['primitive'].pf.calendar
    model = prepared.rt['model']
    source_identity = verify_model_independence(model, cal)
    forbid_solves(model, cal)
    P0 = copy.deepcopy(reference['parameters'])
    original_h0 = float(P0.H0[0])
    original_eta = float(P0.xi_supply[0])
    require(abs(original_eta-.63) <= 1e-12 and
            abs(float(np.asarray(reference['evaluation'].supply_by_loc).sum())-BASELINE_STOCK) <= 1e-9,
            'Baseline supply reference differs')
    q0 = float(np.asarray(reference['solution'].p_eq)[0])
    P_up = copy.deepcopy(P0)
    P_up.H0 = 1.1*np.asarray(P0.H0)
    H_up = float(P_up.H0[0]*(P_up.user_cost_rate*q0/P_up.r_bar[0])**P_up.xi_supply[0])
    up_rule = cal.HousingSupplyRule('static-elastic', q0, H_up, original_eta)
    ref_params = copy.deepcopy(manifest['full_parameter_table'])
    require(len(ref_params) == 31, 'Missing 31 reference parameter rows')
    for row in ref_params:
        row['status'] = 'Inherited reference value; no household reoptimization'
        if row['parameter'] == 'H0':
            row['estimate'] = float(P_up.H0[0])
            row['status'] = 'Experimental absolute supply intercept +10%'
    require(any(r['parameter']=='H0' and float(r['estimate'])==float(P_up.H0[0])
                for r in ref_params), 'H0 parameter row not updated')
    baseline = replay_case('baseline_h0_plus10', reference, P_up, up_rule,
        dict(mode='static-elastic', H0_multiplier=1.1, elasticity=original_eta,
             baseline_absolute_stock=BASELINE_STOCK), prepared, args.output,
        manifest['full_target_table'], ref_params, REFERENCE_SHA, original_h0, original_eta)
    require(time.monotonic()-started < 180, 'Replay budget exceeded after baseline')
    write(args.output/'progress.json', dict(phase='baseline_complete', model_solves=0,
        elapsed_seconds=time.monotonic()-started))

    P_ge = copy.deepcopy(ge['parameters'])
    require(float(P_ge.H0[0]) == original_h0 and
            float(P_ge.xi_supply[0]) == original_eta,
            'Credit endpoint source supply primitives differ')
    qge = float(np.asarray(ge['solution'].p_eq)[0])
    fixed_rule = cal.HousingSupplyRule('fixed-stock', qge, BASELINE_STOCK, 0.)
    ge_params = rows(p['credit_ge_parameters_path'])
    require(len(ge_params)==31, 'Missing 31 GE parameter rows')
    for row in ge_params:
        if row['parameter']=='H0':
            row['status']='Serialized reference H0 retained; inactive under fixed-stock override'
        elif row['parameter']=='housing_supply_elasticity':
            row['status']='Serialized .63 retained; inactive under fixed-stock override; effective elasticity 0'
    fixed = replay_case('credit_fixed_baseline_stock', ge, P_ge, fixed_rule,
        dict(mode='fixed-stock', physical_stock=BASELINE_STOCK,
             serialized_H0_inactive=original_h0,
             serialized_supply_elasticity_inactive=original_eta,
             effective_supply_elasticity=0.), prepared, args.output,
        rows(p['credit_ge_fit_path']), ge_params, GE_SHA, original_h0, original_eta,
        adapter)
    require(time.monotonic()-started < 270, 'Replay budget exceeded after credit case')
    write(args.output/'progress.json', dict(phase='credit_complete', model_solves=0,
        elapsed_seconds=time.monotonic()-started))
    ge_receipt = read(p['credit_ge_receipt_path'])
    require(ge_receipt['checkpoint']['sha256']==GE_SHA and
            abs(float(ge_receipt['closure']['normalized_housing_demand'])-fixed['normalized_demand']) <= 1e-12,
            'Selected GE household solution differs')
    table(args.output/'comparison.csv', [
        dict(case='frozen_baseline', price=q0,
             normalized_demand=float(np.asarray(reference['evaluation'].demand_by_loc).sum()),
             absolute_supply=BASELINE_STOCK, population=1., supply_rule='original elastic'),
        dict(case='baseline_h0_plus10', price=baseline['price'],
             normalized_demand=baseline['normalized_demand'], absolute_supply=baseline['absolute_supply'],
             population=baseline['population'], supply_rule='H0 +10%, eta=.63'),
        dict(case='credit_fixed_baseline_stock', price=fixed['price'],
             normalized_demand=fixed['normalized_demand'], absolute_supply=fixed['absolute_supply'],
             population=fixed['population'], supply_rule='fixed stock at baseline Hs'),
        dict(case='credit_original_elastic_stock', price=qge,
             normalized_demand=float(ge_receipt['closure']['normalized_housing_demand']),
             absolute_supply=float(ge_receipt['closure']['absolute_housing_supply']),
             population=float(ge_receipt['closure']['population_scale']),
             supply_rule='original elastic, selected exact repeat')])
    runtime.verify_sources(dict(contract, objective=manifest['objective']))
    verify_plan(p, __file__)
    require(time.monotonic()-started < 300, 'Five-minute replay budget exceeded')
    write(args.output/'verification.json', dict(status='passed', model_solves=0,
        reference_label=LABEL, cases=[baseline, fixed],
        original_credit_ge_receipt_sha256=sha(p['credit_ge_receipt_path']),
        source_plot_directories=p['source_plot_directories'],
        source_fit_tables=[p['reference_fit_path'], p['credit_ge_fit_path']],
        source_identity=source_identity,
        reference_manifest_sha256=MANIFEST_SHA, plan_sha256=sha(args.plan),
        elapsed_seconds=time.monotonic()-started))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        if '--output' in sys.argv:
            destination = Path(sys.argv[sys.argv.index('--output')+1])
            failure = destination/'failure.json' if destination.is_dir() else destination.with_suffix('.failure.json')
            write(failure, dict(status='failed', model_solves=0, error=str(exc),
                traceback=traceback.format_exc()))
        raise
