#!/usr/bin/env python3
"""Ten bounded prescribed-price lifecycle solves on one frozen wealth grid."""
import argparse
import copy
import csv
import gzip
import json
import math
import os
from pathlib import Path
import pickle
import signal
import subprocess
import sys
import time
import traceback

import run_credit as credit
import run_fixed_price as base
import natural_credit as adapter

LABEL = base.LABEL
FACTORS = (.98, .99, 1., 1.01, 1.02)
CASES = [('grid_control', 1., 'reference'), ('credit', 1., 'credit')] + [
    (('%s_%03d' % (regime, round(1000 * factor))), factor, regime)
    for factor in FACTORS if factor != 1. for regime in ('reference', 'credit')]
MODE = adapter.MODE


def verify_plan(path):
    base.require(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    p = base.read(path)
    base.require(p['schema'] == 'block0506_fixed_grid_price_elasticity_v1' and
                 p['reference_label'] == LABEL and p['reference_manifest_sha256'] == base.MANIFEST_SHA == base.sha(base.MANIFEST),
                 'Wrong reference/manifest')
    base.require(p['driver_sha256'] == base.sha(__file__) and p['credit_mode'] == MODE,
                 'Driver or credit mode differs')
    base.require(p['maximum_lifecycle_solves'] == 10 and p['case_seconds'] == 300 and
                 p['total_seconds'] == 2400 and p['threads'] == 1 and p['memory_gib'] == 16 and
                 p['renormalize_fertility'] is False and p['factors'] == list(FACTORS), 'Budget/economic contract differs')
    base.require(p['economic_changes'] == adapter.CHANGES and p['psi_child_fixed'] == 0.1355551166583114,
                 'Undeclared economic change')
    pre = base.read(p['preflight_receipt_path'])
    base.require(pre['status'] == 'ready' and pre['model_solves'] == 0 and pre['grid_nodes'] == 262 and
                 pre['fixed_common_grid_all_prices'] and all(r['support_ok'] for r in pre['rows']),
                 'Five-price preflight not ready')
    pins = {}
    for row in p['files']:
        target = str(Path(row['path']).resolve())
        base.require(target not in pins and base.sha(target) == row['sha256'], 'Missing/changed source/input: ' + target)
        pins[target] = row['sha256']
    required = [__file__, p['preflight_receipt_path'], p['grid_path']]
    required += [str(Path(__file__).parent / n) for n in ('run_credit.py', 'run_fixed_price.py', 'natural_credit.py', 'preflight.py')]
    required += [p['old_' + regime + '_' + extension + '_path'] for regime in ('reference', 'credit')
                 for extension in ('receipt', 'checkpoint', 'fit', 'parameters')]
    base.require(all(str(Path(x).resolve()) in pins for x in required), 'Required source or exact-replay input unpinned')
    for f in Path(__file__).parent.glob('*.py'):
        base.require(str(f.resolve()) in pins, 'Unpinned overlay source')
        base.require(pre['source_sha256'].get(f.name) == base.sha(f),
                     'Python source differs from passing zero-solve preflight: ' + f.name)
    base.require(pre['grid_sha256'] == base.sha(p['grid_path']) ==
                 'a3b81c459290554f634d0c28409ab846364c0ecf7270bd164f975f092901cbe2' and
                 pre['source_sha256']['natural_credit.py'] == base.sha(Path(__file__).parent / 'natural_credit.py'),
                 'Preflight used other grid/adapter')
    base.require(base.sha(p['old_reference_checkpoint_path']) == p['old_reference_checkpoint_sha256'] ==
                 'f1c91fabca457fccfab22dbca33fae7edc750e4171b0f9796924516ed40c417e' and
                 base.sha(p['old_credit_checkpoint_path']) == p['old_credit_checkpoint_sha256'] ==
                 '7b0b92f0e2ed33cfc80a3c63c2c15d6417b14ab88db3930284e440e94e2e515c',
                 'Prior exact-replay checkpoint identity differs')
    return p


def progress(out, phase, **items):
    base.write(out / 'progress.json', dict(reference_label=LABEL, phase=phase, time_epoch=time.time(), **items))


def exact_prior(packet, fits, parameters, p, regime, out):
    import e5f_current_transition_runtime as native
    old_receipt = base.read(p['old_' + regime + '_receipt_path'])
    base.require(old_receipt['status'] == 'passed' and old_receipt['case'] ==
                 ('grid_control' if regime == 'reference' else 'credit') and
                 old_receipt['checkpoint']['sha256'] == base.sha(p['old_' + regime + '_checkpoint_path']),
                 'Prior control checkpoint/receipt differs')
    with gzip.open(p['old_' + regime + '_checkpoint_path'], 'rb') as stream:
        prior = pickle.load(stream)
    result = native.compare_arrays(prior, packet)
    base.write(out / 'exact_replay_arrays.json', result)
    bad = [k for k, row in result['arrays'].items() if row.get('status') != 'compared' or
           not row.get('exact') or not row.get('finite')]
    base.require(not bad, 'q0 exact numeric replay failed: ' + ', '.join(bad[:10]))
    old_fit = base.csv_rows(p['old_' + regime + '_fit_path'])
    base.require(len(old_fit) == len(fits) == 14, 'Prior 14-row fit missing')
    for actual, old in zip(fits, old_fit):
        for key in ('moment', 'target', 'model', 'gap', 'weight', 'loss_contribution'):
            base.require(actual[key] == old[key] if key == 'moment' else
                         (actual[key] == old[key] or actual[key] != '' and old[key] != '' and
                          float(actual[key]) == float(old[key])), 'q0 target replay differs: ' + actual['moment'])
    old_parameters = base.csv_rows(p['old_' + regime + '_parameters_path'])
    base.require(len(parameters) == len(old_parameters) == 31 and all(
        a['parameter'] == b['parameter'] and float(a['estimate']) == float(b['estimate'])
        for a, b in zip(parameters, old_parameters)), 'q0 parameter replay differs')
    return dict(status='passed', arrays_exact=result['array_count'], target_rows_exact=14,
                parameter_estimates_exact=31)


def child(args, p):
    import numpy as np
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    base.require(base.read(out.parent / 'launch.json')['plan_sha256'] == base.sha(args.plan), 'Launch contract missing')
    index = [x[0] for x in CASES].index(args.case)
    if index:
        previous = base.read(out.parent / 'latest_completed.json')['completed']
        base.require([r['case'] for r in previous] == [x[0] for x in CASES[:index]], 'Prior case order differs')
        for row in previous:
            receipt = out.parent / row['case'] / 'receipt.json'
            base.require(base.sha(receipt) == row['receipt_sha256'] and base.read(receipt)['status'] == 'passed',
                         'Prior case receipt changed')
    progress(out, 'authenticate')
    manifest, contract, objective, runtime, prepared, reference = base.authenticate(out)
    base.require(manifest['checkpoint']['sha256'] == p['reference_checkpoint_sha256'], 'Wrong reference checkpoint')
    P = copy.deepcopy(reference['parameters'])
    old_grid = np.asarray(reference['b_grid'])
    initial_parameters = base.actual_parameters(prepared, P, old_grid)
    base.require(len(initial_parameters) == 31 and all(
        float(initial_parameters[r['parameter']]) == float(r['estimate']) for r in manifest['full_parameter_table']),
        'Reference parameter table changed')
    P.native_inherited_distribution_evidence_dir = str(out / 'inherited_state_diagnostics')
    original_public = base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    reference_pin = base.serialized(vars(reference['parameters']))
    inherited_pin = base.serialized(reference['stationary_g_pre'])
    grid, impact_pre, numerical, embedding = credit.embed_credit_grid(base, P, reference, p['grid_path'])
    base.require(len(grid) == 262 and base.serialized(grid) == base.serialized(base.read(p['grid_path'])['grid']),
                 'Common grid differs')
    base.write(out / 'grid_embedding.json', embedding)
    q0 = float(np.asarray(reference['solution'].p_eq)[0])
    price = np.asarray([q0 * args.factor])
    install = None
    if args.regime == 'credit':
        install = adapter.install(prepared=prepared, reference=reference, P=P, grid=grid,
            price=float(price[0]), output=out, overlay_root=Path(__file__).parent, credit_mode=MODE)
        base.require(install['status'] == 'installed' and install['economic_changes'] == p['economic_changes'],
                     'Credit adapter installation changed')
    actual_public = base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')})
    changed = {k: dict(reference=original_public.get(k), experiment=actual_public.get(k))
               for k in set(original_public) | set(actual_public) if original_public.get(k) != actual_public.get(k)}
    expected = dict(numerical)
    if args.regime == 'credit':
        expected['native_due_stayer_credit'] = False
    base.require(set(changed) == set(expected) and all(actual_public[k] == value for k, value in expected.items()),
                 'Unplanned economic/numerical parameter change')
    base.require(float(P.psi_child) == p['psi_child_fixed'] and
                 base.actual_parameters(prepared, P, grid) == dict(initial_parameters, wealth_grid_nodes=262) and
                 base.serialized(vars(reference['parameters'])) == reference_pin and
                 base.serialized(reference['stationary_g_pre']) == inherited_pin,
                 'Reference preferences/earnings/fiscal/entry/supply changed')
    rt = prepared.rt
    model, cal = rt['model'], rt['primitive'].pf.calendar
    shared = model.precompute_shared(P, grid)
    progress(out, 'one_lifecycle_solve', price=float(price[0]), deadline_epoch=args.deadline)
    base.require(time.time() < args.deadline, 'Case deadline before solve')
    start = time.monotonic()
    sol = model.solve_markov_income_at_prices(price, P, grid, SD=shared, verbose=False, fast_stats=False)
    duration = time.monotonic() - start
    base.require(time.time() < args.deadline, 'Case deadline during solve')
    P._fert2_probs = sol.fert2_probs.copy()
    policy = cal.policy_from_solution(sol, price, P, grid, shared)
    pre, reconstruction = cal.reconstruct_stationary_pre_fertility(sol, policy, P, grid, shared)
    runtime.require_abs_gate(reconstruction['stationary_post_fertility_nesting_l1'], 5e-9, 'Cohort reconstruction')
    runtime.require_abs_gate(reconstruction['stationary_feasibility_projection_mass'], 0., 'Cohort projection')
    supply = cal.HousingSupplyRule('static-elastic', float(price[0]),
        float(P.H0[0] * (P.user_cost_rate * price[0] / P.r_bar[0]) ** P.xi_supply[0]), float(P.xi_supply[0]))
    ev = cal.evaluate_period(price, pre, P, grid, shared, cal.SolveCounter(),
                             supply_rule=supply, supplied_policy=policy)
    packet = dict(parameters=P, b_grid=grid, shared=shared, solution=sol, evaluation=ev,
                  stationary_g_pre=pre, supply_rule=supply, demographic_seed=reference.get('demographic_seed'))
    gates = (credit.credit_audits(base, packet, prepared, out, adapter, True) if args.regime == 'credit'
             else base.gates(packet, prepared, out, stationary=True))
    fertility = {kind: rt['observe_initial_fertility'](ev, P, age_projection=kind)
                 for kind in ('uniform_birth_time', 'constant_post_cell')}
    housing = rt['observe_initial_housing_wealth'](ev, P, grid, shared, diagnostic_enabled=True,
        age_projection='uniform_within_age_cell', diagnostic_allow_family_proxies=True,
        include_wealth=True, include_birth_response=True)
    recent = rt['observe_recent_parent_flow'](ev, P, diagnostic_enabled=True,
        snapshot=rt['SNAPSHOT'], age_projection=rt['AGE_PROJECTION'], diagnostic_allow_residence_proxy=True,
        input_provenance=dict(case_id=args.case, reference_checkpoint_sha256=manifest['checkpoint']['sha256']))
    completed_fertility = float(rt['chain'].extract_moments(sol, P)['tfr'])
    fits = runtime.score_targets(objective, fertility, housing, recent['model_value'], completed_fertility)
    base.require(len(fits) == 14, 'Full fit table missing')
    for row in fits:
        if row['moment'] == 'initial_normalization':
            row['role'] = 'reference replacement benchmark; no normalization in elasticity experiment'
    parameters = copy.deepcopy(manifest['full_parameter_table'])
    base.require(len(parameters) == 31, 'Full parameter table missing')
    for row in parameters:
        row['status'] = 'Fixed reference economic value; no fertility normalization'
        if row['parameter'] == 'wealth_grid_nodes':
            row['estimate'] = 262
            row['status'] = 'Numerical 262-node grid common to every price/regime; original atoms exact'
        if args.regime == 'credit' and row['parameter'] == 'financed_share':
            row['status'] = 'Reference value for provenance; artificial LTV inactive in solvency-only regime'
    base.table(out / 'target_fit.csv', fits)
    base.table(out / 'parameters.csv', parameters)
    base.write(out / 'observers.json', cal.jsonable(dict(fertility=fertility,
        housing_wealth=housing, recent_parent=recent)))
    cohort = credit.aggregates(base, ev, P, grid, model)
    impact_P = copy.deepcopy(P)
    impact_shared = model.precompute_shared(impact_P, grid)
    impact = cal.evaluate_period(price, impact_pre, impact_P, grid, impact_shared,
        cal.SolveCounter(), supply_rule=supply, supplied_policy=policy)
    impact_out = out / 'baseline_state_impact'
    impact_out.mkdir()
    impact_gates = (credit.credit_audits(base, dict(packet, parameters=impact_P, shared=impact_shared,
        evaluation=impact, stationary_g_pre=impact_pre), prepared, impact_out, adapter, False)
        if args.regime == 'credit' else base.gates(dict(packet, parameters=impact_P, shared=impact_shared,
        evaluation=impact, stationary_g_pre=impact_pre), prepared, impact_out, stationary=False))
    impact_summary = credit.aggregates(base, impact, impact_P, grid, model)
    base.write(impact_out / 'summary.json', impact_summary)
    progress(out, 'standard_17_plot_rendering')
    rt['audit'].standard_diagnostics(packet, out, validate_production_young=False)
    plots = sorted(f.name for f in (out / 'standard_diagnostics').glob('*.png'))
    base.require(plots == sorted(manifest['standard_diagnostic_names']), 'Standard plot set differs')
    replay = exact_prior(packet, fits, parameters, p, args.regime, out) if args.factor == 1. else None
    base.require(base.serialized({k: v for k, v in vars(P).items() if not k.startswith('_')}) == actual_public and
                 base.serialized(vars(reference['parameters'])) == reference_pin and
                 base.serialized(reference['stationary_g_pre']) == inherited_pin and
                 base.serialized(grid) == base.serialized(base.read(p['grid_path'])['grid']) and
                 base.serialized(impact_pre) == base.serialized(credit.embed_credit_grid(
                     base, copy.deepcopy(reference['parameters']), reference, p['grid_path'])[1]),
                 'Frozen input or parameters changed during solve/reporting')
    runtime.verify_sources(dict(contract, objective=manifest['objective']))
    base.require(verify_plan(args.plan) == p and time.time() < args.deadline, 'Plan/source/deadline changed')
    with gzip.open(out / 'conditional_cohort_state.pkl.gz', 'wb', compresslevel=1) as stream:
        pickle.dump(packet, stream, protocol=5)
    fit_by_name = {row['moment']: float(row['model']) for row in fits if row['model'] != ''}
    receipt = dict(status='passed', reference_label=LABEL, case=args.case, regime=args.regime,
        price_factor=args.factor, price=float(price[0]), rent=float(P.user_cost_rate * price[0]),
        lifecycle_solves=1, lifecycle_solve_seconds=duration, plan_sha256=base.sha(args.plan),
        reference_manifest_sha256=base.MANIFEST_SHA, reference_checkpoint=manifest['checkpoint'],
        fixed_psi=float(P.psi_child), normalization_performed=False,
        completed_fertility=completed_fertility, mean_first_birth_age=fit_by_name['nchs_mean_age'],
        economic_changes=([dict(object='prescribed_house_price_and_mapped_rent',
                                factor=args.factor, classification='proposed_price_shock')]
                          if args.factor != 1. else []) +
                         ([dict(object='borrowing_constraints', detail=p['economic_changes'][0],
                                classification='proposed_credit_regime')]
                          if args.regime == 'credit' else []),
        parameter_changes=changed, numerical_parameter_changes={k: changed[k] for k in numerical},
        credit_parameter_changes={k: changed[k] for k in changed if k not in numerical},
        grid_embedding=embedding, installation=install, exact_q0_replay=replay,
        cohort_summary=cohort, impact_summary=impact_summary, cohort_gates=gates,
        impact_gates=impact_gates, standard_plot_count=len(plots),
        checkpoint=dict(path=str(out / 'conditional_cohort_state.pkl.gz'),
                        sha256=base.sha(out / 'conditional_cohort_state.pkl.gz')),
        market_clearing_certified=False, demographic_renewal_certified=False,
        interpretation='Prescribed house price and mapped rent. Impact uses one frozen inherited distribution; cohort recomputes normalized-entry composition. No GE or fertility normalization.')
    base.write(out / 'receipt.json', cal.jsonable(receipt))
    progress(out, 'complete', receipt_sha256=base.sha(out / 'receipt.json'))


def make_comparison(out, records):
    receipts = {(r['regime'], r['factor']): base.read(out / r['case'] / 'receipt.json') for r in records}
    names = ('births_per_household', 'first_births', 'second_births', 'third_bin_entries',
             'ownership_rate', 'rooms_per_household')
    rows = []
    for regime in ('reference', 'credit'):
        for scope in ('impact', 'cohort'):
            for factor in FACTORS:
                receipt = receipts[(regime, factor)]
                summary = receipt[scope + '_summary']
                for name in names:
                    rows.append(dict(regime=regime, scope=scope, price_factor=factor,
                                     outcome=name, value=(summary[name]/summary['household_mass']
                                                          if name in ('first_births', 'second_births', 'third_bin_entries')
                                                          else summary[name]),
                                     unit=('births per household' if name in ('births_per_household',
                                           'first_births', 'second_births', 'third_bin_entries') else
                                           'rooms per household' if name == 'rooms_per_household' else 'share'),
                                     prescribed_price=receipt['price'], mapped_rent=receipt['rent']))
                if scope == 'cohort':
                    rows.extend(dict(regime=regime, scope=scope, price_factor=factor,
                                     outcome=name, value=receipt[name], unit=unit,
                                     prescribed_price=receipt['price'], mapped_rent=receipt['rent'])
                        for name, unit in (('completed_fertility', 'children'), ('mean_first_birth_age', 'years')))
    path = out / 'comparison.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    by_key = {(r['regime'], r['scope'], r['outcome'], r['price_factor']): r['value'] for r in rows}
    elasticities = []
    for regime in ('reference', 'credit'):
        for scope in ('impact', 'cohort'):
            for outcome in names + (('completed_fertility', 'mean_first_birth_age') if scope == 'cohort' else ()):
                values = {f: by_key[(regime, scope, outcome, f)] for f in FACTORS}
                for step in (.01, .02):
                    lo, mid, hi = values[1-step], values[1.], values[1+step]
                    def ratio(a, b, fa, fb):
                        return (math.log(b)-math.log(a))/(math.log(fb)-math.log(fa)) if a > 0 and b > 0 else None
                    elasticities.append(dict(regime=regime, scope=scope, outcome=outcome, step=step,
                        central_log_elasticity=ratio(lo, hi, 1-step, 1+step),
                        lower_one_sided_log_elasticity=ratio(lo, mid, 1-step, 1.),
                        upper_one_sided_log_elasticity=ratio(mid, hi, 1., 1+step)))
    with (out / 'elasticities.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(elasticities[0]))
        writer.writeheader(); writer.writerows(elasticities)
    return dict(comparison_sha256=base.sha(path), elasticities_sha256=base.sha(out / 'elasticities.csv'),
                outcome_rows=len(rows), elasticity_rows=len(elasticities))


def controller(args, p):
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    started = time.time(); deadline = started + p['total_seconds']
    base.write(out / 'launch.json', dict(plan=p, plan_sha256=base.sha(args.plan),
        started_epoch=started, deadline_epoch=deadline, slurm_job=os.environ['SLURM_JOB_ID']))
    records = []
    base.write(out / 'latest_completed.json', dict(completed=[], lifecycle_solves=0))
    base.write(out / 'best_so_far.json', dict(status='no_case_completed',
        criterion='smallest absolute prescribed-price housing excess; diagnostic only, no GE certification'))
    for name, factor, regime in CASES:
        base.require(verify_plan(args.plan) == p and time.time() < deadline, 'Plan changed or total deadline reached')
        case_deadline = min(deadline, time.time() + p['case_seconds'])
        command = [sys.executable, str(Path(__file__).resolve()), '--plan', str(args.plan),
                   '--output', str(out / name), '--case', name, '--factor', repr(factor),
                   '--regime', regime, '--deadline', repr(case_deadline)]
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   NUMEXPR_NUM_THREADS='1', NUMBA_NUM_THREADS='1', MPLBACKEND='Agg', PYTHONDONTWRITEBYTECODE='1')
        with (out / (name + '.log')).open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                       env=env, start_new_session=True)
            try:
                while process.poll() is None:
                    progress(out, 'case_running', case=name, factor=factor, completed_cases=len(records),
                             elapsed_seconds=time.time()-started, case_deadline_epoch=case_deadline)
                    if time.time() >= case_deadline:
                        raise TimeoutError('Case/global deadline: ' + name)
                    time.sleep(2)
                base.require(process.returncode == 0, 'Case failed; no retry: ' + name)
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGTERM)
                    try: process.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL); process.wait()
        path = out / name / 'receipt.json'; receipt = base.read(path)
        base.require(receipt['status'] == 'passed' and receipt['lifecycle_solves'] == 1 and
                     receipt['plan_sha256'] == base.sha(args.plan) and receipt['standard_plot_count'] == 17,
                     'Incomplete case receipt')
        row = dict(case=name, factor=factor, regime=regime, receipt_sha256=base.sha(path),
                   abs_market_residual=abs(receipt['cohort_gates']['relative_market_residual']))
        records.append(row)
        base.write(out / 'latest_completed.json', dict(completed=records, lifecycle_solves=len(records),
            remaining_cases=10-len(records)))
        best = min(records, key=lambda x: x['abs_market_residual'])
        base.write(out / 'best_so_far.json', dict(status='case_completed', **best,
            criterion='smallest absolute prescribed-price housing excess; diagnostic only, no GE certification'))
    comparison = make_comparison(out, records)
    base.require(time.time() < deadline, 'Total deadline exceeded during comparison')
    base.write(out / 'completed.json', dict(status='passed', reference_label=LABEL,
        lifecycle_solves=10, elapsed_seconds=time.time()-started, normalization_performed=False,
        market_clearing_certified=False, plan_sha256=base.sha(args.plan), **comparison))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', choices=[x[0] for x in CASES])
    parser.add_argument('--factor', type=float)
    parser.add_argument('--regime', choices=('reference', 'credit'))
    parser.add_argument('--deadline', type=float)
    args = parser.parse_args()
    args.plan = args.plan.resolve(); args.output = args.output.resolve()
    p = verify_plan(args.plan)
    base.require(not args.output.exists(), 'Versioned output exists')
    try:
        if args.case:
            expected = next((f, r) for name, f, r in CASES if name == args.case)
            base.require((args.factor, args.regime) == expected and args.deadline is not None and
                         time.time() < args.deadline <= time.time()+p['case_seconds'], 'Unauthorized child case')
            child(args, p)
        else:
            base.require(args.factor is None and args.regime is None and args.deadline is None,
                         'Controller arguments differ')
            controller(args, p)
    except BaseException as exc:
        if args.output.exists():
            base.write(args.output / 'failure.json', dict(status='failed', error_type=type(exc).__name__,
                error=str(exc), traceback=traceback.format_exc(), retries=0, time_epoch=time.time()))
        raise


if __name__ == '__main__':
    main()
