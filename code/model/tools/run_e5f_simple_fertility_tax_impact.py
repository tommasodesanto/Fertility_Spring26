#!/usr/bin/env python3
"""One-date balanced-budget tax1/tax2 impacts for the selected fertility nests.

Both endpoints rebate all contemporaneous property-tax revenue equally per
household decision unit and inherit the same calibrated population. No future
population advance, recalibration, kernel replacement or figures are performed.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import math
import os
from pathlib import Path
import pickle
import threading
import time
import traceback

import numpy as np
import run_e5f_simple_fertility_policy_impact as common
import run_e5f_post2023_rebated_property_tax_smoke as tax

adapter, audit, solver = common.adapter, common.audit, common.solver
calendar = common.policy.calendar
CASES = ('tax1-equal-rebate', 'tax2-equal-rebate')
SCHEMA = 'e5f_simple_fertility_tax_impact_v1'


def verify_sources(contract):
    if str(Path(__file__).resolve()) not in contract['source_hashes']:
        raise RuntimeError('Tax impact wrapper must be source-pinned')
    common.verify_sources(contract)


def load_contract(path, expected):
    adapter.verify(path, expected)
    c = adapter.read_json(path)
    if c['schema'] != SCHEMA or c['case_seconds'] != 1800:
        raise RuntimeError('Wrong tax impact schema or time budget')
    if time.time() > float(c['launch_deadline_epoch']):
        raise RuntimeError('Tax impact launch deadline exceeded')
    for key in ('selected_summary', 'output_root', 'baseline_smoke_receipt'):
        if not Path(c[key]).is_absolute():
            raise RuntimeError(f'Expected absolute contract path: {key}')
    verify_sources(c)
    return c


def verify_receipt(path, contract, *, expected_hash=None, same_contract=None):
    if expected_hash is not None:
        adapter.verify(path, expected_hash)
    receipt = adapter.read_json(path)
    if receipt['status'] != 'complete':
        raise RuntimeError(f'Incomplete prerequisite: {path}')
    if same_contract is not None and receipt['contract_sha256'] != same_contract:
        raise RuntimeError('Tax1 prerequisite uses a different contract')
    if receipt['selected_summary_sha256'] != contract['selected_summary_sha256']:
        raise RuntimeError('Prerequisite uses a different selected calibration')
    if receipt['selected_checkpoint_sha256'] != contract['checkpoint_sha256']:
        raise RuntimeError('Prerequisite uses a different inherited checkpoint')
    for relative, digest in receipt['artifact_sha256'].items():
        adapter.verify(path.parent / relative, digest)
    return receipt


def fiscal_checks(evaluation, parameters, *, reference=None):
    ledger = tax.fiscal_ledger(evaluation, parameters)
    if any(not math.isfinite(float(value)) for value in ledger.values()):
        raise RuntimeError('Nonfinite fiscal ledger')
    if reference is not None:
        for key, value in ledger.items():
            if not math.isclose(value, reference[key], rel_tol=0, abs_tol=2e-10):
                raise RuntimeError(f'Restored parameters disagree with solved fiscal ledger: {key}')
    supply = float(evaluation.supply_by_loc[0])
    signed_market = float(evaluation.demand_by_loc[0] - supply) / max(abs(supply), 1e-14)
    scale = max(abs(ledger['property_tax_revenue']), abs(ledger['equal_transfer_outlays']), .1)
    residual = np.array([signed_market, ledger['government_budget_residual'] / scale])
    if not np.isfinite(residual).all() or np.max(np.abs(residual)) > 1e-4:
        raise RuntimeError(f'Coupled market/fiscal root gate failed: {residual.tolist()}')
    if abs(ledger['government_budget_residual']) > 2.5e-5:
        raise RuntimeError('Absolute government budget gate failed')
    if not math.isfinite(evaluation.relative_market_residual) or evaluation.relative_market_residual > 2e-4:
        raise RuntimeError('Market gate failed')
    return ledger, residual


def evaluation_gates(current, inherited):
    e = current['evaluation']
    health = calendar.distribution_health(dict(pre=e.g_pre, post=e.g_post_fertility, current=e.g_current))
    if health['nonfinite_distribution_count'] or health['min_distribution_mass'] < -1e-14:
        raise RuntimeError('Distribution health gate failed')
    mass = float(inherited.g_pre.sum())
    gaps = {key: abs(float(getattr(e, key).sum()) - mass)
            for key in ('g_pre', 'g_post_fertility', 'g_current')}
    if any(not math.isfinite(value) or value > 2e-10 for value in gaps.values()):
        raise RuntimeError(f'Impact population denominator changed: {gaps}')
    if not math.isfinite(e.feasibility_projection_mass) or e.feasibility_projection_mass > 1e-6:
        raise RuntimeError('Feasibility projection gate failed')
    return dict(distribution_health=health, mass_gaps=gaps,
                joint_probabilities=common.probability_checks(e))


def lifecycle_quantities(evaluation, P, folder):
    """Whole-node ownership and actual occupied rooms; no ACS age realignment."""
    g, rental_h = evaluation.g_current, evaluation.policy.hR_pol[:, 0]
    def group(ages, children=slice(None), parity=slice(None)):
        selected = g[:, :, :, ages, :, parity, children]
        mass = float(selected.sum())
        owners = float(selected[:, 1:].sum())
        rental = selected[:, 0]
        housing = float(np.sum(rental * rental_h[:, :, ages, :, parity, children]))
        housing += sum(float(selected[:, t].sum()) * float(P.H_own[t-1])
                       for t in range(1, selected.shape[1]))
        return dict(household_mass=mass, owner_rate=owners/mass if mass else None,
                    mean_rooms=housing/mass if mass else None)
    rows = [dict(age_node=float(P.age_start + j*P.da), **group(slice(j, j+1)))
            for j in range(P.J)]
    common.policy.baseline.write_csv(folder / 'lifecycle.csv', rows)
    lo, hi = solver.age_to_index(P, 25), solver.age_to_index(P, 34)
    young = group(slice(lo, hi+1))
    if young['owner_rate'] is None:
        raise RuntimeError('Young ownership has no household support')
    groups = []
    for age_label, ages in (('all_ages', slice(None)), ('young_whole_nodes_25_34', slice(lo, hi+1))):
        for label, children, parity in (
                ('dependent_children', slice(1, None), slice(None)),
                ('no_dependent_children', slice(0, 1), slice(None)),
                ('childless_parity_zero', slice(None), slice(0, 1))):
            groups.append(dict(age_group=age_label, family_group=label,
                               **group(ages, children, parity)))
    common.policy.baseline.write_csv(folder / 'family_group_quantities.csv', groups)
    audit.save_json(folder / 'age_measurement.json', dict(
        definition='unchanged model.age_to_index(P,25) through model.age_to_index(P,34), inclusive',
        age_nodes=[rows[j]['age_node'] for j in range(lo, hi+1)],
        no_annual_age_acs_realignment=True,
        family_definition='no dependent children includes empty nesters; parity-zero childlessness separately reported'))
    return dict(young_ownership_whole_nodes=float(young['owner_rate']),
                young_mean_rooms_whole_nodes=float(young['mean_rooms']))


def solve_endpoint(packet, prepared, inherited, name, folder, progress):
    case = tax.CASES[name]
    P = copy.deepcopy(packet['parameters'])
    state = copy.deepcopy(inherited)
    state.initial_policy = None
    common.assert_model(P, prepared.supply_rule)
    counter = calendar.SolveCounter()
    progress['phase'] = 'coupled_price_transfer_root'
    solved = tax.solve_joint_rebated_period(
        state=state, P=P, b_grid=prepared.b_grid, supply_rule=prepared.supply_rule,
        counter=counter, case=case, tolerance=1e-4)
    # The root caches evaluations while mutating P for trial points. Restore
    # the selected point before using parameters with its returned evaluation.
    tax.set_tax_policy(P, case, solved.transfer)
    common.assert_model(P, prepared.supply_rule)
    ledger, residual = fiscal_checks(solved.evaluation, P, reference=solved.ledger)
    if not np.allclose(residual, solved.residual, rtol=0, atol=2e-10):
        raise RuntimeError('Independent coupled-root residual mismatch')
    audit.save_json(folder / 'coupled_root.json', dict(
        status='root_passed_replay_pending', price=float(solved.price), transfer=float(solved.transfer),
        ledger=ledger, normalized_residuals=residual.tolist(),
        joint_iterations=int(solved.joint_iterations), joint_model_evaluations=int(solved.joint_model_evaluations)))
    progress['phase'] = 'fresh_fixed_price_transfer_replay'
    replay_P = copy.deepcopy(P)
    shared = solver.precompute_shared(replay_P, prepared.b_grid)
    replay = calendar.evaluate_period(
        np.array([solved.price]), state.g_pre.copy(), replay_P, prepared.b_grid,
        shared, counter, supply_rule=prepared.supply_rule)
    reproduction = common.baseline_comparison(replay, solved.evaluation, policies=True)
    common.assert_model(replay_P, prepared.supply_rule)
    replay_ledger, replay_residual = fiscal_checks(replay, replay_P, reference=ledger)
    current = dict(parameters=replay_P, b_grid=prepared.b_grid, evaluation=replay, shared=shared,
                   supply_rule=prepared.supply_rule, state=state, calendar_year=2023,
                   policy_case=name, scope='one_date_rebated_tax_impact_common_inherited_population')
    checkpoint = folder / 'dated_state.pkl.gz'
    with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
        pickle.dump(current, stream, protocol=5)
    progress['phase'] = 'replayed_endpoint_audits'
    gates = evaluation_gates(current, inherited)
    budget = audit.budget_audit(current, folder)
    arrays = audit.policy_array_audit(current, folder)
    if not math.isfinite(float(budget['budget_excess_mass'])) or budget['budget_excess_mass'] > 2e-10:
        raise RuntimeError('Occupied budget gate failed')
    if arrays['occupied_negative_steps'] != 0:
        raise RuntimeError('Occupied value monotonicity gate failed')
    for bounds in arrays['probabilities'].values():
        if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
            raise RuntimeError('Marginal probability gate failed')
    values = common.quantities(replay, replay_P)
    values.update(lifecycle_quantities(replay, replay_P, folder))
    if any(not math.isfinite(float(value)) for value in values.values()):
        raise RuntimeError('Nonfinite impact quantity')
    audit.save_json(folder / 'quantities.json', values)
    audit.save_json(folder / 'fiscal_ledger.json', replay_ledger)
    result = dict(status='complete', policy=name, annual_property_tax_rate=case.annual_tax_rate,
                  property_tax_rate_period_units=float(replay_P.tau_H),
                  quantities=values, ledger=replay_ledger, normalized_root_residuals=replay_residual.tolist(),
                  root_tolerance=1e-4, fiscal_absolute_tolerance=2.5e-5, market_tolerance=2e-4,
                  fiscal_denominator='current household decision-unit mass, including renters and owners',
                  fiscal_revenue_base='all occupied rental and owner housing services at current asset price',
                  joint_iterations=int(solved.joint_iterations), joint_model_evaluations=int(solved.joint_model_evaluations),
                  bellman_solves=int(counter.bellman), fresh_fixed_price_reproduction=reproduction,
                  age_measurement=adapter.read_json(folder / 'age_measurement.json'),
                  family_group_quantities=adapter.read_csv(folder / 'family_group_quantities.csv'),
                  model_flags={key: bool(getattr(replay_P, key)) for key in (
                      'joint_nested_choice', 'fertility_nest_choice', 'two_shock_choice', 'exhaustive_saving_control')},
                  budget=budget, policy_arrays=arrays, gates=gates, checkpoint_sha256=adapter.digest(checkpoint))
    audit.save_json(folder / 'endpoint_receipt.json', result)
    progress['completed_evaluations'] = int(solved.joint_model_evaluations) + 1
    return result


def execute(args):
    c = load_contract(args.contract, args.contract_sha256)
    root = Path(c['output_root'])
    folder = root / args.case
    folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    progress = dict(case=args.case, phase='verifying_inputs', completed_evaluations=0)
    stop = threading.Event()
    def heartbeat():
        while not stop.is_set():
            elapsed = time.monotonic() - started
            audit.save_json(folder / 'heartbeat.json', dict(progress, epoch=time.time(), elapsed_seconds=elapsed))
            if elapsed >= 1800:
                audit.save_json(folder / 'failure.json', dict(error='Declared 1800-second case cap exceeded', **progress))
                os._exit(124)
            stop.wait(30)
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        audit.save_json(folder / 'latest_completed_case.json', dict(status='no_endpoint_completed'))
        audit.save_json(folder / 'best_so_far.json', dict(scope='fixed selected calibration; no search', loss=common.SELECTED_LOSS))
        baseline = verify_receipt(Path(c['baseline_smoke_receipt']), c,
                                  expected_hash=c['baseline_smoke_receipt_sha256'])
        if baseline['case'] != 'smoke':
            raise RuntimeError('Expected the completed baseline exact-loop smoke')
        if args.case == 'tax2-equal-rebate':
            prerequisite = verify_receipt(root / 'tax1-equal-rebate' / 'receipt.json', c,
                                          same_contract=args.contract_sha256)
            if prerequisite['case'] != 'tax1-equal-rebate':
                raise RuntimeError('Tax1 exact-loop smoke prerequisite has wrong case')
        packet, prepared, inherited = common.selected_packet(c)
        result = solve_endpoint(packet, prepared, inherited, args.case, folder, progress)
        audit.save_json(folder / 'latest_completed_case.json', result)
        audit.save_json(folder / 'best_so_far.json', dict(
            scope='fixed selected calibration; no search', loss=common.SELECTED_LOSS,
            completed_endpoint=args.case, selected_summary_sha256=c['selected_summary_sha256']))
        verify_sources(c)
        files = sorted(p for p in folder.rglob('*') if p.is_file() and p.name not in ('receipt.json', 'heartbeat.json'))
        receipt = dict(status='complete', case=args.case, contract_sha256=args.contract_sha256,
                       selected_summary_sha256=c['selected_summary_sha256'], selected_checkpoint_sha256=c['checkpoint_sha256'],
                       scientific_bundle_sha256=c['code_bundle_sha256'],
                       baseline_smoke_receipt_sha256=c['baseline_smoke_receipt_sha256'],
                       elapsed_seconds=time.monotonic()-started, endpoint=result,
                       artifact_sha256={str(p.relative_to(folder)): adapter.digest(p) for p in files},
                       primary_comparison='rebated annual2percent versus rebated annual1percent, common inherited2023 population',
                       scope='impact only; no forward transition or population closure',
                       figures='suppressed_at_author_request', production_promoted=False)
        audit.save_json(folder / 'receipt.json', receipt)
        progress['phase'] = 'complete'
    except BaseException as error:
        progress['phase'] = 'failed'
        audit.save_json(folder / 'failure.json', dict(error=str(error), type=type(error).__name__, traceback=traceback.format_exc()))
        raise
    finally:
        stop.set()
        thread.join(timeout=2)
        audit.save_json(folder / 'heartbeat.json', dict(progress, epoch=time.time(), elapsed_seconds=time.monotonic()-started))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--case', choices=CASES, required=True)
    execute(parser.parse_args())
