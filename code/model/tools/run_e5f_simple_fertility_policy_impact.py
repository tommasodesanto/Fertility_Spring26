#!/usr/bin/env python3
"""Pinned one-date policy impacts from the selected fertility-nest population.

No historical refit, future population advance, saving-kernel substitution,
production promotion, or figures. Every process has an explicit 30-minute cap.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import math
import os
from pathlib import Path
import pickle
import sys
import threading
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code/model'), str(ROOT / 'code/model/tools')]
import numpy as np
import run_e5f_bounded_calibration_refinement as adapter
import run_e5f_joint_nested_finalize as finalize
import run_e5f_independent_numerical_audit as audit
import run_e5f_post2023_policy_mechanisms as policy
import run_e5f_transition_calibration as calibration
from run_e5f_global_saving_quantification import quantities
from intergen_eqscale_seq_optimized import solver

CASES = ('smoke', 'baseline', 'dependent-child-ltv95', 'supply-plus-20')
SELECTED_LOSS = 23.791955301663187
SCHEMA = 'e5f_simple_fertility_policy_impact_v1'


def verify_sources(contract):
    hashes = contract['source_hashes']
    if str(Path(__file__).resolve()) not in hashes:
        raise RuntimeError('The impact driver must be source-pinned')
    for path, digest in hashes.items():
        if not Path(path).is_absolute():
            raise RuntimeError('Source paths must be absolute')
        adapter.verify(path, digest)
    if calibration.code_fingerprint_contract(solver)['bundle_sha256'] != contract['code_bundle_sha256']:
        raise RuntimeError('Scientific source drift')


def load_contract(path, expected):
    adapter.verify(path, expected)
    contract = adapter.read_json(path)
    if contract['schema'] != SCHEMA or contract['case_seconds'] != 1800:
        raise RuntimeError('Unrecognized impact contract or time budget')
    if time.time() > float(contract['launch_deadline_epoch']):
        raise RuntimeError('Impact launch deadline exceeded')
    for name in ('selected_summary', 'output_root'):
        if not Path(contract[name]).is_absolute():
            raise RuntimeError(f'{name} must be absolute')
    return contract


def selected_packet(contract):
    verify_sources(contract)
    path = Path(contract['selected_summary'])
    adapter.verify(path, contract['selected_summary_sha256'])
    receipt = adapter.read_json(path.parent / 'case_receipt.json')
    if receipt['status'] != 'complete':
        raise RuntimeError('Selected calibration receipt is incomplete')
    plan = adapter.load_plan(path.parent.parent / 'plan.json', receipt['plan_sha256'])
    if plan.get('choice_model') != 'fertility_nest':
        raise RuntimeError('Selected plan has the wrong choice specification')
    case = next(x for x in plan['cases'] if x['id'] == receipt['case_id'])
    summary, _, _ = adapter.validate_result(path.parent, plan, case)
    if not math.isclose(summary['best_candidate']['transition_loss'], SELECTED_LOSS,
                        rel_tol=0, abs_tol=1e-12):
        raise RuntimeError('Selected calibration differs from the frozen overnight choice')
    if summary['code_fingerprints']['bundle_sha256'] != contract['code_bundle_sha256']:
        raise RuntimeError('Selected scientific bundle differs')
    if receipt['artifact_sha256'].get('dated_state.pkl.gz') != contract['checkpoint_sha256']:
        raise RuntimeError('Selected checkpoint differs from contracted checkpoint')
    for relative, digest in receipt['artifact_sha256'].items():
        adapter.verify(path.parent / relative, digest)
    finalize.configure_policy_model()
    packet, prepared, state, _ = finalize.prepare(path, summary)
    assert_model(packet['parameters'], packet['supply_rule'])
    state.initial_policy = None
    return packet, prepared, state


def assert_model(parameters, supply_rule):
    expected = dict(joint_nested_choice=True, fertility_nest_choice=True,
                    two_shock_choice=False, exhaustive_saving_control=False)
    for key, value in expected.items():
        if bool(getattr(parameters, key, False)) != value:
            raise RuntimeError(f'Unexpected model flag: {key}')
    if not math.isclose(float(parameters.tenure_choice_kappa), .005, rel_tol=0, abs_tol=1e-15):
        raise RuntimeError('Housing taste scale changed')
    if not math.isclose(float(supply_rule.elasticity), .63, rel_tol=0, abs_tol=1e-15):
        raise RuntimeError('Housing supply elasticity changed')
    for key in ('kappa_fert', 'kappa_fert_continuation'):
        if not math.isfinite(float(getattr(parameters, key))) or float(getattr(parameters, key)) < .005:
            raise RuntimeError('Invalid fertility-nest outer scale')


def baseline_comparison(evaluation, original, *, policies):
    gaps = {}
    for key in ('g_pre', 'g_post_fertility', 'g_current'):
        gaps[key + '_l1'] = float(np.abs(getattr(evaluation, key) - getattr(original, key)).sum())
    gaps['births_absolute'] = abs(float(evaluation.births) - float(original.births))
    gaps['price_max_absolute'] = float(np.max(np.abs(evaluation.policy.price - original.policy.price)))
    if any(not math.isfinite(x) or x > 2e-10 for x in gaps.values()):
        raise RuntimeError(f'Baseline does not reproduce the selected state: {gaps}')
    if policies:
        for key in ('V', 'c_pol', 'hR_pol', 'bp_pol', 'tenure_choice', 'tenure_probs',
                    'loc_probs', 'fert_probs', 'fert2_probs'):
            left, right = getattr(evaluation.policy, key), getattr(original.policy, key)
            if not np.allclose(left, right, rtol=0, atol=2e-10, equal_nan=True):
                raise RuntimeError(f'Fresh fixed-price baseline policy differs: {key}')
        for key in ('probabilities', 'wait_probabilities', 'failure_probabilities', 'products'):
            left = getattr(evaluation.policy.joint_choice, key)
            right = getattr(original.policy.joint_choice, key)
            if not np.allclose(left, right, rtol=0, atol=2e-10, equal_nan=True):
                raise RuntimeError(f'Fresh fixed-price joint policy differs: {key}')
        gaps['all_policy_arrays_within_2e_10'] = True
    return gaps


def probability_checks(evaluation):
    joint = evaluation.policy.joint_choice
    if joint is None:
        raise RuntimeError('Solved policy lacks its joint-choice object')
    result = {}
    for key in ('probabilities', 'wait_probabilities', 'failure_probabilities'):
        value = getattr(joint, key)
        if not np.isfinite(value).all() or np.any(value < 0) or np.any(value > 1):
            raise RuntimeError(f'Invalid joint probabilities: {key}')
        result[key] = dict(minimum=float(value.min()), maximum=float(value.max()))
    occupied = evaluation.g_pre > 0
    sums = joint.probabilities.sum(axis=(-2, -1))
    gap = float(np.max(np.abs(sums[occupied] - 1))) if np.any(occupied) else 0.
    if gap > 2e-11:
        raise RuntimeError(f'Occupied joint menu normalization failed: {gap}')
    result['occupied_menu_sum_max_gap'] = gap
    result['scope'] = ('All joint probability ranges and occupied wait/attempt marginal normalization; '
                       'natural conception/product factorization also checks occupied mass internally.')
    return result


def evaluate_case(packet, prepared, state, name, folder, *, clear):
    folder.mkdir(parents=True, exist_ok=False)
    parameters = copy.deepcopy(packet['parameters'])
    current_state = copy.deepcopy(state)
    current_state.initial_policy = None
    spec = policy.POLICIES[name]
    policy.apply_policy(parameters, spec)
    rule = policy.policy_supply_rule(prepared.supply_rule, spec)
    assert_model(parameters, rule)
    counter = policy.calendar.SolveCounter()
    started = time.monotonic()
    if clear:
        evaluation, shared, fallback = policy.baseline.evaluate_state(
            current_state, parameters, prepared.b_grid, rule, counter, 2e-4, 60)
    else:
        shared = solver.precompute_shared(parameters, prepared.b_grid)
        evaluation = policy.calendar.evaluate_period(
            packet['evaluation'].policy.price.copy(), current_state.g_pre,
            parameters, prepared.b_grid, shared, counter, rule)
        fallback = False
    current = dict(parameters=parameters, b_grid=prepared.b_grid, evaluation=evaluation,
                   shared=shared, supply_rule=rule, state=current_state, calendar_year=2023,
                   policy_case=name, scope='one_date_impact_common_inherited_population')
    checkpoint = folder / 'dated_state.pkl.gz'
    with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
        pickle.dump(current, stream, protocol=5)
    health = policy.calendar.distribution_health(dict(
        pre=evaluation.g_pre, post=evaluation.g_post_fertility, current=evaluation.g_current))
    if health['nonfinite_distribution_count'] or health['min_distribution_mass'] < -1e-14:
        raise RuntimeError('Distribution health gate failed')
    inherited_mass = float(state.g_pre.sum())
    mass_gaps = {key: abs(float(getattr(evaluation, key).sum()) - inherited_mass)
                 for key in ('g_pre', 'g_post_fertility', 'g_current')}
    if any(not math.isfinite(x) or x > 2e-10 for x in mass_gaps.values()):
        raise RuntimeError(f'Common impact population changed: {mass_gaps}')
    if not math.isfinite(evaluation.feasibility_projection_mass) or evaluation.feasibility_projection_mass > 1e-6:
        raise RuntimeError('Feasibility projection gate failed')
    if clear and (not math.isfinite(evaluation.relative_market_residual) or evaluation.relative_market_residual > 2e-4):
        raise RuntimeError('Housing market gate failed')
    joint = probability_checks(evaluation)
    budget = audit.budget_audit(current, folder)
    arrays = audit.policy_array_audit(current, folder)
    if not math.isfinite(budget['budget_excess_mass']) or budget['budget_excess_mass'] > 2e-10 or arrays['occupied_negative_steps'] != 0:
        raise RuntimeError('Occupied budget/value gate failed')
    for bounds in arrays['probabilities'].values():
        if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
            raise RuntimeError('Marginal policy probability range failed')
    comparison = baseline_comparison(evaluation, packet['evaluation'], policies=not clear) if name == 'baseline' else None
    values = quantities(evaluation, parameters)
    if any(not math.isfinite(float(x)) for x in values.values()):
        raise RuntimeError('Nonfinite impact quantity')
    audit.save_json(folder / 'quantities.json', values)
    result = dict(status='complete', policy=name, market_cleared=clear, quantities=values,
                  baseline_comparison=comparison, budget=budget, policy_arrays=arrays,
                  joint_probabilities=joint, distribution_health=health, mass_gaps=mass_gaps,
                  checkpoint_sha256=adapter.digest(checkpoint), grid_fallback=bool(fallback),
                  bellman_solves=counter.bellman, elapsed_seconds=time.monotonic()-started)
    audit.save_json(folder / 'evaluation_receipt.json', result)
    return result


def execute(args):
    contract = load_contract(args.contract, args.contract_sha256)
    root = Path(contract['output_root'])
    folder = root / args.case
    folder.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    state = dict(phase='verifying_inputs', case=args.case, completed_evaluations=0)
    stop = threading.Event()
    def heartbeat():
        while not stop.is_set():
            elapsed = time.monotonic() - start
            audit.save_json(folder / 'heartbeat.json', dict(state, elapsed_seconds=elapsed, epoch=time.time()))
            if elapsed >= contract['case_seconds']:
                audit.save_json(folder / 'failure.json', dict(error='Declared 1800-second case cap exceeded', **state))
                os._exit(124)
            stop.wait(30)
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        audit.save_json(folder / 'latest_completed_case.json', dict(status='no_evaluation_completed'))
        audit.save_json(folder / 'best_so_far.json', dict(scope='fixed selected calibration; no search', loss=SELECTED_LOSS))
        packet, prepared, inherited = selected_packet(contract)
        if args.case != 'smoke':
            receipt = adapter.read_json(root / 'smoke' / 'receipt.json')
            if receipt['status'] != 'complete' or receipt['contract_sha256'] != args.contract_sha256:
                raise RuntimeError('Matching complete smoke receipt required')
            for relative, digest in receipt['artifact_sha256'].items():
                adapter.verify(root / 'smoke' / relative, digest)
        evaluations = []
        steps = [('baseline', False), ('baseline', True)] if args.case == 'smoke' else [(args.case, True)]
        for index, (name, clear) in enumerate(steps):
            state['phase'] = name + ('_market' if clear else '_fixed_price')
            result = evaluate_case(packet, prepared, inherited, name,
                                   folder / ('market' if clear else 'fixed_price'), clear=clear)
            evaluations.append(result)
            state['completed_evaluations'] = index + 1
            audit.save_json(folder / 'latest_completed_case.json', result)
            audit.save_json(folder / 'best_so_far.json', dict(
                scope='fixed selected calibration; no search', loss=SELECTED_LOSS,
                completed_evaluations=index+1, selected_summary_sha256=contract['selected_summary_sha256']))
        verify_sources(contract)
        files = sorted(p for p in folder.rglob('*') if p.is_file() and p.name not in ('heartbeat.json', 'receipt.json'))
        receipt = dict(status='complete', case=args.case, contract_sha256=args.contract_sha256,
                       selected_summary_sha256=contract['selected_summary_sha256'],
                       selected_checkpoint_sha256=contract['checkpoint_sha256'],
                       scientific_bundle_sha256=contract['code_bundle_sha256'],
                       elapsed_seconds=time.monotonic()-start, evaluations=evaluations,
                       artifact_sha256={str(p.relative_to(folder)): adapter.digest(p) for p in files},
                       scope='2023 impact only; same inherited population; no future transition or population closure',
                       figures='suppressed_at_author_request', production_promoted=False)
        audit.save_json(folder / 'receipt.json', receipt)
        state['phase'] = 'complete'
    except BaseException as error:
        state['phase'] = 'failed'
        audit.save_json(folder / 'failure.json', dict(error=str(error), type=type(error).__name__, traceback=traceback.format_exc()))
        raise
    finally:
        stop.set()
        thread.join(timeout=2)
        audit.save_json(folder / 'heartbeat.json', dict(state, elapsed_seconds=time.monotonic()-start, epoch=time.time()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--case', choices=CASES, required=True)
    execute(parser.parse_args())
