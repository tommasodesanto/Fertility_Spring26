#!/usr/bin/env python3
"""Isolated, one-date joint-choice diagnostic; never a calibration launcher."""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import sys
import threading
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code/model'), str(ROOT / 'code/model/tools')]
from e5f_joint_nested_choice import choose, plan_values, scatter_joint_block
import run_e5f_independent_numerical_audit as audit

model, calendar, transition = audit.model, audit.calendar, audit.transition


def write(path, payload):
    audit.save_json(Path(path), payload)


def table(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError('Empty table')
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def array_sha(a):
    a = np.ascontiguousarray(a)
    return hashlib.sha256(str((a.shape, a.dtype.str)).encode() + a.tobytes()).hexdigest()


def install_operators():
    calendar.apply_fertility = transition.apply_sequential_fertility
    calendar.advance_calendar_distribution = transition.advance_sequential_calendar_distribution
    calendar.distribution_rows = transition.independent_child_distribution_rows


def verify_contract(args):
    if audit.digest(args.contract) != args.contract_sha256:
        raise RuntimeError('Experiment contract hash changed')
    contract = json.loads(args.contract.read_text())
    for rel, expected in contract['driver_hashes'].items():
        if audit.digest(ROOT / rel) != expected:
            raise RuntimeError(f'Experimental source changed: {rel}')
    if audit.digest(args.checkpoint) != contract['checkpoint_sha256']:
        raise RuntimeError('Input checkpoint hash differs')
    if audit.baseline.calibration.code_fingerprint_contract(model)['bundle_sha256'] != contract['scientific_bundle_sha256']:
        raise RuntimeError('Scientific bundle differs')
    return contract


def capture_reference(packet, out):
    """Observe the original Bellman, with its own V as frozen continuation."""
    P = copy.deepcopy(packet['parameters'])
    reference = packet['evaluation']
    bg = packet['b_grid']
    if P.I != 1 or not P.sequential_births or not model.independent_child_maturation_active(P):
        raise ValueError('Experiment requires one-market independent-child sequential model')
    if model.readiness_gate_active(P):
        raise ValueError('Readiness extension is outside this experiment')
    if P.tenure_choice_kappa <= 0 or not bool(getattr(P, 'use_tenure_kernel', True)):
        raise ValueError('Reference must use the original compiled tenure logit')
    shape = reference.g_pre.shape
    q = np.empty(shape + (2,)); products = np.empty(shape + (2,), dtype=np.int16)
    original_tenure = model.tenure_logit_kernel
    original_bellman = model.solve_bellman_full_markov_income
    calls = 0; max_reconstruction_gap = 0.0

    def observed_tenure(Vd, *rest):
        nonlocal calls, max_reconstruction_gap
        j = P.J - 1 - calls // shape[4]; z = calls % shape[4]
        if j < 0:
            raise RuntimeError('Unexpected extra tenure call')
        result = original_tenure(Vd, *rest)
        # The deterministic kernel shares every original wealth/feasibility map.
        deterministic = model.tenure_choice_kernel(Vd, *rest[:-1])
        temporary = Vd.copy(); temporary[:, 1:] = -1e10
        rental, rent_product = model.tenure_choice_kernel(temporary, *rest[:-1])
        temporary = Vd.copy(); temporary[:, 0] = -1e10
        owner, own_product = model.tenure_choice_kernel(temporary, *rest[:-1])
        gap = float(np.max(np.abs(np.maximum(rental, owner) - deterministic[0])))
        max_reconstruction_gap = max(max_reconstruction_gap, gap)
        if gap != 0:
            raise RuntimeError('Tenure-specific maxima fail original deterministic maximum')
        shift = float(P.E_loc[0] - P.mu_stay)
        q[:, :, :, j, z, :, :, 0] = rental + shift
        q[:, :, :, j, z, :, :, 1] = owner + shift
        products[:, :, :, j, z, :, :, 0] = rent_product
        products[:, :, :, j, z, :, :, 1] = own_product
        calls += 1
        if z == shape[4] - 1:
            audit.progress(out, 'capture_age_complete', age_index=int(j), blocks=calls)
        return result

    def frozen_bellman(*positional, **keywords):
        if 'continuation_V' in keywords:
            raise RuntimeError('Unexpected preexisting continuation override')
        return original_bellman(*positional, continuation_V=reference.policy.V, **keywords)

    model.tenure_logit_kernel = observed_tenure
    model.solve_bellman_full_markov_income = frozen_bellman
    try:
        shared = model.precompute_shared(P, bg)
        replay = calendar.evaluate_period(reference.policy.price.copy(),
            packet['state'].g_pre.copy(), P, bg, shared, calendar.SolveCounter(), packet['supply_rule'])
    finally:
        model.tenure_logit_kernel = original_tenure
        model.solve_bellman_full_markov_income = original_bellman
    if calls != P.J * shape[4]:
        raise RuntimeError('Missing age/income capture blocks')
    checks = {}
    for field in ['V', 'c_pol', 'hR_pol', 'bp_pol', 'tenure_choice', 'tenure_probs', 'loc_probs', 'fert_probs', 'fert_value']:
        a, b = getattr(reference.policy, field), getattr(replay.policy, field)
        if not np.array_equal(a, b):
            raise RuntimeError(f'Frozen-continuation reference does not reproduce {field}: {np.max(np.abs(a-b))}')
        checks[field] = array_sha(a)
    for field in ['g_pre', 'g_post_fertility', 'g_current']:
        a, b = getattr(reference, field), getattr(replay, field)
        if not np.array_equal(a, b):
            raise RuntimeError(f'Reference distribution does not reproduce {field}')
        checks[field] = array_sha(a)
    if reference.births != replay.births or not np.array_equal(reference.demand_by_loc, replay.demand_by_loc):
        raise RuntimeError('Reference aggregates differ')
    write(out / 'reference_replay.json', {'status': 'exact', 'arrays': checks,
        'blocks': calls, 'maximum_deterministic_reconstruction_gap': max_reconstruction_gap,
        'q_sha256': array_sha(q), 'product_sha256': array_sha(products)})
    # Write the existing seventeen-graph reference packet unchanged.
    audit.standard_diagnostics(packet, out / 'reference', validate_production_young=True)
    return q, products


def choice_diagnostics(packet, current, post, attempt, owner_choice, value, ages, row, out):
    """Supplemental one-date screens, never a replacement equilibrium packet."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    P, origin = packet['parameters'], packet['evaluation'].g_pre
    selectors = [('poorest_wealth_node', 0, 0), ('richest_wealth_node', 0, -1),
        ('youngest_age', 3, 0), ('oldest_age', 3, -1),
        ('inherited_renter', 1, 0), ('childless', 5, 0), ('top_parity', 5, -1),
        ('no_child_at_home', 6, 0), ('highest_child_count', 6, -1)]
    bounds = []
    for label, axis, index in selectors:
        selection = [slice(None)] * origin.ndim; selection[axis] = index
        selection = tuple(selection); mass = origin[selection]
        den = float(mass.sum())
        bounds.append(dict(boundary=label, households=den,
            population_share=den / origin.sum(),
            attempt_rate=float(np.sum(mass * attempt[selection]) / max(den, 1e-30)),
            owner_choice_rate=float(np.sum(mass * owner_choice[selection]) / max(den, 1e-30))))
    mass = origin[:, 1:]; den = float(mass.sum())
    bounds.append(dict(boundary='inherited_owner', households=den, population_share=den/origin.sum(),
        attempt_rate=float(np.sum(mass*attempt[:, 1:])/max(den,1e-30)),
        owner_choice_rate=float(np.sum(mass*owner_choice[:, 1:])/max(den,1e-30))))
    table(out / 'boundary_exposure.csv', bounds)
    finite = np.isfinite(value[:-1]) & np.isfinite(value[1:])
    differences = np.zeros_like(value[:-1])
    np.subtract(value[1:], value[:-1], out=differences, where=finite)
    exposed = (origin[:-1] > 1e-12) & finite
    drops = exposed & (differences < -1e-7)
    row['occupied_value_decreases'] = int(drops.sum())
    row['value_decrease_mass'] = float(origin[:-1][drops].sum())
    row['minimum_exposed_value_increment'] = float(differences[exposed].min()) if exposed.any() else None
    # The original within-product policies remain optimal against frozen future
    # values. Audit their budgets under the newly selected destination mass.
    modified = dict(packet); evaluation = copy.copy(packet['evaluation'])
    evaluation.g_current, evaluation.g_post_fertility = current, post
    modified['evaluation'] = evaluation
    budget = audit.budget_audit(modified, out)
    row['budget_excess_mass'] = budget['budget_excess_mass']
    row['maximum_occupied_budget_excess'] = budget['maximum_occupied_excess']
    # Preserve the benchmark's documented reporting-floor exception, but do
    # not allow offsets: audit positive mass changes at each destination state.
    added = np.maximum(current - packet['evaluation'].g_current, 0)
    if added.sum() > 0:
        evaluation.g_current = added
        incremental_budget = audit.budget_audit(modified, out / 'added_budget_exposure')
        added_bad_mass = incremental_budget['budget_excess_mass']
    else:
        added_bad_mass = 0.0
    row['added_budget_excess_mass'] = float(added_bad_mass)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    age = [r['age'] for r in ages]
    axes[0, 0].plot(age, [r['attempt_rate'] for r in ages], label='Attempt')
    axes[0, 0].plot(age, [r['ownership'] for r in ages], label='Ownership')
    axes[0, 0].set(xlabel='Age node', ylabel='Fraction', ylim=(-.02, 1.02)); axes[0, 0].legend()
    axes[0, 1].plot(age, [r['explicit_birth_rate'] for r in ages])
    axes[0, 1].set(xlabel='Age node', ylabel='Births per household in age cell')
    wealth_mass = origin.sum(axis=tuple(range(1, origin.ndim)))
    for values, label in [(attempt, 'Attempt'), (owner_choice, 'Ownership')]:
        rate = np.divide((origin*values).sum(axis=tuple(range(1, origin.ndim))), wealth_mass,
            out=np.full_like(wealth_mass, np.nan), where=wealth_mass > 1e-12)
        axes[1, 0].plot(packet['b_grid'], rate, label=label)
    axes[1, 0].set(xlabel='Inherited liquid wealth', ylabel='Population-weighted choice fraction', ylim=(-.02, 1.02))
    axes[1, 0].legend()
    labels = ['Rent/wait', 'Rent/attempt', 'Own/wait', 'Own/attempt']
    shares = [row[f'joint_share_t{t}_a{a}'] for t in range(2) for a in range(2)]
    axes[1, 1].bar(labels, shares); axes[1, 1].tick_params(axis='x', labelrotation=20)
    axes[1, 1].set(ylabel='Joint plan share', ylim=(0,1))
    fig.suptitle(f"Supplemental one-date choices: {row['rule']}; kappa={row['outer_scale']:g}, lambda={row['dissimilarity']:g}\nFixed prices, inherited population and baseline future values")
    fig.savefig(out / 'supplemental_choices.png', dpi=160); plt.close(fig)
    if added_bad_mass > 2e-10:
        write(out/'failed_budget_case.json', row)
        raise RuntimeError('Newly occupied conditional policies violate budget gate')


def evaluate_choice(packet, q, products, outer, lam, rule, out):
    started = time.perf_counter()
    P, reference = packet['parameters'], packet['evaluation']
    origin = reference.g_pre
    current = np.zeros_like(origin); post = np.zeros_like(origin)
    attempt = np.zeros_like(origin); owner_choice = np.zeros_like(origin)
    value = np.full_like(origin, -np.inf)
    births = 0.0; dead_mass = 0.0; max_prob_error = 0.0
    joint_mass = np.zeros((2, 2)); conception = model.get_fecundity_by_age(P)
    maps = reference.policy.maps
    for j in range(P.J):
        fertile = P.A_f_start <= j + 1 <= P.A_f_end
        for z in range(origin.shape[4]):
            for nn in range(P.n_parity):
                for cs in range(nn + 1):
                    source = origin[:, :, 0, j, z, nn, cs]
                    if source.sum() == 0:
                        continue
                    available = fertile and nn < P.n_parity - 1
                    dest = (nn + 1, cs + 1) if available else None
                    q0 = q[:, :, 0, j, z, nn, cs]
                    c0 = products[:, :, 0, j, z, nn, cs]
                    q1 = q[:, :, 0, j, z, dest[0], dest[1]] if available else None
                    c1 = products[:, :, 0, j, z, dest[0], dest[1]] if available else None
                    pi = float(conception[j]) if available else 0.0
                    Q = plan_values(q0, q1, pi, P.first_birth_fixed_cost if nn == 0 else 0, available)
                    v, prob = choose(Q, outer, lam, rule)
                    positive = source > 0
                    sums = prob.sum(axis=(-2, -1))
                    dead_mass += float(source[sums == 0].sum())
                    max_prob_error = max(max_prob_error, float(np.max(np.abs(sums[positive] - 1))))
                    if np.any(prob < 0) or not np.isfinite(prob).all():
                        raise RuntimeError('Invalid choice probability')
                    joint_mass += np.sum(source[..., None, None] * prob, axis=(0, 1))
                    attempt[:, :, 0, j, z, nn, cs] = prob[..., 1].sum(axis=-1)
                    owner_choice[:, :, 0, j, z, nn, cs] = prob[..., 1, :].sum(axis=-1)
                    value[:, :, 0, j, z, nn, cs] = v
                    gc, gp, born = scatter_joint_block(source, prob, c0, c1, pi,
                        dest, (nn, cs), maps.tmx_idx[0], maps.tmx_wt[0])
                    current[:, :, 0, j, z] += gc
                    post[:, :, 0, j, z] += gp
                    births += born
    total = float(origin.sum())
    error = max(abs(current.sum() - total), abs(post.sum() - total))
    if dead_mass > 1e-12 or max_prob_error > 1e-12 or error > 2e-10:
        raise RuntimeError(f'Joint forward gate: dead={dead_mass}, probabilities={max_prob_error}, mass={error}')
    # Exact stock-flow accounting: total explicit parity increment equals births.
    parity = np.arange(P.n_parity).reshape((1, 1, 1, 1, 1, P.n_parity, 1))
    parity_increment = float(np.sum((post - origin) * parity))
    if abs(parity_increment - births) > 2e-10:
        raise RuntimeError('Birth flow differs from parity-stock increment')
    demand = calendar.housing_demand_by_location(current, reference.policy.hR_pol, P)
    stock = float(reference.supply_by_loc.sum())
    adjustment = transition.calendar_topcode_birth_accounting(origin, post, births, P)
    row = dict(rule=rule, outer_scale=outer, dissimilarity=lam, inner_scale=outer * lam,
        households=total, explicit_births_per_household=births / total,
        births_per_household=float(adjustment['topcode_adjusted_birth_children']) / total,
        ownership=float(current[:, 1:].sum() / total), rooms=float(demand.sum() / total),
        signed_fixed_price_market_gap=float((demand.sum() - stock) / stock),
        mass_error=float(error), max_probability_error=max_prob_error, infeasible_mass=dead_mass,
        birth_parity_identity_error=abs(parity_increment - births),
        elapsed_seconds=time.perf_counter() - started)
    for tenure in range(2):
        for action in range(2):
            row[f'joint_share_t{tenure}_a{action}'] = float(joint_mass[tenure, action] / total)
    ages = []
    for j in range(P.J):
        pre = origin[:, :, :, j]; mass = float(pre.sum())
        pi = float(conception[j]) if P.A_f_start <= j + 1 <= P.A_f_end else 0
        ages.append(dict(age_index=j, age=float(P.age_start + P.da*j),
            households=mass, attempt_rate=float(np.sum(pre * attempt[:, :, :, j]) / max(mass, 1e-30)),
            explicit_birth_rate=float(pi * np.sum(pre * attempt[:, :, :, j]) / max(mass, 1e-30)),
            ownership=float(current[:, 1:, :, j].sum() / max(mass, 1e-30))))
    out.mkdir(parents=True)
    table(out / 'by_age.csv', ages)
    choice_diagnostics(packet, current, post, attempt, owner_choice, value, ages, row, out)
    arrays_path = out / 'choice_state.npz'
    np.savez_compressed(arrays_path, g_current=current, g_post=post,
                        attempt=attempt, owner_choice=owner_choice, value=value)
    row['choice_state_sha256'] = audit.digest(arrays_path)
    row['elapsed_seconds'] = time.perf_counter() - started
    write(out / 'summary.json', row)
    return row


def run(args):
    contract = verify_contract(args)
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise FileExistsError(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    transition.configure_sequential_model(); install_operators()
    packet = audit.load_checkpoint(args.checkpoint)
    if packet['input_hashes']['code_bundle_sha256'] != contract['scientific_bundle_sha256']:
        raise RuntimeError('Saved source differs')
    if packet['input_hashes']['target_fingerprint'] != contract['reference_target_fingerprint']:
        raise RuntimeError('Reference target contract differs')
    P, e = packet['parameters'], packet['evaluation']
    total = float(e.g_pre.sum())
    birth_accounting = transition.calendar_topcode_birth_accounting(e.g_pre, e.g_post_fertility, e.births, P)
    write(args.outdir / 'reference_aggregates.json', dict(households=total,
        explicit_births_per_household=float(e.births)/total,
        births_per_household=float(birth_accounting['topcode_adjusted_birth_children'])/total,
        ownership=float(e.g_current[:, 1:].sum())/total,
        rooms=float(e.demand_by_loc.sum())/total,
        relative_market_residual=e.relative_market_residual,
        age_start=P.age_start, age_step=P.da))
    audit.budget_audit(packet, args.outdir / 'reference_budget')
    started = time.perf_counter(); stop = threading.Event()
    def heartbeat():
        while not stop.wait(60):
            write(args.outdir / 'runtime_heartbeat.json', {'elapsed_seconds': time.perf_counter()-started,
                'utc': time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())})
    threading.Thread(target=heartbeat, daemon=True).start()
    try:
        write(args.outdir / 'best_so_far.json', {'status': 'diagnostic_only_no_calibration_selection'})
        write(args.outdir / 'latest_completed_case.json', {'status': 'no_case_complete'})
        if args.stage == 'smoke':
            captures = []
            for repeat in range(2):
                location = args.outdir / f'capture_{repeat}'; location.mkdir()
                q, products = capture_reference(packet, location)
                captures.append((array_sha(q), array_sha(products)))
            if captures[0] != captures[1]:
                raise RuntimeError('Independent choice captures differ')
            np.savez_compressed(args.outdir / 'captured_values.npz', q=q, products=products)
            cases = [(.5, .5), (.5, 1.0)]
        else:
            if not args.smoke or audit.digest(args.smoke / 'summary.json') != args.smoke_sha256:
                raise RuntimeError('Exact smoke receipt missing or changed')
            smoke = json.loads((args.smoke / 'summary.json').read_text())
            if smoke['status'] != 'complete' or smoke['contract_sha256'] != args.contract_sha256:
                raise RuntimeError('Smoke not complete for this contract')
            if audit.digest(args.smoke / 'captured_values.npz') != smoke['capture_sha256']:
                raise RuntimeError('Captured values changed')
            with np.load(args.smoke / 'captured_values.npz') as data:
                q, products = data['q'], data['products']
            cases = [(float(a), float(b)) for a,b in contract['cases']]
        rows = []
        for index, (outer, lam) in enumerate(cases):
            for rule in ['joint', 'sequential']:
                if time.perf_counter() - started > contract[f'{args.stage}_seconds']:
                    raise TimeoutError('Stage budget exhausted')
                row = evaluate_choice(packet, q, products, outer, lam, rule,
                                      args.outdir / f'case_{index:02d}_{rule}')
                rows.append(row); table(args.outdir / 'cases.csv', rows)
                write(args.outdir / 'latest_completed_case.json', row)
                audit.progress(args.outdir, 'case_complete', case=index, rule=rule, **{
                    k: row[k] for k in ['ownership', 'births_per_household', 'signed_fixed_price_market_gap']})
            if lam == 1:
                a,b=rows[-2:]
                for key in ['ownership', 'births_per_household', 'rooms']:
                    if abs(a[key]-b[key]) > 1e-12:
                        raise RuntimeError('Flat joint logit equality fails')
        verify_contract(args)
        write(args.outdir / 'summary.json', {'status':'complete', 'stage':args.stage,
            'contract_sha256':args.contract_sha256, 'completed_cases':len(rows),
            'elapsed_seconds':time.perf_counter()-started,
            'capture_sha256':audit.digest(args.outdir/'captured_values.npz') if args.stage=='smoke' else smoke['capture_sha256'],
            'interpretation':'fixed price, inherited population and continuation; not equilibrium or calibration',
            'production_changed':False})
    finally:
        stop.set()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract',type=Path,required=True); parser.add_argument('--contract-sha256',required=True)
    parser.add_argument('--checkpoint',type=Path,required=True); parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--stage',choices=['smoke','panel'],required=True)
    parser.add_argument('--smoke',type=Path); parser.add_argument('--smoke-sha256')
    args=parser.parse_args()
    try:
        run(args)
    except Exception as exc:
        args.outdir.mkdir(parents=True,exist_ok=True)
        write(args.outdir/'failure.json',{'type':type(exc).__name__,'message':str(exc)})
        raise
