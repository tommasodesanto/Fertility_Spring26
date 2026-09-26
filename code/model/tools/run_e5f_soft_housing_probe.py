#!/usr/bin/env python3
"""Torch-only, one softened-housing equilibrium against an authenticated control.

All preferences, prices' supply rule, earnings, entry and measurement contracts
are inherited. Only parent housing services change. Fertility utility is not
renormalized: normalized entry is maintained and reproduction gaps are reported.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import threading
import time

PIN = 'c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a'
BASE = Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2')


def write(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def table(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader(); writer.writerows(rows)


def incentives(packet, out):
    """Same pre-birth weighting and finite-logit gap as the saved-point audit."""
    import numpy as np
    from diagnose_saved_utility_birth_incentives import summarize
    from intergen_eqscale_seq_optimized.parameters import get_fecundity_by_age
    P, sol = packet['parameters'], packet['solution']
    g, policy = packet['stationary_g_pre'], packet['evaluation'].policy
    fec = get_fecundity_by_age(P)
    groups = {'first_birth': [], 'later_births': [], 'all_births': [], 'first_birth_nonpositive_liquid_wealth': []}
    age_rows = []
    post = g.copy()
    for j in range(int(P.J)):
        if not P.A_f_start <= j+1 <= P.A_f_end: continue
        for n in range(int(P.n_parity)-1):
            records = []
            for m in range(n+1):
                w = g[:, :, :, j, :, n, m]
                pair = policy.fert_probs[:, :, :, j, :, :2] if n == 0 else policy.fert2_probs[:, :, :, j, :, :, n-1, m]
                p0, p1 = pair[..., 0], pair[..., 1]
                occupied = w > 0
                unavailable = (p0 == 0) & (p1 == 0)
                assert float(w[unavailable].sum()) < 1e-12
                np.testing.assert_allclose((p0+p1)[occupied & ~unavailable], 1., atol=1e-12, rtol=0)
                born = w*p1*fec[j]
                post[:, :, :, j, :, n, m] -= born
                post[:, :, :, j, :, n+1, m+1] += born
                kappa = P.kappa_fert if n == 0 else P.kappa_fert_continuation
                record = (w[occupied], p0[occupied], p1[occupied], np.full(occupied.sum(), fec[j]), np.full(occupied.sum(), kappa))
                records.append(record)
                groups['first_birth' if n == 0 else 'later_births'].append(record)
                groups['all_births'].append(record)
                if n == 0:
                    low = occupied & (packet['b_grid'][:, None, None, None] <= 0.)
                    groups['first_birth_nonpositive_liquid_wealth'].append((w[low], p0[low], p1[low], np.full(low.sum(),fec[j]), np.full(low.sum(),kappa)))
            row = summarize(f'child_{n+1}_age_{P.age_start+j*P.da}', records)
            row.update(age=float(P.age_start+j*P.da), birth_number=n+1)
            age_rows.append(row)
    assert float(np.abs(post-sol.g_beginning_distribution).sum()) < 1e-10
    rows = [summarize(k, v) for k, v in groups.items()]
    table(out/'incentives.csv', rows); table(out/'incentives_by_age.csv', age_rows)
    return rows


def report(packet, runtime, tax, objective, reference, out, case):
    import numpy as np
    P, sol, evaluation = (packet[k] for k in ('parameters', 'solution', 'evaluation'))
    grid, shared = packet['b_grid'], packet['shared']
    prices = np.asarray(sol.p_eq)
    primitive, model = runtime['primitive'], runtime['model']
    gates = primitive.pf.transition.operator_gates(sol, evaluation.policy, packet['stationary_g_pre'], P, grid, shared)
    reconstructed, reconstruction = primitive.pf.calendar.reconstruct_stationary_pre_fertility(sol,evaluation.policy,P,grid,shared)
    np.testing.assert_allclose(reconstructed,packet['stationary_g_pre'],atol=1e-12,rtol=0)
    gates.update(reconstruction)
    for name in ('one_step_constant_path_nesting_l1', 'mature_flow_abs_error', 'birth_flow_abs_error', 'topcode_adjusted_birth_flow_abs_error'):
        assert abs(float(gates[name])) < 5e-9, (name, gates[name])
    assert abs(gates['zero_entry_mass_accounting_residual']) < 2e-8
    assert gates['stationary_feasibility_projection_mass'] < 1e-6
    assert evaluation.relative_market_residual < 2e-4
    budget = primitive.dated_budget(evaluation, P, shared, grid, float(P.user_cost_rate*prices[0]))
    assert budget['budget_excess_mass'] <= 2e-10
    purchase = runtime['accounting'].audit_purchase_accounting(evaluation, P, shared, grid, model)
    fiscal = runtime['certify_initial_pension'](evaluation.g_current, P, marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    arrays = runtime['audit'].policy_array_audit(packet, out)
    assert not arrays['occupied_negative_steps']
    assert not any(x['nonfinite'] or x['minimum'] < 0 or x['maximum'] > 1 for x in arrays['probabilities'].values())
    early = dict(fertility={p: runtime['observe_initial_fertility'](evaluation, P, age_projection=p) for p in ('uniform_birth_time', 'constant_post_cell')},
        housing_wealth=runtime['observe_initial_housing_wealth'](evaluation, P, grid, shared, diagnostic_enabled=True, age_projection='uniform_within_age_cell', diagnostic_allow_family_proxies=True, include_wealth=True, include_birth_response=True))
    recent = runtime['observe_recent_parent_flow'](evaluation, P, diagnostic_enabled=True, snapshot=runtime['SNAPSHOT'], age_projection=runtime['AGE_PROJECTION'], diagnostic_allow_residence_proxy=True, input_provenance={'case_id': case})
    rows = tax.target_rows(objective, early, recent['model_value'], float(runtime['chain'].extract_moments(sol,P)['tfr']))
    assert len(rows) == 13
    if case == 'control':
        inherited = {r['moment']: r for r in csv.DictReader((reference/'target_fit.csv').open())}
        for row in rows:
            for key in ('target','model','gap','weight','loss_contribution'):
                if row[key] != '':
                    assert float(row[key]) == float(inherited[row['moment']][key]), (row['moment'],key)
    table(out/'target_fit.csv', rows)
    parameters = list(csv.DictReader((reference/'parameters.csv').open()))
    actual = tax.actual_parameters(P)
    for row in parameters:
        key = row['parameter']
        if key in actual: assert float(row['estimate']) == float(actual[key]), key
        if key == 'pension_period': row['estimate'] = float(P.pension)
        if key == 'psi_child': row['status'] = 'held at reference; no fertility renormalization'
        elif row['status'] == 'experimental free coordinate': row['status'] = 'held at reference estimate; no recalibration'
    parameters.append(dict(parameter='softness_fraction', estimate=0. if case=='control' else .1, lower='', upper='', status='fixed experimental restriction; not estimated', near_bound=''))
    table(out/'parameters.csv', parameters)
    incentive_rows = incentives(packet, out)
    # Actual post-tenure renter occupancy, not conditional renter policies.
    mass = np.asarray(evaluation.g_current)[:,0,:,:,:,:,1:]
    housing = np.asarray(evaluation.policy.hR_pol)[:,0,:,:,:,:,1:]
    hb = float(np.max(shared.hb_flat))
    parent_rental = dict(mass=float(mass.sum()),
        fraction_below_old_floor=float(mass[housing < hb].sum()/mass.sum()),
        mean_housing=float((mass*housing).sum()/mass.sum()),old_floor=hb)
    runtime['audit'].standard_diagnostics(packet, out, validate_production_young=False)
    assert len(list((out/'standard_diagnostics').glob('*.png'))) == 17
    result = dict(case=case, loss=sum(float(r['loss_contribution']) for r in rows if r['loss_contribution']!=''), prices=prices.tolist(),
        target_fit=rows, incentives=incentive_rows, parent_rental=parent_rental, housing_observer=early['housing_wealth'], recent_parent_observer=recent,
        market_residual=float(evaluation.relative_market_residual), fiscal=fiscal, budget=budget, purchase=purchase, operator_gates=gates,
        fixed_psi=float(P.psi_child), birth_replacement_relative_gap=float(sol.adult_entry_stationary_relative_gap),
        closure='Normalized entry/age composition, housing and pension clearing; changed reproduction gap reported, not closed by refitting psi',
        caveats=['Frozen old ACS housing target retained for comparison, not adopted AHS target.', 'First-birth housing observer is a stationary proxy, not the matched empirical event study.', 'Gap inversion excludes current taste draw but retains future shock-inclusive values.'])
    write(out/'receipt.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['inspect', 'run'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Torch allocation required')
    import numpy as np
    import collect_e5f_utility_comparison as collector
    args.output.mkdir(parents=True, exist_ok=False)
    os.environ[collector.ENV_PIN] = PIN
    contract, runner = collector._load_contract(BASE/'launch_v1/contract.json')
    _, _, _, tax, _, objective, runtime, _ = runner.setup(contract, 'floor_linear', args.output/'runtime')
    selection = collector.read_json(BASE/'results/run_001/floor_linear/selected.json')
    selected = collector.scientific_checkpoint(selection['original_case_output'], contract, 'floor_linear', runtime, tax)
    packet = selected['packet']; P = packet['parameters']; model = runtime['model']
    flags = {k: getattr(P,k,None) for k in ('exhaustive_saving_control','bellman_mode','interp_method','use_full_kernel','howard_iter','joint_nested_choice','preference_spec','owner_h_bar_scale','Nb','J','Nz','n_parity','n_child_states')}
    preflight = dict(flags=flags, pins=selected['pins'], contract_sha256=PIN, fixed_psi=P.psi_child,
        grid_shape=list(packet['solution'].V.shape), chosen_reference_solve_seconds=selected['receipt']['chosen_solve_seconds'],
        experiment='Only parent housing services softened with zero-anchored softplus delta=.1*hP; childless services unchanged.',
        fixed='Earnings, wealth, entry, mortality, bequests, credit, shocks, child benefits, fiscal rule, supply, housing products, targets and weights.',
        budget='One control Bellman replay; one softened stationary housing equilibrium; no search, no normalization, 25 minute allocation including reporting.')
    write(args.output/'preflight.json', preflight)
    if args.stage == 'inspect': print(json.dumps(preflight)); return
    import e5f_soft_housing_adapter as adapter
    started = time.monotonic()
    state = {'phase': 'self_tests', 'complete': False}
    def beat():
        while not state['complete']:
            write(args.output/'heartbeat.json', dict(**state, elapsed_seconds=time.monotonic()-started))
            time.sleep(30)
    threading.Thread(target=beat, daemon=True).start()
    tests = adapter.run_self_tests()
    write(args.output/'tests.json', tests)
    grid = packet['b_grid']
    # Verify the reconstructed native runtime on the actual control's prices.
    state['phase'] = 'control_bellman_replay'
    Pc = copy.deepcopy(P)
    native = model.solve_markov_income_at_prices(np.asarray(packet['solution'].p_eq), Pc, grid, SD=model.precompute_shared(Pc,grid))
    comparison = {}
    for key in ('V','c_pol','hR_pol','bp_pol','fert_probs','fert2_probs','g'):
        comparison[key] = float(np.max(np.abs(np.asarray(getattr(native,key))-np.asarray(getattr(packet['solution'],key)))))
        np.testing.assert_allclose(getattr(native,key),getattr(packet['solution'],key), atol=1e-10, rtol=1e-12, err_msg='control replay '+key)
    write(args.output/'control_replay.json', comparison)
    control_out = args.output/'control'; control_out.mkdir()
    state['phase'] = 'control_report'
    control = report(packet,runtime,tax,objective,selected['case'],control_out,'control')
    write(args.output/'latest_completed.json', control)
    del native
    state['phase'] = 'soft_equilibrium'
    installed = adapter.install(model, args.output/'adapter', softness=.1)
    write(args.output/'adapter_receipt.json', installed)
    Ps = copy.deepcopy(P)
    sol, Ps, prices, fiscal = runtime['solve_balanced_initial_equilibrium'](model=model, parameters=Ps, b_grid=grid,
        initial_prices=np.asarray(packet['solution'].p_eq).copy(), payroll_tax=float(P.tau_pay), marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    shared = model.precompute_shared(Ps, grid); Ps._fert2_probs = sol.fert2_probs.copy()
    calendar = runtime['primitive'].pf.calendar
    policy = calendar.policy_from_solution(sol,prices,Ps,grid,shared)
    pre, reconstruction = calendar.reconstruct_stationary_pre_fertility(sol,policy,Ps,grid,shared)
    assert reconstruction['stationary_post_fertility_nesting_l1'] < 5e-9
    supply = calendar.HousingSupplyRule('static-elastic', float(prices[0]), float(Ps.H0[0]*(Ps.user_cost_rate*prices[0]/Ps.r_bar[0])**Ps.xi_supply[0]),float(Ps.xi_supply[0]))
    evaluation = calendar.evaluate_period(prices,pre,Ps,grid,shared,calendar.SolveCounter(),supply_rule=supply,supplied_policy=policy)
    soft_packet = dict(parameters=Ps,b_grid=grid,evaluation=evaluation,shared=shared,supply_rule=supply,solution=sol,stationary_g_pre=pre,ancestry_contract_sha256=PIN,demographic_seed=packet.get('demographic_seed'))
    soft_out = args.output/'soft'; soft_out.mkdir()
    with gzip.open(soft_out/'initial_state.pkl.gz','wb',compresslevel=1) as stream: pickle.dump(soft_packet,stream,protocol=5)
    state['phase'] = 'soft_report'
    soft = report(soft_packet,runtime,tax,objective,selected['case'],soft_out,'soft')
    write(args.output/'latest_completed.json', soft)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    hp = float(np.max(shared.hb_flat)); h = np.linspace(.001,4,500)
    axes[0].plot(h[h>hp], h[h>hp]-hp,label='Hard requirement')
    axes[0].plot(h,[adapter.soft_services(float(x),hp,.1)[0] for x in h],label='Soft requirement')
    axes[0].axvline(hp,color='gray',ls=':',lw=1)
    axes[0].set(xlabel='Housing (rooms)',ylabel='Effective housing services'); axes[0].legend()
    for j,(case,result) in enumerate((('Hard',control),('Soft',soft))):
        rooms = next(r['model'] for r in result['target_fit'] if r['moment']=='first_birth_rooms')
        axes[1].bar(j,rooms); axes[1].text(j,rooms,f'{rooms:.3f}',ha='center',va='bottom')
        shares = [100*r['birth_share_positive_interior'] for r in result['incentives'][:2]]
        axes[2].bar(np.arange(2)+(j-.5)*.34,shares,.34,label=case)
    axes[1].set(xticks=[0,1],xticklabels=['Hard','Soft'],ylabel='First-birth housing change (rooms)',title='Inherited stationary proxy')
    axes[2].set(xticks=[0,1],xticklabels=['First births','Later births'],ylabel='Births with positive attempt-minus-wait gap (%)'); axes[2].legend()
    fig.suptitle('Fixed preferences; softened housing requirement only',fontsize=12)
    fig.savefig(args.output/'comparison.png',dpi=160); plt.close(fig)
    write(args.output/'complete.json', dict(status='completed_experimental_comparison',elapsed_seconds=time.monotonic()-started, job=os.environ['SLURM_JOB_ID'], control=control,soft=soft,preflight=preflight))
    state['complete'] = True
    print(json.dumps(dict(status='complete',elapsed_seconds=time.monotonic()-started)))


if __name__ == '__main__': main()
