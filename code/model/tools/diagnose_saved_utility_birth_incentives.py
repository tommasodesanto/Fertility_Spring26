#!/usr/bin/env python3
"""Torch-only, no-solve incentive audit of a pinned utility-comparison checkpoint.

Load through the comparison's authenticated collector. Use saved pre-fertility
mass, and validate its timing by reproducing the saved post-fertility population
and birth total. Invert both stored logit probabilities (not 1-p, which loses
precision). Zero probabilities are reported separately, never inverted.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np


def summarize(label, records):
    w, p0, p1, pi, kappa = (np.concatenate([r[i] for r in records]) for i in range(5))
    births = w * p1 * pi
    interior = (p0 > 0) & (p1 > 0)
    negative, positive = interior & (p1 < p0), interior & (p1 > p0)
    gap = np.log(p1[interior]) - np.log(p0[interior])
    weights = w[interior]
    order = np.argsort(gap)
    cumulative = np.cumsum(weights[order])
    quantiles = [None] * 3
    if cumulative.size and cumulative[-1] > 0:
        quantiles = np.interp(np.array([.1, .5, .9]) * cumulative[-1], cumulative, gap[order]).tolist()
    mass, birth_total = float(w.sum()), float(births.sum())
    def risk(mask):
        return float(w[mask].sum() / mass) if mass else None
    def birth_share(mask):
        return float(births[mask].sum() / birth_total) if birth_total else None
    return dict(group=label, risk_mass=mass, expected_births=birth_total,
        risk_share_negative_interior=risk(negative), risk_share_positive_interior=risk(positive),
        risk_share_tied_interior=risk(interior & (p1 == p0)),
        risk_share_zero_try=risk((p1 == 0) & (p0 > 0)),
        risk_share_zero_wait=risk((p0 == 0) & (p1 > 0)),
        risk_share_unavailable=risk((p0 == 0) & (p1 == 0)),
        birth_share_negative_interior=birth_share(negative),
        birth_share_positive_interior=birth_share(positive),
        birth_share_boundary=birth_share(~interior),
        gap_over_kappa_p10=quantiles[0], gap_over_kappa_median=quantiles[1],
        gap_over_kappa_p90=quantiles[2])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    import collect_e5f_utility_comparison as collector
    contract, runner = collector._load_contract(args.contract)
    arm = 'floor_linear'
    _, _, _, tax, _, _, runtime, _ = runner.setup(contract, arm, args.output / 'runtime')
    selection = collector.read_json(args.run_root / arm / 'selected.json')
    assert selection['arm'] == arm
    assert selection['contract_sha256'] == os.environ[collector.ENV_PIN]
    selected = collector.scientific_checkpoint(selection['original_case_output'], contract, arm, runtime, tax)
    packet, receipt = selected['packet'], selected['receipt']
    P, solution = packet['parameters'], packet['solution']
    policy = packet['evaluation'].policy
    from intergen_eqscale_seq_optimized.parameters import (
        get_fecundity_by_age, independent_child_maturation_active,
        readiness_gate_active, readiness_settled_state,
    )
    assert bool(P.sequential_births) and not bool(getattr(P, 'joint_nested_choice', False))
    assert independent_child_maturation_active(P)
    assert int(P.A_f_end) < int(P.J)
    g = np.asarray(packet['stationary_g_pre'])
    assert np.all(np.isfinite(g)) and np.all(g >= 0)
    post = g.copy()
    fec = np.asarray(get_fecundity_by_age(P))
    assert np.all((fec >= 0) & (fec <= 1))
    first, continuation = np.asarray(policy.fert_probs), np.asarray(policy.fert2_probs)
    np.testing.assert_array_equal(first, solution.fert_probs)
    # Zero projection means the authenticated pre-choice population is unchanged.
    assert receipt['operator_gates']['stationary_feasibility_projection_mass'] == 0
    groups = {'first_birth': [], 'later_births': [], 'all_births': []}
    age_rows = []
    unavailable_mass = 0.0
    settled = readiness_settled_state(P)
    for j in range(P.J):
        if not P.A_f_start <= j + 1 <= P.A_f_end:
            continue
        age = P.age_start + j * P.da
        for n in range(P.n_parity - 1):
            child_states = [settled] if n == 0 else range(n + 1)
            margin = 'first_birth' if n == 0 else 'later_births'
            scale = float(P.kappa_fert if n == 0 else (
                P.kappa_fert_continuation if P.kappa_fert_continuation is not None else P.kappa_fert))
            assert np.isfinite(scale) and scale > 0
            local_records = []
            for m in child_states:
                w = g[:, :, :, j, :, n, m]
                pair = first[:, :, :, j, :, :2] if n == 0 else continuation[:, :, :, j, :, :, n - 1, m]
                p0, p1 = pair[..., 0], pair[..., 1]
                occupied = w > 0
                assert np.all(np.isfinite(pair[occupied]))
                assert np.all((pair[occupied] >= 0) & (pair[occupied] <= 1 + 1e-12))
                # The native solver zeros BOTH probabilities when both branch
                # values hit its dead-state cutoff. Account for that mass rather
                # than misclassifying it as logit underflow or a forced action.
                unavailable = (p0 == 0) & (p1 == 0)
                unavailable_mass += float(w[unavailable].sum())
                np.testing.assert_allclose((p0 + p1)[occupied & ~unavailable], 1, atol=1e-12, rtol=0)
                born = w * p1 * fec[j]
                post[:, :, :, j, :, n, m] -= born
                destination = 1 if n == 0 else m + 1
                post[:, :, :, j, :, n + 1, destination] += born
                record = (w[occupied], p0[occupied], p1[occupied],
                    np.full(occupied.sum(), fec[j]), np.full(occupied.sum(), scale))
                local_records.append(record)
                groups[margin].append(record)
                groups['all_births'].append(record)
            row = summarize(f'birth_{n + 1}_age_{age:g}', local_records)
            row.update(age=float(age), birth_number=n + 1, conception_probability=float(fec[j]), shock_scale=scale)
            age_rows.append(row)
    summaries = [summarize(name, records) for name, records in groups.items()]
    assert unavailable_mass <= 1e-12, unavailable_mass
    post_error = float(np.abs(post - solution.g_beginning_distribution).sum())
    birth_error = abs(summaries[-1]['expected_births'] - receipt['operator_gates']['stationary_birth_flow'])
    assert post_error < 1e-10, post_error
    assert birth_error < 1e-10, birth_error
    for filename, rows in [('summary.csv', summaries), ('by_age_birth_number.csv', age_rows)]:
        with (args.output / filename).open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    config_names = ['sequential_births', 'joint_nested_choice', 'child_state_mode',
        'age_start', 'da', 'J', 'A_f_start', 'A_f_end', 'kappa_fert',
        'kappa_fert_continuation', 'first_birth_fixed_cost', 'psi_child',
        'fecundity_omega1', 'fecundity_omega2', 'utility_comparison_arm',
        'utility_child_benefit_exponent', 'transfer_floor_Gn']
    result = dict(status='verified_saved_policy_diagnostic', new_solves=0,
        case=str(selected['case']), pins=selected['pins'],
        source_manifest_sha256=receipt['source_manifest_sha256'],
        contract_sha256=os.environ[collector.ENV_PIN],
        script_sha256=collector.sha256(__file__), job_id=os.environ.get('SLURM_JOB_ID'),
        config={name: getattr(P, name, None) for name in config_names},
        readiness_gate=readiness_gate_active(P), fecundity_by_age=fec.tolist(),
        post_fertility_population_l1_error=post_error, birth_flow_error=birth_error,
        unavailable_choice_mass=unavailable_mass,
        summaries=summaries,
        interpretation='Attempt-minus-wait gap excluding the current taste draw, but retaining solved future shock-inclusive values. Negative means prefer waiting now; not a never-have-children comparison or a zero-shock equilibrium.',
        boundary_rule='Zero try/wait probabilities reported separately: no inverse-logit magnitude assigned; numerical underflow and infeasible alternatives are not distinguished.',
        normalization='Aggregate mean fertility was imposed in the original point; this diagnostic does not independently validate tastes or identify shock dominance.',
        calibration_tables=dict(targets=str(selected['case'] / 'target_fit.csv'), parameters=str(selected['case'] / 'parameters.csv')))
    (args.output / 'receipt.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for n in (1, 2, 3):
        rows = [r for r in age_rows if r['birth_number'] == n]
        label = f'Child {n}'
        for ax, key in zip(axes, ['risk_share_positive_interior', 'birth_share_positive_interior']):
            ax.plot([r['age'] for r in rows], [100*r[key] if r[key] is not None else np.nan for r in rows], marker='o', label=label)
    axes[0].set_title('At-risk households with positive gap')
    axes[1].set_title('Births from states with positive gap')
    for ax in axes:
        ax.set_xlabel('Parent age'); ax.set_ylabel('Percent'); ax.set_ylim(-2, 102); ax.legend(); ax.grid(alpha=.2)
    fig.suptitle('Saved floor/linear point: trying now versus waiting\nCurrent taste draw excluded; future taste uncertainty retained', fontsize=11)
    fig.savefig(args.output / 'incentives.png', dpi=150)
    plt.close(fig)
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    main()
