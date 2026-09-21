"""Review the saved, fixed-contract September 21 panel; never run the model."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'output/model/native_financing_diagnostic_20260919/specification_followup/quantification_v1/sensitivity'
OUT = BASE / 'lead_analysis'


def read(path):
    return json.loads(path.read_text())


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def save_csv(name, rows):
    with (OUT / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    summary_path = BASE / 'collected/production/summary.json'
    summary = read(summary_path)
    plan = read(BASE / 'collected/prior_plan.json')
    assert hashlib.sha256((BASE / 'collected/prior_plan.json').read_bytes()).hexdigest() == '6f3665488e63db302c14ae810dc5c68f44d403ee85b6dfd6b48cb88fa1334ccd'
    rows = summary['rows']
    assert summary['status'] == 'verified_selection' and len(rows) == 28
    assert all(row['status'] == 'valid' for row in rows), 'Partial panel needs a separate explicit review'
    anchor = next(row for row in rows if row['case_id'] == 'smoke_anchor_01')
    selected = summary['selected']
    full_targets = []
    for row in rows:
        receipt = row['receipt']
        fit = receipt['target_fit']
        assert len(fit) == 13 and len(receipt['parameters']) == 17
        assert receipt['contract_sha256'] == plan['objective_canonical_sha256']
        assert receipt['candidate_payload_fingerprint'] == plan['candidate_payload_fingerprint']
        assert receipt['source_fingerprints'] == anchor['receipt']['source_fingerprints']
        signature = [(x['restriction_id'], x['target'], x.get('actual_weight'), x.get('scored')) for x in fit]
        assert fingerprint(signature) == plan['target_signature']
        assert receipt['_summary']['status'] == 'verified_scored_candidate'
        assert receipt['normalization']['absolute_gap'] <= 5e-4
        loss = 0.
        for moment in fit:
            gap = float(moment['model']) - float(moment['target'])
            assert np.isclose(gap, moment['gap'], atol=1e-12, rtol=0)
            if moment['scored']:
                contribution = gap * gap * moment['actual_weight']
                assert np.isclose(contribution, moment['loss_contribution'], atol=1e-10, rtol=1e-12)
                loss += contribution
            full_targets.append(dict(case_id=row['case_id'], **{k: moment[k] for k in ('restriction_id', 'label', 'target', 'model', 'gap', 'actual_weight', 'loss_contribution', 'scored')}))
        assert np.isclose(loss, row['objective'], atol=1e-9, rtol=1e-12)
        for param in receipt['parameters']:
            if param['structural_coordinate']:
                name = param['parameter']; lower, upper = plan['parameter_bounds'][name]
                assert lower <= param['estimate'] <= upper
                assert abs(param['estimate'] - row['parameters'][name]) <= 1e-12
    for row in rows:
        reference = anchor if row['group'] == 'anchor' else selected if row['group'] == 'final' else None
        if reference:
            assert row['numeric_signature'] == reference['numeric_signature']
            assert row['array_signature'] == reference['array_signature']
    smoke = read(BASE / 'collected/smoke/smoke_summary.json')
    assert anchor['array_signature'] == smoke['reference_array_signature']
    assert selected['objective'] == min(r['objective'] for r in rows)
    assert sum(r['stationary_solves'] for r in rows) == summary['stationary_solves_actual'] == 168
    verified_gates = []
    for name in ('smoke_anchor_01', 'smoke_anchor_02', 'selected', 'final_rep_01', 'final_rep_02'):
        gate = read(BASE / 'collected' / name / 'receipts/verified_evaluation_receipt.json')
        assert gate['status'] == 'verified' and gate['numerical_gates_verified'] is True
        assert gate['source_fingerprints'] == selected['receipt']['source_fingerprints']
        verified_gates.append(name)
    graph_manifest = read(BASE / 'collected/graph_hash_verification.json')
    graph_sets = {}
    for item in graph_manifest['hashes']:
        graph = BASE / 'collected' / item['set'] / 'standard_diagnostics' / item['file']
        digest = hashlib.sha256(graph.read_bytes()).hexdigest()
        assert digest == item['sha256']
        graph_sets.setdefault(item['set'], {})[item['file']] = digest
    assert len(graph_sets['selected']) == 17
    assert graph_sets['selected'] == graph_sets['final_rep_01'] == graph_sets['final_rep_02']
    save_csv('all_target_fit.csv', full_targets)
    selected_targets = [{k: x[k] for k in ('restriction_id','label','target','model','gap','actual_weight','loss_contribution','scored')} for x in selected['receipt']['target_fit']]
    save_csv('selected_target_fit.csv', selected_targets)
    params = []
    for row in selected['receipt']['parameters']:
        name = row['parameter']; bounds = plan['parameter_bounds'].get(name)
        lower, upper = bounds if bounds else ('', '')
        distance = min(row['estimate'] - lower, upper - row['estimate']) / (upper - lower) if bounds else None
        side = ('lower' if row['estimate'] - lower < upper - row['estimate'] else 'upper') if bounds and distance <= .01 else ''
        params.append(dict(parameter=name, estimate=row['estimate'], actual_lower=lower, actual_upper=upper,
                           near_bound=side if bounds else 'not searched', fractional_bound_distance=distance,
                           restriction_or_status=row['status'], structural_coordinate=row['structural_coordinate']))
    save_csv('selected_parameters.csv', params)
    scored = [r for r in anchor['receipt']['target_fit'] if r['scored']]
    names = [r['restriction_id'] for r in scored]
    weights = np.sqrt([r['actual_weight'] for r in scored])
    def moment_vector(row):
        moments = {r['restriction_id']: r['model'] for r in row['receipt']['target_fit']}
        return np.array([moments[n] for n in names])
    base_vector = moment_vector(anchor)
    by_key = {(r.get('coordinate'), r.get('group'), r.get('sign')): r for r in rows if r['group'] in ('full', 'half')}
    def derivative(name, group):
        minus, plus = by_key.get((name, group, -1)), by_key.get((name, group, 1))
        if minus and plus:
            return (moment_vector(plus) - moment_vector(minus)) / (plus['step'] + minus['step']), 'central'
        row = minus or plus
        if row:
            return (moment_vector(row) - base_vector) / (row['sign'] * row['step']), 'one-sided'
        return None, 'missing'
    coordinates = list(plan['parameter_bounds'])
    columns, half_columns, diagnostics, response_rows = [], [], [], []
    for name in coordinates:
        full, method = derivative(name, 'full')
        assert full is not None
        step = next(r['step'] for r in rows if r.get('coordinate') == name and r['group'] == 'full')
        half, half_method = derivative(name, 'half')
        weighted = weights * full * step
        weighted_half = weights * half * step if half is not None else weighted
        columns.append(weighted); half_columns.append(weighted_half)
        relative = float(np.linalg.norm(weighted_half-weighted) / max(np.linalg.norm(weighted), 1e-30)) if half is not None else None
        cosine = float(np.dot(weighted_half, weighted) / max(np.linalg.norm(weighted_half)*np.linalg.norm(weighted), 1e-30)) if half is not None else None
        diagnostics.append(dict(parameter=name, method=method, full_step=step, weighted_response_norm=float(np.linalg.norm(weighted)), half_step_available=half is not None, half_relative_difference=relative, full_half_cosine=cosine))
        for i, target in enumerate(names):
            response_rows.append(dict(parameter=name, target=target, method=method, full_step=step, derivative=float(full[i]), weighted_response_per_full_step=float(weighted[i]), half_derivative=float(half[i]) if half is not None else None))
    J = np.column_stack(columns); Jh = np.column_stack(half_columns)
    singular = np.linalg.svd(J, compute_uv=False)
    singular_half = np.linalg.svd(Jh, compute_uv=False)
    save_csv('paired_response_map.csv', response_rows)
    save_csv('step_stability.csv', diagnostics)
    receipt = dict(status='reviewed_complete_saved_panel', summary_sha256=hashlib.sha256(summary_path.read_bytes()).hexdigest(), objective=plan['objective_canonical_sha256'], target_fingerprint=plan['target_signature'], objective_evaluations=28, verified_target_rows=len(full_targets), verified_loss_contributions=28*12, stationary_solves=168, rejected_or_incomplete=0, exact_anchor_repeats=True, exact_selected_repeats=True, original_reference_arrays_equal=True, selected_loss=selected['objective'], anchor_loss=anchor['objective'], improvement_percent=100*(1-selected['objective']/anchor['objective']), bounds=plan['parameter_bounds'], selected_array_signatures=selected['array_signature'], selected_source_fingerprints=selected['receipt']['source_fingerprints'], scaled_singular_values=singular.tolist(), scaled_condition=float(singular[0]/singular[-1]), half_substitution_condition=float(singular_half[0]/singular_half[-1]), coordinate_scaling='delta parameter divided by planned full step; rows multiplied by sqrt(actual weight); no log transform', response_checks=diagnostics, limitations='Finite local descriptive map, not an identification proof. Beta and h_P one-sided at upper bounds. Four coordinates lack half steps. No missing or rejected probes. Variation across finite differences may reflect curvature or numerical/grid effects; exact repetition alone does not settle accuracy.')
    receipt.update(source_fingerprints_identical_all_28=True, numerical_gate_receipts_checked=verified_gates,
                   standard_graphs_rehashed=51, standard_graphs_equal_across_three_sets=True,
                   source_verification_scope='Certified pinned-run receipts checked for equality; no fresh morning rehash of the remote source tree.')
    (OUT/'review_receipt.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 7))
    lim = float(np.max(np.abs(J)))
    image = ax.imshow(J, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
    ax.set_xticks(range(len(coordinates)), coordinates, rotation=45, ha='right')
    ax.set_yticks(range(len(names)), names)
    ax.set_title('Supplemental: local weighted moment responses\nOne planned full step; central except inward beta and h_P')
    fig.colorbar(image, ax=ax, label='Moment change × square root of working weight')
    fig.tight_layout(); fig.savefig(OUT/'supplemental_response_map.png', dpi=150); plt.close(fig)
    print(json.dumps({k: receipt[k] for k in ('status','selected_loss','improvement_percent','scaled_condition','half_substitution_condition','response_checks')}, indent=2))


if __name__ == '__main__':
    main()
