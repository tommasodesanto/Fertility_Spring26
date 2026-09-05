#!/usr/bin/env python3
"""Analyze saved empirical support receipts and unchanged model calibration fits."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PRIMARY = ROOT / 'code/data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1'
REFINEMENT = ROOT / 'output/model/e5f_morning_refinement_20260905a'
DEFAULT = ROOT / 'output/model/e5f_first_birth_measurement_review_20260905a'


def rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def write_rows(path, data):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0])); w.writeheader(); w.writerows(data)


def empirical(output):
    smoke = json.loads((output / 'smoke/run_receipt.json').read_text())
    full = json.loads((output / 'full/run_receipt.json').read_text())
    verification = json.loads((output / 'post_run_verification.json').read_text())
    if smoke['status'] != 'pass' or not verification['primary_outputs_unchanged'] or not verification['source_builder_unchanged']:
        raise RuntimeError('Sample replay or primary-source verification has not passed')
    support_path = output / 'input_cohort_event_support.csv'
    assert hashlib.sha256(support_path.read_bytes()).hexdigest() == verification['input_support_sha256']
    assert verification['input_support_replays_identical']
    data = rows(support_path)
    treated = [r for r in data if r['never_treated'] == '0']
    assert sum(int(r['observations']) for r in data) == 49872
    base_cohorts = {int(float(r['first_child_year'])) for r in treated if float(r['K']) == -2}
    means = []
    weights = {}
    for k in (-2, -1, 0, 1, 2, 3, 4):
        sub = [r for r in treated if float(r['K']) == k]
        W = sum(float(r['weight_sum']) for r in sub)
        by_cohort = {}
        for r in sub:
            g = int(float(r['first_child_year']))
            by_cohort[g] = by_cohort.get(g, 0.) + float(r['weight_sum']) / W
        weights[k] = by_cohort
        means.append({'event_time': k, 'observations': sum(int(r['observations']) for r in sub),
                      'mean_birth_year': sum(float(r['first_child_year']) * float(r['weight_sum']) for r in sub) / W,
                      'post1997_weight_share': sum(float(r['weight_sum']) for r in sub if r['post1997'] == '1') / W,
                      'birth2007plus_weight_share': sum(float(r['weight_sum']) for r in sub if float(r['first_child_year']) >= 2007) / W,
                      'no_kminus2_cohort_weight_share': sum(float(r['weight_sum']) for r in sub if int(float(r['first_child_year'])) not in base_cohorts) / W,
                      'scope': 'pre_regression_input_sample_not_final_estimation_sample'})
    write_rows(output / 'empirical_support_summary.csv', means)
    primary = rows(PRIMARY / 'target_receipt.csv')[0]
    saved_events = rows(PRIMARY / 'event_study_estimates.csv')
    event = {int(float(r['relative_time'])): float(r['b']) for r in saved_events}
    event_se = {int(float(r['relative_time'])): r['se'] for r in saved_events}
    target = float(primary['estimate'])
    contrasts = [
        {'contrast': 'primary_-1_to_3', 'estimate': target, 'standard_error': primary['standard_error'], 'status': 'unchanged_primary'},
        {'contrast': 'alternative_-2_to_2', 'estimate': event[2], 'standard_error': event_se[2], 'status': 'diagnostic_point_from_existing_event_table'},
        {'contrast': 'pooled_pre_post_two_year_cells', 'estimate': (event[2]+event[3]-event[-1])/2, 'standard_error': 'not_available', 'status': 'diagnostic_point_from_existing_event_table'},
        {'contrast': 'post_birth_0_to_4_different_estimand', 'estimate': event[4]-event[0], 'standard_error': 'not_available', 'status': 'different_estimand_not_a_replacement'}]
    write_rows(output / 'event_window_diagnostics.csv', contrasts)
    union = set(weights[-1]) | set(weights[3])
    tv = sum(abs(weights[-1].get(g, 0.) - weights[3].get(g, 0.)) for g in union) / 2
    return {'target': target, 'target_se': float(primary['standard_error']),
            'input_cohort_weight_total_variation': tv,
            'treated_calendar_range': [min(float(r['first_child_year'])+float(r['K']) for r in treated), max(float(r['first_child_year'])+float(r['K']) for r in treated)],
            'support': means, 'event_window_diagnostics': contrasts,
            'analysis_scope': 'input support and saved primary estimates; no new regression coefficients or covariance',
            'full_replay_status': full['status'], 'full_replay_error': full.get('error'),
            'full_replay_elapsed_seconds': full['elapsed_seconds'],
            'common_cohort_weight_effect': 'not_quantified',
            'standard_errors_of_window_differences': 'not_available',
            'post_run_verification': verification}


def normalization_check(output):
    receipt = json.loads((output / 'design_check/smoke/run_receipt.json').read_text())
    assert receipt['status'] == 'pass'
    assert all(receipt['sample_design_assertions'].values())
    data = rows(output / 'input_cohort_event_support.csv')
    treated = [r for r in data if r['never_treated'] == '0']
    groups = {int(float(r['first_child_year'])) for r in treated}
    baseline = {int(float(r['first_child_year'])) for r in treated if float(r['K']) == -2}
    weights = {}
    for k in (-1, 3):
        sub = [r for r in treated if float(r['K']) == k]
        total = sum(float(r['weight_sum']) for r in sub)
        weights[k] = {g: sum(float(r['weight_sum']) for r in sub if int(float(r['first_child_year'])) == g)/total for g in groups}
    detail = [{'cohort': g, 'weight_minus1': weights[-1][g], 'weight_plus3': weights[3][g],
               'no_observed_kminus2': g not in baseline,
               'null_loading': weights[3][g]-weights[-1][g] if g not in baseline else 0.}
              for g in sorted(groups)]
    write_rows(output / 'normalization_null_directions.csv', detail)
    # Use a witness observed at BOTH endpoints, avoiding an endpoint-only example.
    witnesses = [r for r in detail if r['no_observed_kminus2'] and r['weight_minus1'] > 0 and r['weight_plus3'] > 0]
    witness = max(witnesses, key=lambda r: abs(r['null_loading']))
    g = witness['cohort']
    # Shift every supported coefficient of cohort g by one room and every
    # person fixed effect in that cohort by minus one. This leaves each row
    # unchanged because its indicators partition the observations exactly.
    fitted_value_shifts = [int(int(float(r['first_child_year'])) == g) * (int(float(r['K']) != -2)-1) for r in treated]
    assert max(abs(x) for x in fitted_value_shifts) == 0
    assert abs(witness['null_loading']) > 1e-6
    return {'status': 'pass', 'cohorts_without_observed_reference': len(groups-baseline),
            'witness': witness, 'maximum_absolute_fitted_value_shift': max(abs(x) for x in fitted_value_shifts),
            'interpretation': 'normalization dependence of nominal common-reference contrast; not a corrected target or estimated bias',
            'sample_design_receipt': 'design_check/smoke/run_receipt.json'}


def model(output):
    data = []
    fits = []; parameters = []
    for stage in ('round1', 'round2'):
        grouped = {}
        for r in rows(REFINEMENT / stage / 'report/all_target_fits.csv'):
            fits.append({'stage': stage, **r})
            z = grouped.setdefault(r['case_id'], {'stage': stage, 'case_id': r['case_id'], 'label': r['label'], 'loss': 0.})
            z[r['moment']] = float(r['model']); z['loss'] += float(r['loss_contribution'])
        for r in rows(REFINEMENT / stage / 'report/all_parameters.csv'):
            parameters.append({'stage': stage, **r})
        data.extend(grouped.values())
    assert len(data) == 35 and len(fits) == 35*12
    fingerprints = set()
    for r in data:
        p = REFINEMENT / r['stage'] / f"task_{int(r['case_id']):03d}" / 'summary.json'
        fingerprints.add(json.loads(p.read_text())['target_fingerprint'])
    assert fingerprints == {'3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497'}
    selected = next(r for r in data if r['stage'] == 'round2' and r['case_id'] == '9')
    branches = rows(REFINEMENT / 'round2/task_009/cases/task_001/dated_first_birth_housing_did.csv')
    terminal = branches[-1]
    response = float(terminal['treated_mean_housing']) - float(terminal['control_mean_housing'])
    assert abs(response-selected['housing_increment_0to1']) < 1e-12
    write_rows(output / 'observed_same_target_candidates.csv', data)
    write_rows(output / 'all_target_fits.csv', fits)
    write_rows(output / 'all_parameters.csv', parameters)
    write_rows(output / 'selected_dated_housing_branches.csv', branches)
    return {'evaluated_coordinate_and_joint_cases': len(data), 'target_fingerprint': fingerprints.pop(),
            'selected': selected, 'minimum_observed_rooms': min(r['aggregate_mean_occupied_rooms_18_85'] for r in data),
            'maximum_observed_birth_response': max(r['housing_increment_0to1'] for r in data),
            'destination_treated_rooms': float(terminal['treated_mean_housing']),
            'destination_control_rooms': float(terminal['control_mean_housing']),
            'destination_continuation_birth_rate': float(terminal['treated_continuation_births']) / float(terminal['destination_mass']),
            'no_new_equilibrium_solved': True}


def plots(output, empirical_result, model_result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    events = rows(PRIMARY / 'event_study_estimates.csv')
    x = np.array([float(r['relative_time']) for r in events]); y = np.array([float(r['b']) for r in events]); se = np.array([float(r['se']) for r in events])
    fig, ax = plt.subplots(figsize=(8.7, 4.))
    ax.errorbar(x, y, yerr=1.96*se, fmt='o-', color='#204b6b', markersize=4, capsize=3, linewidth=1)
    ax.axhline(0, color='0.55', linewidth=.7); ax.axvline(0, color='0.55', linewidth=.7, linestyle='--')
    ax.scatter([-1,3], [y[np.where(x==-1)[0][0]],y[np.where(x==3)[0][0]]],color='#b84f35',s=55,zorder=5)
    ax.set(xlabel='Years relative to first birth (tails pooled)', ylabel='Rooms relative to event year -2',title='Full empirical event-study path; the primary contrast remains -1 to +3',xticks=range(-7,12,2))
    fig.tight_layout(); fig.savefig(output/'empirical_event_profile.png',dpi=190); plt.close(fig)
    data = rows(output/'observed_same_target_candidates.csv')
    fig, ax = plt.subplots(figsize=(8.7,4.))
    for stage, marker, color in [('round1','o','#7597ad'),('round2','s','#204b6b')]:
        z=[r for r in data if r['stage']==stage]
        ax.scatter([float(r['aggregate_mean_occupied_rooms_18_85']) for r in z], [float(r['housing_increment_0to1']) for r in z], marker=marker,c=color,s=35,label=stage)
    ax.scatter([5.779970481941968],[.7202462623815278], marker='*',s=140,c='#b84f35',label='Target pair')
    z=model_result['selected']; ax.scatter([z['aggregate_mean_occupied_rooms_18_85']],[z['housing_increment_0to1']],facecolors='none',edgecolors='#b84f35',s=150,linewidths=1.5,label='Selected review candidate')
    ax.set(xlabel='Mean occupied rooms',ylabel='First-birth housing response (rooms)',title='Observed 35-case neighborhood; this is not a global feasibility frontier')
    ax.legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(output/'observed_housing_fit.png',dpi=190);plt.close(fig)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=DEFAULT);args=p.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    e=empirical(output);e['normalization_check']=normalization_check(output);m=model(output);plots(output,e,m)
    payload={'empirical':e,'model':m}
    (output/'measurement_summary.json').write_text(json.dumps(payload,indent=2))
    print(json.dumps({'empirical':{k:v for k,v in e.items() if k not in ('source_read_only_replay','support')},'model':m},indent=2))
