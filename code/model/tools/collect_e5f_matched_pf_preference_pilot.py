"""Verify locally collected timing-pilot receipts and compare dated birth counts.

Read-only with respect to numerical outputs; never launches or reruns a model.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from e5f_matched_pf_birth_path import compare_birth_path, read_rows

FINGERPRINT = '3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497'
PARENT_SUMMARY = '6025e0c3734bd90e2f210172e109d2413fc79a192daf03a32f55d669f2b1144c'
SHAPES = (-0.5, 0., 0.5)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def write_csv(path, rows):
    if rows:
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def verify_case(directory, shape, phase, source, parent, parent_fit):
    s, c = load(directory/'summary.json'), load(directory/'contract.json')
    trial = directory/'evaluation_001'
    e = load(trial/'summary.json')
    require(not (directory/'failure.json').exists(), 'Case also has a failure receipt')
    require(s['status'] == 'passed_conditional_preference_pilot'
            and s['shape'] == shape and s['evaluations'] == 1, 'Case identity/status')
    require(not s['production_promoted'] and not s['calibrated_history'], 'Diagnostic scope')
    require(c['shape_coefficient'] == shape and c['conditional_only'] is True, 'Case contract')
    years = [2007 + 4*i for i in range(6 if phase == 'smoke' else 100)]
    require(c['path_date_count'] == len(years), 'Contract horizon')
    require(e['mapping_valid'] and e['target_fingerprint'] == FINGERPRINT, 'Mapping or targets')
    require(e['contract_sha256'] == c['contract_sha256'], 'Contract identity')
    for name, expected in c['source_sha256'].items():
        require(digest(source/name) == expected, 'Source mismatch: ' + name)
    for name, expected in e['artifact_sha256'].items():
        require(digest(trial/name) == expected, 'Artifact mismatch: ' + name)
    for name, gate in e['gates'].items():
        require(gate['passed'] and abs(gate['value']) <= gate['tolerance'], 'Gate: ' + name)
    path = read_rows(trial/'transition_path.csv')
    require([int(r['calendar_year']) for r in path] == years == e['years'], 'Dates')
    residuals = []
    for row, expected, psi in zip(path, e['residual'], c['psi_path']):
        value = (float(row['housing_demand'])-float(row['housing_supply']))/float(row['housing_supply'])
        require(math.isfinite(value) and abs(value-expected) < 1e-12, 'Market residual')
        require(float(row['psi_child']) == psi, 'Preference-path application')
        require(float(row['feasibility_frontier_projection_mass']) <= 1e-6, 'Feasibility')
        require(int(row['nonfinite_distribution_count']) == 0, 'Nonfinite distribution')
        residuals.append(value)
    require(abs(max(map(abs, residuals))-e['maximum_market_residual']) < 1e-12, 'Market maximum')
    psi = c['psi_path']
    require(psi[0] == 0.2900515293650047 and psi[4] == -0.03133626685497737
            and all(v == psi[4] for v in psi[4:]), 'Fixed preference endpoints')
    require(all(a >= b for a,b in zip(psi,psi[1:])), 'Preference monotonicity')
    observations = load(trial/'observed_dates.json')
    require([r['calendar_year'] for r in observations] == years, 'Budget observer dates')
    for row, residual in zip(observations, residuals):
        budget = row['budget']
        require(budget['budget_excess_mass'] == 0
                and budget['maximum_occupied_excess'] <= budget['budget_tolerance'], 'Household budget')
        require(budget['actual_rent'] > 0 and abs(row['market_residual']-residual) < 1e-12, 'Dated observer')
    fits = read_rows(trial/'target_fit.csv')
    require(len(fits) == 12, 'Full target table')
    for row, reference in zip(fits, parent_fit):
        require(all(row[k] == reference[k] for k in ('moment','target','weight')), 'Target/weight drift')
        gap = float(row['model'])-float(row['target'])
        loss = gap*gap*float(row['weight'])
        require(abs(gap-float(row['gap'])) < 1e-12, 'Fit gap')
        require(math.isclose(loss,float(row['loss_contribution']),rel_tol=1e-12,abs_tol=1e-12), 'Loss row')
    require(math.isclose(sum(float(r['loss_contribution']) for r in fits),e['loss'],rel_tol=1e-12), 'Loss sum')
    require(e['artifact_sha256']['parameters.csv'] == parent['artifact_sha256']['parameters.csv'], 'Parameter drift')
    if phase == 'main' and shape == 0:
        require(s['mode'] == 'replay' and s['residual_replay_gap'] <= 2e-10, 'Baseline replay')
        for name in ('target_fit.csv','parameters.csv','measurement.json','transition_path.csv'):
            require(e['artifact_sha256'][name] == parent['artifact_sha256'][name], 'Exact baseline artifact replay')
    return dict(shape=shape, phase=phase, source_files_verified=len(c['source_sha256']),
        artifact_files_verified=len(e['artifact_sha256']), gates_verified=len(e['gates']),
        dates_verified=len(years), household_budget_excess_mass=0,
        maximum_market_residual=e['maximum_market_residual'], market_tolerance=2e-4,
        market_residual_passes=e['maximum_market_residual'] <= 2e-4,
        inherited_objective=e['loss'], mode=s['mode'], seconds=s['elapsed_seconds'],
        residual_replay_gap=s['residual_replay_gap'], terminal_distance=e['terminal_distance'],
        historical_equilibrium_certified=False, production_promoted=False), path, fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot-root', required=True, type=Path)
    parser.add_argument('--source-root', required=True, type=Path)
    parser.add_argument('--parent-evaluation', required=True, type=Path)
    args = parser.parse_args()
    root, parent_dir = args.pilot_root, args.parent_evaluation
    require(digest(parent_dir/'summary.json') == PARENT_SUMMARY, 'Pinned parent summary')
    parent = load(parent_dir/'summary.json')
    for name in ('target_fit.csv','parameters.csv','transition_path.csv'):
        require(digest(parent_dir/name) == parent['artifact_sha256'][name], 'Parent artifact: '+name)
    parent_fit = read_rows(parent_dir/'target_fit.csv')
    empirical = read_rows(root/'fertility_data/empirical_blocks.csv')
    baseline = compare_birth_path(read_rows(parent_dir/'transition_path.csv'), empirical)
    output = root/'computation'
    output.mkdir(exist_ok=True)
    verified, pending, births, all_fits = [], [], [], []
    for phase in ('smoke','main'):
        for i, shape in enumerate(SHAPES):
            case = output/phase/f'case_{i}'
            if not (case/'summary.json').exists():
                pending.append(dict(phase=phase,case=i,reason='No locally collected completed summary'))
                continue
            receipt, path, fits = verify_case(case,shape,phase,args.source_root,parent,parent_fit)
            verified.append(receipt)
            if phase == 'main':
                comparison = compare_birth_path(path,empirical,
                    anchor_first_block_births=baseline['first_block_births'])
                write_json(case/'birth_path_comparison.json',comparison)
                receipt['birth_shape_mse'] = comparison['shape_mean_squared_gap']
                receipt['first_block_birth_change'] = comparison['first_block_change_from_anchor']
                births.extend(dict(shape=shape,**row) for row in comparison['rows'])
                all_fits.extend(dict(shape=shape,**row) for row in fits)
    write_json(output/'verified_receipts.json',dict(verified=verified,pending=pending,
        source_commit='e399c90e',production_promoted=False,
        comparison_scope='Inherited parameters; changed preference paths use fixed parent prices'))
    write_json(output/'baseline_birth_comparison.json',baseline)
    write_csv(output/'all_main_birth_comparisons.csv',births)
    write_csv(output/'all_main_target_fits.csv',all_fits)
    parameters = read_rows(parent_dir/'parameters.csv')
    require(len(parameters)==15 and sum(r['is_free_parameter']=='True' for r in parameters)==11,'Full parameters')
    write_csv(output/'all_inherited_parameters.csv',[
        dict(**r,held_fixed_in_this_pilot=True) for r in parameters])
    print(json.dumps(dict(verified_cases=len(verified),pending_cases=pending,
        main_cases=[r for r in verified if r['phase']=='main']),indent=2))


if __name__ == '__main__':
    main()
