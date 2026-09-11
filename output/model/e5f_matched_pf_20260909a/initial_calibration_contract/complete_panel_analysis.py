"""Score the complete saved 19-case panel using externally pinned working weights.

Pure numerical diagnostic: no model import/solve, target adoption, optimization,
proposal selection, parameter estimation or inferential claim. Requires explicit
expected byte hashes supplied independently by the caller. Original observation
flags remain unchanged. Writes only the five complete_panel_* result files here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
NAMES = ('beta_annual', 'kappa_fert', 'kappa_fert_continuation', 'chi', 'H0',
         'theta0', 'theta1', 'first_birth_fixed_cost', 'h_P')
IDS = ('initial_normalization', 'cps_childlessness', 'cps_exactly_one',
       'nchs_mean_age', 'nchs_share30', 'wealth_earnings', 'bequest_wealth',
       'old_dispersion', 'mean_rooms', 'ownership_30_55', 'first_birth_rooms',
       'family_rooms', 'recent_parent_ownership')
RECENT = BASE / 'initial_fit_readout/recent_parent_probe/completed_17362130'
RECENT_MOMENT = 'ownership_current_birth_from_empty_dependent_home_minus_current_empty_home_30_55'
PROJECTION = 'uniform_birth_time'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    require(not isinstance(value, bool), 'Boolean supplied as a number')
    value = float(value)
    require(math.isfinite(value), 'Nonfinite numeric input')
    return value


def same(a, b, label, tolerance=2e-12):
    require(math.isclose(number(a), number(b), rel_tol=tolerance,
                         abs_tol=tolerance), label)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def indexed(rows, field):
    result = {row[field]: row for row in rows}
    require(len(result) == len(rows), 'Duplicate ' + field)
    return result


def lookup(data, path):
    for key in path.split('.'):
        data = data[key]
    return number(data)


def coordinate(name, value):
    value = number(value)
    require(value > 0 and (name != 'beta_annual' or value < 1),
            'Parameter outside the logarithmic coordinate domain')
    return math.log(-math.log(value)) if name == 'beta_annual' else math.log(value)


def write_csv(path, rows):
    require(bool(rows), 'Cannot write empty table')
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('weights', 'lead-decision', 'panel', 'recent-panel'):
        parser.add_argument('--expected-' + key + '-sha256', required=True)
    args = parser.parse_args()
    inputs = {}

    def read(path, expected=None, kind='json'):
        path = Path(path).resolve()
        data = path.read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        if expected is not None:
            require(actual == expected, 'Expected input SHA256 mismatch: ' + str(path))
        require(str(path) not in inputs or inputs[str(path)] == actual,
                'Input changed during read: ' + str(path))
        inputs[str(path)] = actual
        if kind == 'csv':
            return list(csv.DictReader(data.decode().splitlines()))
        return json.loads(data)

    weights = read(HERE / 'working_weights.csv', args.expected_weights_sha256, 'csv')
    decision = read(HERE / 'lead_methodological_decision.json', args.expected_lead_decision_sha256)
    panel = read(BASE / 'initial_sensitivity_panel/raw_case_summary.json', args.expected_panel_sha256)
    recent = read(RECENT / 'recent_parent_panel.json', args.expected_recent_panel_sha256)
    panel_receipt = read(BASE / 'initial_sensitivity_panel/panel_receipt.json')
    recent_receipt = read(RECENT / 'recent_parent_panel_receipt.json')
    plan = read(BASE / 'initial_sensitivity_panel/cases.json')
    old_fit = indexed(read(BASE / 'initial_fit_readout/target_fit.csv', kind='csv'), 'restriction_id')
    old_derivatives = read(BASE / 'initial_sensitivity_panel/local_derivatives.csv', kind='csv')
    old_analysis = read(BASE / 'initial_sensitivity_panel/local_analysis.json')
    parameter_reference = read(BASE / 'initial_fit_readout/parameters.csv', kind='csv')
    require(decision['primary_cps_projection'] == PROJECTION, 'Wrong primary CPS projection')
    require(decision['scored_restrictions'] == 12 and decision['structural_coordinates'] == 9,
            'Lead decision changed the task dimensions')
    require(decision['recent_parent_model_moment'] == RECENT_MOMENT, 'Wrong recent-parent analogue')
    require([r['restriction_id'] for r in weights] == list(IDS), 'Incomplete/reordered weight contract')
    require(set(old_fit) == set(IDS), 'Incomplete original target mapping')
    decision_rows = indexed(decision['rows'], 'restriction_id')
    require(set(decision_rows) == set(IDS), 'Incomplete lead decision rows')
    require(panel['status'] == 'complete' and panel['case_count'] == 19, 'Incomplete raw panel')
    require(panel_receipt['raw_case_summary_sha256'] == args.expected_panel_sha256,
            'Panel receipt does not bind this raw panel')
    require(recent['case_count'] == 19 and recent['named_analogue_accepted_by_lead'] is True,
            'Recent-parent collection or lead declaration incomplete')
    require(recent_receipt['status'] == 'complete_verified_collection'
            and recent_receipt['panel_cases'] == 19, 'Recent-parent receipt incomplete')
    require(recent['actual_parameters_source_sha256'] == args.expected_panel_sha256,
            'Recent-parent observations refer to a different parameter panel')
    require(old_analysis['raw_panel_sha256'] == args.expected_panel_sha256,
            'Old derivative analysis belongs to another panel')
    cases = indexed(panel['cases'], 'case_id')
    observations = indexed(recent['rows'], 'case_id')
    planned = indexed(plan, 'case_id')
    expected_cases = {'baseline'} | {n + suffix for n in NAMES for suffix in ('_minus', '_plus')}
    require(set(cases) == set(observations) == set(planned) == expected_cases,
            'Exactly baseline plus both directions of nine coordinates required')
    require(len(parameter_reference) == 17, 'Complete 17-parameter reference required')
    reference_parameters = indexed(parameter_reference, 'parameter')

    for row in weights:
        rid = row['restriction_id']
        same(row['target'], old_fit[rid]['target'], 'Approved target changed: ' + rid)
        same(row['target'], decision_rows[rid]['target'], 'Lead/CSV target mismatch: ' + rid)
        require(row['role'] == decision_rows[rid]['role'], 'Role mismatch: ' + rid)
        if rid == 'initial_normalization':
            require(row['working_weight'] == row['working_scale'] == '', 'Normalization must be unscored')
            same(row['target'], 2.1, 'Normalization changed')
        else:
            scale, weight = number(row['working_scale']), number(row['working_weight'])
            require(scale > 0 and weight > 0, 'Every scored row needs a positive finite weight')
            same(weight, 1 / scale ** 2, 'Weight/scale inconsistency: ' + rid)
            same(weight, decision_rows[rid]['working_weight'], 'Lead/CSV weight mismatch: ' + rid)

    fit_rows, parameter_rows, case_rows = [], [], []
    values = {}
    case_order = [r['case_id'] for r in sorted(panel['cases'], key=lambda r: r['index'])]
    for case_id in case_order:
        case, obs = cases[case_id], observations[case_id]
        require(case['numerical_gates_verified'] is True and case['graph_count'] == 17,
                'Original numerical or diagnostic gates failed: ' + case_id)
        require(case['contract_sha256'] == planned[case_id]['contract_sha256'], 'Case contract mismatch')
        require(obs['checkpoint_sha256'] == case['checkpoint']['sha256'], 'Checkpoint mismatch')
        require(obs['actual_parameters_source_sha256'] == args.expected_panel_sha256, 'Parameter source mismatch')
        require(obs['diagnostic_name'] == RECENT_MOMENT, 'Wrong recent-parent observation')
        for name in NAMES:
            same(obs[name], case['parameters'][name], 'Recent-parent parameter mismatch: ' + name)
        packet = read(obs['observation_path'], obs['observation_sha256'])
        require(packet['available'] is True and packet['moment'] == RECENT_MOMENT,
                'Recent-parent packet missing/incorrect')
        require(packet['production_eligible'] is False and packet['target_contract_activated'] is False,
                'Original observation flags must remain preserved')
        same(packet['model_value'], obs['model_value'], 'Recent-parent collection value mismatch')
        require(packet['metadata']['policy_input_provenance']['checkpoint_sha256']
                == case['checkpoint']['sha256'], 'Packet checkpoint mismatch')
        rates = []
        for group_name in ('selected_birth', 'current_empty'):
            group = packet['groups'][group_name]
            denominator, numerator = number(group['denominator']), number(group['owner_numerator'])
            require(denominator > 0 and 0 <= numerator <= denominator, 'Invalid group masses')
            same(numerator / denominator, group['ownership_rate'], 'Group ratio mismatch')
            rates.append(number(group['ownership_rate']))
        same(rates[0] - rates[1], obs['model_value'], 'Recent-parent group contrast mismatch')
        directory = BASE / 'initial_sensitivity_panel/collected' / case_id / 'repetition_01'
        early = read(directory / 'early_measurement.json')
        require(early == case['early_measurement'], 'Collected early packet differs from raw panel')
        params = read(directory / 'parameters.csv', kind='csv')
        require(len(params) == 17 and set(indexed(params, 'parameter')) == set(reference_parameters),
                'Incomplete per-case parameters')
        for row in params:
            name = row['parameter']
            value = number(row['estimate'])
            same(value, case['parameters'][name], 'Raw/CSV parameter disagreement')
            reference = reference_parameters[name]
            for key in ('lower', 'upper', 'transform'):
                require(row[key] == reference[key], 'Parameter restriction changed: ' + name)
            near = row['near_bound'].lower() == 'true'
            require(row['near_bound'].lower() in ('true', 'false'), 'Invalid bound flag')
            if name in NAMES:
                lower, upper = number(row['lower']), number(row['upper'])
                require(lower <= value <= upper, 'Parameter outside bounds')
                require(near == (min(value - lower, upper - value) <= .01 * (upper - lower)),
                        'Incorrect original one-percent near-bound flag')
            parameter_rows.append(dict(case_id=case_id, **row,
                analysis_status='saved diagnostic input; not a selected calibrated estimate',
                checkpoint_sha256=case['checkpoint']['sha256']))
        normalized = number(case['normalization']['completed_fertility'])
        require(abs(normalized - 2.1) <= 5e-4, 'Separate normalization tolerance failed')
        losses, moment_values = [], []
        for weight_row in weights:
            rid = weight_row['restriction_id']
            path = (RECENT_MOMENT if rid == 'recent_parent_ownership' else old_fit[rid]['model_observation'])
            if rid == 'initial_normalization':
                model = normalized
            elif rid == 'recent_parent_ownership':
                model = number(obs['model_value'])
            else:
                require('constant_post_cell' not in path, 'Wrong primary projection mapping')
                model = lookup(early, path)
            target = number(weight_row['target'])
            gap = model - target
            scored = rid != 'initial_normalization'
            weight = number(weight_row['working_weight']) if scored else None
            scale = number(weight_row['working_scale']) if scored else None
            contribution = weight * gap ** 2 if scored else None
            if scored:
                losses.append((rid, contribution))
                moment_values.append(model)
            fit_rows.append(dict(case_id=case_id, **weight_row, model=model, gap=gap,
                scored=scored, standardized_gap=gap / scale if scored else None,
                loss_contribution=contribution, model_observation=path, cps_projection=PROJECTION,
                model_approximation='maintained synchronized dependent-residence proxy; not exact ACS'
                    if rid == 'recent_parent_ownership' else old_fit[rid].get('model_measurement_caveat', ''),
                checkpoint_sha256=case['checkpoint']['sha256']))
        require(len(losses) == 12, 'Missing scored restriction')
        total = math.fsum(loss for _, loss in losses)
        require(math.isfinite(total) and total >= 0, 'Invalid complete loss')
        dominant = sorted(losses, key=lambda item: item[1], reverse=True)
        summary = dict(case_id=case_id, complete_12_moment_loss=total,
            normalization=normalized, normalization_gap=normalized - 2.1,
            scored_moments=12, target_restrictions=13, parameter_rows=17,
            recent_parent_model=number(obs['model_value']),
            recent_parent_contribution=dict(losses)['recent_parent_ownership'],
            numerical_gates_verified=True, checkpoint_sha256=case['checkpoint']['sha256'],
            dominant_restriction=dominant[0][0], dominant_contribution=dominant[0][1],
            dominant_loss_share=dominant[0][1] / total if total else 0.,
            second_restriction=dominant[1][0], second_contribution=dominant[1][1],
            third_restriction=dominant[2][0], third_contribution=dominant[2][1],
            calibrated_estimate=False, inference_claim=False)
        case_rows.append(summary)
        values[case_id] = np.asarray(moment_values, dtype=float)

    scored_ids = list(IDS[1:])
    target = np.array([number(r['target']) for r in weights[1:]])
    sqrt_weight = np.sqrt([number(r['working_weight']) for r in weights[1:]])
    old_lookup = {(r['restriction'], r['parameter']): r for r in old_derivatives
                  if r['projection'] == PROJECTION}
    require(len(old_lookup) == 99, 'Original eleven-by-nine derivative reference incomplete')
    recent_derivative = indexed(recent['raw_derivatives'], 'parameter')
    J, Jminus, Jplus = [np.empty((12, 9), dtype=float) for _ in range(3)]
    derivatives, asymmetry = [], []
    for column, name in enumerate(NAMES):
        cm, cp, cb = [cases[name + suffix] for suffix in ('_minus', '_plus')] + [cases['baseline']]
        for changed_case in (cm, cp):
            require([n for n in NAMES if changed_case['parameters'][n] != cb['parameters'][n]] == [name],
                    'Direction changes more than its declared structural coordinate')
        um, up, ub = [coordinate(name, c['parameters'][name]) for c in (cm, cp, cb)]
        require(min(um, up) < ub < max(um, up), 'Actual log-coordinate changes do not bracket baseline')
        mm, mp, mb = [values[c['case_id']] for c in (cm, cp, cb)]
        central = (mp - mm) / (up - um)
        side_minus, side_plus = (mb - mm) / (ub - um), (mp - mb) / (up - ub)
        J[:, column], Jminus[:, column], Jplus[:, column] = central, side_minus, side_plus
        raw_asymmetry = float(np.linalg.norm(side_plus - side_minus) / max(np.linalg.norm(central), 1e-12))
        weighted_asymmetry = float(np.linalg.norm(sqrt_weight * (side_plus - side_minus))
                                   / max(np.linalg.norm(sqrt_weight * central), 1e-12))
        asymmetry.append(dict(parameter=name, raw_relative_column_asymmetry=raw_asymmetry,
            weighted_relative_column_asymmetry=weighted_asymmetry,
            weighted_column_norm=float(np.linalg.norm(sqrt_weight * central)),
            beta_case_labels_follow_beta_not_coordinate=name == 'beta_annual'))
        for index, rid in enumerate(scored_ids):
            if rid != 'recent_parent_ownership':
                old = old_lookup[(rid, name)]
                for actual, key in ((central[index], 'centered_derivative'),
                                    (side_minus[index], 'minus_derivative'), (side_plus[index], 'plus_derivative')):
                    same(actual, old[key], 'Old eleven-row derivative reproduction failed', 1e-10)
            else:
                old = recent_derivative[name]
                for actual, key in ((central[index], 'central_derivative'),
                                    (side_minus[index], 'minus_side_derivative'), (side_plus[index], 'plus_side_derivative')):
                    same(actual, old[key], 'Recent-parent derivative reproduction failed', 1e-10)
            derivatives.append(dict(restriction_id=rid, parameter=name, cps_projection=PROJECTION,
                coordinate='log(-log(beta_annual))' if name == 'beta_annual' else 'log(parameter)',
                parameter_minus=cm['parameters'][name], parameter_base=cb['parameters'][name],
                parameter_plus=cp['parameters'][name], coordinate_minus=um, coordinate_base=ub,
                coordinate_plus=up, actual_minus_coordinate_step=um - ub, actual_plus_coordinate_step=up - ub,
                model_minus=mm[index], model_base=mb[index], model_plus=mp[index], target=target[index],
                central_derivative=central[index], minus_side_derivative=side_minus[index],
                plus_side_derivative=side_plus[index], signed_slope_asymmetry=side_plus[index] - side_minus[index],
                relative_slope_asymmetry=abs(side_plus[index] - side_minus[index]) / max(abs(central[index]), 1e-12),
                working_weight=sqrt_weight[index] ** 2,
                weighted_central_derivative=sqrt_weight[index] * central[index],
                weighted_minus_side_derivative=sqrt_weight[index] * side_minus[index],
                weighted_plus_side_derivative=sqrt_weight[index] * side_plus[index]))
    A = sqrt_weight[:, None] * J
    require(np.isfinite(A).all(), 'Nonfinite weighted Jacobian')
    singular_values = np.linalg.svd(A, compute_uv=False)
    rank_tolerance = max(A.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.count_nonzero(singular_values > rank_tolerance))
    condition = float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0 else None
    require(len(fit_rows) == 247 and len(parameter_rows) == 323 and len(derivatives) == 108,
            'Incomplete output tables')
    # Recheck bytes immediately before writing; callers pin substantive inputs.
    for path, expected in inputs.items():
        require(digest(path) == expected, 'Input changed during analysis: ' + path)
    result = dict(schema='e5f_complete_panel_working_md_diagnostic_v1',
        status='complete_19_case_12_moment_numerical_diagnostic',
        objective_name=decision['name'], objective_definition=decision['objective'],
        working_weight_interpretation=decision['weight_interpretation'], cps_projection=PROJECTION,
        recent_parent_approximation=decision['recent_parent_approximation'],
        calibrated_estimate=False, calibrated_benchmark_promoted=False, production_policy_eligible=False,
        optimization_performed=False, inference_claim=False, candidates_selected=False,
        target_restrictions=13, scored_moments=12, structural_coordinates=9, case_count=19,
        normalization=dict(target=2.1, absolute_tolerance=5e-4, scored=False),
        row_order=scored_ids, column_order=list(NAMES), case_order=case_order,
        cases=case_rows, working_weight_rows=weights,
        raw_jacobian=J.tolist(), weighted_jacobian=A.tolist(),
        raw_minus_side_jacobian=Jminus.tolist(), raw_plus_side_jacobian=Jplus.tolist(),
        weighted_svd=dict(singular_values=singular_values.tolist(), numerical_rank=rank,
            condition_number=condition, rank_absolute_tolerance=float(rank_tolerance),
            rank_rule='max(matrix_shape) * float64 machine epsilon * largest singular value',
            coordinate_definition='log(-log(beta)) for beta; log actual parameter otherwise',
            interpretation='Local numerical conditioning only; no strong/global identification or inference claim'),
        column_asymmetry=asymmetry,
        asymmetry_definition='Euclidean norm of plus-minus one-sided slope difference divided by central-column norm; separately raw and sqrt(weight)-scaled; denominator floor 1e-12',
        derivative_definition='Two-sided secant in actual transformed coordinates, plus one-sided secants through baseline; beta minus/plus names follow beta level and reverse coordinate order',
        verification=dict(original_11_by_9_derivatives_reproduced=True,
            recent_parent_9_derivatives_reproduced=True, all_case_checkpoint_links_verified=True,
            all_collected_early_packets_equal_raw_panel=True, all_17_parameter_rows_verified=True,
            all_source_packets_preserved=True, all_13_targets_preserved=True,
            normalization_excluded_from_objective=True, source_files_unchanged_at_completion=True),
        externally_required_input_pins=dict(weights=args.expected_weights_sha256,
            lead_decision=args.expected_lead_decision_sha256, panel=args.expected_panel_sha256,
            recent_panel=args.expected_recent_panel_sha256),
        input_file_sha256=inputs, analysis_script_sha256=digest(__file__),
        outputs=dict(summary_csv='complete_panel_analysis.csv', fit_csv='complete_panel_fit.csv',
            parameters_csv='complete_panel_parameters.csv', derivatives_csv='complete_panel_derivatives.csv'))
    write_csv(HERE / 'complete_panel_analysis.csv', case_rows)
    write_csv(HERE / 'complete_panel_fit.csv', fit_rows)
    write_csv(HERE / 'complete_panel_parameters.csv', parameter_rows)
    write_csv(HERE / 'complete_panel_derivatives.csv', derivatives)
    (HERE / 'complete_panel_analysis.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(status=result['status'], baseline=case_rows[0], weighted_svd=result['weighted_svd'],
                         rows=dict(fit=len(fit_rows), parameters=len(parameter_rows), derivatives=len(derivatives))), indent=2))


if __name__ == '__main__':
    main()
