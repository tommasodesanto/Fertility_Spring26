"""Compare completed experiment and exact canonical baseline; never solves."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[4]
OUTPUT = PROJECT / 'output/model/experiments/birth_count_choice/current_params_v1'
BASELINE = PROJECT / 'output/model/local_solution/cases/20261003T175652812716Z_b1c72f13'
SNAPSHOT = PROJECT / 'code/model/production/reference_inputs/bundle.json'


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def read_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def read_case(path):
    case = path.resolve(strict=True)
    metadata = json.loads((case / 'metadata.json').read_text())
    if metadata.get('status') != 'complete' or metadata.get('closure', {}).get('status') != 'converged':
        raise RuntimeError('Refusing incomplete or unconverged case: ' + str(case))
    contract = json.loads((case / 'input_contract.json').read_text())
    fits, parameters = read_rows(case / 'target_fit.csv'), read_rows(case / 'parameters.csv')
    if len(fits) != 14 or len(parameters) != 31:
        raise RuntimeError('Completed case lacks all 14 fit and 31 parameter rows')
    if len(list((case / 'standard_diagnostics').glob('*.png'))) != 17:
        raise RuntimeError('Completed case lacks the standard 17 diagnostics')
    if not (case / 'native_result.npz').is_file():
        raise RuntimeError('Completed case lacks its native result archive')
    target = [{key: row[key] for key in ('moment', 'target', 'weight', 'role')} for row in fits]
    target_pin, weight_pin = canonical(target), canonical(dict(base_contract=target, multipliers={}))
    pins = json.loads(SNAPSHOT.read_text())
    provenance = contract.get('parameter_file', {}).get('provenance')
    metadata['_target_provenance_check'] = ('saved parameter-file provenance plus recomputed complete tables' if provenance is not None
                                              else 'baseline predates parameter files; recomputed complete saved tables against immutable canonical snapshot')
    for key, actual in (('target_fingerprint', target_pin), ('weight_fingerprint', weight_pin)):
        if actual != pins[key] or (provenance is not None and actual != provenance[key]):
            raise RuntimeError('Complete target/weight fingerprint differs: ' + key)
    return case, metadata, fits, parameters, target_pin, weight_pin


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def display(value):
    if value == '' or value is None: return '—'
    if isinstance(value, bool): return str(value)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value).replace('|', '\\|')
    return f'{number:.3e}' if number and abs(number) < .001 else f'{number:.3f}'


def markdown_table(headers, rows):
    return ['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join('---' for _ in headers) + ' |'] + [
        '| ' + ' | '.join(display(value) for value in row) + ' |' for row in rows]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(argv)
    baseline, old, old_fit, old_parameters, target_pin, weight_pin = read_case(BASELINE)
    case, new, new_fit, new_parameters, new_target, new_weight = read_case(OUTPUT / 'latest')
    if (target_pin, weight_pin) != (new_target, new_weight):
        raise RuntimeError('Cannot compare differing target/weight contracts')
    if old['parameters'] != new['parameters'] or len(new['parameters']) != 10:
        raise RuntimeError('Experiment changed the ten supplied parameter coordinates')
    if [row['moment'] for row in old_fit] != [row['moment'] for row in new_fit]:
        raise RuntimeError('Target order differs')
    comparisons = []
    for a, b in zip(old_fit, new_fit):
        if any(a[key] != b[key] for key in ('moment','target','weight','role')):
            raise RuntimeError('Empirical target contract differs')
        comparisons.append(dict(moment=a['moment'], target=a['target'], baseline_model=a['model'],
            experiment_model=b['model'], baseline_gap=a['gap'], experiment_gap=b['gap'],
            weight=a['weight'], baseline_loss_contribution=a['loss_contribution'],
            experiment_loss_contribution=b['loss_contribution'], role=a['role']))
    old_by_name = {row['parameter']: row for row in old_parameters}
    if set(old_by_name) != {row['parameter'] for row in new_parameters}:
        raise RuntimeError('Parameter table identities differ')
    parameter_comparisons = []
    for row in new_parameters:
        earlier = old_by_name[row['parameter']]
        if any(earlier[key] != row[key] for key in ('lower','upper')):
            raise RuntimeError('Parameter restrictions differ: ' + row['parameter'])
        parameter_comparisons.append(dict(parameter=row['parameter'], baseline_estimate=earlier['estimate'],
            experiment_estimate=row['estimate'], lower=row['lower'], upper=row['upper'],
            baseline_near_bound=earlier['near_bound'], experiment_near_bound=row['near_bound'],
            baseline_status=earlier['status'], experiment_status=row['status']))
    old_loss = sum(float(row['loss_contribution']) for row in old_fit if row['role']=='scored')
    new_loss = sum(float(row['loss_contribution']) for row in new_fit if row['role']=='scored')
    lines = ['# Joint birth-count experiment: complete fit comparison', '',
        f'Baseline loss: **{old_loss:.6f}**. Experimental loss: **{new_loss:.6f}**.',
        'Both cases passed the inherited stationary-GE acceptance and exact-repeat gates.',
        'The ten supplied parameter coordinates are unchanged; this is a specification experiment without recalibration.', '',
        f'Baseline: `{baseline}`', f'Experiment: `{case}`',
        f'Target fingerprint: `{target_pin}`', f'Weight fingerprint: `{weight_pin}`',
        f'Baseline provenance check: {old["_target_provenance_check"]}.',
        f'Experiment provenance check: {new["_target_provenance_check"]}.', '',
        'The economic change is the joint intended-count menu with existing-age Binomial success probabilities. '
        'Birth-order flows count children; recent-parent ownership weights successful households once. '
        'Existing target definitions, within-period interpolation, entry, post-interest timing and soft credit are retained.', '',
        '## All 14 target rows', '',
        'Gap means model minus target. Blank normalization weights and zero-weight validation rows remain in the table. '
        'CSV values retain full precision.', '']
    lines += markdown_table(['Moment','Target','Baseline','Experiment','Baseline gap','Experiment gap','Weight','Baseline loss','Experiment loss'],
        [[r[key] for key in ('moment','target','baseline_model','experiment_model','baseline_gap','experiment_gap','weight',
                             'baseline_loss_contribution','experiment_loss_contribution')] for r in comparisons])
    lines += ['', '## All 31 parameter records', '',
        'Restrictions are the inherited reference bounds; they are advisory for this fixed-input GE. '
        'Near-bound indicators use the inherited one-percent-of-range screen.', '']
    lines += markdown_table(['Parameter','Baseline','Experiment','Lower','Upper','Near bound: baseline','Near bound: experiment'],
        [[r[key] for key in ('parameter','baseline_estimate','experiment_estimate','lower','upper',
                             'baseline_near_bound','experiment_near_bound')] for r in parameter_comparisons])
    lines += ['', '## Convergence and accounting', '']
    convergence = []
    for key in ('status','closure_mode','renewal_residual','absolute_housing_residual','actual_paygo_residual',
                'population_scale','fixed_h0_population_scale','implied_H0_at_population_one',
                'adjusted_births_per_normalized_household','actual_entry_per_normalized_household',
                'normalized_housing_demand','physical_housing_supply','standard_plot_count','target_fit_rows','parameter_rows'):
        convergence.append([key, old['closure'].get(key), new['closure'].get(key)])
    convergence.insert(1,['price', old['price'], new['price']])
    for key in ('renewal_adjusted_distribution_l1','mass_residual','entry_gap'):
        convergence.append(['native_population_step.'+key, old['closure'].get('native_population_step',{}).get(key),
                            new['closure'].get('native_population_step',{}).get(key)])
    lines += markdown_table(['Object','Baseline','Experiment'], convergence)
    lines += ['', '## Inspection and interpretation', '',
        f'The standard packet is at `{case / "standard_diagnostics"}`; cached policy and aggregate plots are in that case.',
        'The descriptive arrays attempt_hazard_by_age, first_birth_hazard_by_age, fert_by_age and first_birth_age_distribution use pre-birth exposure rather than the older post-birth exposure. '
        'The active target observers are separate and retain their definitions.',
        'The all-zero dead-menu identity for less than 1e-12 occupied household mass preserves the inherited numerical '
        'convention; it changes no economic parameter or target. The triggering failed witness was 2.97e-38 household mass.',
        'This first converged experimental steady state does not establish an accepted calibration, identification, '
        'grid adequacy, transition readiness or paper adoption.', '',
        '[Full fit CSV](comparison.csv) · [Full parameter CSV](parameters_comparison.csv)', '']
    write_csv(OUTPUT/'comparison.csv', comparisons)
    write_csv(OUTPUT/'parameters_comparison.csv', parameter_comparisons)
    (OUTPUT/'RESULTS.md').write_text('\n'.join(lines))
    print(OUTPUT/'RESULTS.md')


if __name__ == '__main__': main()
