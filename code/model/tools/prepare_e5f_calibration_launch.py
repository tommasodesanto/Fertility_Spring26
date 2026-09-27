#!/usr/bin/env python3
"""Freeze the reviewed single-specification calibration contract on Torch.

This standard-library-only builder never imports the model or launches work.
Stage every source/runtime edit first. An existing output is never overwritten;
the controller separately requires a verified exact-loop smoke receipt.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

ROOT = Path('/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a')
WORK = ROOT / 'calibration_code_integration_20260927_v1'
COMMON = ('H0', 'beta_annual', 'chi', 'first_birth_fixed_cost', 'kappa_fert',
          'kappa_fert_continuation', 'theta0')
FREE = COMMON + ('delta_alpha_jump', 'child_benefit_curvature')
EARLY_VALUE = 0.8095276384290021
EARLY_MAPPING = 'fertility.uniform_birth_time.moments.mean_children_ever_born_capped3_age25'
ROOM_MAPPING = 'housing_wealth.moments.aggregate_mean_occupied_rooms_ahs_uncapped_18_85'
TOOL_NAMES = {
    'driver': 'run_e5f_utility_overnight_calibration.py',
    'calibration_runtime': 'e5f_calibration_runtime.py',
    'estate_audit': 'e5f_overnight_estate_audit.py',
    'recovery_runner': 'run_e5f_utility_comparison.py',
    'recovery_search': 'run_e5f_utility_comparison_search.py',
    'recovery_policy': 'e5f_utility_recovery_policy_v1.py',
}


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def pin(path):
    path = Path(path).resolve(strict=True)
    return {'path': str(path), 'sha256': sha(path)}


def checked(record):
    if pin(record['path'])['sha256'] != record['sha256']:
        raise ValueError('Inherited artifact changed: ' + record['path'])
    return read(record['path'])


def write(path, value):
    with Path(path).open('x') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def inventory(root):
    """Pin every staged file except Python interpreter caches; reject symlinks."""
    files = {}
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if '__pycache__' in relative.parts or path.suffix in ('.pyc', '.pyo'):
            continue
        if path.is_symlink():
            raise ValueError('Staged source/runtime must not contain symlinks: ' + str(path))
        if path.is_file():
            files[str(relative)] = sha(path)
    if not files:
        raise ValueError('Empty staged inventory: ' + str(root))
    return {'root': str(root), 'files': files, 'file_count': len(files),
            'exclusions': ['__pycache__ directories', '*.pyc', '*.pyo']}


def initial_point(path):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    values = {row['parameter']: float(row['estimate']) for row in rows}
    if len(values) != len(rows):
        raise ValueError('Duplicate diagnostic parameter')
    point = {name: values[name] for name in COMMON}
    point.update(delta_alpha_jump=.2, child_benefit_curvature=.14)
    if not all(math.isfinite(value) for value in point.values()):
        raise ValueError('Nonfinite initial parameter')
    return point


def new_provenance(old, early, ahs, early_path, ahs_path):
    provenance = copy.deepcopy(old)
    provenance.update(displayed_rows=14, positive_weight_rows=13,
                      free_parameter_counts={'single_share': 9},
                      status='reviewable proposed calibration contract; production blocked by author',
                      identification='13 scored moments for 9 free parameters; rank not certified')
    rooms = next(row for row in provenance['target_rows'] if row['id'] == 'mean_rooms')
    rooms['superseded_reference'] = copy.deepcopy(rooms)
    rooms.update(target=ahs['mean_rooms'], label='Mean occupied rooms, ages 18–85 (AHS 2007)',
        definition='Weighted mean literal ROOMS in the occupied national 2007 AHS sample; no cap at nine.',
        sample=ahs['sample'], model_observation=ROOM_MAPPING,
        source={'builder': 'code/data/ahs_supply_snapshot/build_ahs_2007_room_target.py',
                'path': str(ahs_path), 'record_id': 'mean_rooms',
                'contract_id': 'ahs_2007_literal_rooms_18_85',
                'raw_source_url': ahs['source_url'], 'raw_source_sha256': ahs['source_sha256']},
        standard_error=ahs['standard_error'], uncertainty_status=ahs['variance_method'],
        status='author-adopted AHS mean; inherited objective weight retained',
        weight_status='inherited working weight; not inverse AHS variance',
        mapping_warning='Uniform model exposure within four-year age cells; literal AHS ROOMS has public-use topcode 21. Model rooms are uncapped. Housing supply normalization remains an empirical normalization.')
    cell = early['estimates']['25']
    provenance['target_rows'].append({
        'id': 'early_fertility', 'label': 'Children ever born, capped at three, age 25',
        'target': cell['mean_children_ever_born_capped3'], 'weight': 100.0,
        'definition': 'Pooled supplement-weighted E[min(FREVER,3)] among women of completed age 25.',
        'sample': early['sample'], 'model_observation': EARLY_MAPPING,
        'source': {'builder': 'code/data/cps_fertility/build_early_fertility_target.py',
                   'path': str(early_path), 'record_id': 'estimates/25/mean_children_ever_born_capped3',
                   'contract_id': 'cps_2004_2006_exact_age25_capped3_20260926'},
        'sample_records': cell['n'], 'standard_error': cell['capped3_bootstrap_se'],
        'uncertainty_status': early['uncertainty'],
        'status': 'empirical target prepared for proposed first calibration',
        'weight_status': 'proposed working weight 100; author decision outstanding; not inverse bootstrap variance',
        'mapping_warning': 'Observed completed age 25 corresponds to [25,26). Fixed uniform birth-time interpolation in [22,26) uses post weight .875 and literal ever-born weights [0,1,2,3]. Household reproductive-member exposure is not certified female exposure. Bootstrap is not CPS design-consistent.',
        'model_projection': {'age_interval': [25, 26], 'model_cell': [22, 26],
                             'post_weight': .875, 'ever_born_weights': [0, 1, 2, 3],
                             'age_projection': 'uniform_birth_time'},
        'candidate_estimates': {'uncapped_mean': cell['mean_children_ever_born']},
    })
    return provenance


def build(args):
    work = args.work_root.resolve(strict=True)
    # Never run source hashing or scientific preparation on the author's Mac.
    if sys.platform != 'linux' or not work.is_relative_to(ROOT):
        raise RuntimeError('This preparation command is restricted to the Torch scratch project')
    output = (args.output or work / 'launch_v1').resolve()
    if not output.is_relative_to(work) or output.exists():
        raise ValueError('Output must be a new directory inside the staged work root')
    source, runtime_tools = work / 'source', work / 'tools'
    base = read(args.base_contract)
    old_objective = checked(base['parent_objective'])
    old_provenance = checked(base['files']['target_provenance'])
    early_path = args.early_target.resolve(strict=True)
    ahs_path = args.ahs_target.resolve(strict=True)
    early, ahs = read(early_path), read(ahs_path)
    if early['estimates']['25']['mean_children_ever_born_capped3'] != EARLY_VALUE:
        raise ValueError('Early fertility target differs from the reviewed empirical receipt')
    if abs(ahs['mean_rooms'] - 5.7294342401) > 1e-10:
        raise ValueError('AHS target differs from the reviewed empirical receipt')
    prior_rows = {r['restriction_id']: r for r in old_objective['target_rows']}
    if len(prior_rows) != 13 or len(old_provenance['target_rows']) != 13:
        raise ValueError('Expected exactly thirteen inherited display rows')
    for row in old_provenance['target_rows']:
        inherited = prior_rows[row['id']]
        if inherited['target'] != row['target'] or inherited['actual_weight'] != row['weight']:
            raise ValueError('Inherited target/provenance mismatch: ' + row['id'])
    provenance = new_provenance(old_provenance, early, ahs, early_path, ahs_path)
    by_id = {row['id']: row for row in provenance['target_rows']}
    objective = copy.deepcopy(old_objective)
    # The intact historical objective is copied separately; historical experiment
    # settings cannot masquerade as current launch settings in the derived file.
    for key in ('paired_contract', 'experimental_changes', 'source_provenance'):
        objective.pop(key, None)
    objective.update(contract_id='single_share_calibration_20260927_v1',
        objective_name='single_share_calibration_20260927',
        status='reviewable proposed minimum-distance objective; production blocked by author',
        calibrated_smm=False, benchmark_certified=False,
        identification_note=provenance['identification'], cps_projection='uniform_birth_time')
    restrictions = {r['parameter']: r for r in old_objective['parameter_restrictions']}
    objective['parameter_restrictions'] = [copy.deepcopy(restrictions[name]) for name in COMMON]
    objective['parameter_restrictions'].extend([
        {'parameter': 'delta_alpha_jump', 'lower': 0., 'upper': .25, 'transform': 'softzero'},
        {'parameter': 'child_benefit_curvature', 'lower': 0., 'upper': .8, 'transform': 'softzero'}])
    for name in ('mean_rooms', 'early_fertility'):
        p = by_id[name]
        # Replace all row provenance for the replaced measurement; never leave an
        # ACS authoritative record attached to the new AHS target.
        row = dict(restriction_id=name, label=p['label'], target=p['target'], actual_weight=p['weight'],
            model_observation=p['model_observation'], definition=p['definition'], sample=p['sample'],
            empirical_builder=p['source']['builder'], empirical_source_path=p['source']['path'],
            empirical_provenance_contract=p['source']['path'], empirical_provenance_contract_id=p['source']['contract_id'],
            empirical_record_id=p['source']['record_id'], empirical_standard_error=p['standard_error'],
            uncertainty_status=p['uncertainty_status'], warning=p['mapping_warning'],
            weight_status=p['weight_status'], role='scored_working_restriction', calibrated_smm=False)
        if name in prior_rows:
            index = next(i for i, old in enumerate(objective['target_rows']) if old['restriction_id'] == name)
            objective['target_rows'][index] = row
        else:
            objective['target_rows'].append(row)
    objective['mapping_verification'] = [r for r in objective.get('mapping_verification', [])
                                          if r.get('restriction_id') not in ('mean_rooms', 'early_fertility')]
    objective['mapping_verification'].extend(
        {'restriction_id': name, 'mapping': by_id[name]['model_observation'],
         'status': 'explicit registry mapping; exact-loop numerical validation required'}
        for name in ('mean_rooms', 'early_fertility'))
    objective.setdefault('weight_rationales', {}).update(
        mean_rooms='Retain inherited working weight while replacing the empirical mean and observer with literal AHS rooms; not inverse AHS variance.',
        early_fertility='Proposed working weight 100; author decision outstanding; not inverse bootstrap variance.')
    point = initial_point(args.initial_parameters)
    for row in objective['parameter_restrictions']:
        if not row['lower'] <= point[row['parameter']] <= row['upper']:
            raise ValueError('Initial parameter outside inherited bounds: ' + row['parameter'])
    if len(objective['target_rows']) != 14 or sum(r['actual_weight'] is not None for r in objective['target_rows']) != 13:
        raise ValueError('Wrong scored/display row count')
    source_inventory = inventory(source)
    tools_inventory = inventory(runtime_tools)
    required_source = ('code/model/intergen_eqscale_seq_optimized/solver.py',
                       'code/model/tools/e5f_initial_fertility_observer.py',
                       'code/model/tools/e5f_initial_housing_observer.py',
                       'code/model/tools/e5f_stationary_paygo.py')
    for relative in required_source:
        if relative not in source_inventory['files']:
            raise ValueError('Required staged source missing: ' + relative)
    files = {key: pin(runtime_tools / filename) for key, filename in TOOL_NAMES.items()}
    files.update({'runtime_file:' + relative: {'path': str(runtime_tools / relative), 'sha256': digest}
                  for relative, digest in tools_inventory['files'].items()})
    files['generator'] = pin(__file__)
    output.mkdir()
    reference = output / 'reference'
    reference.mkdir()
    for name, path in [('objective.json', base['parent_objective']['path']),
                       ('target_provenance.json', base['files']['target_provenance']['path']),
                       ('base_contract.json', args.base_contract), ('early_fertility_target.json', early_path),
                       ('ahs_2007_room_target.json', ahs_path), ('initial_parameters.csv', args.initial_parameters)]:
        destination = reference / name
        shutil.copyfile(path, destination)
        files['reference:' + name] = pin(destination)
    write(output / 'source_inventory.json', source_inventory)
    write(output / 'runtime_inventory.json', tools_inventory)
    write(output / 'target_provenance.json', provenance)
    objective['current_target_provenance'] = pin(output / 'target_provenance.json')
    objective['inherited_objective'] = pin(reference / 'objective.json')
    write(output / 'objective.json', objective)
    for name in ('source_inventory', 'runtime_inventory', 'target_provenance'):
        files[name] = pin(output / (name + '.json'))
    economic_changes = [
        {'object': 'utility', 'status': 'author-adopted for first calibration',
         'change': 'Same native compensated housing-share specification as matched diagnostic; first-child loading estimated, later-child loading fixed zero; no child housing floor.'},
        {'object': 'child_benefit_curvature', 'status': 'estimation proposed; bounds await author decision',
         'change': 'Estimate curvature in proposed bounds [0,.8]; starting value .14 matches the diagnostic receipt, which reports exponent .86 despite its historical folder suffix.'},
        {'object': 'mean_rooms', 'status': 'author-adopted', 'change': 'Replace ACS capped-nine mean with national AHS 2007 literal-room mean; retain old working weight.'},
        {'object': 'early_fertility', 'status': 'proposed objective addition; weight awaits author decision', 'change': 'Add exact-age-25 children-ever-born mean capped at three; proposed working weight 100.'},
        {'object': 'pensions', 'status': 'author-adopted', 'change': 'Fixed annual pension/income ratio from empirical receipt; payroll tax derived from stationary demographics.'},
        {'object': 'estate_entry_funding', 'status': 'provisional working closure', 'change': 'Positive net estates fund positive entrant assets; residual exits through a nonutility sink. Net-estate valuation and physical/financial settlement remain provisional.'},
        {'object': 'income_entry', 'status': 'unchanged inherited inputs', 'change': 'B15 earnings, inherited entrant-wealth marginal and income/wealth rank coupling retained; no zero-entry substitution.'},
    ]
    contract = dict(schema='e5f_single_share_overnight_v1', status='reviewed_smoke',
        approval={'authority': 'Author halted calibration launch pending instructions; preparation only',
                  'production_authorized': False,
                  'condition': 'Author instructions and approval, then external exact-loop smoke verification required before search'},
        production_blockers=['Author has halted launch pending instructions',
                             'Early-fertility weight 100 remains proposed',
                             'Curvature bounds [0,.8] remain proposed',
                             'External exact-loop smoke receipt and lead review are required'],
        source_root=str(source), source_manifest=pin(output / 'source_inventory.json'),
        runtime_tools=str(runtime_tools), files=files, base_contract=pin(args.base_contract),
        objective=pin(output / 'objective.json'), objective_canonical_sha256=canonical(objective),
        target_weight_fingerprint=canonical(objective['target_rows']),
        initial_point=point, seed=20260927,
        fixed={'sigma': 2., 'delta_alpha': 0., 'reference_rent': base['reference_rent'],
               'tenure_choice_kappa': .005, 'q_period': .08243216,
               'pension_ratio': base['pension_ratio'],
               'other_settings': 'Inherited selected checkpoint; no other default changes'},
        normalization={'initial_psi': .0645126953125, 'initial_step': .05,
                       'warm_price': True, 'maximum_stationary_solves': 23},
        budget={'total_seconds': 28800, 'search_seconds': 23400, 'repeat_seconds': 3600,
                'export_seconds': 1800, 'workers': 10, 'points_per_round': 10, 'rounds': 30,
                'objective_cap_seconds': 3100, 'smoke_seconds': 3700},
        proposal_widths=dict(zip(FREE, (.3, .001, .03, .03, .015, .03, .01, .025, .08))),
        proposal_width_status='Provisional lead search scales in raw parameter units; inherited seven bounds unchanged',
        thread_policy='No change to inherited thread counts; controller records environment',
        economic_changes=economic_changes,
        pending_observer_mismatches=[{'moment': row['id'], 'warning': row['mapping_warning']}
                                    for row in provenance['target_rows'] if row.get('mapping_warning')],
        deferred_review_blocks=['tenure choice scale .005', 'interest rate 2%, credit, DUE and payment-to-income restrictions',
            'rental menu and cap', 'income and entrant-wealth compatibility',
            'estate valuation, SCF child-directed measurement and older-family-income proxy',
            'weights, identification and bounds',
            'numerical/economic tests: finer owner housing grid and conception schedule'],
        identification={'free_parameters': 9, 'displayed_rows': 14, 'scored_rows': 13,
                        'normalizations': 1, 'local_rank_certified': False},
        scope='Stationary first calibration only; no policy transition or grid/conception refinement')
    write(output / 'contract.json', contract)
    return {'contract': pin(output / 'contract.json'), 'source_file_count': source_inventory['file_count'],
            'runtime_file_count': tools_inventory['file_count'], 'displayed_rows': 14, 'scored_rows': 13,
            'free_parameters': 9, 'launch_performed': False,
            'required_next_step': 'Author instructions; this preparation does not authorize a calibration launch'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-root', type=Path, default=WORK)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--base-contract', type=Path, default=ROOT / 'utility_four_arm_preparation_20260925_v2/launch_v1/contract.json')
    parser.add_argument('--early-target', type=Path, default=ROOT / 'early_fertility_target_20260926/output/early_fertility_target.json')
    parser.add_argument('--ahs-target', type=Path, required=True)
    parser.add_argument('--initial-parameters', type=Path, default=ROOT / 'housing_fertility_cost_diagnostic_20260926/combined_001/matched_shares_020_1.0/parameters.csv')
    print(json.dumps(build(parser.parse_args()), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
