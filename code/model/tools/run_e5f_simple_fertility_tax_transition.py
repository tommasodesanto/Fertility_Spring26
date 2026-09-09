#!/usr/bin/env python3
"""Bounded closed-household rebated-tax diagnostic, not a population forecast.

Each date reuses the checked impact solve, then advances the maintained queue.
No calibrated parameter, supply anchor, numerical gate or model kernel changes.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import math
import os
from pathlib import Path
import pickle
import threading
import time
import traceback

import numpy as np
import run_e5f_simple_fertility_tax_impact as impact

common, adapter, audit = impact.common, impact.adapter, impact.audit
baseline, transition = common.policy.baseline, common.policy.transition
CASES = impact.CASES
SCHEMA = 'e5f_simple_fertility_tax_transition_v1'
YEARS = list(range(2023, 2064, 4))
SCOPE = ('Temporary-equilibrium household-unit diagnostic; closed M=0, rho=1, '
         'inherited birth queue and births/2.1; not production or a resident-population forecast')


def verify_sources(c):
    for source in (__file__, baseline.__file__, transition.__file__, impact.calendar.__file__):
        if str(Path(source).resolve()) not in c['source_hashes']:
            raise RuntimeError(f'Transition helper must be source-pinned: {source}')
    impact.verify_sources(c)


def load_contract(path, expected, *, launching=True):
    adapter.verify(path, expected)
    c = adapter.read_json(path)
    closure = dict(outside_flow=0.0, retention=1.0,
                   birth_to_entry_conversion=1.0/2.1, next_bridge_year=None)
    if (c['schema'] != SCHEMA or c['years'] != YEARS or c['closure'] != closure
            or c['date_seconds'] != 1800 or c['stage_seconds'] != 14400):
        raise RuntimeError('Unexpected transition specification or budget')
    if launching and time.time() > c['launch_deadline_epoch']:
        raise RuntimeError('Transition launch deadline expired')
    for key in ('selected_summary', 'output_root'):
        if not Path(c[key]).is_absolute():
            raise RuntimeError(f'Nonabsolute contract path: {key}')
    if set(c['impact_receipts']) != set(CASES):
        raise RuntimeError('Both verified tax impact receipts are required')
    verify_sources(c)
    if transition.EFFECTIVE_BIRTH_TO_HOUSEHOLD_CONVERSION != 1.0/2.1:
        raise RuntimeError('Birth-to-entry conversion changed')
    lag = impact.tax.renewal_lag_gate()
    if lag['birth_to_population_effect_years'] != 20:
        raise RuntimeError('Expected twenty-year population lag')
    return c


def impact_receipt(c, case):
    spec = c['impact_receipts'][case]
    path = Path(spec['path'])
    if not path.is_absolute():
        raise RuntimeError('Impact receipt path must be absolute')
    r = impact.verify_receipt(path, c, expected_hash=spec['sha256'])
    if r['case'] != case:
        raise RuntimeError('Impact case mismatch')
    return path, r


def load_packet(path):
    with gzip.open(path, 'rb') as stream:
        return pickle.load(stream)


def save_packet(path, value):
    with gzip.open(path, 'wb', compresslevel=1) as stream:
        pickle.dump(value, stream, protocol=5)


def numeric_comparison(left, right, label):
    """Compare the complete saved quantity dictionary with an absolute gate."""
    if set(left) != set(right):
        raise RuntimeError(f'{label}: quantity keys differ')
    gaps = {}
    for key in left:
        a, b = left[key], right[key]
        if a is None or b is None:
            if a is not None or b is not None:
                raise RuntimeError(f'{label}: missing quantity {key}')
            continue
        gap = abs(float(a)-float(b))
        if not math.isfinite(gap) or gap > 2e-10:
            raise RuntimeError(f'{label}: {key} differs by {gap}')
        gaps[key] = gap
    return gaps


def compare_policy_tree(left, right, label='policy'):
    """Check every stored policy/map/joint field, including auxiliary arrays."""
    if isinstance(left, np.ndarray):
        if not isinstance(right, np.ndarray) or left.shape != right.shape or not np.allclose(
                left, right, rtol=0, atol=2e-10, equal_nan=True):
            raise RuntimeError(f'{label}: stored array differs')
    elif hasattr(left, '__dict__'):
        if not hasattr(right, '__dict__'):
            raise RuntimeError(f'{label}: object representation differs')
        compare_policy_tree(vars(left), vars(right), label)
    elif isinstance(left, dict):
        if not isinstance(right, dict) or set(left) != set(right):
            raise RuntimeError(f'{label}: field set differs')
        for key in left:
            compare_policy_tree(left[key], right[key], f'{label}.{key}')
    elif isinstance(left, (list, tuple)):
        if type(left) != type(right) or len(left) != len(right):
            raise RuntimeError(f'{label}: sequence differs')
        for index, (a, b) in enumerate(zip(left, right)):
            compare_policy_tree(a, b, f'{label}[{index}]')
    elif isinstance(left, (float, np.floating)):
        if not np.isclose(left, right, rtol=0, atol=2e-10, equal_nan=True):
            raise RuntimeError(f'{label}: scalar differs')
    elif left != right:
        raise RuntimeError(f'{label}: value differs')


def compare_date(folder, reference):
    current = load_packet(folder/'dated_state.pkl.gz')
    original = load_packet(reference/'dated_state.pkl.gz')
    gaps = common.baseline_comparison(current['evaluation'], original['evaluation'], policies=True)
    compare_policy_tree(current['evaluation'].policy, original['evaluation'].policy)
    gaps['complete_policy_tree'] = True
    gaps['quantities'] = numeric_comparison(
        adapter.read_json(folder/'quantities.json'),
        adapter.read_json(reference/'quantities.json'), 'date replay')
    gaps['fiscal_ledger'] = numeric_comparison(
        adapter.read_json(folder/'fiscal_ledger.json'),
        adapter.read_json(reference/'fiscal_ledger.json'), 'fiscal replay')
    for filename, labels in (('family_group_quantities.csv', ('age_group', 'family_group')),
                             ('lifecycle.csv', ())):
        left, right = adapter.read_csv(folder/filename), adapter.read_csv(reference/filename)
        if len(left) != len(right):
            raise RuntimeError(f'{filename}: row count differs')
        for a, b in zip(left, right):
            if any(a[key] != b[key] for key in labels):
                raise RuntimeError(f'{filename}: group labels differ')
            numeric_comparison({k:(None if v == '' else v) for k,v in a.items() if k not in labels},
                               {k:(None if v == '' else v) for k,v in b.items() if k not in labels}, filename)
    gaps['group_and_lifecycle_quantities'] = True
    return gaps


def verify_stage(root, stage, case, c, contract_sha):
    folder = root/stage/case
    path = folder/'receipt.json'
    expected = (folder/'receipt.sha256').read_text().strip()
    adapter.verify(path, expected)
    receipt = impact.verify_receipt(path, c, same_contract=contract_sha)
    if (receipt['stage'] != stage or receipt['case'] != case or
            receipt['years'] != (YEARS[:2] if stage == 'smoke' else YEARS)):
        raise RuntimeError('Wrong transition stage, case or calendar')
    return receipt, expected


def verify_advance(state, next_state, row):
    for field, births_key in (('scheduled_entries', 'birth_children_topcode_adjusted'),
                              ('scheduled_raw_entries', 'birth_children')):
        old, new = getattr(state, field), getattr(next_state, field)
        if len(old) != 4 or len(new) != 4 or not np.isfinite(old+new).all() or min(old+new) < 0:
            raise RuntimeError(f'Invalid four-slot queue: {field}')
        expected = list(old[1:]) + [float(row[births_key])/2.1]
        if not np.allclose(new, expected, rtol=0, atol=2e-10):
            raise RuntimeError(f'Queue advance identity failed: {field}')
    if not math.isclose(row['entrant_flow_next'], state.scheduled_entries[0], rel_tol=0, abs_tol=2e-10):
        raise RuntimeError('Next entrants differ from inherited mature queue slot')


def execute(args):
    c = load_contract(args.contract, args.contract_sha256)
    root = Path(c['output_root'])
    prerequisites = {}
    for case in CASES:
        impact_receipt(c, case)
        if args.stage == 'full':
            _, digest = verify_stage(root, 'smoke', case, c, args.contract_sha256)
            prerequisites[case] = digest
    folder = root/args.stage/args.case
    folder.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    progress = dict(case=args.case, stage=args.stage, phase='loading_selected', completed_dates=0)
    clock = dict(date_started=None)
    stop = threading.Event()

    def heartbeat():
        while not stop.is_set():
            now = time.monotonic()
            elapsed = now-started
            date_elapsed = 0 if clock['date_started'] is None else now-clock['date_started']
            audit.save_json(folder/'heartbeat.json', dict(progress, epoch=time.time(),
                            elapsed_seconds=elapsed, date_elapsed_seconds=date_elapsed))
            if elapsed >= c['stage_seconds'] or date_elapsed >= c['date_seconds']:
                audit.save_json(folder/'failure.json', dict(progress, error='Declared time budget exceeded'))
                os._exit(124)
            stop.wait(30)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        audit.save_json(folder/'latest_completed_case.json', dict(status='none_completed'))
        audit.save_json(folder/'best_so_far.json', dict(scope='fixed selected calibration; no search',
                       loss=common.SELECTED_LOSS, selected_summary_sha256=c['selected_summary_sha256']))
        packet, prepared, state = common.selected_packet(c)
        mass_2023 = float(state.g_pre.sum())
        initial = dict(distribution_sha256=baseline.array_sha256(state.g_pre),
                       queue=list(state.scheduled_entries), raw_queue=list(state.scheduled_raw_entries))
        reference_path, _ = impact_receipt(c, args.case)
        years = YEARS[:2] if args.stage == 'smoke' else YEARS
        rows, comparisons = [], {}
        for index, year in enumerate(years):
            clock['date_started'] = time.monotonic()
            progress.update(calendar_year=year, phase='solving_date')
            date_dir = folder/f'date_{year}'
            date_dir.mkdir()
            common.assert_model(packet['parameters'], prepared.supply_rule)
            result = impact.solve_endpoint(packet, prepared, state, args.case, date_dir, progress)
            current = load_packet(date_dir/'dated_state.pkl.gz')
            current.update(calendar_year=year, scope=SCOPE)
            save_packet(date_dir/'dated_state.pkl.gz', current)
            result.update(calendar_year=year, scope=SCOPE,
                          checkpoint_sha256=adapter.digest(date_dir/'dated_state.pkl.gz'))
            audit.save_json(date_dir/'endpoint_receipt.json', result)
            common.assert_model(current['parameters'], prepared.supply_rule)
            checked = {}
            if year == 2023:
                checked['verified_impact'] = compare_date(date_dir, reference_path.parent)
            if args.stage == 'full' and index < 2:
                checked['smoke'] = compare_date(date_dir, root/'smoke'/args.case/f'date_{year}')
            row, next_state = baseline.advance_from_evaluation(
                label=args.case, period_from_2007=baseline.TRANSITION_PERIODS+index,
                evaluation=current['evaluation'], state=state, P=current['parameters'],
                b_grid=prepared.b_grid, shared=current['shared'], supply_rule=prepared.supply_rule,
                outside_flow=0., retention=1., initial_mass_2007=prepared.initial_mass_2007,
                mass_2023=mass_2023, next_bridge_year=None, grid_fallback=False)
            if not math.isfinite(row['mass_accounting_residual']) or abs(row['mass_accounting_residual']) > 2e-10:
                raise RuntimeError('Transition mass accounting exceeds 2e-10')
            if row['calendar_year'] != year:
                raise RuntimeError('Transition calendar mismatch')
            verify_advance(state, next_state, row)
            next_state.initial_policy = None
            after = dict(calendar_year=year+4, state=next_state, scope=SCOPE)
            save_packet(date_dir/'next_state.pkl.gz', after)
            row.update(result['quantities'])
            row.update(result['ledger'])
            row.update(policy_case=args.case, household_mass=row['adult_population'],
                       total_adjusted_births=row['birth_children_topcode_adjusted'],
                       births_per_household=row['topcode_adjusted_births_per_adult'],
                       annual_property_tax_rate=impact.tax.CASES[args.case].annual_tax_rate,
                       scope=SCOPE)
            audit.save_json(date_dir/'advance.json', dict(
                calendar_year=year, next_calendar_year=year+4, row=row,
                inherited_queue=list(state.scheduled_entries),
                inherited_raw_queue=list(state.scheduled_raw_entries),
                next_queue=list(next_state.scheduled_entries),
                next_raw_queue=list(next_state.scheduled_raw_entries),
                next_household_mass=float(next_state.g_pre.sum()),
                next_distribution_sha256=baseline.array_sha256(next_state.g_pre)))
            if args.stage == 'full' and index < 2:
                reference_next = load_packet(root/'smoke'/args.case/f'date_{year}'/'next_state.pkl.gz')['state']
                gap = float(np.abs(next_state.g_pre-reference_next.g_pre).sum())
                if gap > 2e-10 or not math.isfinite(gap):
                    raise RuntimeError('Smoke next-state distribution mismatch')
                checked['next_state_l1'] = gap
                numeric_comparison(dict(enumerate(next_state.scheduled_entries)),
                                   dict(enumerate(reference_next.scheduled_entries)), 'smoke queue')
                numeric_comparison(dict(enumerate(next_state.scheduled_raw_entries)),
                                   dict(enumerate(reference_next.scheduled_raw_entries)), 'smoke raw queue')
            comparisons[str(year)] = checked
            audit.save_json(date_dir/'reproduction.json', checked)
            rows.append(row)
            baseline.write_csv(folder/'path.csv', rows)
            audit.save_json(folder/'latest_completed_case.json', dict(status='running', latest=row, completed_dates=len(rows)))
            state = next_state
            progress.update(completed_dates=len(rows), phase='date_complete')
            clock['date_started'] = None
        audit.save_json(folder/'latest_completed_case.json', dict(status='complete', latest=rows[-1],
                       completed_dates=len(rows), scope=SCOPE))
        verify_sources(c)
        files = sorted(p for p in folder.rglob('*') if p.is_file() and p.name not in ('receipt.json', 'receipt.sha256', 'heartbeat.json'))
        receipt = dict(status='complete', stage=args.stage, case=args.case, years=years,
                       contract_sha256=args.contract_sha256,
                       selected_summary_sha256=c['selected_summary_sha256'],
                       selected_checkpoint_sha256=c['checkpoint_sha256'],
                       scientific_bundle_sha256=c['code_bundle_sha256'],
                       impact_receipts=c['impact_receipts'], smoke_receipt_sha256=prerequisites,
                       closure=c['closure'], common_initial_state=initial, reproduction=comparisons,
                       elapsed_seconds=time.monotonic()-started, scope=SCOPE, production_promoted=False,
                       artifact_sha256={str(p.relative_to(folder)):adapter.digest(p) for p in files})
        audit.save_json(folder/'receipt.json', receipt)
        (folder/'receipt.sha256').write_text(adapter.digest(folder/'receipt.json')+'\n')
        progress['phase'] = 'complete'
    except BaseException as error:
        progress['phase'] = 'failed'
        audit.save_json(folder/'failure.json', dict(error=str(error), type=type(error).__name__, traceback=traceback.format_exc()))
        raise
    finally:
        stop.set()
        thread.join(timeout=2)
        audit.save_json(folder/'heartbeat.json', dict(progress, epoch=time.time(), elapsed_seconds=time.monotonic()-started))


def collect(args):
    c = load_contract(args.contract, args.contract_sha256, launching=False)
    root = Path(c['output_root'])
    receipt_hashes, paths = {}, {}
    initial_states = []
    for case in CASES:
        impact_receipt(c, case)
        smoke, smoke_hash = verify_stage(root, 'smoke', case, c, args.contract_sha256)
        full, full_hash = verify_stage(root, 'full', case, c, args.contract_sha256)
        initial_states.extend((smoke['common_initial_state'], full['common_initial_state']))
        receipt_hashes[case] = dict(smoke=smoke_hash, full=full_hash)
        if full['smoke_receipt_sha256'][case] != smoke_hash:
            raise RuntimeError('Full stage used a different smoke receipt')
        paths[case] = adapter.read_csv(root/'full'/case/'path.csv')
        if len(paths[case]) != len(YEARS):
            raise RuntimeError('Collected path has missing dates')
        for year in YEARS[:2]:
            compare_date(root/'full'/case/f'date_{year}', root/'smoke'/case/f'date_{year}')
        reference, _ = impact_receipt(c, case)
        compare_date(root/'full'/case/'date_2023', reference.parent)
    for case in CASES:
        full, _ = verify_stage(root, 'full', case, c, args.contract_sha256)
        if full['smoke_receipt_sha256'] != {k:v['smoke'] for k,v in receipt_hashes.items()}:
            raise RuntimeError('Full run did not pin both current smoke receipts')
    initial0 = adapter.read_json(root/'full'/CASES[0]/'receipt.json')['common_initial_state']
    initial1 = adapter.read_json(root/'full'/CASES[1]/'receipt.json')['common_initial_state']
    if initial0 != initial1:
        raise RuntimeError('Tax paths did not inherit the same population and queues')
    if any(value != initial_states[0] for value in initial_states):
        raise RuntimeError('Smoke/full branches differ in inherited population or queues')
    effects = []
    for year, base, reform in zip(YEARS, paths[CASES[0]], paths[CASES[1]]):
        if int(base['calendar_year']) != year or int(reform['calendar_year']) != year:
            raise RuntimeError('Collected calendars differ')
        mass_gap = abs(float(base['household_mass'])-float(reform['household_mass']))
        if year < 2043 and mass_gap > 2e-10:
            raise RuntimeError(f'Household mass responds before birth-entry lag at {year}: {mass_gap}')
        row = dict(calendar_year=year, household_mass_absolute_difference=mass_gap)
        for key in ('births_per_household', 'total_adjusted_births', 'household_mass',
                    'asset_price', 'young_mean_rooms_whole_nodes', 'housing_demand_per_adult'):
            denominator = float(base[key])
            if denominator == 0:
                raise RuntimeError(f'Zero comparison denominator: {key}')
            row[key+'_pct_change'] = 100*(float(reform[key])/denominator-1)
        for key in ('owner_rate', 'young_ownership_whole_nodes'):
            row[key+'_pp_change'] = 100*(float(reform[key])-float(base[key]))
        families = []
        for case in CASES:
            groups = adapter.read_csv(root/'full'/case/f'date_{year}'/'family_group_quantities.csv')
            families.append(next(x for x in groups if x['age_group']=='young_whole_nodes_25_34'
                                 and x['family_group']=='dependent_children'))
        row['young_dependent_owner_rate_pp_change'] = 100*(float(families[1]['owner_rate'])-float(families[0]['owner_rate']))
        row['young_dependent_rooms_pct_change'] = 100*(float(families[1]['mean_rooms'])/float(families[0]['mean_rooms'])-1)
        effects.append(row)
    baseline.write_csv(root/'comparison.csv', effects)
    lines = ['# Rebated property-tax transition diagnostic', '', SCOPE+'.', '',
             'Annual 2% versus 1%, both rebated equally to household decision units. '
             'The same calibrated supply rule and inherited 2023 state are retained.', '',
             '| Year | Births/HH % | Total births % | Households % | Price % | Young parent ownership pp | Young parent rooms % |',
             '|---|---:|---:|---:|---:|---:|---:|']
    for row in effects:
        lines.append('| '+str(row['calendar_year'])+' | '+' | '.join(f'{row[key]:.4f}' for key in (
            'births_per_household_pct_change', 'total_adjusted_births_pct_change',
            'household_mass_pct_change', 'asset_price_pct_change',
            'young_dependent_owner_rate_pp_change', 'young_dependent_rooms_pct_change'))+' |')
    lines.extend(['', 'Births are top-code adjusted; each path also records raw explicit births. '
                  'Young parents are household units in the unchanged whole-node 25–34 group with dependent children. '
                  'Group means can change through composition. Household mass must remain identical before 2043 under the maintained entry lag. '
                  'The historical-to-future household-entry handoff and age alignment remain limitations.', '',
                  'Every date has a fiscal ledger, full policy checkpoint, lifecycle and family quantities, '
                  'entry/queue accounting, and numerical gates. No new calibration or welfare claim.'])
    (root/'READOUT.md').write_text('\n'.join(lines)+'\n')
    audit.save_json(root/'comparison_receipt.json', dict(
        status='complete', contract_sha256=args.contract_sha256, receipt_sha256=receipt_hashes,
        years=YEARS, pre2043_household_mass_gate=2e-10, scope=SCOPE, production_promoted=False,
        comparison_sha256=adapter.digest(root/'comparison.csv'), readout_sha256=adapter.digest(root/'READOUT.md')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--case', choices=CASES)
    parser.add_argument('--stage', choices=('smoke', 'full'))
    parser.add_argument('--collect', action='store_true')
    args = parser.parse_args()
    if args.collect:
        if args.case or args.stage:
            parser.error('--collect cannot be combined with --case or --stage')
        collect(args)
    else:
        if not args.case or not args.stage:
            parser.error('--case and --stage are required for a run')
        execute(args)
