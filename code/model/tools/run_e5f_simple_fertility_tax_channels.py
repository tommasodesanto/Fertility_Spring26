#!/usr/bin/env python3
"""Eight fixed-component tax/price/rebate cells and their Shapley decomposition.

Mixed cells are conditional household responses, not market-clearing or
balanced-budget policies. All cells inherit the same selected 2023 population.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import itertools
import math
import os
from pathlib import Path
import pickle
import threading
import time
import traceback

import numpy as np
import run_e5f_simple_fertility_policy_impact as common
import run_e5f_simple_fertility_tax_impact as taxcheck
import run_e5f_post2023_rebated_property_tax_smoke as tax
from run_e5f_rebated_tax_shapley_diagnosis import shapley_decomposition
from intergen_eqscale_seq_optimized import joint_nested

adapter, audit, calendar, solver = common.adapter, common.audit, common.policy.calendar, common.solver
SCHEMA = 'e5f_simple_fertility_tax_channels_v1'
ENDPOINTS = ('tax1-equal-rebate', 'tax2-equal-rebate')
BITS = tuple(itertools.product((0, 1), repeat=3))
METRICS = ('births_per_household', 'rooms_per_household', 'ownership',
           'young_ownership_whole_nodes', 'young_mean_rooms_whole_nodes')


def verify_sources(c):
    if str(Path(__file__).resolve()) not in c['source_hashes']:
        raise RuntimeError('Channels wrapper must be source-pinned')
    common.verify_sources(c)


def load_contract(path, expected, *, collecting=False):
    adapter.verify(path, expected)
    c = adapter.read_json(path)
    if c['schema'] != SCHEMA or c['case_seconds'] != 1800:
        raise RuntimeError('Wrong channels schema or time budget')
    if not collecting and time.time() > float(c['launch_deadline_epoch']):
        raise RuntimeError('Channel-cell launch deadline exceeded')
    for key in ('selected_summary', 'output_root'):
        if not Path(c[key]).is_absolute():
            raise RuntimeError(f'Expected absolute path: {key}')
    if set(c['endpoint_receipts']) != set(ENDPOINTS) or set(c['endpoint_receipt_sha256']) != set(ENDPOINTS):
        raise RuntimeError('Exactly the two rebated endpoints are required')
    verify_sources(c)
    return c


def endpoint_receipts(c):
    result = []
    for name in ENDPOINTS:
        path = Path(c['endpoint_receipts'][name])
        if not path.is_absolute():
            raise RuntimeError('Endpoint receipt path must be absolute')
        receipt = taxcheck.verify_receipt(path, c, expected_hash=c['endpoint_receipt_sha256'][name])
        if receipt['case'] != name or receipt['scientific_bundle_sha256'] != c['code_bundle_sha256']:
            raise RuntimeError('Endpoint case or scientific source mismatch')
        endpoint = receipt['endpoint']
        result.append(dict(path=path, receipt=receipt,
                           annual_tax=endpoint['annual_property_tax_rate'],
                           price=endpoint['quantities']['asset_price'],
                           transfer=endpoint['ledger']['equal_transfer_period_units']))
    return result


def verify_cell(path, c, expected_contract, expected_id):
    receipt = taxcheck.verify_receipt(path, c, same_contract=expected_contract)
    if receipt['cell'] != expected_id or receipt['bits'] != list(BITS[expected_id]):
        raise RuntimeError('Channel cell identity mismatch')
    if receipt['scientific_bundle_sha256'] != c['code_bundle_sha256']:
        raise RuntimeError('Channel cell scientific bundle mismatch')
    return receipt


def birth_by_age_parity(e, P, folder):
    rows, total, post_gap = [], 0., 0.
    for j in range(P.J):
        post, effective, births, attempts, risk = joint_nested.factor_age(
            e.g_pre[:, :, :, j], e.policy.joint_choice, P, j)
        post_gap += float(np.abs(post - e.g_post_fertility[:, :, :, j]).sum())
        total += float(births.sum())
        for n in range(P.n_parity):
            den = float(risk[n])
            rows.append(dict(age_node=float(P.age_start + j*P.da), parity=n,
                             at_risk_households=den, explicit_births=float(births[n]),
                             birth_attempts=float(attempts[n]),
                             birth_rate=float(births[n])/den if den else None,
                             attempt_rate=float(attempts[n])/den if den else None))
        del post, effective
    if abs(total-float(e.births)) > 2e-10 or post_gap > 2e-10:
        raise RuntimeError('Age/parity birth reporting does not reproduce the owned joint kernel')
    common.policy.baseline.write_csv(folder/'births_by_age_parity.csv', rows)
    return dict(explicit_births=total, explicit_birth_gap=abs(total-float(e.births)), post_fertility_l1=post_gap,
                scope='Explicit modeled births and attempts from owned joint kernel; no topcode adjustment at age/parity level')


def run_cell(c, contract_hash, cell, folder, progress):
    ends = endpoint_receipts(c)
    root = Path(c['output_root'])
    if cell not in (0, 7):
        for endpoint_id in (0, 7):
            verify_cell(root/f'cell_{endpoint_id:03d}'/'receipt.json', c, contract_hash, endpoint_id)
    packet, prepared, inherited = common.selected_packet(c)
    tax_bit, price_bit, rebate_bit = BITS[cell]
    P = copy.deepcopy(packet['parameters'])
    case = tax.TaxCase(name=f'channel_cell_{cell}', label='fixed component diagnostic',
                       annual_tax_rate=ends[tax_bit]['annual_tax'], rebate_revenue=True)
    transfer, price = ends[rebate_bit]['transfer'], ends[price_bit]['price']
    tax.set_tax_policy(P, case, transfer)
    common.assert_model(P, prepared.supply_rule)
    counter = calendar.SolveCounter()
    progress['phase'] = 'fixed_component_household_solve'
    shared = solver.precompute_shared(P, prepared.b_grid)
    e = calendar.evaluate_period(np.array([price]), inherited.g_pre.copy(), P,
                                prepared.b_grid, shared, counter, supply_rule=prepared.supply_rule)
    current = dict(parameters=P, b_grid=prepared.b_grid, evaluation=e, shared=shared,
                   supply_rule=prepared.supply_rule, state=inherited, calendar_year=2023,
                   diagnostic_cell=cell, component_bits=BITS[cell],
                   scope='fixed components; mixed cells are not equilibria')
    checkpoint = folder/'dated_state.pkl.gz'
    with gzip.open(checkpoint, 'wb', compresslevel=1) as stream:
        pickle.dump(current, stream, protocol=5)
    progress['phase'] = 'household_accounting_audits'
    gates = taxcheck.evaluation_gates(current, inherited)
    budget, arrays = audit.budget_audit(current, folder), audit.policy_array_audit(current, folder)
    if not math.isfinite(float(budget['budget_excess_mass'])) or budget['budget_excess_mass'] > 2e-10:
        raise RuntimeError('Occupied budget gate failed')
    if arrays['occupied_negative_steps']:
        raise RuntimeError('Occupied value monotonicity gate failed')
    for bounds in arrays['probabilities'].values():
        if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
            raise RuntimeError('Marginal probability range gate failed')
    ledger = tax.fiscal_ledger(e, P)
    if any(not math.isfinite(float(value)) for value in ledger.values()):
        raise RuntimeError('Nonfinite fiscal ledger')
    values = common.quantities(e, P)
    values.update(taxcheck.lifecycle_quantities(e, P, folder))
    if any(not math.isfinite(float(value)) for value in values.values()):
        raise RuntimeError('Nonfinite cell quantity')
    birth_report = birth_by_age_parity(e, P, folder)
    reproduction = None
    if cell in (0, 7):
        end = ends[cell//7]
        original = audit.load_checkpoint(end['path'].parent/'dated_state.pkl.gz')
        reproduction = common.baseline_comparison(e, original['evaluation'], policies=True)
        for key, value in end['receipt']['endpoint']['quantities'].items():
            if not math.isclose(values[key], value, rel_tol=0, abs_tol=2e-10):
                raise RuntimeError(f'Endpoint quantity reproduction failed: {key}')
        taxcheck.fiscal_checks(e, P, reference=end['receipt']['endpoint']['ledger'])
        del original
    audit.save_json(folder/'quantities.json', values)
    audit.save_json(folder/'fiscal_ledger.json', ledger)
    result = dict(cell=cell, bits=list(BITS[cell]), annual_property_tax_rate=case.annual_tax_rate,
                  asset_price=price, equal_transfer_period_units=transfer, quantities=values,
                  family_group_quantities=adapter.read_csv(folder/'family_group_quantities.csv'),
                  age_measurement=adapter.read_json(folder/'age_measurement.json'),
                  ledger=ledger, gates=gates, budget=budget, policy_arrays=arrays,
                  birth_reporting=birth_report, endpoint_reproduction=reproduction,
                  bellman_solves=int(counter.bellman), checkpoint_sha256=adapter.digest(checkpoint),
                  equilibrium_gates_applied=cell in (0,7),
                  scope='Fixed tax/price/rebate components: mixed cells intentionally need not clear housing or balance government budgets')
    audit.save_json(folder/'cell_result.json', result)
    progress['completed_evaluations'] = 1
    return result


def collect(c, contract_hash, folder):
    endpoint_receipts(c)
    receipts = [verify_cell(Path(c['output_root'])/f'cell_{i:03d}'/'receipt.json', c, contract_hash, i)
                for i in range(8)]
    results = [receipt['result'] for receipt in receipts]
    metric_values = {key: {BITS[i]: row['quantities'][key] for i,row in enumerate(results)} for key in METRICS}
    groups = [{(r['age_group'],r['family_group']): r for r in row['family_group_quantities']} for row in results]
    if any(set(g) != set(groups[0]) for g in groups):
        raise RuntimeError('Family group definitions differ across cells')
    omitted = []
    for group in sorted(groups[0]):
        for metric in ('owner_rate', 'mean_rooms'):
            name = ':'.join((*group,metric))
            values = [g[group][metric] for g in groups]
            if any(value in (None, '') for value in values):
                omitted.append(dict(metric=name, reason='no household support in at least one cell'))
                continue
            metric_values[name] = {BITS[i]: float(value) for i,value in enumerate(values)}
    rows = []
    for metric, values in metric_values.items():
        if any(not math.isfinite(float(x)) for x in values.values()):
            raise RuntimeError(f'Nonfinite decomposition input: {metric}')
        contributions = shapley_decomposition(values)
        base, reform = float(values[(0,0,0)]), float(values[(1,1,1)])
        own = 'ownership' in metric or metric.endswith('owner_rate')
        for component, value in contributions.items():
            rows.append(dict(metric=metric, component=component, baseline=base, reform=reform,
                             contribution_level=value, contribution_percent_of_baseline=100*value/base if base else None,
                             contribution_percentage_points=100*value if own else None,
                             total_change_level=reform-base, add_up_gap=sum(contributions.values())-(reform-base)))
    common.policy.baseline.write_csv(folder/'shapley_components.csv', rows)
    summary = dict(status='complete', contract_sha256=contract_hash, cell_count=8,
                   primary_comparison='rebated annual2percent versus rebated annual1percent',
                   components=['tax_rate','asset_price','equal_rebate'], decomposition=rows,
                   omitted_metrics=omitted, decomposition_add_up_tolerance=1e-12,
                   selected_summary_sha256=c['selected_summary_sha256'],
                   scope='Exact eight-cell Shapley allocation of the impact change; mixed cells are conditional responses, not equilibrium policies',
                   production_promoted=False, figures='suppressed_at_author_request')
    audit.save_json(folder/'summary.json', summary)
    lines = ['# Property-tax impact channels', '',
             'Comparison: annual 2% versus 1%, both with equal household rebates, from the same inherited population.', '',
             'Eight fixed tax, price and rebate combinations are evaluated. Mixed combinations are conditional household responses, not equilibria.', '',
             '| Outcome | Tax | Asset price | Equal rebate | Total |',
             '|---|---:|---:|---:|---:|']
    for metric in METRICS:
        selected = [r for r in rows if r['metric'] == metric]
        own = 'ownership' in metric
        key = 'contribution_percentage_points' if own else 'contribution_percent_of_baseline'
        suffix = ' pp' if own else '%'
        values = {r['component']: r[key] for r in selected}
        numbers = [values[name] for name in ('tax_rate','asset_price','equal_rebate')]
        lines.append('| '+metric+' | '+' | '.join(f'{x:+.6f}{suffix}' for x in (*numbers,sum(numbers)))+' |')
    lines += ['', 'Young ownership and rooms use unchanged whole-node model age aggregation; no annual-age ACS alignment is implied.',
              'Family groups, lifecycle quantities and explicit age/parity birth-attempt rates are saved in each cell. No new calibration or future transition was run.', '']
    (folder/'READOUT.md').write_text('\n'.join(lines))
    return summary


def execute(args):
    c = load_contract(args.contract, args.contract_sha256, collecting=args.collect)
    label = 'report' if args.collect else f'cell_{args.cell:03d}'
    folder = Path(c['output_root'])/label
    folder.mkdir(parents=True, exist_ok=False)
    started, stop = time.monotonic(), threading.Event()
    progress = dict(phase='verifying_inputs', task=label, completed_evaluations=0)
    def heartbeat():
        while not stop.is_set():
            elapsed = time.monotonic()-started
            audit.save_json(folder/'heartbeat.json', dict(progress, epoch=time.time(), elapsed_seconds=elapsed))
            if elapsed >= 1800:
                audit.save_json(folder/'failure.json', dict(error='Declared 1800-second cap exceeded', **progress))
                os._exit(124)
            stop.wait(30)
    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        audit.save_json(folder/'latest_completed_case.json', dict(status='none_completed'))
        audit.save_json(folder/'best_so_far.json', dict(scope='fixed selected calibration; no search', loss=common.SELECTED_LOSS))
        result = collect(c,args.contract_sha256,folder) if args.collect else run_cell(c,args.contract_sha256,args.cell,folder,progress)
        audit.save_json(folder/'latest_completed_case.json', result)
        verify_sources(c)
        files = sorted(p for p in folder.rglob('*') if p.is_file() and p.name not in ('receipt.json','heartbeat.json'))
        receipt = dict(status='complete', contract_sha256=args.contract_sha256,
                       selected_summary_sha256=c['selected_summary_sha256'], selected_checkpoint_sha256=c['checkpoint_sha256'],
                       scientific_bundle_sha256=c['code_bundle_sha256'], elapsed_seconds=time.monotonic()-started,
                       artifact_sha256={str(p.relative_to(folder)):adapter.digest(p) for p in files},
                       result=result, figures='suppressed_at_author_request', production_promoted=False)
        if not args.collect:
            receipt.update(cell=args.cell,bits=list(BITS[args.cell]))
        audit.save_json(folder/'receipt.json',receipt)
        progress['phase']='complete'
    except BaseException as error:
        progress['phase']='failed'
        audit.save_json(folder/'failure.json',dict(error=str(error),type=type(error).__name__,traceback=traceback.format_exc()))
        raise
    finally:
        stop.set()
        thread.join(timeout=2)
        audit.save_json(folder/'heartbeat.json',dict(progress,epoch=time.time(),elapsed_seconds=time.monotonic()-started))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract',type=Path,required=True)
    parser.add_argument('--contract-sha256',required=True)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--cell',type=int,choices=range(8))
    mode.add_argument('--collect',action='store_true')
    execute(parser.parse_args())
