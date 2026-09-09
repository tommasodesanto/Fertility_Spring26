#!/usr/bin/env python3
"""Render the unchanged 17 diagnostics from one verified saved tax-path date.

This driver reads solved states. It never solves household or market problems,
changes numerical kernels, or writes into the numerical transition output tree.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import time
import traceback

import numpy as np
import run_e5f_simple_fertility_tax_transition as transition
from intergen_eqscale_seq_optimized import diagnostics

common, adapter, audit = transition.common, transition.adapter, transition.audit
BASE_GRAPHS = {
    'ownership_by_age.png', 'fertility_by_age.png', 'housing_market.png',
    'market_clearing_by_market.png', 'market_clearing_residuals.png',
    'tenure_services.png', 'owner_rungs.png', 'housing_prices.png',
    'income_state_outcomes.png', 'ownership_by_age_income_state.png',
    'liquid_wealth_by_age_income_state.png', 'housing_by_age_income_state.png',
    'fertility_policy_by_age_income_state.png',
}
REQUIRED_ARTIFACTS = (
    'dated_state.pkl.gz', 'endpoint_receipt.json', 'advance.json', 'reproduction.json',
    'quantities.json', 'fiscal_ledger.json', 'budget_summary.json',
    'policy_array_summary.json', 'age_measurement.json', 'family_group_quantities.csv', 'lifecycle.csv',
)


def finite_gate(value, maximum, label):
    number = float(value)
    if not math.isfinite(number) or abs(number) > maximum:
        raise RuntimeError(f'{label} gate failed: {value}')


def expected_graphs(P):
    result = set(BASE_GRAPHS)
    for age in np.asarray(getattr(P, 'diagnostic_policy_ages', [30., 42.]), dtype=float).reshape(-1):
        j = diagnostics.age_to_index(P, float(age))
        node = int(round(P.age_start + j*P.da))
        result.add(f'wealth_dist_childless_renter_age{node}.png')
        result.add(f'policy_childless_renter_age{node}.png')
    if len(result) != 17:
        raise RuntimeError('Saved model does not imply the established seventeen-graph set')
    return result


def verify_date(c, contract_hash, case, year, stage):
    """Verify the completed stage receipt, then only artifacts used for this date."""
    folder = Path(c['output_root'])/stage/case
    path = folder/'receipt.json'
    receipt_hash = (folder/'receipt.sha256').read_text().strip()
    adapter.verify(path, receipt_hash)
    receipt = adapter.read_json(path)
    years = transition.YEARS[:2] if stage == 'smoke' else transition.YEARS
    if (receipt['status'] != 'complete' or receipt['stage'] != stage or receipt['case'] != case
            or receipt['years'] != years or year not in years
            or receipt['contract_sha256'] != contract_hash):
        raise RuntimeError('Completed matching transition stage/date required')
    if (receipt['selected_summary_sha256'] != c['selected_summary_sha256']
            or receipt['selected_checkpoint_sha256'] != c['checkpoint_sha256']
            or receipt['scientific_bundle_sha256'] != c['code_bundle_sha256']
            or receipt['closure'] != c['closure']):
        raise RuntimeError('Transition stage provenance differs from its frozen contract')
    adapter.verify(c['selected_summary'], c['selected_summary_sha256'])
    dated = folder/f'date_{year}'
    verified = {}
    for name in REQUIRED_ARTIFACTS:
        relative = f'date_{year}/{name}'
        if relative not in receipt['artifact_sha256']:
            raise RuntimeError(f'Date artifact absent from completed receipt: {relative}')
        expected = receipt['artifact_sha256'][relative]
        adapter.verify(dated/name, expected)
        verified[name] = expected
    endpoint = adapter.read_json(dated/'endpoint_receipt.json')
    advance = adapter.read_json(dated/'advance.json')
    reproduction = adapter.read_json(dated/'reproduction.json')
    if (endpoint['status'] != 'complete' or endpoint['policy'] != case
            or endpoint['calendar_year'] != year or advance['calendar_year'] != year
            or advance['next_calendar_year'] != year+4
            or advance['row']['calendar_year'] != year):
        raise RuntimeError('Dated validity or calendar metadata missing')
    if endpoint['checkpoint_sha256'] != verified['dated_state.pkl.gz']:
        raise RuntimeError('Endpoint checkpoint digest differs from completed-stage digest')
    if reproduction != receipt['reproduction'][str(year)]:
        raise RuntimeError('Date reproduction receipt differs from completed stage')
    if year == 2023 and 'verified_impact' not in reproduction:
        raise RuntimeError('Initial date lacks verified impact reproduction')
    if stage == 'full' and year in transition.YEARS[:2] and 'smoke' not in reproduction:
        raise RuntimeError('Full path lacks required smoke reproduction')
    finite_gate(advance['row']['mass_accounting_residual'], 2e-10, 'cohort mass')
    finite_gate(endpoint['budget']['budget_excess_mass'], 2e-10, 'household budget')
    if endpoint['policy_arrays']['occupied_negative_steps'] != 0:
        raise RuntimeError('Occupied value audit failed')
    for value in endpoint['normalized_root_residuals']:
        finite_gate(value, 1e-4, 'coupled root')
    finite_gate(endpoint['ledger']['government_budget_residual'], 2.5e-5, 'fiscal balance')
    finite_gate(endpoint['quantities']['market_residual'], 2e-4, 'housing market')
    finite_gate(endpoint['quantities']['feasibility_projection'], 1e-6, 'feasibility projection')
    for value in endpoint['gates']['mass_gaps'].values():
        finite_gate(value, 2e-10, 'dated mass')
    for bounds in endpoint['policy_arrays']['probabilities'].values():
        if bounds['nonfinite'] or bounds['minimum'] < 0 or bounds['maximum'] > 1:
            raise RuntimeError('Saved marginal probability audit failed')
    return dated, endpoint, dict(stage_receipt_sha256=receipt_hash, verified_date_artifact_sha256=verified)


def execute(args):
    c = transition.load_contract(args.transition_contract, args.transition_contract_sha256, launching=False)
    output_root = args.output_root.resolve()
    if not args.output_root.is_absolute():
        raise RuntimeError('Graph output root must be absolute')
    numerical_root = Path(c['output_root']).resolve()
    if output_root == numerical_root or output_root.is_relative_to(numerical_root):
        raise RuntimeError('Graphs must be written outside the numerical output tree')
    dated, endpoint, provenance = verify_date(c,args.transition_contract_sha256,args.case,args.year,args.stage)
    folder = output_root/args.case/f'date_{args.year}'
    folder.mkdir(parents=True,exist_ok=False)
    started = time.monotonic()
    source_paths = (Path(__file__).resolve(),Path(audit.__file__).resolve(),Path(diagnostics.__file__).resolve())
    graphics_hashes = {str(path):adapter.digest(path) for path in source_paths}
    try:
        packet = transition.load_packet(dated/'dated_state.pkl.gz')
        if packet['calendar_year'] != args.year or packet['policy_case'] != args.case:
            raise RuntimeError('Saved checkpoint date or case differs')
        common.assert_model(packet['parameters'],packet['supply_rule'])
        names = expected_graphs(packet['parameters'])
        audit.standard_diagnostics(packet,folder,validate_production_young=False)
        graph_folder = folder/'standard_diagnostics'
        graphs = sorted(graph_folder.glob('*.png'))
        if {p.name for p in graphs} != names:
            raise RuntimeError('Stable seventeen-graph filename set is incomplete or changed')
        units = adapter.read_json(folder/'market_quantity_units.json')
        transition.verify_sources(c)
        for path,digest in graphics_hashes.items():
            adapter.verify(path,digest)
        files = sorted(p for p in folder.rglob('*') if p.is_file())
        manifest = dict(status='complete',case=args.case,calendar_year=args.year,stage=args.stage,
                        transition_contract_sha256=args.transition_contract_sha256,
                        selected_summary_sha256=c['selected_summary_sha256'],
                        selected_checkpoint_sha256=c['checkpoint_sha256'],
                        scientific_bundle_sha256=c['code_bundle_sha256'],
                        input_checkpoint=str(dated/'dated_state.pkl.gz'),
                        input_checkpoint_sha256=endpoint['checkpoint_sha256'],
                        input_provenance=provenance,graphics_source_sha256=graphics_hashes,
                        standard_png_count=len(graphs),standard_filenames=sorted(names),
                        graph_sha256={p.name:adapter.digest(p) for p in graphs},
                        artifact_sha256={str(p.relative_to(folder)):adapter.digest(p) for p in files},
                        units=units,age_measurement=endpoint['age_measurement'],
                        elapsed_seconds=time.monotonic()-started,
                        scope='Read-only rendering of verified temporary-equilibrium household-unit diagnostic; not a resident-population forecast',
                        numerical_solves=0,kernel_changes=False,production_promoted=False)
        audit.save_json(folder/'graph_manifest.json',manifest)
    except BaseException as error:
        audit.save_json(folder/'failure.json',dict(error=str(error),type=type(error).__name__,traceback=traceback.format_exc()))
        raise


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--transition-contract',type=Path,required=True)
    parser.add_argument('--transition-contract-sha256',required=True)
    parser.add_argument('--case',choices=transition.CASES,required=True)
    parser.add_argument('--year',type=int,choices=transition.YEARS,required=True)
    parser.add_argument('--stage',choices=('smoke','full'),default='full')
    parser.add_argument('--output-root',type=Path,required=True)
    execute(parser.parse_args())
