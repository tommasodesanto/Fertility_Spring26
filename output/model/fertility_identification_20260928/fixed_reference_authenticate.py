#!/usr/bin/env python3
"""Read-only Torch authentication; writes a small, exclusive reference manifest.

2007 stationary reference — block0506, September 28 verified export.
No model solves, rendering, calibration, checkpoint copies, or active-code edits.
Run inside the existing immutable snapshot's container mount, under Slurm.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26')
PHYSICAL = Path('/scratch/td2248/projects/fertility_night_calibration_20260928_v1/project')
B = ROOT / 'output/model/fertility_identification_20260928'
EXPECTED = '68323aadd2c9ad221742842ace9ab108e40437303f0d34da00e7cd83b89f5abf'
LABEL = '2007 stationary reference — block0506, September 28 verified export'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def physical(path):
    return str(PHYSICAL / Path(path).relative_to(ROOT))


def json_value(value):
    if hasattr(value, 'shape') and getattr(value, 'size', 0) > 1024:
        h = hashlib.sha256()
        for start in range(0, value.size, 65536):
            h.update(value.flat[start:start + 65536].tobytes())
        return dict(serialized_array=True, shape=list(value.shape), dtype=str(value.dtype),
                    size=int(value.size), sha256_c_order_bytes=h.hexdigest())
    if hasattr(value, 'tolist'):
        return json_value(value.tolist())
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {'nonfinite_float': repr(value)}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f'Unserialized parameter type: {type(value)}')


def main():
    check(sys.platform == 'linux' and os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Torch Slurm only')
    out = B / 'fixed_reference_manifest.json'
    check(not out.exists(), 'Never overwrite the frozen reference manifest')
    check(sha(B / 'contract_v1/contract.json') == EXPECTED, 'Contract pin')
    c = read(B / 'contract_v1/contract.json')
    pins = list(c['files'].values()) + [c['source_manifest'], c['source_contract'], c['native_ancestry_contract']]
    for pin in pins:
        check(sha(pin['path']) == pin['sha256'], 'Source pin: ' + pin['path'])
    inventory = read(c['source_manifest']['path'])['files']
    for name, digest in inventory.items():
        check(sha(ROOT / name) == digest, 'Scientific source: ' + name)
    spec = c['lanes']['primary']
    check(sha(spec['objective']['path']) == spec['objective']['sha256'], 'Objective pin')
    objective = read(spec['objective']['path'])
    fingerprint = hashlib.sha256(json.dumps(objective['target_rows'], sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    check(fingerprint == spec['target_weight_fingerprint'], 'Target/weight fingerprint')
    export = B / 'resume_v1/selected_export/primary'
    hashes = read(export / 'artifact_hashes.json')
    for name, digest in hashes.items():
        check(sha(export / name) == digest, 'Export artifact: ' + name)
    receipt = read(export / 'receipt.json')
    er = read(export / 'export_receipt.json')
    fits, params = rows(export / 'target_fit.csv'), rows(export / 'parameters.csv')
    check(len(fits) == 14 and len(params) == 31, 'Complete tables')
    check(receipt['point'] == c['initial_point'] == er['selected']['point'], 'Frozen coordinates')
    check(receipt['source_manifest_sha256'] == c['source_manifest']['sha256'], 'Receipt source identity')
    check(receipt['target_system_sha256'] == spec['objective']['sha256'], 'Receipt objective identity')
    check(receipt['target_weight_fingerprint'] == fingerprint, 'Receipt target identity')
    target_rows = {r['restriction_id']: r for r in objective['target_rows']}
    for r in fits:
        target = target_rows[r['moment']]
        check(float(r['target']) == target['target'], 'Target: ' + r['moment'])
        check(abs(float(r['gap']) - (float(r['model']) - float(r['target']))) < 1e-12, 'Gap arithmetic')
        if r['weight']:
            check(float(r['weight']) == target['actual_weight'], 'Weight: ' + r['moment'])
            check(abs(float(r['loss_contribution']) - float(r['weight']) * float(r['gap'])**2) < 1e-12, 'Loss arithmetic')
    loss = math.fsum(float(r['loss_contribution']) for r in fits if r['loss_contribution'])
    check(abs(loss - 19.581310760138322) < 1e-12, 'Common-primary loss')
    anchor = Path(c['anchor']['case_path'])
    for key in ('receipt', 'target_fit', 'parameters'):
        pin = c['anchor'][key]
        check(sha(pin['path']) == pin['sha256'], 'Overnight anchor ' + key)
    check(rows(anchor / 'parameters.csv') == params, 'All 31 anchor parameters')
    anchor_fits = {r['moment']: r for r in rows(anchor / 'target_fit.csv')}
    for r in fits:
        for key in ('target', 'model', 'gap'):
            check(r[key] == anchor_fits[r['moment']][key], 'Anchor physical row ' + r['moment'])
    checkpoints = []
    for record in [er['selected']] + er['repeats']:
        p = Path(record['case_path'])
        rr = read(p / 'receipt.json')
        check(sha(p / 'receipt.json') == record['receipt_sha256'], 'Case receipt')
        digest = sha(p / 'initial_state.pkl.gz')
        check(digest == record['checkpoint_sha256'] == rr['case_checkpoint_sha256'], 'Case checkpoint')
        check(rows(p / 'parameters.csv') == params and rows(p / 'target_fit.csv') == fits, 'Exact repeated tables')
        checkpoints.append(dict(case=record['case'], physical_path=physical(p / 'initial_state.pkl.gz'), sha256=digest))
    anchor_receipt = read(anchor / 'receipt.json')
    anchor_digest = sha(anchor / 'initial_state.pkl.gz')
    check(anchor_digest == anchor_receipt['case_checkpoint_sha256'], 'Original block0506 checkpoint')
    check(sha(export / 'initial_state.pkl.gz') == receipt['case_checkpoint_sha256'] == er['repeats'][0]['checkpoint_sha256'], 'Export checkpoint alias')
    names = c['standard_diagnostic_names']
    check(len(names) == 17, 'Standard plot count')
    for name in names:
        values = [sha(Path(r['case_path']) / 'standard_diagnostics' / name) for r in er['repeats']]
        check(values[0] == values[1] == hashes['standard_diagnostics/' + name], 'Repeated plot ' + name)
    print('Source, checkpoint, 14/31 tables and 17 plot authentication PASS; zero solves', flush=True)
    # Register the exact frozen runtime classes, without solving, before unpickling.
    sys.path[:0] = [str(ROOT / 'code/model/tools'), str(ROOT / 'tmp/e5f_overnight_local_20260927/portable/tools_v4')]
    import e5f_evening_calibration_runtime as runtime
    prepared = runtime.setup(dict(c, objective=spec['objective']), objective,
                             B / ('fixed_reference_auth_runtime_' + os.environ['SLURM_JOB_ID']))
    with gzip.open(export / 'initial_state.pkl.gz', 'rb') as stream:
        packet = pickle.load(stream)
    P = packet['parameters']
    actual = dict(prepared.tax.actual_parameters(P), delta_alpha_jump=P.delta_alpha_jump,
                  child_benefit_curvature=P.child_benefit_curvature, tenure_choice_kappa=P.tenure_choice_kappa)
    actual.update(psi_child=P.psi_child, child_benefit_CRRA_coefficient=(1-P.child_benefit_curvature)*P.psi_child,
        theta1=P.theta1, sigma=P.sigma, alpha_cons=P.alpha_cons, delta_alpha=P.delta_alpha, h_P=0.,
        utility_reference_rent=P.utility_reference_rent, q_annual=(1+P.q)**(1/P.period_years)-1,
        financed_share=P.phi[0], housing_supply_elasticity=P.xi_supply[0], payroll_tax=P.tau_pay,
        pension_period=P.pension, annual_depreciation=prepared.ancestor.ANNUAL_DEP,
        period_depreciation=P.delta, annual_property_tax=prepared.ancestor.ANNUAL_PROPERTY_TAX,
        period_property_tax=P.tau_H, selling_cost=P.psi, rental_cap=P.hR_max,
        wealth_grid_nodes=len(packet['b_grid']), income_states=len(P.z_grid))
    for r in params:
        check(float(actual[r['parameter']]) == float(r['estimate']), 'Serialized parameter ' + r['parameter'])
        if r['lower']:
            lo, hi, value = float(r['lower']), float(r['upper']), float(r['estimate'])
            check(lo <= value <= hi, 'Parameter bound')
            check((min(value-lo, hi-value) <= .01*(hi-lo)) == (r['near_bound'] == 'True'), 'Near-bound flag')
    all_parameters = json_value(vars(P))
    manifest = dict(schema='fixed_calibration_reference_v1', label=LABEL,
        created_utc=datetime.now(timezone.utc).isoformat(), slurm_job=os.environ['SLURM_JOB_ID'],
        status='authenticated_existing_reference_zero_solves', model_solves=0, rendering=0,
        preservation='No reference checkpoint, source or selected-export artifact modified or copied. Never overwrite this manifest; explicit author switch requires a new identity and comparison.',
        local_export=str(export), torch_project_root=str(PHYSICAL), container_project_root=str(ROOT),
        torch_export=physical(export), checkpoint=checkpoints[1], other_case_checkpoints=checkpoints,
        original_block0506_checkpoint=dict(physical_path=physical(anchor / 'initial_state.pkl.gz'), sha256=anchor_digest),
        contract=dict(path=str(B / 'contract_v1/contract.json'), sha256=EXPECTED),
        source_manifest=c['source_manifest'], source_contract=c['source_contract'], native_ancestry_contract=c['native_ancestry_contract'],
        source_files_verified=len(inventory), controller_and_ancestry_pins_verified=len(pins),
        scientific_identity=receipt['scientific_identity'], scientific_candidate_id=receipt['scientific_candidate_id'],
        objective=spec['objective'], target_weight_fingerprint=fingerprint,
        common_primary_loss=er['selected']['primary_rescore'], recomputed_loss=loss,
        full_target_table=fits, full_parameter_table=params, actual_serialized_parameters=all_parameters,
        parameter_class=f'{type(P).__module__}.{type(P).__name__}',
        full_parameters_note='Every serialized instance field authenticated; arrays over 1024 entries pinned by dimensions, dtype and SHA256 of C-order bytes. Other values retained verbatim. All 31 reported values checked against the checkpoint; no constructor defaults substituted.',
        artifact_hashes=hashes, standard_diagnostic_names=names,
        repeats='Two exact 14-row/31-row table repeats and all 17 PNG hash matches verified again.',
        economic_contract=dict(reference_year=2007, deferred_transition_year=2023, calibration_normalization_target=2.1,
            frozen_psi_child=P.psi_child, counterfactual_renormalization=False,
            preserve='All preferences, earnings, entry distribution, timing, credit and closure primitives remain frozen unless an experiment explicitly names its changes.',
            equilibrium_endogenous='Prices, choices and distributions may adjust only under the specified experiment and closure; payroll/pension and population closure must be reconciled before production runs.',
            warning='Old transition and frictionless-credit results are not authenticated counterfactuals for this reference.'),
        inherited_gates={k: receipt[k] for k in ('adult_entry_gate','fiscal','household_budget','purchase_accounting','policy_array_gates','market_residual')},
        model_observer_warnings=receipt['model_observer_warnings'],
        limitations=['Numerical identity and existing gates are verified, not global optimality or economic validity.',
            'Standard policy panels include buyer-conditional choices and full-grid states; weight realized behavior by occupied pre-choice mass.',
            'High-wealth ownership decline, age30 housing downturn and retirement profiles remain unresolved.',
            'Legacy diagnostic summary moments use different definitions from the authoritative target-fit observers.',
            'Estate recipient/counterparty/physical settlement and empirical observer approximations remain provisional.',
            'Torch scratch retention is external to this manifest; checkpoint remains at its existing path.'])
    encoded = json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + '\n'
    check(len(encoded.encode()) < 250000, 'Manifest must remain compact')
    with out.open('x') as stream:
        stream.write(encoded)
    print(json.dumps(dict(status='PASS', manifest=str(out), sha256=sha(out), size_bytes=out.stat().st_size,
                         parameter_fields=len(all_parameters), source_files=len(inventory), model_solves=0)))


if __name__ == '__main__':
    main()
