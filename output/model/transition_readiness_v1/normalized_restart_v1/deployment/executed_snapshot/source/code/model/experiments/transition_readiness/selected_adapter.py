#!/usr/bin/env python3
"""Authenticate a selected calibration before any native setup or checkpoint load.

This seam accepts an engine only when its full economic contract equals the
authenticated selected contract. The historical block0506 loader is not such an
engine for the current floor/credit economy.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


class ReadinessBlocked(ValueError):
    pass


def require(test, message):
    if not test:
        raise ReadinessBlocked(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def pinned(item):
    require(isinstance(item, dict) and set(item) == {'path', 'sha256'}, 'Expected exact path/SHA-256 artifact')
    path = Path(item['path'])
    require(path.is_absolute() and path.is_file(), 'Missing absolute pinned artifact: ' + str(path))
    require(isinstance(item['sha256'], str) and len(item['sha256']) == 64 and sha(path) == item['sha256'],
            'Artifact hash mismatch: ' + str(path))
    return path


def read(item):
    return json.loads(pinned(item).read_text())


CONTRACT_KEYS = {'preferences', 'earnings', 'entry', 'grid', 'credit', 'sale_repayment',
                 'mortality_repayment', 'fiscal', 'population', 'housing', 'geography',
                 'estate', 'targets', 'timing'}
CLASSIFICATIONS = {'estimated', 'empirically_normalized', 'externally_fixed', 'outstanding'}


def authenticate(manifest_pin):
    """Return verified metadata; no imports, setup, numerical work or unpickling."""
    m = read(manifest_pin)
    require(m.get('schema') == 'selected_calibration_transition_adapter_v1', 'Unsupported selected manifest schema')
    require(m.get('author_selected') is True and m.get('provisional') is False,
            'Selected calibration must be author-selected and nonprovisional')
    contract = m.get('economic_contract', {})
    require(set(contract) == CONTRACT_KEYS, 'Complete economic contract required')
    for name, row in contract.items():
        require(isinstance(row, dict) and set(row) == {'identity', 'classification', 'evidence'}, 'Invalid contract row: ' + name)
        require(isinstance(row['identity'], str) and row['identity'] and row['classification'] in CLASSIFICATIONS,
                'Unclassified closure object: ' + name)
        require(row['classification'] != 'outstanding', 'Production closure outstanding: ' + name)
        pinned(row['evidence'])
    sources = read(m['source_manifest'])
    require(isinstance(sources, dict) and sources.get('pins'), 'Source manifest must contain nonempty pins')
    for item in sources['pins']:
        pinned(item)
    target = read(m['target_contract'])
    require(target.get('fingerprint') == contract['targets']['identity'], 'Target fingerprint differs from selected contract')
    params = read(m['effective_parameters'])
    require(isinstance(params, dict) and bool(params), 'Complete effective parameters required')
    state = read(m['state_identity'])
    require(state.get('checkpoint_sha256') == m['checkpoint']['sha256'], 'State receipt belongs to a different checkpoint')
    require(state.get('source_manifest_sha256') == m['source_manifest']['sha256'], 'State engine source identity differs')
    require(state.get('effective_parameters_sha256') == m['effective_parameters']['sha256'], 'State parameters identity differs')
    require(state.get('initial_queues_verified') is True and state.get('initial_population_verified') is True,
            'Initial population and dated birth-entry queues require native verification')
    for key in ('entry', 'grid'):
        require(state.get(key + '_identity') == contract[key]['identity'], 'Selected ' + key + ' identity differs')
    checkpoint = pinned(m['checkpoint'])
    repeat = read(m['selected_price_repeat'])
    require(repeat.get('passed') is True and repeat.get('fresh_native_repeat') is True, 'Fresh selected-price repeat is required')
    for key in ('checkpoint_sha256', 'source_manifest_sha256', 'target_fingerprint', 'effective_parameters_sha256', 'state_identity_sha256'):
        expected = {'checkpoint_sha256': m['checkpoint']['sha256'], 'source_manifest_sha256': m['source_manifest']['sha256'],
                    'target_fingerprint': target['fingerprint'], 'effective_parameters_sha256': m['effective_parameters']['sha256'],
                    'state_identity_sha256': m['state_identity']['sha256']}[key]
        require(repeat.get(key) == expected, 'Selected repeat identity differs: ' + key)
    gates = repeat.get('gates')
    require(isinstance(gates, dict) and bool(gates) and all(v is True for v in gates.values()), 'Selected repeat gates failed or absent')
    require(repeat.get('economic_contract') == contract, 'Repeat economic contract differs')
    # Concrete receipts emitted by utility_floor_round2_v1/runner.py and its
    # phase_b_pilot.py; an arbitrary single boolean is not a native repeat proof.
    native_repeat = read(m['native_selected_repeat'])
    require(native_repeat.get('status') == 'exact_full_ge_repeat_passed' and native_repeat.get('target_rows') == 14
            and native_repeat.get('parameter_rows') == 31 and len(native_repeat.get('standard_plot_hashes', {})) == 17,
            'Native exact repeat must retain all 14 targets, 31 parameters and 17 plots')
    native_ge = read(m['native_ge_receipt'])
    require(native_ge.get('status') == 'passed' and isinstance(native_ge.get('exact_repeat'), dict),
            'Native selected GE and fresh selected-price repeat have not passed')
    for key in ('renewal_residual', 'population_scale'):
        require(native_ge.get('selected', {}).get(key) == native_ge['exact_repeat'].get(key)
                and native_ge['exact_repeat'].get(key) is not None, 'Native GE repeat differs: ' + key)
    for row in (native_ge['selected'], native_ge['exact_repeat']):
        require(abs(float(row.get('renewal_residual', math.inf))) <= 1e-6,
                'Native demographic renewal gate fails')
        require(abs(float(row.get('actual_paygo_residual', math.inf))) <= 1e-6, 'Native actual PAYGO gate fails')
        scale = float(row.get('population_scale', 0))
        require(math.isfinite(scale) and scale > 0 and row.get('outside_entry') == 0
                and row.get('birth_to_entry_conversion') == 1/2.1, 'Native population closure differs')
        step = row.get('native_population_step', {})
        require(abs(float(step.get('mass_residual', math.inf))) <= max(1., scale)*2e-8,
                'Native population step mass gate fails')
    return {'manifest': m, 'contract': contract, 'checkpoint': checkpoint, 'parameters': params,
            'source_manifest': sources, 'target_contract': target, 'state_identity': state}


def setup_authenticated(manifest_pin, *, engine_contract, setup):
    """Callable integration seam: only invoke setup after every authentication gate."""
    authenticated = authenticate(manifest_pin)
    require(engine_contract == authenticated['contract'],
            'Transition engine floor/credit/entry/closure contract differs; legacy substitution is forbidden')
    require(callable(setup), 'A real native setup callable is required')
    return setup(authenticated)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--manifest-sha256', required=True)
    p.add_argument('--production', action='store_true')
    args = p.parse_args()
    try:
        result = authenticate({'path': str(args.manifest), 'sha256': args.manifest_sha256})
        require(not args.production, 'Production blocked: this CLI supplies authentication only; native dated-engine validation and full-horizon certification are still required')
        print(json.dumps({'status': 'authenticated_metadata', 'model_calls': 0, 'production_ready': False,
                          'checkpoint_sha256': result['manifest']['checkpoint']['sha256']}))
    except (ReadinessBlocked, KeyError, OSError, json.JSONDecodeError) as exc:
        p.exit(2, 'BLOCKED: ' + str(exc) + '\n')


if __name__ == '__main__':
    main()
