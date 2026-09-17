"""Observe a pinned revised initial checkpoint without solving a new economy."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import pickle
import time

import run_e5f_matched_pf_smoke as primitive
from e5f_initial_fertility_observer import observe_initial_fertility
from e5f_initial_housing_observer import observe_initial_housing_wealth
from e5f_parenthood_utility import validate_parenthood_utility
from e5f_stationary_paygo import certify_initial_pension


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', type=Path, required=True)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    primitive.verify(args.contract, args.contract_sha256)
    contract = json.loads(args.contract.read_text())
    if contract.get('schema') != 'e5f_initial_measurement_probe_v1' or contract.get('calibrated_smm') is not False:
        raise ValueError('Explicit diagnostic-only contract required')
    root = Path(__file__).resolve().parents[3]
    required = {str(p.relative_to(root)) for p in (root/'code/model').rglob('*.py')}
    if not required.issubset(contract['source_sha256']):
        raise ValueError('Incomplete source pins')
    for name, pin in contract['source_sha256'].items():
        path = (root/name).resolve()
        if not path.is_relative_to(root):
            raise ValueError('Source path outside snapshot')
        primitive.verify(path, pin)
    primitive.verify(contract['checkpoint'], contract['checkpoint_sha256'])
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    write = primitive.pf.calendar.write_json_atomic
    write(args.output/'contract.json', contract)
    write(args.output/'heartbeat.json', dict(phase='loading', elapsed_seconds=0.))
    primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = primitive.pf.transition.advance_sequential_calendar_distribution
    with gzip.open(contract['checkpoint'], 'rb') as stream:
        packet = pickle.load(stream)
    P, evaluation = packet['parameters'], packet['evaluation']
    validate_parenthood_utility(P)
    if not P.exhaustive_saving_control or P.joint_nested_choice:
        raise ValueError('Active sequential optimizer required')
    fiscal = certify_initial_pension(evaluation.g_current, P, marginal_tolerance=1e-9, fiscal_tolerance=1e-6)
    fertility = {}
    for projection in ('uniform_birth_time', 'constant_post_cell'):
        fertility[projection] = observe_initial_fertility(evaluation, P, age_projection=projection)
        write(args.output/f'fertility_{projection}.json', fertility[projection])
    write(args.output/'heartbeat.json', dict(phase='housing_and_matched_birth_branch', elapsed_seconds=time.monotonic()-started))
    housing = observe_initial_housing_wealth(evaluation, P, packet['b_grid'], packet['shared'],
        diagnostic_enabled=True, age_projection='uniform_within_age_cell',
        diagnostic_allow_family_proxies=True, include_wealth=True, include_birth_response=True)
    write(args.output/'housing_wealth.json', housing)
    result = dict(schema='e5f_initial_measurement_probe_v1', status='complete_diagnostic',
        calibrated_smm=False, empirical_target_contract_activated=False,
        checkpoint_sha256=contract['checkpoint_sha256'], source_files_verified=len(contract['source_sha256']),
        fiscal=fiscal, fertility=fertility, housing_wealth=housing,
        elapsed_seconds=time.monotonic()-started, household_or_equilibrium_solves=0)
    write(args.output/'summary.json', result)
    write(args.output/'heartbeat.json', dict(phase='complete', elapsed_seconds=time.monotonic()-started))


if __name__ == '__main__':
    main()
