"""Pinned read-only checkpoint/ledger comparison; no optimization or model solve."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify(path, expected):
    if len(expected) != 64 or sha256(path) != expected:
        raise ValueError(f'Pinned SHA256 mismatch: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract', required=True, type=Path)
    parser.add_argument('--contract-sha256', required=True)
    parser.add_argument('--source-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    verify(args.contract, args.contract_sha256)
    c = json.loads(args.contract.read_text())
    verify(__file__, c['runner_sha256'])
    if (c['schema'] != 'e5f_historical_paygo_prefix_checkpoint_comparison_v1'
            or c['years'] != [2007, 2011, 2015, 2019, 2023]
            or c['source_commit'] != '527ab397218e18611ba6930e62283eefc4225688'
            or not 0 < c['fiscal_tolerance'] <= 1e-6
            or not 0 < c['initial_factor_tolerance'] <= 1e-12
            or not 0 < c['age_mass_tolerance'] <= 1e-12
            or c['resources'] != {'cpus': 1, 'memory_gb': 8, 'slurm_minutes': 5}
            or c['seconds'] != 300 or c['production_eligible']
            or c['full_future_marginals_available'] or c['root_coordinates_reduced']):
        raise ValueError('Diagnostic source, timing or unchanged numerical contract violated')
    source = args.source_root.resolve()
    manifest = c['source_sha256']
    files = {str(p.relative_to(source)) for p in (source/'code/model').rglob('*')
             if p.is_file() and '__pycache__' not in p.parts}
    if (files != set(manifest) or len(files) != c['source_inventory_count']
            or sum(name.endswith('.py') for name in files) != c['source_python_count']):
        raise ValueError('Complete code/model source inventory differs from frozen commit')
    for relative, pin in manifest.items():
        path = (source/relative).resolve()
        if not path.is_relative_to(source) or Path(relative).is_absolute():
            raise ValueError('Unsafe source manifest path')
        verify(path, pin)
    for name in ('initial_checkpoint', 'actual_root', 'preflight', 'observed_ages', 'history_contract'):
        verify(c[name]['path'], c[name]['sha256'])
    actual = json.loads(Path(c['actual_root']['path']).read_text())
    preflight = json.loads(Path(c['preflight']['path']).read_text())
    historical = json.loads(Path(c['history_contract']['path']).read_text())
    if (historical['initial_checkpoint'] != c['initial_checkpoint']
            or preflight['status'] != 'passed' or not actual['fresh_path_matches_final']
            or not actual['final']['mapping_valid'] or actual['years'][:5] != c['years']):
        raise ValueError('Historical source linkage or final mapping/replay check failed')
    args.output.mkdir(parents=True, exist_ok=False)

    def save(name, value):
        path = args.output/name
        temporary = path.with_suffix(path.suffix+'.tmp')
        temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
        os.replace(temporary, path)

    save('heartbeat.json', {'phase': 'pins_verified_loading_checkpoint', 'elapsed_seconds': time.monotonic()-start})
    sys.path[:0] = [str(source/'code/model/tools'), str(source/'code/model')]
    import numpy as np
    from e5f_historical_paygo_prefix import predict_historical_paygo_prefix, compare_historical_paygo_prefix
    from e5f_social_security import fiscal_accounts
    from run_e5f_open_population_transition import reweight_distribution_to_acs_2007_ages

    # Trusted pickle is opened only after its original checkpoint hash passes.
    with gzip.open(c['initial_checkpoint']['path'], 'rb') as stream:
        packet = pickle.load(stream)
    P, stationary = packet['parameters'], np.asarray(packet['stationary_g_pre'])
    if (P.I != 1 or P.J != 17 or P.Nb != 120 or stationary.ndim != 7
            or stationary.shape[0] != 120 or P.tau_pay != .179
            or not P.exhaustive_saving_control or not np.isfinite(stationary).all()
            or np.any(stationary < 0.)):
        raise ValueError('Original full-grid sequential checkpoint contract failed')
    bridge = preflight['initial_bridge']['age_reweight']
    ages = P.age_start + np.arange(P.J)*P.da
    np.testing.assert_array_equal(ages, np.asarray(bridge['model_age_grid']))
    np.testing.assert_allclose(stationary.sum(axis=(0, 1, 2, 4, 5, 6)),
        bridge['stationary_age_mass'], rtol=0., atol=c['age_mass_tolerance'])
    g2007 = reweight_distribution_to_acs_2007_ages(stationary, bridge)
    M2007 = g2007.sum(axis=(0, 1, 2, 5, 6))
    # Full 2007 distribution checks the factor units independently of the
    # predictor's singleton-axis representation, at the original period benefit.
    full_accounts = fiscal_accounts(g2007, P)
    preflight_accounts = preflight['initial_bridge']['announcement_state_accounts_at_old_pension']
    initial_factor_gaps = {}
    for name in ('payroll_tax_base_period', 'retiree_benefit_exposure', 'payroll_tax_revenue',
                 'pension_outlays', 'pension_period_units', 'implied_balanced_pension_period'):
        v, target = float(full_accounts[name]), float(preflight_accounts[name])
        gap = abs(v-target)/max(abs(v), abs(target))
        initial_factor_gaps[name] = gap
        if not np.isfinite(gap) or gap > c['initial_factor_tolerance']:
            raise ValueError(f'Full 2007 factor/unit reproduction failed: {name}')
    with Path(c['observed_ages']['path']).open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    expected_keys = {(year, float(age)) for year in c['years'] for age in ages}
    lookup = {(int(row['year']), float(row['model_age'])): row for row in rows}
    if len(rows) != len(lookup) or set(lookup) != expected_keys:
        raise ValueError('Observed age CSV must contain exactly five dates and every model age once')
    totals = np.array([[float(lookup[year, float(age)]['observed_target_mass']) for age in ages]
                       for year in c['years']])
    observed_actual = np.array([[float(lookup[year, float(age)]['actual_after_bridge_mass']) for age in ages]
                                for year in c['years']])
    np.testing.assert_allclose(observed_actual, totals, rtol=0., atol=c['age_mass_tolerance'])
    np.testing.assert_allclose(M2007.sum(axis=1), totals[0], rtol=0., atol=c['age_mass_tolerance'])
    prediction = predict_historical_paygo_prefix(initial_age_income_mass=M2007, parameters=P,
        observed_age_masses=totals, years=c['years'])
    comparison = compare_historical_paygo_prefix(prediction=prediction,
        actual_ledgers=actual['final']['payload']['fiscal_accounts'][:5],
        actual_years=actual['years'][:5], fiscal_tolerance=c['fiscal_tolerance'])
    marginal = prediction.pop('age_income_marginals')
    prediction['pensions_period'] = prediction['pensions_period'].tolist()
    np.savez_compressed(args.output/'predicted_marginals.npz', years=np.asarray(c['years']),
                        age_income_marginals=marginal, initial_age_income_mass=M2007)
    save('prediction.json', prediction)
    save('comparison.json', comparison)
    summary = dict(status='passed_diagnostic_comparison' if comparison['all_comparisons_pass'] else 'failed_diagnostic_comparison',
        contract_sha256=args.contract_sha256, runner_sha256=sha256(__file__),
        source_commit=c['source_commit'], source_files_verified=len(files),
        source_python_files_verified=c['source_python_count'],
        input_pins={k:c[k] for k in ('initial_checkpoint','actual_root','preflight','observed_ages','history_contract')},
        years=c['years'], predicted_pensions_period=prediction['pensions_period'],
        initial_full_distribution_shape=list(g2007.shape), initial_full_fiscal_accounts=full_accounts,
        initial_preflight_relative_factor_gaps=initial_factor_gaps,
        maximum_predicted_pension_scaled_budget_residual=max(abs(row['predicted_pension_scaled_budget_residual_on_actual']) for row in comparison['rows']),
        maximum_relative_ledger_factor_gap=max(max(row['relative_ledger_factor_gaps'].values()) for row in comparison['rows']),
        actual_trial_budgets_pass=comparison['actual_trial_budgets_pass'],
        full_2007_prechoice_distribution_reconstructed=True,
        future_full_age_income_marginals_verified=False,
        limitation='2011--2023 full marginals were not saved; ledger agreement verifies only their fiscal aggregates. No 2027 prediction or reduced root.',
        historical_equilibrium_certified=False, production_eligible=False, root_coordinates_reduced=False,
        bellman_solves=0, ge_solves=0, elapsed_seconds=time.monotonic()-start)
    save('summary.json', summary)
    save('heartbeat.json', {'phase': 'completed', 'elapsed_seconds': summary['elapsed_seconds']})
    print(json.dumps({k:summary[k] for k in ('status','predicted_pensions_period','maximum_predicted_pension_scaled_budget_residual','maximum_relative_ledger_factor_gap','elapsed_seconds')}), flush=True)
    return 0 if comparison['all_comparisons_pass'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
