"""Run the predeclared one-case full-history regression using existing checks."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
import run_e5f_bounded_calibration_refinement as checks
import collect_e5f_transition_calibration as collector

base = Path(__file__).resolve().parent
contract = json.loads((base / 'history_contract.json').read_text())
argv = contract['argv']
center = Path(argv[argv.index('--panel-center-json') + 1])
assert checks.digest(center) == contract['expected_center_sha256']
assert checks.digest(contract['reference_summary']) == contract['reference_summary_sha256']
out = Path(argv[argv.index('--outdir') + 1])
assert not out.exists()
subprocess.run([sys.executable] + argv, check=True)
s = checks.read_json(out / 'summary.json')
b = s['best_candidate']
expected = SimpleNamespace(
    expected_source_sha256=checks.SOURCE,
    expected_target_set='e5_fullhistory_roomsfix_h1_20260817',
    expected_target_fingerprint=checks.TARGET,
    expected_code_bundle_sha256=argv[argv.index('--expected-code-bundle-sha256') + 1],
    expected_model_profile='e5f-income-entry', expected_panel_seed=2026090506,
    expected_center_sha256=contract['expected_center_sha256'],
    expected_panel_design='mixed', expected_local_radius=.005,
    expected_housing_supply_elasticity=.63, expected_tenure_choice_kappa=.005)
collector.validate_expected_contract(s, s['panel_design'], expected)
collector.validate_renewal_accounting(s)
collector.validate_calibration_scope(s)
assert s['target_count'] == 12 and s['transition_free_parameter_count'] == 11
moments = collector.validate_target_fit(checks.read_csv(out / 'target_fit_long.csv'), s, b)
collector.validate_dated_housing_ledger(out, s, b, moments)
collector.population_bridge_contract(s['population_bridge'])
for value, limit in [(b['max_market_residual'], 2e-4),
    (b['max_mass_residual'], 2e-10), (abs(b['population_target_gap']), 2e-10),
    (s['stationary_measurement_max_abs_gap'], 2e-8),
    (abs(b['terminal_synthetic_childless_consistency_gap']), 2e-10)]:
    assert value <= limit
receipt = checks.compare_reference(out, Path(contract['reference_summary']))
assert checks.read_csv(out / 'parameter_table.csv') == checks.read_csv(
    Path(contract['reference_summary']).parent / 'parameter_table.csv')
receipt.update(status='pass', all_parameters_and_bounds_exact=True,
    code_bundle_sha256=expected.expected_code_bundle_sha256,
    summary_sha256=checks.digest(out / 'summary.json'),
    contract_sha256=checks.digest(base / 'history_contract.json'),
    checker_sha256=checks.digest(__file__), no_calibration_search=True)
checks.write_json(base / 'history_receipt.json', receipt)
print('EXACT_HISTORY_PASS', json.dumps(receipt), flush=True)
