"""Observe one admitted historical date by replaying its saved forecast coordinates.

No root search, calibration, policy run, or later-vintage substitution is performed.
The full original forecast is reproduced, but only its first date is exported as
realized history. The retained finite-horizon/no-rebate limitations remain explicit.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import pickle
import sys
import time

import numpy as np


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path, expected):
    if sha(path) != expected:
        raise ValueError('Changed pinned input: ' + str(path))


def load(path):
    with gzip.open(path, 'rb') as stream:
        return pickle.load(stream)


def save(path, value):
    def serial(item):
        if isinstance(item, np.ndarray):
            return item.tolist()
        if isinstance(item, np.generic):
            return item.item()
        raise TypeError(type(item).__name__)
    Path(path).write_text(json.dumps(value, default=serial, allow_nan=False, indent=2) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--forecast-stage', type=Path, required=True)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--plan-sha256', required=True)
    ap.add_argument('--start-year', type=int, choices=(2007, 2011, 2015), required=True)
    ap.add_argument('--admitted-fit', type=Path, required=True)
    ap.add_argument('--admitted-fit-sha256', required=True)
    ap.add_argument('--inherited-checkpoint', type=Path)
    ap.add_argument('--inherited-sha256')
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    verify(args.plan, args.plan_sha256)
    verify(args.admitted_fit, args.admitted_fit_sha256)
    plan = read(args.plan)
    if plan.get('stationary_restart_2019'):
        raise ValueError('A stationary patch is not carried historical evidence')
    if (args.start_year == 2007) != (args.inherited_checkpoint is None):
        raise ValueError('Only 2007 starts from the approved initial household state')
    if bool(args.inherited_checkpoint) != bool(args.inherited_sha256):
        raise ValueError('Carried checkpoint and its SHA256 must be supplied together')
    admitted_history = read(args.admitted_fit)
    rows = [r for r in admitted_history if r['year'] == args.start_year]
    if len(rows) != 1:
        raise ValueError('The requested year must have one admitted fitted forecast')
    admitted = rows[0]
    if Path(admitted['folder']).resolve() != args.forecast_stage.resolve():
        raise ValueError('Forecast folder is not the admitted historical vintage')
    if admitted['error_abs'] > plan['fertility_fit_tolerance']:
        raise ValueError('Requested forecast was not admitted by the fertility gate')
    for path, digest in plan['file_sha256'].items():
        verify(path, digest)

    # The shared profile helper has historical sys.path setup at module scope.
    # Restore our path before loading the exact plan-pinned economic modules.
    saved_path = list(sys.path)
    from collect_e5f_patch_readout import profile
    sys.path[:] = saved_path
    root = Path(plan['source_root'])
    sys.path[:0] = [str(args.plan.parent), str(root / 'code/model/tools'), str(root / 'code/model')]
    import run_e5f_candidate_terminal as td
    import run_e5f_transition_calibration as fertility
    import e5f_successive_surprises as surprise
    import e5f_balanced_terminal as terminal_module
    import run_e5f_matched_pf_smoke as primitive
    from e5f_approved_initial_state import build_approved_initial_state
    from e5f_initial_fertility_observer import observe_initial_fertility
    surprise_path = str(Path(surprise.__file__).resolve())
    if surprise_path not in plan['file_sha256']:
        raise ValueError('Forecast module is not among the original plan source pins')
    verify(surprise_path, plan['file_sha256'][surprise_path])
    primitive.pf.transition.configure_sequential_model()
    primitive.pf.calendar.apply_fertility = primitive.pf.transition.apply_sequential_fertility
    primitive.pf.calendar.advance_calendar_distribution = primitive.pf.transition.advance_sequential_calendar_distribution
    contract = plan['terminal_template']
    td.validate_contract(contract)
    td.verify_sources(contract)
    for key in ('initial_checkpoint', 'initial_summary', 'initial_contract'):
        td.verify(contract[key]['path'], contract[key]['sha256'])
    seed = load(contract['initial_checkpoint']['path'])
    initial_summary = read(contract['initial_summary']['path'])
    td.validate_initial_receipts(contract, seed, initial_summary, read(contract['initial_contract']['path']))
    old = build_approved_initial_state(packet=seed, normalization=initial_summary['normalization'],
        outside_origin_entry_share=plan['outside_origin_entry_share'],
        preference_change_2023=0., fertility_tolerance=5e-4)
    demographics = seed['demographic_seed']
    if args.start_year == 2007:
        inherited = surprise.InheritedState(2007, old.initial_state)
    else:
        previous = [r for r in admitted_history if r['year'] == args.start_year - 4]
        if len(previous) != 1:
            raise ValueError('Missing admitted predecessor for the inherited state')
        expected_checkpoint = (Path(previous[0]['folder']).parents[1]
            / f'realized_state_{args.start_year}.pkl.gz')
        if args.inherited_checkpoint.resolve() != expected_checkpoint.resolve():
            raise ValueError('Inherited checkpoint is not the admitted predecessor output')
        verify(args.inherited_checkpoint, args.inherited_sha256)
        inherited = load(args.inherited_checkpoint)
        if inherited.year != args.start_year:
            raise ValueError('Carried household state has the wrong observation date')
        # For the first resumed date, retain the original controller's exact pin.
        resumed = plan.get('resume_fitted_prefix', {}).get('checkpoint')
        if args.start_year == 2011 and resumed != dict(path=str(args.inherited_checkpoint), sha256=args.inherited_sha256):
            raise ValueError('2011 checkpoint differs from the original resume plan')
    del seed
    native = args.forecast_stage / 'vintage' / str(args.start_year)
    receipt_path = native / 'root_receipt.json'
    receipt = read(receipt_path)
    if (not receipt.get('finite_horizon_market_fiscal_converged')
            or receipt['start_year'] != args.start_year or receipt['psi'] != admitted['psi']):
        raise ValueError('Native receipt does not certify the admitted finite forecast')
    terminal_dir = args.forecast_stage.parent / 'terminal'
    terminal_path = terminal_dir / 'terminal_state.pkl.gz'
    terminal_receipt_path = terminal_dir / 'root_receipt.json'
    terminal_packet = load(terminal_path)
    terminal_receipt = read(terminal_receipt_path)
    payload = terminal_receipt['final']['payload']
    terminal = terminal_module.BalancedTerminalEndpoint(terminal_packet['parameters'],
        terminal_packet['b_grid'], terminal_packet['policy'], terminal_packet['endpoint'],
        terminal_packet['social_security'], payload['diagnostics'], payload['household_gates'])
    source_paths = [args.plan, args.admitted_fit, receipt_path, native / 'expected_transition.csv',
        native / 'realized_period.json', terminal_path, terminal_receipt_path,
        Path(__file__), Path(observe_initial_fertility.__code__.co_filename), Path(profile.__code__.co_filename)]
    if args.inherited_checkpoint:
        source_paths.append(args.inherited_checkpoint)
    source_hashes = {str(p): sha(p) for p in source_paths}
    args.out.mkdir(parents=True, exist_ok=False)
    observed = {}
    dates = []

    def observer(index, evaluation, parameters, grid, shared):
        if index != len(dates):
            raise RuntimeError('Forecast callback indices are nonconsecutive')
        year = args.start_year + 4 * index
        dates.append(year)
        save(args.out / 'latest_date.json', dict(calendar_year=year))
        if index == 0:
            if parameters.period_years != 4 or parameters.da != 4:
                raise ValueError('The saved observation clock requires four-year periods')
            observed.update(calendar_year=year, profile=profile(evaluation, parameters, grid),
                fertility=fertility.period_fertility_diagnostics(evaluation, parameters),
                fertility_stock_timing=observe_initial_fertility(evaluation, parameters,
                    age_projection='uniform_birth_time'), top_bin_weight=float(parameters.tfr_top_bin_weight))

    started = time.monotonic()
    final = receipt['final']
    path = surprise.evaluate_forecast(inherited=inherited, old_state=old, demographics=demographics,
        prices=final['prices'], pensions=final['fiscal_values'], psi=receipt['psi'],
        terminal=terminal, observer=observer)
    with (native / 'expected_transition.csv').open() as stream:
        prior = list(csv.DictReader(stream))
    if len(path.rows) != len(prior) or dates != [int(r['calendar_year']) for r in prior]:
        raise RuntimeError('Saved forecast dates or row count changed')
    keys = ('asset_price', 'renter_price', 'housing_demand', 'housing_supply', 'owner_rate',
        'birth_children_topcode_adjusted', 'pension_period_units', 'payroll_tax_revenue', 'pension_outlays')
    discrepancies = [abs(float(a[k]) - float(b[k])) for a, b in zip(prior, path.rows) for k in keys]
    maximum = max(discrepancies)
    if not np.isfinite(discrepancies).all() or maximum > 2e-10:
        raise RuntimeError('Original fixed-coordinate path did not reproduce: ' + str(maximum))
    totals = observed['profile']['totals']
    if abs(totals['rooms'] - path.rows[0]['housing_demand']) > 2e-10:
        raise RuntimeError('Profile rooms differ from the realized housing market')
    if not 0 <= totals['capped_rooms'] <= totals['rooms'] + 2e-10:
        raise RuntimeError('Invalid capped housing total')
    if abs(observed['fertility']['period_tfr_topcode_adjusted'] - admitted['model']) > 2e-10:
        raise RuntimeError('Current-date fertility differs from the admitted fit')
    for name, digest in source_hashes.items():
        verify(name, digest)
    verification = dict(status='PASS', calendar_year=args.start_year, replay_maximum_abs=maximum,
        seconds=time.monotonic() - started, forecast_dates=dates, source_sha256=source_hashes,
        inherited_checkpoint_sha256=args.inherited_sha256,
        initial_checkpoint_sha256=contract['initial_checkpoint']['sha256'],
        finite_converged=True, horizon_verified=False, production_eligible=False,
        fiscal_regime='retained no-rebate diagnostic', admitted_fit=admitted,
        interpretation='Only the first date is realized history; later early-vintage expectations are not exported as history')
    save(args.out / 'historical_stock_observation.json', dict(observed, verification=verification))
    save(args.out / 'verification.json', verification)
    print(json.dumps(dict(status='PASS', year=args.start_year, replay_maximum_abs=maximum)))


if __name__ == '__main__':
    main()
