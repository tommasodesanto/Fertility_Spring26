"""One isolated stationary GE with the joint intended-birth menu; no recalibration."""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

PROJECT = Path(__file__).resolve().parents[4]
OUTPUT = PROJECT / 'output/model/experiments/birth_count_choice/current_params_v1'


def selected_inputs():
    from model.parameter_files import load_parameter_file
    from model.inputs import load_inputs
    config = load_parameter_file(PROJECT / 'code/model/parameters/best_params.py')
    P, grid = load_inputs(config['parameters'], config['external_inputs'], config['native_overrides'])
    if P.H0.tolist() != [6.40569359569417] or config['price_guess'] != 0.7760569760205563:
        raise RuntimeError('Canonical current-parameter anchor changed; review before experiment')
    return config, P, grid


def preflight():
    """Authenticate all source/target contracts and adapters without any solve."""
    import json
    import tempfile
    import time
    import numpy as np
    from model.calibration import _contract, effective_input_fingerprint
    from model.reporting import build_context, production_model_facade
    config, P, grid = selected_inputs()
    before = effective_input_fingerprint(P, grid)
    P.birth_count_choice_enabled = True
    target, target_pin, weight_pin, _ = _contract()
    if target_pin != config['provenance']['target_fingerprint'] or weight_pin != config['provenance']['weight_fingerprint']:
        raise RuntimeError('Parameter file target/weight provenance drift')
    probe = dict(vars(P)); probe.pop('birth_count_choice_enabled')
    from types import SimpleNamespace
    if effective_input_fingerprint(SimpleNamespace(**probe), grid) != before:
        raise RuntimeError('Experiment changed an input other than count-choice flag')
    with tempfile.TemporaryDirectory(prefix='birth_count_preflight_') as directory:
        context = build_context(P, grid, directory, price_start=config['price_guess'],
                                deadline=time.time() + 120, max_lifecycle=32, closure='fixed_h0')
        adapters = context['birth_count_observer_adapters']
        if len(context['manifest']['full_parameter_table']) != 31 or len(target) != 14:
            raise RuntimeError('Full fit/parameter table dimensions changed')
        facade = production_model_facade()
        from model.engine import solver
        if facade.solve_markov_income_at_prices is not solver.solve_markov_income_at_prices:
            raise RuntimeError('Facade accidentally binds a different engine')
        if '/birth_count_choice/model/engine/' not in facade.solve_markov_income_at_prices.__code__.co_filename:
            raise RuntimeError('Facade does not bind the copied engine')
    return dict(status='passed_zero_solves', lifecycle_solves=0,
                canonical_input_fields=len(probe), copied_grid_nodes=len(grid),
                input_fingerprint_before_flag=before, changed_input_fields=['birth_count_choice_enabled'],
                target_rows=len(target), parameter_rows=31,
                target_fingerprint=target_pin, weight_fingerprint=weight_pin,
                budget_seconds=900, max_lifecycle=32, output_root=str(OUTPUT),
                observer_functions_adapted=[r['function'] for r in adapters])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight', action='store_true', help='Authenticate inputs and reporting; zero lifecycle solves')
    args = parser.parse_args(argv)
    # --help exits before importing numerical or reporting dependencies.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    if args.preflight:
        import json
        print(json.dumps(preflight(), indent=2))
        return
    from model.workflow import run_stationary
    config, _, _ = selected_inputs()
    result, case = run_stationary(config['parameters'], config['external_inputs'], config['native_overrides'],
        price_guess=config['price_guess'], budget_seconds=900, max_lifecycle=32,
        closure='fixed_h0', output_root=OUTPUT,
        parameter_file_metadata={key: config[key] for key in ('config_source','config_sha256','config_text','provenance')})
    from model.calibration import residual_from_report
    residual, fits, params = residual_from_report(result.report_directory)
    print(f'Converged experimental stationary GE: {case}')
    print(f'price={result.price:.12g}; loss={float(residual @ residual):.12g}; {len(fits)} fit rows; {len(params)} parameter rows')
    return case


if __name__ == '__main__':
    main()
