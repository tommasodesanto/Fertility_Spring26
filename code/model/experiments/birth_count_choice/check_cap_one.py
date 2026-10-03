"""Explicit one-core cap-one fixed-price regression; this command does solve once."""
from __future__ import annotations
import argparse
from pathlib import Path
import os
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, help='Canonical completed baseline case; default canonical latest')
    parser.add_argument('--out', type=Path, required=True, help='Fresh isolated verification directory')
    args = parser.parse_args(argv)
    for key in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        os.environ[key] = '1'
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import json
    import signal
    import numpy as np
    import numba
    from run import PROJECT, selected_inputs
    from model.storage import load_case
    from model.equilibrium import solve_at_price
    from model.calibration import effective_input_fingerprint
    baseline, directory = load_case(args.baseline or PROJECT / 'output/model/local_solution/latest')
    config, P, grid = selected_inputs()
    if baseline.parameters != config['parameters']:
        raise RuntimeError('Cached baseline ten-coordinate inputs differ')
    if not np.array_equal(grid, baseline.b_grid):
        raise RuntimeError('Cached baseline grid differs')
    # Require a compatible canonical full input contract, not just its label.
    contract = json.loads((Path(directory) / 'input_contract.json').read_text())
    if contract['external_inputs'] != config['external_inputs'] or contract['native_overrides'] != config['native_overrides']:
        raise RuntimeError('Cached baseline fixed inputs differ')
    args.out.mkdir(parents=True, exist_ok=False)
    input_pin = effective_input_fingerprint(P, grid)
    P.birth_count_choice_enabled = True; P.birth_count_choice_cap = 1
    numba.set_num_threads(1)
    def timeout(*unused): raise TimeoutError('Cap-one 300-second lifecycle budget reached')
    prior = signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL,300)
    try:
        outcome = solve_at_price(P,grid,baseline.price)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0); signal.signal(signal.SIGALRM,prior)
    compared = {}
    for name, old in vars(baseline.solution).items():
        if not isinstance(old,np.ndarray) or old.dtype.hasobject: continue
        new = getattr(outcome['solution'],name,None)
        if new is None or np.shape(new) != old.shape:
            raise RuntimeError('Cap-one baseline array absent/shape differs: ' + name)
        max_error = float(np.max(np.abs(np.asarray(new,dtype=float)-np.asarray(old,dtype=float)))) if old.size else 0.
        equal = bool(np.allclose(new,old,atol=2e-12,rtol=2e-12,equal_nan=True))
        compared[name] = dict(matches_tolerance=equal,max_abs_error=max_error)
    descriptive = ['attempt_hazard_by_age','first_birth_hazard_by_age',
                   'fert_by_age','first_birth_age_distribution']
    required_core = ['g','g_beginning_distribution','g_stay_distribution',
                     'aggregate_housing_demand','rental_demand_by_market','owner_demand_by_size']
    core_checks = {}
    for name in required_core:
        old_value = getattr(baseline.solution, name, None)
        new_value = getattr(outcome['solution'], name, None)
        if old_value is None or new_value is None:
            raise RuntimeError('Required cap-one core comparison absent: '+name)
        error = float(np.max(np.abs(np.asarray(new_value)-np.asarray(old_value))))
        core_checks[name] = dict(max_abs_error=error, passes_5e9_tolerance=error <= 5e-9)
    probability_checks = {}
    for name in ('tenure_probs','loc_probs','fert_probs','fert2_probs'):
        value = np.asarray(getattr(outcome['solution'], name))
        probability_checks[name] = dict(minimum=float(value.min()),maximum=float(value.max()),
            finite=bool(np.isfinite(value).all()),
            within_roundoff_bounds=bool(np.isfinite(value).all() and value.min() >= -2e-12 and value.max() <= 1+2e-12))
    core_passed = all(row['passes_5e9_tolerance'] for row in core_checks.values())
    probability_passed = all(row['within_roundoff_bounds'] for row in probability_checks.values())
    receipt = dict(status='core_nesting_passed_raw_differences_disclosed' if core_passed and probability_passed else 'failed',
                   core_nesting_checks=core_checks,probability_roundoff_checks=probability_checks,
                   raw_array_comparison_count=len(compared),full_raw_array_equivalence_claimed=False,
                   revised_descriptive_arrays=descriptive,
                   revised_definition='pre-birth rather than post-birth exposure; active target observers separate',
                   derived_descriptive_scalars=['mean_age_first_birth','share_first_births_age30plus',
                                               'childless_chosen_45','childless_clock_45'],
                   lifecycle_solves=1, fixed_price_partial_equilibrium=True,
                   baseline_case=str(directory),price=float(baseline.price),
                   canonical_input_fingerprint_before_flags=input_pin,
                   only_new_flags={'birth_count_choice_enabled':True,'birth_count_choice_cap':1},
                   comparison_tolerance={'atol':2e-12,'rtol':2e-12}, arrays=compared)
    (args.out/'cap_one_comparison.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='arrays'},indent=2))
    if receipt['status'] == 'failed': raise RuntimeError('Cap-one core nesting/probability check fails; inspect receipt')


if __name__ == '__main__': main()
