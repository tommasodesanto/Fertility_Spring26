"""One bounded partial-equilibrium exercise at the authenticated soft candidate.

Only house price and its linked user-cost rent change. N=1, housing supply
coefficient, earnings, entry, preferences, fiscal parameters and credit are held
fixed. The normalized lifecycle distribution is not a birth-renewal equilibrium.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import signal
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
PACKETS = ROOT / 'output/model/fixed_reference_economics_20260928'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selection', required=True, type=Path)
    ap.add_argument('--reference-report', required=True, type=Path)
    ap.add_argument('--reference-arrays', required=True, type=Path)
    ap.add_argument('--price-multiplier', required=True, type=float)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    if not .5 <= args.price_multiplier <= 1.5:
        raise ValueError('Diagnostic price multiplier must be in [0.5,1.5]')
    if args.out.exists():
        raise ValueError('Refusing an existing output directory')
    args.out.mkdir(parents=True)
    started = time.time()

    def alarm(*_):
        raise TimeoutError('10-minute single-solve budget')

    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, 600)
    sys.path.insert(0, str(PACKETS / 'normalized_calibration_v2'))
    import run_psi as v2
    # Exact executed package used by the original soft evaluator.
    origin = ROOT / 'output/model/publication_refactor_20260929/small_credit_replication_v1/arms/indexed/source'
    sys.path.insert(0, str(origin))
    from small_credit_lab import credit
    from small_credit_lab.engine import solver, diagnostics
    try:
        v2.native.verify_sources()
        saved = json.loads(args.selection.read_text())
        source = ROOT / saved['source']
        assert hashlib.sha256(source.read_bytes()).hexdigest() == saved['source_sha256']
        point = saved['selected']['parameters']
        assert json.loads(source.read_text())[saved['source_key']]['best'] == saved['selected']
        assert saved['selected']['weight_fingerprint'] == v2.weight_fingerprint({})
        assert v2.native.target_identity(saved['selected']['target_fit']) == v2.CONFIG['base_target_contract']
        _, bounds, _ = v2.inputs.seed_and_bounds('floor_s0')
        bounds = {k: tuple(v) for k, v in bounds.items()}
        assert bounds['h_P'] == (.1, 2.3)
        bounds.update(h_P=(.1, 2.6), psi_child=tuple(v2.CONFIG['psi_bounds']))
        v2.inputs.LANES['floor_s0'].update(seed=dict(point), bounds=bounds, free_coordinates=list(point))
        P, grid = v2.inputs.proposal('floor_s0')
        P, _ = v2.inputs.entry(P, grid, 'nonnegative_mean')
        P = v2.native.utility_checks(P, grid, 'floor_s0', args.out)
        credit.bind_engine_credit(P, 'corrected', 0.)
        closure = json.loads((args.reference_report / 'closure.json').read_text())
        assert closure['normalized_population'] == 1. and P.N_target == 1.
        P.H0 = np.array([closure['H0_derived']])
        price = np.array([closure['price'] * args.price_multiplier])
        P.native_inherited_distribution_evidence_dir = str(args.out / 'inherited_state_failures')
        v2.write(args.out / 'start.json', dict(status='running', deadline_epoch=started+600,
            single_core=True, maximum_lifecycle_solves=1, price=price.tolist(),
            experimental_change=__doc__, fixed_parameters=point,
            reference_arrays_sha256=hashlib.sha256(args.reference_arrays.read_bytes()).hexdigest()))
        sd = solver.precompute_shared(P, grid)
        sol = solver.solve_markov_income_at_prices(price, P, grid, SD=sd, fast_stats=False)
        solve_seconds = time.time()-started
        assert abs(float(sol.total_mass)-1.) < 5e-9
        arrays = {k:v for k,v in vars(sol).items() if isinstance(v,np.ndarray) and v.dtype != object}
        arrays['parameters.H0'] = P.H0
        arrays['normalization.population'] = np.array([1.])
        # A no-change exercise must reproduce every saved policy and distribution
        # array exactly before the price response is trusted.
        comparison = None
        if args.price_multiplier == 1.:
            with np.load(args.reference_arrays, allow_pickle=False) as ref:
                names = ('V','c_pol','bp_pol','hR_pol','c_pol_stay','bp_pol_stay',
                         'tenure_probs','fert_probs','fert2_probs','g','g_beginning_distribution')
                comparison = {k:float(np.max(np.abs(arrays[k]-ref[k]))) for k in names}
            assert max(comparison.values()) < 1e-10, comparison
        for key in ('tenure_probs','fert_probs','fert2_probs'):
            assert np.isfinite(arrays[key]).all() and arrays[key].min() >= 0 and arrays[key].max() <= 1
        path = args.out / 'solution_arrays.npz'
        np.savez_compressed(path, **arrays)
        diagnostics.write_diagnostics(sol, P, args.out / 'standard_diagnostics')
        expected = {p.name for p in (args.reference_report / 'standard_diagnostics').glob('*.png')}
        produced = {p.name for p in (args.out / 'standard_diagnostics').glob('*.png')}
        assert len(expected) == 17 and produced == expected, (expected-produced,produced-expected)
        descriptor = dict(id='soft_price_'+format(args.price_multiplier,'.3g').replace('.','p'),
            label=f'Soft timing · fixed price {100*(args.price_multiplier-1):+.0f}%',
            arrays=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            price=float(price[0]),timing='transaction_inside')
        v2.write(args.out / 'explorer_case.json', descriptor)
        v2.write(args.out / 'completed.json', dict(status='partial_equilibrium_diagnostic',
            price=float(price[0]), price_multiplier=args.price_multiplier, fixed_H0=float(P.H0[0]),
            lifecycle_solves=1, solve_and_setup_seconds=solve_seconds, elapsed_seconds=time.time()-started,
            control_maximum_array_errors=comparison, standard_plot_count=17,
            normalized_population=float(sol.total_mass),
            housing_demand=float(sol.aggregate_housing_demand),housing_supply=float(sol.aggregate_housing_supply),
            housing_excess=float(sol.aggregate_housing_excess),
            equilibrium_claim=False, experimental_not_adopted=True,
            distribution='Normalized lifecycle distribution with fixed entrant law; birth renewal not imposed',
            rent_rule='r=user_cost_rate*p; owner price and rent change together',
            accounting_validation='Source pins, no-change control, probability/mass checks; no full native observer certificate for PE',
            explorer_case=descriptor))
    except BaseException as exc:
        v2.write(args.out / 'failure.json',dict(type=type(exc).__name__,message=str(exc),elapsed_seconds=time.time()-started))
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
