"""Origination-only mortgage vs revolving collateral at the 14.402 base (Oct 8 2026, New York).

Fixed price 0.77941391535061, rebate T held at 0.1921730548652237 (as output/model/rental_menu_precaution_20261007).
Engine: EXPERIMENT COPY tmp/origination_only_20261008/root = byte copy of tmp/rental_menu_precaution_20261007/root plus one
default-off switch P.stayer_no_new_borrowing (kernels.full_owner_block_kernel arg due_no_new_borrowing): in the DUE stayer
branch the floor min(b, -phi p H) becomes min(b, 0). Purchase, resize, sale, renter, death and net-estate floors unchanged.

Cells:
  R80  revolving collateral (current rule), phi 0.80   = benchmark (must reproduce the saved 14.402 case)
  R95  revolving collateral, phi 0.95 (purchase AND stayer ceiling)  = cell A6_LTV95 of the rental-menu packet
  O80  origination-only (stayers b' >= min(b, 0)), phi 0.80 at origination
  O95  origination-only, phi 0.95 at origination
  R80F..O95F  same four cells with the production sale-screen fix B (screen removed)
  R80I, R95I  revolving, with the sale-to-rent screen R b + S >= 0 replaced by R b + S + y >= 0 (DIAGNOSTIC, copied engine)
"""
import copy, hashlib, json, sys, time, traceback
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ENGINE = REPO / 'tmp/origination_only_20261008/root/code/model/experiments/birth_count_choice'
RUN = REPO / 'tmp/rebate_entry_impl_20261006/local_mac/runs/round3_production_chain_0/run'
sys.path.insert(0, str(ENGINE))
LOSS = 14.402423836558775
PRICE_EXPECTED = 0.77941391535061
T_EXPECTED = 0.1921730548652237
CELLS = ['R80', 'R95', 'O80', 'O95', 'R80I', 'R95I', 'R80F', 'R95F', 'O80F', 'O95F']


def write(path, obj):
    def enc(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.floating, np.integer)): return o.item()
        if isinstance(o, (set, tuple)): return list(o)
        return str(o)
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, default=enc) + '\n')


def point():
    c = json.loads((RUN / 'completed.json').read_text()); s = c['selected']
    assert c['status'] == 'selected_numerically_verified' and abs(c['native_loss'] - LOSS) < 1e-12
    price, T = float(s['price']), float(s['closure']['property_tax_rebate']['transfer'])
    assert abs(price - PRICE_EXPECTED) < 1e-12 and abs(T - T_EXPECTED) < 1e-12, (price, T)
    return dict(s['parameters']), float(s['H0_derived']), price, T


def load_base(extra_external=None):
    from model.inputs import load_inputs
    from model.estate_contract import apply_experiment_flags, experiment_flags
    pt, h0, price, T = point()
    P, grid = load_inputs(parameters=pt, external_inputs=dict({'H0': [h0]}, **(extra_external or {})), native_overrides={}, entry_law='levels_rank')
    apply_experiment_flags(P, experiment_flags(1))
    assert P.property_tax_rebate_closure == 'balanced'
    P.property_tax_lump_sum_transfer = T
    return P, grid, price, T


def build(cell):
    P0, grid, price, T = load_base()
    if cell[1:3] == '95':
        P, grid2, _, _ = load_base({'phi': [0.95] * 4}); assert np.array_equal(grid, grid2)
    else:
        P = copy.deepcopy(P0)
    assert bool(P.native_purchase_income) and bool(P.native_due_stayer_credit) and not bool(P.use_pti_constraint)
    assert float(P.unsecured_credit_limit) == 0.0
    if cell[0] == 'O':
        P.stayer_no_new_borrowing = True
    if cell.endswith('F'):
        P.sale_screen_removed = True   # production fix B: sale-to-rent checked by the renter budget only
    if cell.endswith('I'):
        P.sale_screen_includes_income = True   # DIAGNOSTIC: sale-to-rent screen R b + S + y >= 0
    from model.inputs import validate_inputs
    validate_inputs(P, grid)
    spec = dict(cell=cell, price=price, rebate_T_fixed=T, base_loss=LOSS, phi=float(np.asarray(P.phi)[0]),
                stayer_rule="b' >= max(min(b, 0), death floor)" if cell[0] == 'O' else "b' >= max(min(b, -phi p H), death floor)",
                unchanged=['purchase: R b + y >= (1-phi) p H and b\' >= -phi p H', 'resize: fresh origination at -phi p H\', old debt settled from sale',
                           'renters b\' >= 0', 'no amortization requirement', 'fixed price and rent, rebate T held, preferences and parameters as 14.402'])
    return P, grid, price, T, spec


def solve(cell):
    from model.equilibrium import solve_at_price
    from model import fiscal_closure as fiscal
    dest = HERE / 'cells' / cell; dest.mkdir(parents=True, exist_ok=False)
    P, grid, price, T, spec = build(cell)
    write(dest / 'spec.json', spec)
    t0 = time.monotonic()
    out = solve_at_price(P, grid, price); sol, Q, sd = out['solution'], out['P'], out['shared']
    revenue, grants, mass, residual = fiscal.budget_terms(sol)
    write(dest / 'fiscal.json', dict(mode='fixed', transfer=T, revenue=revenue, grant_outlays=grants, household_mass=mass, residual=residual))
    arrays = {k: v for k, v in vars(sol).items() if isinstance(v, np.ndarray) and not v.dtype.hasobject}
    arrays.update({'shared.' + k: v for k, v in vars(sd).items() if isinstance(v, np.ndarray) and not v.dtype.hasobject})
    for k in ['_bp_pol_stay', '_g_stay_distribution']:
        if isinstance(getattr(Q, k, None), np.ndarray): arrays['P.' + k] = getattr(Q, k)
    np.savez_compressed(dest / 'solution_arrays.npz', **arrays)
    write(dest / 'scalar_stats.json', {k: float(v) for k, v in vars(sol).items() if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool)})
    secs = time.monotonic() - t0
    write(dest / 'solve_completed.json', dict(cell=cell, price=price, seconds=secs, childless_rate=float(getattr(sol, 'childless_rate', np.nan)),
                                              own_rate=float(getattr(sol, 'own_rate', np.nan)), mean_age_first_birth=float(getattr(sol, 'mean_age_first_birth', np.nan)),
                                              executed_stayer_no_new_borrowing=bool(getattr(Q, 'stayer_no_new_borrowing', False))))
    print(f'{cell} solved in {secs:.1f} s childless={float(getattr(sol, "childless_rate", np.nan)):.5f} own={float(getattr(sol, "own_rate", np.nan)):.5f}', flush=True)


if __name__ == '__main__':
    for c in sys.argv[1].split(','):
        assert c in CELLS, c
        try: solve(c)
        except Exception as exc:
            d = HERE / 'failed'; d.mkdir(exist_ok=True)
            write(d / f'{c}.json', dict(error=repr(exc), traceback=traceback.format_exc())); print('FAILED', c, repr(exc)[:400], flush=True)
