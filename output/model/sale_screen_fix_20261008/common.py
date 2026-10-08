"""Case definitions at the Mac round-3 best 14.402 (chain 0; copy of refit_best_15p07_credit/common.py with RUN and loss changed).

Same case names and definitions as chain 11 (baseline, phi095, phi100, price110, Shapley arms M/F/O..., family-space
worlds), plus ge_phi095 / ge_phi100. Differences from chain 11 are only those of the engine: the balanced property-tax
rebate T is re-solved in every run (so it is endogenous in every arm, including M), and entry wealth follows the
PSID levels law, whose common scale is proportional to entry earnings, so the money arm M scales entry wealth by s
automatically (checked in the audit).
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
import os
ENGINE = REPO / 'tmp/sale_screen_fix_20261008/root/code/model/experiments/birth_count_choice'   # Oct 8 sale-screen fix root
LEGACY = os.environ.get('SALE_SCREEN_LEGACY', '0') == '1'
RUN = REPO / 'tmp/rebate_entry_impl_20261006/local_mac/runs/round3_production_chain_0/run'
sys.path.insert(0, str(ENGINE))
LAMBDA = 1.10
MONEY_SCALARS = ['theta1', 'c_min', 'c_bar_0', 'c_bar_n', 'transfer_floor_G0', 'transfer_floor_Gn',
                 'unsecured_credit_limit', 'estate_tax_exemption', 'birth_entry_grant_amount']
ROOM_SCALARS = ['hR_max', 'h_own_min', 'h_own_max', 'h_bar_0', 'h_bar_jump', 'h_bar_n',
                'owner_size_cost_ref', 'hbar_child_rooms']
SHAPLEY = ['M', 'F', 'O', 'MF', 'MO', 'FO', 'MFO']
WORLDS_A = {'cap4': dict(native=dict(hR_max=4.0)), 'cap20': dict(native=dict(hR_max=20.0)),
            'plus1room': dict(param=dict(h_P=-1.0), native=dict(hbar_child_rooms=1.0)),
            'hP20': dict(param=dict(h_P=2.0)), 'hP15': dict(param=dict(h_P=1.5))}
WORLD_CASES = [f'{w}_p{p}' for w in WORLDS_A for p in ('80', '95')]
CASES = ['ge_baseline', 'baseline', 'phi095', 'phi100', 'price110'] + SHAPLEY + WORLD_CASES + ['ge_phi095', 'ge_phi100']


def point():
    c = json.loads((RUN / 'completed.json').read_text()); s = c['selected']
    assert c['status'] == 'selected_numerically_verified' and abs(c['native_loss'] - 14.402423836558775) < 1e-12
    return dict(s['parameters']), float(s['H0_derived']), float(s['price']), float(s['closure']['property_tax_rebate']['transfer'])


def _load(params, external, native):
    from model.inputs import load_inputs
    from model.estate_contract import apply_experiment_flags, experiment_flags
    P, grid = load_inputs(parameters=params, external_inputs=external, native_overrides=native, entry_law='levels_rank')
    apply_experiment_flags(P, experiment_flags(1))
    assert P.property_tax_rebate_closure == 'balanced'
    return P, grid


def build_inputs(case):
    pt, h0, price0, T = point()
    P0, grid0 = _load(pt, {'H0': [h0]}, {})
    params, external, native = dict(pt), {'H0': [h0]}, {}
    alpha = float(P0.alpha_cons); s, g = LAMBDA ** (-(1 - alpha)), LAMBDA ** alpha
    spec = dict(case=case, s=s, g=g, H0=h0, rebate_start=T, price=price0)
    price = price0
    if case in ('phi095', 'ge_phi095'): external['phi'] = [0.95] * 4
    elif case in ('phi100', 'ge_phi100'): external['phi'] = [1.0] * 4
    elif case == 'price110': price = price0 * LAMBDA
    elif case in SHAPLEY:
        if 'M' in case:
            external['w_hat'] = np.asarray(P0.w_hat, float) * s
            for k in MONEY_SCALARS: external[k] = float(getattr(P0, k)) * s
            params['theta0'] = pt['theta0'] * s ** (float(P0.sigma) - 1.0)
        if 'F' in case:
            params['h_P'] = pt['h_P'] * g; native['hbar_child_rooms'] = float(P0.hbar_child_rooms) * g
        if 'O' in case:
            for k in ROOM_SCALARS:
                if k != 'hbar_child_rooms': native[k] = float(getattr(P0, k)) * g
            native['H_own'] = np.asarray(P0.H_own, float) * g
    elif case in WORLD_CASES:
        w, p = case.rsplit('_p', 1); d = WORLDS_A[w]
        native.update(d.get('native', {}))
        if 'param' in d: params['h_P'] = pt['h_P'] - 1.0 if w == 'plus1room' else d['param']['h_P']
        if p == '95': external['phi'] = [0.95] * 4
    elif case not in ('baseline', 'ge_baseline'):
        raise ValueError(case)
    P, grid = _load(params, external, native)
    P.property_tax_lump_sum_transfer = T; P0.property_tax_lump_sum_transfer = T
    P.legacy_pre_income_sale_screen = LEGACY; spec['legacy_pre_income_sale_screen'] = LEGACY
    spec.update(price=price, external={k: np.asarray(v).tolist() for k, v in external.items() if k != 'H0'},
                native={k: np.asarray(v).tolist() for k, v in native.items()},
                params_changed={k: [pt[k], params[k]] for k in params if params[k] != pt[k]})
    return P, grid, price, spec, P0, grid0


def changed_fields(P0, P):
    out = []
    for k in sorted(set(vars(P0)) | set(vars(P))):
        a, b = getattr(P0, k, None), getattr(P, k, None)
        try:
            if not np.array_equal(np.asarray(a), np.asarray(b)): out.append(k)
        except Exception:
            if str(a) != str(b): out.append(k)
    return out


def entry_mean(P, grid):
    C = np.asarray(P.fixed_reference_entry_conditional); zw = np.asarray(P.z_weights, float); zw = zw / zw.sum()
    return float((zw * (np.asarray(grid)[:, None] * C).sum(0)).sum())


def write(path, obj):
    def enc(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, (np.floating, np.integer)): return o.item()
        if isinstance(o, (set, tuple)): return list(o)
        return str(o)
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, default=enc) + '\n')
