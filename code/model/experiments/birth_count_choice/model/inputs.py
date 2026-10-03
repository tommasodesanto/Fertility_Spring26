"""Portable selected inputs; no calibration import and no lifecycle work.

Entry wealth retains the provisional nonnegative five-bin-mean law and its
native linear projection onto the finite 120-node wealth grid.
"""
from __future__ import annotations
import copy, hashlib, json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping
import numpy as np
ROOT = Path(__file__).resolve().parent
SNAPSHOT = Path(__file__).resolve().parents[5] / 'code/model/production/reference_inputs'
DEFAULT_PARAMETERS = {
    'beta_annual': 0.9663191380998087, 'chi': 1.0500762402240174,
    'first_birth_fixed_cost': 0.3045994545418478, 'kappa_fert': 0.11652185618155607,
    'kappa_fert_continuation': 0.40070847699255835, 'theta0': 0.10097400014250629,
    'h_P': 2.593759507364224, 'child_benefit_curvature': 0.0629608477756522,
    'tenure_choice_kappa': 0.014123854930856623, 'psi_child': 0.17892072066041628,
}
DEFAULT_PRICE = 0.7760569760205563
# Native period inputs are authoritative: annual equivalents are reporting only.
DEFAULT_INPUTS = {
    'sigma': 2.0, 'alpha_cons': 0.733, 'theta1': 0.008193084126995582,
    'theta_n': 0.0, 'R_gross': 1.08243216, 'delta': 0.05545379079326218,
    'tau_H': 0.042393443095490375, 'psi': 0.06, 'phi': [0.8]*4,
    'unsecured_credit_limit': 0.0, 'c_min': 0.04, 'owner_size_cost': 0.0,
    'owner_size_cost_power': 2.0, 'owner_size_cost_ref': 6.0,
    'retirement_income_z_scale': 0.0, 'fecundity_omega1': 0.02,
    'fecundity_omega2': 0.134, 'property_tax_lump_sum_transfer': 0.0,
    'H0': [6.40569359569417], 'eta_supply': [1.75], 'xi_supply': [0.63],
    'r_bar': [0.16],
    'income': [[2.650830656801071, 2.650830656801071, 3.4664708588937074,
        3.4664708588937074, 4.078201010463186, 4.078201010463186,
        4.078201010463186, 4.017027995306238, 4.017027995306238,
        4.017027995306238, 3.8131179447830785, 3.8131179447830785,
        0.917784047463731, 0.917784047463731, 0.917784047463731,
        0.917784047463731, 0.917784047463731]],
    'survival_probs': [1.0]*12 + [0.9391263063710125, 0.9184976343249724,
        0.8849521927812863, 0.8300468061015381],
}

def _decode(v, arrays):
    if isinstance(v, list): return [_decode(x, arrays) for x in v]
    if not isinstance(v, dict): return v
    if '__ndarray__' in v: return arrays[v['__ndarray__']].copy()
    if '__npscalar__' in v: return np.dtype(v['__npscalar__']).type(v['value'])
    if '__tuple__' in v: return tuple(_decode(x, arrays) for x in v['__tuple__'])
    if '__dict__' in v: return {k:_decode(x, arrays) for k,x in v['__dict__'].items()}
    if '__float__' in v: return float(v['__float__'])
    raise ValueError('invalid reference input encoding')

def _base():
    meta = json.loads((SNAPSHOT/'bundle.json').read_text())
    if hashlib.sha256((SNAPSHOT/'arrays.npz').read_bytes()).hexdigest() != meta['arrays_sha256']:
        raise RuntimeError('production primitive array snapshot drift')
    with np.load(SNAPSHOT/'arrays.npz', allow_pickle=False) as a:
        P=SimpleNamespace(**{k:_decode(v,a) for k,v in meta['parameters'].items()})
        grid=a['b_grid'].copy()
    return P,grid

def bind_parameters(P, grid, parameters):
    """Apply any supplied ten-coordinate edits to a copy, without search bounds."""
    unknown=set(parameters)-set(DEFAULT_PARAMETERS)
    if unknown: raise ValueError('unknown parameter coordinates: '+', '.join(sorted(unknown)))
    Q=copy.deepcopy(P)
    for key,value in parameters.items():
        value=float(value)
        if not np.isfinite(value): raise ValueError('nonfinite parameter: '+key)
        if key=='beta_annual':
            if not 0<value<1: raise ValueError('beta_annual must lie in (0,1)')
            Q.beta=value**float(Q.period_years);Q.rho=1/Q.beta-1;Q.rho_hat=Q.rho
        elif key=='h_P': Q.hbar_first_child_jump=value
        else: setattr(Q,key,value)
    Q.eps_fert=float(Q.kappa_fert)
    validate_inputs(Q,grid)
    return Q

def _assign(P, values, label):
    for key,value in values.items():
        if key.startswith('_') or not hasattr(P,key): raise ValueError(label+' unknown native field: '+key)
        old=getattr(P,key)
        if isinstance(old,np.ndarray):
            new=np.asarray(value,dtype=old.dtype)
            if new.shape!=old.shape: raise ValueError(label+' wrong shape: '+key)
            setattr(P,key,new.copy())
        elif isinstance(old,bool):
            if not isinstance(value,bool): raise ValueError(label+' Boolean required: '+key)
            setattr(P,key,value)
        elif isinstance(old,(float,int,np.number)):
            if not np.isscalar(value) or not np.isfinite(value): raise ValueError(label+' finite scalar required: '+key)
            if isinstance(old,(int,np.integer)) and float(value)!=int(value): raise ValueError(label+' integer required: '+key)
            setattr(P,key,type(old)(value))
        else: setattr(P,key,copy.deepcopy(value))

def validate_inputs(P,grid):
    grid=np.asarray(grid)
    if grid.ndim!=1 or len(grid)!=int(P.Nb) or not np.isfinite(grid).all() or not np.all(np.diff(grid)>0): raise ValueError('invalid wealth grid')
    if np.asarray(P.z_grid).shape!=(int(P.Nz),) or np.asarray(P.z_weights).shape!=(int(P.Nz),): raise ValueError('invalid income grid shape')
    for key,value in vars(P).items():
        if isinstance(value,np.ndarray) and value.dtype.kind in 'fci' and not np.isfinite(value).all(): raise ValueError('nonfinite native field: '+key)
    for name in ['z_weights','survival_probs']:
        v=np.asarray(getattr(P,name))
        if (v<0).any() or (v>1).any(): raise ValueError('invalid probability: '+name)
    if not np.isclose(P.z_weights.sum(),1,rtol=0,atol=2e-12): raise ValueError('income weights must sum to one')
    if np.asarray(P.Pi_z).shape!=(P.Nz,P.Nz) or (P.Pi_z<0).any() or not np.allclose(P.Pi_z.sum(axis=1),1,rtol=0,atol=2e-12): raise ValueError('invalid income transition matrix')
    if np.asarray(P.income).shape!=(len(P.H0),P.J) or np.asarray(P.survival_probs).shape!=(P.J-1,): raise ValueError('invalid lifecycle shapes')
    if (np.asarray(P.phi)<0).any() or (np.asarray(P.phi)>1).any(): raise ValueError('financed shares must be in [0,1]')
    if not np.array_equal(np.asarray(P.phi),np.full_like(np.asarray(P.phi),P.phi[0])): raise ValueError('nonuniform financed shares unsupported by the active native contract')
    if (np.asarray(P.H0)<=0).any() or (np.asarray(P.r_bar)<=0).any() or (np.asarray(P.xi_supply)<0).any(): raise ValueError('invalid supply primitives')
    if P.R_gross<=0 or P.user_cost_rate<=0 or P.sigma<=0 or not 0<P.alpha_cons<1: raise ValueError('invalid return or utility primitives')
    if P.chi<=0 or P.hbar_first_child_jump<0 or P.hR_max<=0: raise ValueError('invalid housing service domain')
    if not 0<=P.child_benefit_curvature<1 or min(P.psi_child,P.first_birth_fixed_cost,P.kappa_fert,P.kappa_fert_continuation,P.tenure_choice_kappa,P.theta0)<0: raise ValueError('invalid preference domain')
    if not 0<=P.delta<1 or not 0<=P.psi<1 or P.tau_H<0: raise ValueError('invalid housing rates')
    C=np.asarray(P.fixed_reference_entry_conditional)
    if C.shape!=(P.Nb,P.Nz) or (C<0).any() or not np.allclose(C.sum(axis=0),1,rtol=0,atol=2e-12): raise ValueError('invalid conditional entry mass')

def load_inputs(parameters=None, external_inputs=None, native_overrides=None):
    """Load pure primitives; reject unsupported fiscal/entry changes explicitly."""
    base,grid=_base();P=bind_parameters(base,grid,{**DEFAULT_PARAMETERS,**(parameters or {})})
    external=dict(external_inputs or {});native=dict(native_overrides or {})
    overlap=set(external)&set(native)
    for k in overlap:
        if not np.array_equal(external[k],native[k]): raise ValueError('contradictory external/native override: '+k)
    values={**external,**native}
    # Structural or entry-law edits require an explicitly constructed new
    # bundle. Gross earnings/payroll changes have a native fixed-tax PAYGO map.
    protected={'J_R','period_years','z_grid','z_weights','Pi_z','entry_wealth_ratio_nodes','entry_wealth_ratio_weights','fixed_reference_entry_conditional','fixed_reference_entry_grid','Nb','Nz'}
    for key in protected & set(values):
        if not np.array_equal(values[key],getattr(base,key)):
            raise ValueError('unsupported entry-law or structural edit: '+key+'; construct a consistent input bundle')
    fiscal_primitives={'w_hat','income_age_profile','tau_pay'}
    fiscal_changed=any(k in values and not np.array_equal(values[k],getattr(base,k)) for k in fiscal_primitives)
    derived_fiscal={'income','pension','pension_by_loc'}
    if not fiscal_changed:
        for key in derived_fiscal & set(values):
            if not np.array_equal(values[key],getattr(base,key)):
                raise ValueError('derived fiscal input '+key+' requires gross earnings or payroll primitives')
    derived={'q','user_cost_rate','rho','rho_hat','eps_fert'}
    if derived & set(values): raise ValueError('edit authoritative primitives rather than derived fields: '+', '.join(sorted(derived & set(values))))
    maps={'beta_annual':'beta','h_P':'hbar_first_child_jump'}
    for key in set(parameters or {}) & (set(values)|set(maps)):
        attr=maps.get(key,key)
        if attr in values and not np.array_equal(values[attr],getattr(P,attr)): raise ValueError('contradictory parameter/native override: '+key)
    _assign(P,{k:v for k,v in values.items() if k not in derived_fiscal},'inputs')
    if fiscal_changed:
        from .engine.e5f_stationary_paygo import bind_initial_balanced_pension
        P,_=bind_initial_balanced_pension(P,payroll_tax=P.tau_pay)
        for key in derived_fiscal & set(values):
            if not np.array_equal(values[key],getattr(P,key)):
                raise ValueError('supplied '+key+' conflicts with native fixed-payroll balanced pension mapping')
    P.q=float(P.R_gross)-1.;P.user_cost_rate=P.q+float(P.delta)+float(P.tau_H)
    P.rho=1./float(P.beta)-1.;P.rho_hat=P.rho;P.eps_fert=float(P.kappa_fert)
    P.pension_by_loc=np.asarray(P.income)[:,int(P.J_R)].copy();P.pension=float(P.pension_by_loc[0])
    validate_inputs(P,grid)
    return P,grid.copy()
