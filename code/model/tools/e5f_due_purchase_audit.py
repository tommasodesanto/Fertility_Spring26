"""Origin-specific DUE audit, extending the authenticated purchase ledger.

Source closure: this file, e5f_overnight_estate_audit.py and the byte-identical
inherited e5f_earnings_wealth_contract.py must all be pinned by the run manifest.
The inherited file's hash is independently checked below before import. Its
transaction accounting runs unchanged on the complete origin distribution;
only its final buyer-debt mass excludes the separately audited owner stayers.
No policies, probabilities, distributions, or model parameters are changed.
"""
from __future__ import annotations
import copy
import hashlib
import importlib.util
from pathlib import Path
import numpy as np
from e5f_overnight_estate_audit import policy_mass_branches

INHERITED_PURCHASE_SHA256 = '71eebc1d42e28ccd304f921535f0d154f3a994451b5c53bcb8097233ce18669a'
SUPPORTS_NATIVE_DUE_STAYER_CREDIT = True
_INHERITED = None


def _inherited_audit():
    global _INHERITED
    if _INHERITED is None:
        source = Path(__file__).with_name('e5f_earnings_wealth_contract.py')
        if hashlib.sha256(source.read_bytes()).hexdigest() != INHERITED_PURCHASE_SHA256:
            raise RuntimeError('Inherited purchase audit differs from its reviewed source pin')
        spec = importlib.util.spec_from_file_location('_due_authenticated_purchase_contract', source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _INHERITED = module.audit_purchase_accounting
    return _INHERITED


def audit_purchase_accounting(evaluation, P, shared, grid, model):
    inherited = _inherited_audit()
    if not bool(getattr(P, 'native_due_stayer_credit', False)):
        return inherited(evaluation, P, shared, grid, model)
    branches = policy_mass_branches(evaluation, P)
    buyers, stayers = branches
    proxy = copy.copy(evaluation)
    proxy.g_current = buyers[0]
    # All transaction checks still see original g_post_fertility and choices.
    receipt = inherited(proxy, P, shared, grid, model)
    g, saving, _ = stayers
    bg = np.asarray(grid, dtype=float)
    if bg.shape != (g.shape[0],) or not np.isfinite(bg).all():
        raise ValueError('DUE stayer audit requires aligned finite wealth grid')
    if (not 0 <= float(P.psi) < 1 or not np.isfinite(P.psi)
            or int(P.J) != g.shape[3]):
        raise ValueError('Invalid DUE liquidation cost or age support')
    if bool(getattr(P, 'use_age_survival', False)):
        survival = np.asarray(P.survival_probs, dtype=float)
        if (survival.shape != (P.J-1,) or not np.isfinite(survival).all()
                or np.any(survival < 0) or np.any(survival > 1)):
            raise ValueError('Invalid DUE age survival schedule')
    else:
        survival = np.ones(P.J-1)
    violation = estate_violation = largest = 0.
    price = float(evaluation.policy.price[0])
    for j in range(P.J):
        death_possible = j == P.J-1 or survival[j] < 1.
        for ten, house in enumerate(P.H_own, start=1):
            collateral = -np.asarray(shared.phi_choice[0,ten]) * price * float(house)
            # Inherited principal b, not capitalized R*b: interest is serviced
            # in the separate dated budget audit, never rolled into this bound.
            principal_floor = np.minimum(bg[:,None,None,None], collateral[None,None,:,:])
            mass = g[:,ten,0,j]
            final = saving[:,ten,0,j]
            shortfall = principal_floor-final
            violation += float(mass[shortfall > 1e-9].sum())
            if death_possible:
                estate_floor = -(1.-float(P.psi))*price*float(house)
                estate_gap = estate_floor-final
                estate_violation += float(mass[estate_gap > 1e-9].sum())
                shortfall = np.maximum(shortfall,estate_gap)
            occupied = mass > 1e-12
            if np.any(occupied): largest=max(largest,float(shortfall[occupied].max()))
    receipt.update(stayer_mass=float(g.sum()),
        stayer_principal_floor_violation_mass=violation,
        stayer_death_estate_violation_mass=estate_violation,
        maximum_occupied_stayer_debt_shortfall=largest,
        source_contract='origin_specific_due_v1')
    if (max(violation,estate_violation)>receipt['mass_tolerance']
            or largest>receipt['wealth_tolerance']):
        raise RuntimeError(f'DUE stayer accounting gate failed: {receipt}')
    return receipt
