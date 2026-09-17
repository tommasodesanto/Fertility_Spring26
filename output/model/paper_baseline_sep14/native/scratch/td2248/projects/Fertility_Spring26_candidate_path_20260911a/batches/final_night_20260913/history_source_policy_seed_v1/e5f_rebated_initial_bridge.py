"""Reuse the retained historical bridge after verifying a rebated initial state.

The temporary contract substitution changes only the fiscal assertion in the
old builder. It neither changes parameters nor bypasses its normalization,
pension, supply, grid, renewal, or historical-age checks.
"""
from __future__ import annotations
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np


def build_rebated_initial_state(*, packet, normalization,
                               outside_origin_entry_share,
                               preference_change_2023=0., fertility_tolerance=5e-4):
    import e5f_approved_initial_state as retained
    import e5f_balanced_history as balanced
    from e5f_social_security import fiscal_accounts
    _, primitive, _, _ = balanced._runtime()
    P = packet['parameters']
    e = packet['evaluation']
    mass = float(e.g_current.sum())
    transfer = float(P.property_tax_lump_sum_transfer)
    if not np.isfinite(transfer) or transfer < 0 or mass <= 0:
        raise ValueError('Invalid rebated initial transfer or household mass')
    revenue = float(primitive.model.property_tax_revenue_from_distribution(
        e.g_current, e.policy.hR_pol, e.policy.price, P))
    outlays = transfer * mass
    residual = (revenue-outlays)/max(abs(revenue), abs(outlays), 1e-12)
    accounts = fiscal_accounts(e.g_current, P)
    pension_gap = ((accounts['payroll_tax_revenue']-accounts['pension_outlays']) /
                   max(abs(accounts['payroll_tax_revenue']), abs(accounts['pension_outlays']), 1e-12))
    if not np.isfinite([residual, pension_gap]).all() or max(abs(residual), abs(pension_gap)) > 1e-6:
        raise ValueError('Initial property-tax rebate or PAYGO budget fails')
    class RebatedContract:
        def __init__(self, *args):
            self.annual_property_tax = .01
            self.period_property_tax = .04
            self.transfer = transfer
            self.label = 'equal_rebate_verified_initial'
        def validate(self, Q):
            if (float(Q.period_years) != 4. or not np.isclose(Q.tau_H, .04, rtol=0, atol=1e-15)
                    or float(Q.property_tax_lump_sum_transfer) != transfer):
                raise ValueError('Rebated fiscal contract differs from verified parameters')
    with patch.object(retained, 'CalibrationFiscalContract', RebatedContract):
        old = retained.build_approved_initial_state(
            packet=packet, normalization=normalization,
            outside_origin_entry_share=outside_origin_entry_share,
            preference_change_2023=preference_change_2023,
            fertility_tolerance=fertility_tolerance)
    old.diagnostics['initial_rebate_accounts'] = dict(revenue=revenue, outlays=outlays,
        relative_residual=residual, pension_relative_residual=pension_gap)
    old.diagnostics['announcement'] = 'Each successive preference surprise is expected to persist'
    return old
