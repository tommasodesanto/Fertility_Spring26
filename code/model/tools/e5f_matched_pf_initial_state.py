"""Rebuild an arm's normalized old state using the retained calibration recipe.

This module does not fit targets or solve the anticipated historical transition.
It supplies the inherited 2007 household distribution and the original 2007
housing-supply anchor. A PF caller must use the entire returned preference path
in backward induction, including at its first date; the old stationary policy
is a pre-announcement object, not the household's 2007 transition policy.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
import math
from types import SimpleNamespace
from typing import Any, Callable
import time

import numpy as np

import run_e5f_matched_pf_smoke as primitive
import run_e5f_perfect_foresight_transition as pf
import run_e5f_transition_calibration as calibration


@dataclass(frozen=True)
class CalibrationFiscalContract:
    annual_property_tax: float
    period_property_tax: float
    transfer: float
    closure: str

    def validate(self, P):
        # This is the historical calibration fiscal specification. A rebated
        # policy baseline needs its own joint transfer/price solve downstream.
        if (self.closure != 'retained_calibration_unrebated'
                or not math.isclose(self.annual_property_tax, .01, rel_tol=0, abs_tol=1e-15)
                or not math.isclose(self.period_property_tax, .04, rel_tol=0, abs_tol=1e-15)
                or self.transfer != 0.):
            raise ValueError('Old initialization requires explicit retained 1% annual tax, 4% four-year tax, zero-transfer calibration contract')
        if (float(P.period_years) != 4. or not math.isclose(float(P.tau_H), self.period_property_tax, rel_tol=0, abs_tol=1e-15)
                or float(getattr(P, 'property_tax_lump_sum_transfer', 0.)) != self.transfer):
            raise ValueError('Input parameters differ from the explicit calibration fiscal contract')


@dataclass
class NormalizedOldState:
    parameters: Any
    b_grid: np.ndarray
    solution: Any
    policy: Any
    shared: Any
    stationary_g_pre: np.ndarray
    initial_state: pf.PFInitialState
    historical_conditioning: pf.HistoricalConditioning
    supply_rule: Any
    years: np.ndarray
    psi_path: np.ndarray
    diagnostics: dict[str, Any]


def historical_preference_path(old_psi, selected_summary):
    """Retain the estimated 2007-to-2023 change when renormalizing the level."""
    selected_old = float(selected_summary['old_psi_child'])
    selected_new = float(selected_summary['best_candidate']['new_psi_child'])
    change = selected_new - selected_old
    if not all(math.isfinite(x) for x in (old_psi, selected_old, selected_new)):
        raise ValueError('Preference endpoints must be finite')
    path = np.array([pf.transition.preference_shifter_at_date(i, old_psi, old_psi + change, 4)
                     for i in range(5)], dtype=float)
    return np.array([2007, 2011, 2015, 2019, 2023]), path, change


def verify_selected_supply_anchor(rule, selected_summary):
    """Check that the input packet still refers to its original 2007 anchor."""
    saved = selected_summary['housing_supply']
    expected = dict(initial_price=saved['retained_date0_asset_price'],
                    initial_stock=saved['date0_normalized_housing_stock'],
                    elasticity=saved['transition_housing_supply_elasticity'])
    for name, target in expected.items():
        if not math.isclose(float(getattr(rule, name)), float(target), rel_tol=0, abs_tol=1e-12):
            raise ValueError(f'Selected supply anchor mismatch: {name}')
    if rule.mode != 'static-elastic' or not math.isclose(float(rule.elasticity), .63, rel_tol=0, abs_tol=1e-15):
        raise ValueError('Expected original static-elastic rule with externally fixed elasticity .63')
    return expected


def initialize_normalized_old_state(
    *, parameters, b_grid, selected_summary, selected_supply_rule, arm,
    fiscal_contract: CalibrationFiscalContract,
    completed_fertility_tolerance: float, max_stationary_solves: int,
    deadline_monotonic: float, progress: Callable[[dict[str, Any]], None] | None = None,
) -> NormalizedOldState:
    """Normalize completed fertility to 2.1 without reconstructing P from defaults.

    The caller owns source/input hashing, checkpointing, and the hard process
    watchdog. ``deadline_monotonic`` is checked between stationary solves;
    ``max_stationary_solves`` caps calls to the existing scalar normalizer.
    No empirical targets, weights, estimated parameter coordinates, numerical
    tolerances, or inherited household-state primitives are changed here.
    """
    if arm not in ('sequential', 'nested'):
        raise ValueError('Specify sequential or nested arm')
    if (not math.isfinite(completed_fertility_tolerance)
            or not 0 < completed_fertility_tolerance <= 5e-4):
        raise ValueError('Explicit old-fertility tolerance must be no looser than retained 5e-4')
    if type(max_stationary_solves) is not int or not 1 <= max_stationary_solves <= 23:
        raise ValueError('Explicit stationary solve cap must be 1–23')
    if not math.isfinite(deadline_monotonic) or deadline_monotonic <= time.monotonic():
        raise ValueError('A future monotonic deadline is required')
    if progress is not None and not callable(progress):
        raise ValueError('Progress must be callable')
    fiscal_contract.validate(parameters)
    primitive.verify_selected_parameters(parameters, selected_summary, selected_supply_rule)
    original_anchor = verify_selected_supply_anchor(selected_supply_rule, selected_summary)
    grid = np.asarray(b_grid, dtype=float).copy()
    if (grid.ndim != 1 or len(grid) != 120 or not np.isfinite(grid).all()
            or np.any(np.diff(grid) <= 0) or int(parameters.Nb) != len(grid)
            or int(parameters.J) != 17 or int(parameters.I) != 1):
        raise ValueError('Initializer requires the exact 120-node, 17-age, one-market input profile')
    chain, model = pf.transition.configure_sequential_model()
    if not np.array_equal(grid, model.make_grid(parameters)):
        raise ValueError('Input wealth grid differs from the selected parameter grid')
    pf.calendar.apply_fertility = pf.transition.apply_sequential_fertility
    pf.calendar.advance_calendar_distribution = pf.transition.advance_sequential_calendar_distribution
    base = copy.deepcopy(parameters)
    base.joint_nested_choice = arm == 'nested'
    base.fertility_nest_choice = arm == 'nested'
    base.two_shock_choice = False
    base.exhaustive_saving_control = True
    if not bool(getattr(base, 'normalize_transition_mass_roundoff', False)):
        raise ValueError('Expected retained transition roundoff normalization; initializer does not change it')
    if min(float(base.kappa_fert), float(base.kappa_fert_continuation)) < .005:
        raise ValueError('Common parameters violate nested GEV scale restriction')
    # Preserve the original old-GE supply primitives (including xi_supply).
    # The externally fixed .63 is applied to the dated rule after the same
    # normalized-old-state/ACS2007 anchor construction as in calibration.
    retained_elasticity = float(selected_summary['housing_supply']['retained_housing_supply_elasticity'])
    if not np.allclose(np.asarray(base.xi_supply), retained_elasticity, rtol=0, atol=1e-15):
        raise ValueError('Original stationary supply primitive changed')
    if not 0 < float(base.tol_eq) <= 2.5e-5:
        raise ValueError('Original old-GE tolerance must be no looser than 2.5e-5')
    seed_price = float(original_anchor['initial_price'])
    records = []
    def solve_with_explicit_parameters(overrides, verbose=False):
        if set(overrides) != {'psi_child'}:
            raise ValueError('Old normalizer attempted an unauthorized primitive override')
        if len(records) >= max_stationary_solves or time.monotonic() >= deadline_monotonic:
            raise TimeoutError('Old-normalization solve/time budget exhausted')
        P = copy.deepcopy(base)
        P.psi_child = float(overrides['psi_child'])
        record = dict(index=len(records), psi_child=P.psi_child, status='started')
        records.append(record)
        if progress:
            progress(dict(record))
        started = time.monotonic()
        sol, P, price = model.solve_markov_income_equilibrium(np.array([seed_price]), P, grid, verbose=False)
        if not bool(sol.converged) or not bool(sol.timings.get('strict_converged', False)):
            raise RuntimeError('Old-state GE normalization solve failed its unchanged strict gate')
        record.update(status='complete', elapsed_seconds=time.monotonic() - started,
                      asset_price=float(price[0]), market_error=float(sol.timings['best_eq_error']))
        if progress:
            progress(dict(record))
        return sol, P, price
    adapter = SimpleNamespace(run_model_cp_dt=solve_with_explicit_parameters,
                              extract_moments=chain.extract_moments)
    sol, P, price, _, normalization = calibration.solve_old_steady_state(
        adapter, {}, initial_psi=float(selected_summary['old_psi_child']),
        completed_fertility_target=2.1,
        completed_fertility_tolerance=completed_fertility_tolerance, normalize=True)
    fiscal_contract.validate(P)
    shared = model.precompute_shared(P, grid)
    P._fert2_probs = np.asarray(sol.fert2_probs, dtype=float).copy()
    policy = pf.calendar.policy_from_solution(sol, price, P, grid, shared)
    stationary, reconstruction = pf.calendar.reconstruct_stationary_pre_fertility(sol, policy, P, grid, shared)
    operators = pf.transition.operator_gates(sol, policy, stationary, P, grid, shared)
    operators.update(reconstruction)
    keys = ('stationary_post_fertility_nesting_l1', 'one_step_constant_path_nesting_l1',
            'mature_flow_abs_error', 'birth_flow_abs_error', 'topcode_adjusted_birth_flow_abs_error')
    if any(not math.isfinite(float(operators[k])) or abs(float(operators[k])) > 5e-9 for k in keys):
        raise RuntimeError(f'Normalized stationary operator gate failed: {operators}')
    if abs(float(operators['zero_entry_mass_accounting_residual'])) > 2e-8 or float(operators['stationary_feasibility_projection_mass']) > 1e-6:
        raise RuntimeError(f'Normalized stationary mass/feasibility gate failed: {operators}')
    old_renewal = calibration.closure.topcode_consistent_renewal_accounting(sol, P)
    conversion = pf.transition.effective_birth_to_household_conversion(2.1)
    outside_share = float(selected_summary['outside_origin_entry_share'])
    renewal = pf.transition.stationary_renewal_from_births(float(sol.entry_rate),
        float(old_renewal['topcode_adjusted_birth_children']), outside_share, conversion)
    if abs(renewal['queue_B_over_E'] - 1.) > completed_fertility_tolerance or abs(renewal['identity_residual']) > 2e-12:
        raise RuntimeError(f'Normalized old renewal gate failed: {renewal}')
    ages = P.age_start + np.arange(P.J) * P.da
    age_mass = stationary.sum(axis=(0, 1, 2, 4, 5, 6))
    survival = np.asarray(P.survival_probs)[:P.J - 1] if bool(P.use_age_survival) else np.ones(P.J - 1)
    _, structural_gate = calibration.validated_structural_stationary_age_mass(
        age_mass, entry_flow=float(sol.entry_rate), structural_survival=survival)
    age_reweight = pf.transition.acs_2007_age_reweight_diagnostic(age_mass, ages,
        float(sol.entry_rate), periods=4, period_years=4., structural_survival=survival)
    initial_g = pf.transition.reweight_distribution_to_acs_2007_ages(stationary, age_reweight)
    supply, supply_normalization = pf.calendar.normalize_date0_housing_supply(
        initial_g, policy, P, grid, shared, 'static-elastic')
    supply = pf.calendar.HousingSupplyRule(supply.mode, supply.initial_price, supply.initial_stock, .63)
    supply_normalization.update(retained_housing_supply_elasticity=retained_elasticity,
        transition_housing_supply_elasticity=.63, elasticity_status='externally_fixed_profile_not_estimated')
    slots = int(selected_summary['renewal_accounting_contract']['birth_vintage_queue_waiting_slots'])
    if slots != 4:
        raise ValueError('Retained historical queue requires four waiting slots')
    initial = pf.PFInitialState(initial_g, [renewal['queue_mature_flow_B']] * slots,
        [conversion * float(sol.total_births_kfe)] * slots)
    conditioning = pf.HistoricalConditioning(2007, float(initial_g.sum()),
        {1: 2011, 2: 2015, 3: 2019, 4: 2023}, renewal['outside_flow_M'], renewal['retention_rho'])
    conditioning.validate(P, 4, initial, conversion)
    years, psi_path, preference_change = historical_preference_path(float(normalization['psi_child']), selected_summary)
    diagnostics = dict(arm=arm, normalization=normalization, stationary_solves=records,
        choice_flags={k: getattr(P, k) for k in ('joint_nested_choice', 'fertility_nest_choice',
            'two_shock_choice', 'exhaustive_saving_control', 'normalize_transition_mass_roundoff')},
        stationary_operator_gates=operators, structural_age_gate=structural_gate,
        renewal=renewal, supply_normalization=supply_normalization, selected_reference_anchor=original_anchor,
        old_stationary_supply_elasticity=retained_elasticity, dated_supply_elasticity=.63,
        fiscal_contract=vars(fiscal_contract), preference_change_2023=preference_change,
        announcement='All historical preference and equilibrium paths are known at 2007; PF caller must solve date0 using future continuation',
        supply_anchor='Arm-specific normalized old policy after ACS2007 age reweight; never rebased at2023',
        production_policy_fiscal_baseline='Outstanding separate rebated-price/transfer equilibrium; initializer uses retained calibration fiscal contract',
        historical_equilibrium_solved=False)
    return NormalizedOldState(P, grid, sol, policy, shared, stationary, initial, conditioning,
                              supply, years, psi_path, diagnostics)
