"""Complete the historical bridge from a verified revised initial solution.

No household or price solve, supply recalibration, or empirical-target change is
performed here. The caller pins the solution, normalization and empirical inputs.
"""
from __future__ import annotations
import copy
import numpy as np
from e5f_matched_pf_initial_state import NormalizedOldState, CalibrationFiscalContract
from e5f_parenthood_utility import validate_parenthood_utility
from e5f_stationary_paygo import certify_initial_pension
from e5f_social_security import fiscal_accounts
import run_e5f_transition_calibration as calibration
import run_e5f_perfect_foresight_transition as pf


def build_approved_initial_state(*, packet, normalization, outside_origin_entry_share,
                                 preference_change_2023, fertility_tolerance):
    """Reweight the normalized pre-announcement state and retain its supply law.

    The historical preference change is supplied explicitly by a fitting caller;
    it is learned in2007. Queues retain the approved 1/2.1 conversion and four
    waiting vintages. Observed age reweighting never adjusts the supply level.
    """
    if (not np.isfinite(fertility_tolerance) or not 0<fertility_tolerance<=5e-4
            or not np.isfinite(preference_change_2023)
            or not np.isfinite(outside_origin_entry_share) or not 0<outside_origin_entry_share<1):
        raise ValueError('Explicit finite preference, outside-entry share and tight normalization gate required')
    P=copy.deepcopy(packet['parameters'])
    validate_parenthood_utility(P)
    if (P.I!=1 or P.J!=17 or P.Nb!=120 or P.tau_pay!=.179
            or bool(P.joint_nested_choice) or not bool(P.exhaustive_saving_control)
            or not np.array_equal(P.xi_supply,np.array([.63]))):
        raise ValueError('Approved sequential initial economic contract changed')
    fiscal_contract=CalibrationFiscalContract(.01,.04,0.,'retained_calibration_unrebated')
    fiscal_contract.validate(P)
    if (not np.isfinite([normalization['target'],normalization['completed_fertility'],normalization['psi_child']]).all()
            or normalization['target']!=2.1 or abs(normalization['completed_fertility']-2.1)>fertility_tolerance
            or normalization['psi_child']!=float(P.psi_child)):
        raise ValueError('Supplied initial normalization is not verified at this preference intercept')
    sol=packet['solution']; grid=np.asarray(packet['b_grid']).copy()
    from intergen_eqscale_seq_optimized.calibration import extract_moments
    actual_fertility=float(extract_moments(sol,P)['tfr'])
    if (not np.isfinite(actual_fertility) or abs(actual_fertility-2.1)>fertility_tolerance
            or not np.isclose(actual_fertility,normalization['completed_fertility'],rtol=0,atol=1e-12)):
        raise ValueError('Solution fertility differs from the supplied normalization receipt')
    policy=packet['evaluation'].policy; shared=packet['shared']
    pre=np.asarray(packet['stationary_g_pre']).copy()
    pension=certify_initial_pension(pre,P,marginal_tolerance=1e-9,fiscal_tolerance=1e-6)
    supply=packet['supply_rule']
    if supply.mode!='static-elastic' or supply.elasticity!=.63:
        raise ValueError('Explicit approved initial supply law required')
    for price in (float(policy.price[0]),float(policy.price[0])*1.1):
        original=P.H0[0]*(P.user_cost_rate*price/P.r_bar[0])**P.xi_supply[0]
        if not np.isclose(supply.quantity(np.array([price]))[0],original,rtol=1e-12,atol=0):
            raise ValueError('Historical supply law differs from the calibrated initial curve')
    ages=P.age_start+np.arange(P.J)*P.da
    age_mass=pre.sum(axis=(0,1,2,4,5,6))
    survival=np.asarray(P.survival_probs)[:P.J-1] if P.use_age_survival else np.ones(P.J-1)
    _,age_gate=calibration.validated_structural_stationary_age_mass(age_mass,
        entry_flow=float(sol.entry_rate),structural_survival=survival)
    age_reweight=pf.transition.acs_2007_age_reweight_diagnostic(age_mass,ages,
        float(sol.entry_rate),periods=4,period_years=4.,structural_survival=survival)
    g2007=pf.transition.reweight_distribution_to_acs_2007_ages(pre,age_reweight)
    births=calibration.closure.topcode_consistent_renewal_accounting(sol,P)
    conversion=pf.transition.effective_birth_to_household_conversion(2.1)
    renewal=pf.transition.stationary_renewal_from_births(float(sol.entry_rate),
        float(births['topcode_adjusted_birth_children']),outside_origin_entry_share,conversion)
    if abs(renewal['queue_B_over_E']-1)>fertility_tolerance or abs(renewal['identity_residual'])>2e-12:
        raise RuntimeError(f'Approved stationary birth/entry identity failed: {renewal}')
    initial=pf.PFInitialState(g2007,[renewal['queue_mature_flow_B']]*4,
                             [conversion*float(sol.total_births_kfe)]*4)
    conditioning=pf.HistoricalConditioning(2007,float(g2007.sum()),
        {1:2011,2:2015,3:2019,4:2023},renewal['outside_flow_M'],renewal['retention_rho'])
    conditioning.validate(P,4,initial,conversion)
    years=np.array([2007,2011,2015,2019,2023])
    path=float(P.psi_child)+float(preference_change_2023)*np.arange(5)/4.
    diagnostics=dict(schema='e5f_approved_parenthood_initial_state_v1',arm='sequential',
        normalization=copy.deepcopy(normalization),verified_solution_fertility=actual_fertility,stationary_pension=pension,
        structural_age_gate=age_gate,age_reweight=age_reweight,renewal=renewal,
        outside_origin_entry_share=float(outside_origin_entry_share),
        fiscal_contract=vars(fiscal_contract),birth_to_entry_conversion=conversion,
        announcement_state_accounts_at_old_pension=fiscal_accounts(g2007,P),
        announcement_pension_status='Reweighted age distribution requires its own dated balanced pension',
        supply_anchor='Unchanged initial GE supply law; no age-reweight or2023 rebasing',
        old_stationary_supply_elasticity=.63,dated_supply_elasticity=.63,
        preference_change_2023=float(preference_change_2023),
        announcement='Entire historical path known in2007; preference intercept constant after2023',
        historical_equilibrium_verified=False,calibrated_smm=False)
    return NormalizedOldState(P,grid,sol,policy,shared,pre,initial,conditioning,
                              supply,years,path,diagnostics)
