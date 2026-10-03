"""Adopted October 3 post-interest chain 13; verified old-target continuation.

One stationary solution uses these supplied inputs; it does not search for a
better parameter vector. Copy this file or edit toy_params.py for experiments.
"""
# PARAMETERS: ten calibration coordinates. beta_annual is annual; costs and
# housing taste/fertility controls use the retained native model units.
PARAMETERS = {
    "beta_annual": 0.9663191380998087,  # Annual patience; the loader compounds to four years.
    "chi": 1.0500762402240174,  # Owner housing-service premium in preferences.
    "first_birth_fixed_cost": 0.3045994545418478,  # Utility/value cost conditional on a successful first birth.
    "kappa_fert": 0.11652185618155607,  # First-birth choice dispersion.
    "kappa_fert_continuation": 0.40070847699255835,  # Later-birth choice dispersion.
    "theta0": 0.10097400014250629,  # Bequest motive strength.
    "h_P": 2.593759507364224,  # Parenthood housing-floor jump, in rooms.
    "child_benefit_curvature": 0.0629608477756522,  # Curvature of child utility.
    "tenure_choice_kappa": 0.014123854930856623,  # Tenure choice dispersion.
    "psi_child": 0.17892072066041628,  # Level of utility from children.
}

# EXTERNAL_INPUTS: period-native primitives (one period is four years).
# phi is the financed share under the retained soft purchase rule; .8 gives
# a (1-phi)=.2 down-payment threshold, not a strict liquid-cash requirement.
# All four entries must match. unsecured_credit_limit is renter debt capacity in native wealth units.
# R_gross is gross period interest; delta and tau_H are period rates; psi is
# selling cost. H0 is the physical housing supply coefficient, not quantity.
# w_hat/income_age_profile are gross-earnings controls; tau_pay is the fixed
# payroll tax. Disposable income and a balanced pension are derived from them.
# Housing is in rooms; entry wealth/income processes and the finite grid retain
# the authenticated snapshot. Unsupported structural edits fail explicitly.
EXTERNAL_INPUTS = {
    # Preferences, equivalence-scale/bequest controls.
    "sigma": 2.0, "alpha_cons": 0.733, "theta1": 0.008193084126995582,
    "theta_n": 0.0,
    # Finance and four-year interest/depreciation/property-tax rates.
    "R_gross": 1.08243216, "delta": 0.05545379079326218,
    "tau_H": 0.042393443095490375, "psi": 0.06, "phi": [0.8, 0.8, 0.8, 0.8],
    # Credit capacity, in native wealth units.
    "unsecured_credit_limit": 0.0, "c_min": 0.04, "owner_size_cost": 0.0,
    "owner_size_cost_power": 2.0, "owner_size_cost_ref": 6.0,
    "retirement_income_z_scale": 0.0, "fecundity_omega1": 0.02,
    "fecundity_omega2": 0.134, "property_tax_lump_sum_transfer": 0.0,
    # Housing supply: physical coefficient, normalization and elasticity.
    "H0": [6.40569359569417], "eta_supply": [1.75], "xi_supply": [0.63],
    "r_bar": [0.16],
    # Gross earnings plus fixed payroll tax determine disposable income and
    # the balanced initial pension. The entrant wealth law stays fixed.
    "w_hat": [1.0], "income_age_profile": [
        0.720554272517321, 0.720554272517321, 0.9422632794457273,
        0.9422632794457273, 1.1085450346420322, 1.1085450346420322,
        1.1085450346420322, 1.0919168591224018, 1.0919168591224018,
        1.0919168591224018, 1.0364896073903, 1.0364896073903,
        0.720554272517321, 0.720554272517321, 0.720554272517321,
        0.720554272517321, 0.720554272517321],
    "tau_pay": 0.08028070961950022,
    # Four-year survival probabilities between successive age cells.
    "survival_probs": [1.0] * 12 + [0.9391263063710125, 0.9184976343249724,
        0.8849521927812863, 0.8300468061015381],
}

NATIVE_OVERRIDES = {}

PRICE_GUESS = 0.7760569760205563

CLOSURE = "fixed_h0"

BUDGET_SECONDS = 1800

PROVENANCE = {
    "role": "author_adopted_working_continuation",
    "reference": "post-interest soft chain 13; retained old wealth target",
    "verified_receipt": "output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/production_alternative_chain_13/run/completed.json",
    "native_loss": 13.771131463467462,
    'reference_receipt_sha256': '12880d38b335bab36c23f06f1272771d04db485bbe26bbc98c80b3acdb37cbdb',
    'input_snapshot_sha256': '26ab79ace7df6b3429997ef43ba45598988abeafd662f82596ff16c83b35e448',
    'input_arrays_sha256': '3877c55f683248d4b81352c04558b1cb93fb150e330bf805caa73ff72db57396',
    'target_fingerprint': 'db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1',
    'weight_fingerprint': '2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0',
    "not_certified": "global optimum, grid adequacy, transition or paper baseline",
}
