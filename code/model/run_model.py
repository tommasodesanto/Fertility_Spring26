"""Editable local stationary-GE runner.  Importing this file is inert."""
from __future__ import annotations

import argparse
from pathlib import Path

# Author-adopted October 3 post-interest chain 13. Edit the ten selected
# coordinates and direct period-native inputs, then run one stationary GE call.
PARAMETERS = {
    "beta_annual": 0.9663191380998087,
    "chi": 1.0500762402240174, "first_birth_fixed_cost": 0.3045994545418478,
    "kappa_fert": 0.11652185618155607, "kappa_fert_continuation": 0.40070847699255835,
    "theta0": 0.10097400014250629, "h_P": 2.593759507364224,
    "child_benefit_curvature": 0.0629608477756522,
    "tenure_choice_kappa": 0.014123854930856623, "psi_child": 0.17892072066041628,
}
EXTERNAL_INPUTS = {
    "sigma": 2.0, "alpha_cons": 0.733, "theta1": 0.008193084126995582,
    "theta_n": 0.0, "R_gross": 1.08243216, "delta": 0.05545379079326218,
    "tau_H": 0.042393443095490375, "psi": 0.06, "phi": [0.8, 0.8, 0.8, 0.8],
    "unsecured_credit_limit": 0.0, "c_min": 0.04, "owner_size_cost": 0.0,
    "owner_size_cost_power": 2.0, "owner_size_cost_ref": 6.0,
    "retirement_income_z_scale": 0.0, "fecundity_omega1": 0.02,
    "fecundity_omega2": 0.134, "property_tax_lump_sum_transfer": 0.0,
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
    "survival_probs": [1.0] * 12 + [0.9391263063710125, 0.9184976343249724,
        0.8849521927812863, 0.8300468061015381],
}
NATIVE_OVERRIDES = {}  # Advanced direct native-field edits.
PRICE_GUESS = 0.7760569760205563
CLOSURE = "fixed_h0"
BUDGET_SECONDS = 1800
OUTPUT_ROOT = Path(__file__).resolve().parents[2] / "output/model/local_solution"


def main(argv=None):
    # Make the executable safe to inspect: argparse handles --help and rejects
    # unsupported arguments before importing the solver/workflow or reserving a case.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args(argv)
    from production.workflow import run_stationary
    result, case = run_stationary(PARAMETERS, EXTERNAL_INPUTS, NATIVE_OVERRIDES,
                                  price_guess=PRICE_GUESS, budget_seconds=BUDGET_SECONDS,
                                  closure=CLOSURE, output_root=OUTPUT_ROOT)
    print("Completed local stationary GE.")
    status = (result.closure.get("status", "converged")
              if isinstance(result.closure, dict) else "converged")
    print(f"Status: {status}; price={result.price:.12g}")
    if isinstance(result.closure, dict):
        for key in ("renewal_residual", "absolute_housing_residual",
                    "actual_paygo_residual", "population_scale",
                    "fixed_h0_population_scale", "implied_H0_at_population_one"):
            if key in result.closure:
                print(f"{key}: {result.closure[key]}")
        if {"fixed_h0_population_scale", "implied_H0_at_population_one"} <= set(result.closure):
            print("Scale interpretations: fixed-H0 population scale and population-one implied H0.")
    print(f"Latest cached case: {OUTPUT_ROOT / 'latest'}")
    return case


if __name__ == "__main__":
    main()
