"""Editable local stationary-GE runner.  Importing this file is inert."""
from __future__ import annotations

from pathlib import Path

# These literal defaults are the author-adopted October 3 post-interest chain 13.
PARAMETERS = {
    "beta_annual": 0.9663191380998087, "chi": 1.0500762402240174,
    "first_birth_fixed_cost": 0.3045994545418478, "kappa_fert": 0.11652185618155607,
    "kappa_fert_continuation": 0.40070847699255835, "theta0": 0.10097400014250629,
    "h_P": 2.593759507364224, "child_benefit_curvature": 0.0629608477756522,
    "tenure_choice_kappa": 0.014123854930856623, "psi_child": 0.17892072066041628,
}
EXTERNAL_INPUTS = {
    "sigma": 2.0, "alpha_cons": 0.733, "theta1": 0.008193084126995582,
    "R_gross": 1.08243216, "delta": 0.05545379079326218,
    "tau_H": 0.042393443095490375, "psi": 0.06, "phi": [0.8, 0.8, 0.8, 0.8],
    "unsecured_credit_limit": 0.0, "H0": [6.40569359569417], "xi_supply": [0.63],
}
NATIVE_OVERRIDES = {}  # Advanced direct P-field edits.
PRICE_GUESS = 0.7760569760205563
CLOSURE = "fixed_h0"
BUDGET_SECONDS = 1800
OUTPUT_ROOT = Path(__file__).resolve().parents[2] / "output/model/local_solution"


def main():
    from production.workflow import run_stationary
    result, case = run_stationary(PARAMETERS, EXTERNAL_INPUTS, NATIVE_OVERRIDES,
                                  price_guess=PRICE_GUESS, budget_seconds=BUDGET_SECONDS,
                                  closure=CLOSURE, output_root=OUTPUT_ROOT)
    print(f"Completed local stationary GE at price {result.price:.12g}: {case}")
    print(f"Latest cached case: {OUTPUT_ROOT / 'latest'}")
    return case


if __name__ == "__main__":
    main()
