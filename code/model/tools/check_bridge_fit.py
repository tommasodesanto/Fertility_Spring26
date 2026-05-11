"""Solve the model at the slide bridge benchmark theta and compare moments
to the live calibration targets. No geography inversion, no SMM search."""

from __future__ import annotations

import numpy as np

from dt_cp_model.parameters import build_calibration_setup
from dt_cp_model.solver import run_model_cp_dt
from dt_cp_model.theta import apply_theta
from dt_cp_model.objective import extract_moments


BRIDGE_THETA = np.array([
    0.935,    # beta
    0.1226,   # b_entry_fixed
    0.0692,   # psi_child
    0.9000,   # h_bar_jump
    0.9084,   # h_bar_n
    0.0871,   # c_bar_n
    2.3,      # kappa_fert
    1.0464,   # chi
    1.5,      # kappa_loc
    5.0,      # mu_move
    0.53,     # theta0
    0.25,     # theta_n
    2.59,     # h_bar_0
])

# Slide-bench geography wedges and housing grids (extracted from
# repair_candidate_anchor_lowfert_sprint_2026_04_19_job_006_strictscreen_strict_eq_localbaseline.mat)
BRIDGE_OVERRIDES = {
    "E_loc": np.array([0.0, 0.1318]),
    "r_bar": np.array([0.04, 0.0762]),
    "H_own": np.array([5.3725, 5.84, 6.3075, 6.875, 7.625, 8.68]),
    "hR_max": 5.1,
}

# Reference moments straight from the saved MAT moments struct (MATLAB output).
MATLAB_MOMENTS = {
    "tfr": 1.8492,
    "childless_rate": 0.1724,
    "mean_age_first_birth": 34.7381,
    "tfr_gradient": 0.3267,
    "own_rate": 0.8119,            # MATLAB stores own_rate_3055 here
    "own_gradient": 0.2276,
    "own_family_gap": 0.2704,
    "prime_childless_renter_median_rooms": 5.10,
    "prime_childless_owner_median_rooms": 7.625,
    "housing_increment_0to1": 0.4056,
    "housing_increment_1to2": 0.4318,
    "young_liquid_wealth_to_income": 0.4872,
    "center_share_nonparents": 0.4197,
    "center_share_newparents": 0.3598,
    "migration_rate": 0.0253,
    "old_age_own_rate": 0.8537,
    "old_age_parent_childless_gap": 0.0806,
}


def main() -> None:
    setup = build_calibration_setup("benchmark")
    P = apply_theta(setup.P_base, BRIDGE_THETA, setup.names)
    for k, v in BRIDGE_OVERRIDES.items():
        setattr(P, k, v)
    P.max_iter_eq = 200

    sol, P_out, p_eq = run_model_cp_dt(P, verbose=False)
    moments = extract_moments(sol, P_out, p_eq, 0.0, 0.0, 0.0, True)

    targets = setup.targets
    weights = setup.weights

    print(f"{'moment':40s} {'target':>10s} {'matlab':>10s} {'python':>10s} {'py-target':>10s} {'py-matlab':>10s} {'wt*err^2':>10s}")
    print("-" * 105)
    total_loss = 0.0
    for name, tgt in targets.items():
        py = float(getattr(moments, name, np.nan))
        ml = MATLAB_MOMENTS.get(name, np.nan)
        denom = abs(tgt) if abs(tgt) > 1e-6 else 1.0
        dev = (py - tgt) / denom
        contrib = weights[name] * dev * dev
        total_loss += contrib
        print(f"{name:40s} {tgt:10.4f} {ml:10.4f} {py:10.4f} {py - tgt:+10.4f} {py - ml:+10.4f} {contrib:10.4f}")
    print("-" * 105)
    print(f"{'TOTAL SMM LOSS (no geo penalty)':40s} {'':>10s} {'':>10s} {'':>10s} {'':>10s} {'':>10s} {total_loss:10.4f}")
    print()
    print(f"Python p_eq = {p_eq}")
    print(f"Python pop_share = {sol.pop_share}")


if __name__ == "__main__":
    main()
