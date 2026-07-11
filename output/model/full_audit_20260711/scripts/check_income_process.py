#!/usr/bin/env python3
"""Income-process audit checks (no model solve; construction-level only)."""
import json
import math
import numpy as np

from intergen_housing_fertility.local_panel import (
    DEFAULT_INCOME_GRID_5,
    DEFAULT_INCOME_WEIGHTS_5,
    DEFAULT_RHO_Z,
    income_process_overrides,
    income_process_fingerprint,
)
from intergen_housing_fertility.parameters import (
    setup_parameters,
    apply_overrides,
    make_persistent_transition_matrix,
)
from intergen_housing_fertility.production_profile import production_profile_overrides
from intergen_housing_fertility.solver import income_at_state, annual_gross_income_at_state

MATCHED_ANNUAL_RHO = 0.9601845894041878
MATCHED_ANNUAL_INNOVATION_SD = 0.06453733259357768

out = {}

# ---- 1. Rouwenhorst construction --------------------------------------
inc = income_process_overrides(5, "rouwenhorst", MATCHED_ANNUAL_INNOVATION_SD, MATCHED_ANNUAL_RHO)
z, w, Pi = np.asarray(inc["z_grid"]), np.asarray(inc["z_weights"]), np.asarray(inc["Pi_z"])
rho_period = float(inc["income_shock_persistence"])
logz = np.log(z)
log_mean = float(w @ logz)
log_var = float(w @ (logz - log_mean) ** 2)
lev_mean = float(w @ z)
lev_var = float(w @ (z - lev_mean) ** 2)
# stationary dist of Pi vs binomial weights
evals, evecs = np.linalg.eig(Pi.T)
statv = np.real(evecs[:, np.argmin(np.abs(evals - 1.0))])
statv = statv / statv.sum()
# first-order log autocorr
cov1 = float(np.sum(w[:, None] * Pi * (logz[:, None] - log_mean) * (logz[None, :] - log_mean)))
out["rouwenhorst"] = {
    "z_grid": z.tolist(),
    "z_weights": w.tolist(),
    "Pi_row0": Pi[0].tolist(),
    "rho_period_recorded": rho_period,
    "rho_annual^4": MATCHED_ANNUAL_RHO**4,
    "rho_annual_vs_0.85^{1/4}": [MATCHED_ANNUAL_RHO, 0.85 ** 0.25, MATCHED_ANNUAL_RHO - 0.85 ** 0.25],
    "sigma_eps_period_implied": MATCHED_ANNUAL_INNOVATION_SD
    * math.sqrt(sum(MATCHED_ANNUAL_RHO ** (2 * k) for k in range(4))),
    "stationary_log_mean": log_mean,
    "stationary_log_var": log_var,
    "stationary_log_sd": math.sqrt(log_var),
    "stationary_level_mean": lev_mean,
    "stationary_level_sd": math.sqrt(lev_var),
    "binomial_vs_eig_stationary_maxabs": float(np.max(np.abs(w - statv))),
    "invariance_wPi_minus_w_maxabs": float(np.max(np.abs(w @ Pi - w))),
    "row_sums_maxerr": float(np.max(np.abs(Pi.sum(1) - 1))),
    "log_autocorr_1": cov1 / log_var,
    "fingerprint": income_process_fingerprint(inc),
}
# theoretical sigma_log from continuous formula
sig_eps4 = out["rouwenhorst"]["sigma_eps_period_implied"]
sigma_log_theory = sig_eps4 / math.sqrt(1 - (MATCHED_ANNUAL_RHO ** 4) ** 2)
out["rouwenhorst"]["sigma_log_theory"] = sigma_log_theory
out["rouwenhorst"]["sigma_log_theory_sq"] = sigma_log_theory ** 2
# NOTE: discrete grid log-var equals sigma_log^2 * (N-1)/(N-1)?? compute directly:
# For Rouwenhorst the discrete stationary variance of the grid equals sigma_log^2 exactly.

# ---- 2. Old ("current") process ---------------------------------------
inc0 = income_process_overrides(5, "current")
z0, w0, Pi0 = np.asarray(inc0["z_grid"]), np.asarray(inc0["z_weights"]), np.asarray(inc0["Pi_z"])
logz0 = np.log(z0)
lm0 = float(w0 @ logz0)
lv0 = float(w0 @ (logz0 - lm0) ** 2)
levm0 = float(w0 @ z0)
levv0 = float(w0 @ (z0 - levm0) ** 2)
cov10 = float(np.sum(w0[:, None] * Pi0 * (logz0[:, None] - lm0) * (logz0[None, :] - lm0)))
out["old_current"] = {
    "z_grid": z0.tolist(),
    "z_weights": w0.tolist(),
    "stationary_log_mean": lm0,
    "stationary_log_var": lv0,
    "stationary_log_sd": math.sqrt(lv0),
    "stationary_level_mean": levm0,
    "stationary_level_sd": math.sqrt(levv0),
    "log_autocorr_1": cov10 / lv0,
    "rho": DEFAULT_RHO_Z,
}
out["matching_hypothesis"] = {
    "snapshot_claimed_log_var": 0.0533671307381197,
    "old_grid_log_var": lv0,
    "new_process_log_var": log_var,
    "old_minus_new_log_var": lv0 - log_var,
    "rho_annual_equals_0.85^0.25": abs(MATCHED_ANNUAL_RHO - 0.85 ** 0.25) < 1e-12,
    "level_sd_old": math.sqrt(levv0),
    "level_sd_new": math.sqrt(lev_var),
}

# ---- 3. Age profile / pension under the combined spec -----------------
base = {
    "J": 17,
    "Nb": 20,
    "n_house": 5,
}
prof = production_profile_overrides()
overrides = {**base, **prof, **inc, "q": (1.02) ** 4 - 1, "delta": 1 - (1 - 0.011) ** 4}
P = apply_overrides(setup_parameters(), overrides)
if P is not None:
    ages = P.age_start + np.arange(P.J) * P.da
    out["age_profile"] = {
        "ages": ages.tolist(),
        "J_R": int(P.J_R),
        "first_retired_age": float(P.age_start + P.J_R * P.da),
        "income_age_profile": np.asarray(P.income_age_profile).tolist(),
        "work_mean_of_profile": float(np.mean(P.income_age_profile[: P.J_R])),
        "profile_age18_relative": float(P.income_age_profile[0]),
        "raw_value_age18_before_norm": 0.650,
        "expected_norm_ratio_0.650_over_workmean": 0.650 / (np.mean([0.65, 0.65, 0.85, 0.85, 1.0, 1.0, 1.0, 0.985, 0.985, 0.985, 0.935, 0.935])),
        "P_income_row0": np.asarray(P.income[0]).tolist(),
        "pension_per_period": float(P.pension),
        "tau_pay": float(P.tau_pay),
        "period_scale_used": float(P.period_years),
        "income_at_state_j0_zlow": income_at_state(P, 0, 0, float(z[0])),
        "income_at_state_j12_zlow_equals_pension": income_at_state(P, 0, 12, float(z[0])),
        "annual_gross_entry_income_zmid": annual_gross_income_at_state(P, 0, 0, float(z[2])),
        "z_grid_in_P": np.asarray(P.z_grid).tolist(),
        "Pi_z_row0_in_P": np.asarray(P.Pi_z)[0].tolist(),
        "income_shock_persistence_in_P": float(P.income_shock_persistence),
    }

print(json.dumps(out, indent=2))
