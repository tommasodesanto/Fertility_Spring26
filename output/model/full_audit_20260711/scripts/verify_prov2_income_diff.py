"""PROV-2 verifier: quantify the difference between the two income processes
that the merge-order flip selects between, without any model solve.

Process A (dirty tree / overnight, income wins): Rouwenhorst(5) with
rho_annual=0.9601845894041878, sd_annual=0.06453733259357768, PERIOD_YEARS=4.
Process B (HEAD order, profile wins): PRODUCTION_Z_GRID=[0.6..1.4],
weights [.1,.2,.4,.2,.1], mixture persistence 0.85.
"""
import numpy as np

from intergen_housing_fertility.local_panel import income_process_overrides
from intergen_housing_fertility.production_profile import production_profile_overrides


def stats(name, z_grid, z_weights, Pi_z, rho_label):
    z_grid = np.asarray(z_grid, float)
    z_weights = np.asarray(z_weights, float)
    Pi_z = np.asarray(Pi_z, float)
    mean = z_weights @ z_grid
    log_sd = np.sqrt(z_weights @ (np.log(z_grid) - z_weights @ np.log(z_grid)) ** 2)
    # first-order autocorrelation of log z under (weights, Pi)
    lz = np.log(z_grid)
    mu = z_weights @ lz
    var = z_weights @ (lz - mu) ** 2
    cov = 0.0
    for i in range(len(z_grid)):
        for j in range(len(z_grid)):
            cov += z_weights[i] * Pi_z[i, j] * (lz[i] - mu) * (lz[j] - mu)
    rho_eff = cov / var
    print(f"{name}: mean={mean:.4f} log-sd={log_sd:.4f} grid=[{z_grid.min():.3f},{z_grid.max():.3f}] "
          f"period-autocorr(logz)={rho_eff:.4f} ({rho_label})")
    print("  z_grid   :", np.round(z_grid, 4))
    print("  z_weights:", np.round(z_weights, 4))
    return log_sd, rho_eff


inc = income_process_overrides(5, "rouwenhorst", 0.06453733259357768, 0.9601845894041878)
prof = production_profile_overrides()

sd_a, rho_a = stats("A Rouwenhorst (overnight actual)", inc["z_grid"], inc["z_weights"], inc["Pi_z"],
                    f"rho_period declared={inc['income_shock_persistence']:.4f}")
sd_b, rho_b = stats("B production profile (HEAD order)", prof["z_grid"], prof["z_weights"], prof["Pi_z"],
                    f"mixture persistence={prof['income_shock_persistence']}")

print(f"\nlog-sd ratio A/B = {sd_a/sd_b:.3f}; autocorr A vs B = {rho_a:.3f} vs {rho_b:.3f}")
print("key overlap (income dict keys also set by profile):")
overlap = sorted(set(inc) & set(prof))
print(" ", overlap)
