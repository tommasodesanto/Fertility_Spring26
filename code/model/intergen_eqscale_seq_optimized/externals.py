"""Author-decided external inputs with verified provenance (2026-07-23)."""
from __future__ import annotations

from intergen_eqscale_seq_optimized.local_panel import income_process_overrides

# Floden and Linde (2001, RED 4(2), 406-437), Table IV, United States:
# GMM estimates of the persistent hourly-wage process on PSID 1988-92
# heads. rho = 0.9136 (s.e. 0.0090); innovation variance
# sigma_eps^2 = 0.0426 (s.e. 0.0048). The separately identified
# measurement-error variance (0.0421) is excluded, so this is the clean
# persistent innovation.
FL_RHO_ANNUAL = 0.9136
FL_INNOVATION_VARIANCE_ANNUAL = 0.0426

# Heathcote, Storesletten, and Violante (2017, QJE 132(4), 1693-1754),
# Section II and Figure I: post-government income satisfies
# y_tilde = lambda * y**(1 - tau); benchmark US estimate tau = 0.181
# (s.e. 0.002), OLS in logs, PSID 2000-06 with TAXSIM (CBO-based
# robustness value 0.200). The log-linear tax function preserves the
# AR(1) exactly and scales the innovation s.d. by (1 - tau).
HSV_TAU_US = 0.181

FLHSV_ANNUAL_INNOVATION_SD = (
    FL_INNOVATION_VARIANCE_ANNUAL ** 0.5
) * (1.0 - HSV_TAU_US)


def flhsv_income_overrides(n_states: int = 5) -> dict:
    """After-tax household income process: FL persistence, HSV-scaled risk."""
    return income_process_overrides(
        n_states, "rouwenhorst", FLHSV_ANNUAL_INNOVATION_SD, FL_RHO_ANNUAL
    )
