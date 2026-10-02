"""Extracted from code/model/intergen_eqscale_seq_optimized/child_preferences.py (sha256 69544858b0137e09d37ac9e07ee6470ddc2d5849236b2c3e33d270737164f09d).

Mechanical copy by refactor_lab/materialize.py: only reachable top-level
definitions, bodies byte-identical; import edits listed in the receipt.
"""
from __future__ import annotations
import math
import numpy as np


def apply_child_preferences(P, alpha, benefit, material_multiplier):
    """Apply the declared specification before native family-type compression.

    P.psi_child stores b in b*m**(1-kappa), where m is children at home.
    With compensated shares, material utility is CRRA of A(m)*Q/e(m),
    Q=c**alpha(m)*s**(1-alpha(m)). A=K(alpha0,r*)/K(alpha(m),r*) uses
    a fixed reference rent, never the current equilibrium rent.
    Absent both options, no array or floating-point operation is changed.
    """
    curvature = float(getattr(P, "child_benefit_curvature", 0.0))
    compensated = bool(getattr(P, "compensated_child_housing_shares", False))
    if not math.isfinite(curvature) or not 0 <= curvature < 1:
        raise ValueError("Child-benefit curvature must be finite and in [0,1)")
    if curvature == 0 and not compensated:
        return
    if str(getattr(P, "child_state_mode", "")) != "independent_count":
        raise ValueError("Native child preferences require independent children-at-home states")
    if getattr(P, "utility_comparison_arm", None) is not None:
        raise ValueError("Choose native child preferences or the legacy comparison adapter, not both")
    if curvature != 0:
        for n in range(int(P.n_parity)):
            for m in range(1, min(n + 1, int(P.n_child_states))):
                benefit[n, m] = P.psi_child * float(m) ** (1.0 - curvature)
    if not compensated:
        return
    if (str(getattr(P, "preference_spec", "")) != "eqscale"
            or str(getattr(P, "eqscale_form", "")) != "power"
            or bool(getattr(P, "child_room_floor", False))
            or float(getattr(P, "hbar_first_child_jump", 0.0)) != 0
            or float(getattr(P, "hbar_child_rooms", 0.0)) != 0
            or float(P.delta_alpha) != 0):
        raise ValueError("Compensated first-child shares require power equivalence scale, no floor and no later loading")
    base_alpha = float(P.alpha_cons)
    loading = float(P.delta_alpha_jump)
    reference_rent = float(getattr(P, "utility_reference_rent", math.nan))
    if (not math.isfinite(reference_rent) or reference_rent <= 0
            or not math.isfinite(base_alpha) or not math.isfinite(loading)
            or loading < 0 or not .05 <= base_alpha - loading <= base_alpha <= .95):
        raise ValueError("Explicit positive reference rent and unclipped interior housing shares required")
    log_k = alpha * np.log(alpha) + (1 - alpha) * np.log((1 - alpha) / reference_rent)
    log_k0 = (base_alpha * math.log(base_alpha)
              + (1 - base_alpha) * math.log((1 - base_alpha) / reference_rent))
    factors = np.where(alpha == base_alpha, 1.0, np.exp(log_k0 - log_k))
    material_multiplier *= factors ** (1.0 - float(P.sigma))
