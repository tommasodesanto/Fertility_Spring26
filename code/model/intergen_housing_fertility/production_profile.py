"""Named, reproducible production configuration for the July 9 repair cycle."""

from __future__ import annotations

from typing import Any

import numpy as np

from .calibration import PERIOD_YEARS


PRODUCTION_PROFILE_NAME = "intergen_july9_repair_v1"
PRODUCTION_J = 17
PRODUCTION_SEARCH_NB = 120
PRODUCTION_VERIFY_NB = 240
PRODUCTION_H_OWN = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
PRODUCTION_RENTER_CAP = 6.0
PRODUCTION_INCOME_STATES = 5
PRODUCTION_TARGET_SET = "candidate_replacement_post_audit_v1"
PRODUCTION_MAX_ITER_EQ = 25

# These are the source-controlled bounds that generated the reported July runs.
PRODUCTION_SEARCH_BOUNDS: tuple[tuple[str, float, float], ...] = (
    ("beta_annual", 0.940, 0.995),
    ("alpha_cons", 0.400, 0.950),
    ("c_bar_0", 0.020 * PERIOD_YEARS, 0.320 * PERIOD_YEARS),
    ("c_bar_n", 0.050, 1.500),
    ("h_bar_0", 1.000, 6.000),
    ("h_bar_jump", 0.050, 2.500),
    ("h_bar_n", 0.020, 2.000),
    ("psi_child", 0.000, 0.350),
    ("kappa_fert", 1.000, 12.000),
    ("tenure_choice_kappa", 0.000, 0.120),
    ("chi", 0.400, 1.150),
    ("theta0", 0.000, 2.000),
    ("theta_n", 0.000, 1.500),
)

PRELIMINARY_FIXED_BETA_ANNUAL = 0.94
# The local July 9 transcript retained four decimals; the exact remote scratch
# artifact was unavailable when this source profile was created.
PRELIMINARY_FIXED_THETA_N = 0.9811
PRELIMINARY_THETA_N_PRECISION_NOTE = (
    "Four-decimal value transcribed from the July 9 run log; replace only from "
    "a recovered machine-readable task artifact."
)
PRELIMINARY_CHI_RUNGS: tuple[float, ...] = (0.40, 0.55, 0.70, 0.85, 1.00, 1.15)


def production_profile_metadata() -> dict[str, Any]:
    return {
        "name": PRODUCTION_PROFILE_NAME,
        "J": PRODUCTION_J,
        "search_Nb": PRODUCTION_SEARCH_NB,
        "verification_Nb": PRODUCTION_VERIFY_NB,
        "H_own": PRODUCTION_H_OWN.tolist(),
        "renter_cap": PRODUCTION_RENTER_CAP,
        "income_states": PRODUCTION_INCOME_STATES,
        "target_set": PRODUCTION_TARGET_SET,
        "max_iter_eq": PRODUCTION_MAX_ITER_EQ,
        "bounds": [
            {"name": name, "lower": lower, "upper": upper}
            for name, lower, upper in PRODUCTION_SEARCH_BOUNDS
        ],
    }


def production_profile_overrides() -> dict[str, Any]:
    """Return model objects that must not drift with generic run arguments."""

    return {
        "H_own": PRODUCTION_H_OWN.copy(),
        "hR_max": PRODUCTION_RENTER_CAP,
    }


def validate_production_profile(
    profile_name: str,
    *,
    J: int,
    Nb: int,
    n_house: int,
    income_states: int,
    target_set: str,
    stage: str = "search",
) -> None:
    if str(profile_name) != PRODUCTION_PROFILE_NAME:
        raise ValueError(f"unknown production profile: {profile_name}")
    stage_clean = str(stage).strip().lower()
    expected_nb = PRODUCTION_SEARCH_NB if stage_clean == "search" else PRODUCTION_VERIFY_NB
    if stage_clean not in {"search", "verify"}:
        raise ValueError(f"unknown production profile stage: {stage}")
    received = {
        "J": int(J),
        "Nb": int(Nb),
        "n_house": int(n_house),
        "income_states": int(income_states),
        "target_set": str(target_set),
    }
    expected = {
        "J": PRODUCTION_J,
        "Nb": expected_nb,
        "n_house": len(PRODUCTION_H_OWN),
        "income_states": PRODUCTION_INCOME_STATES,
        "target_set": PRODUCTION_TARGET_SET,
    }
    if received != expected:
        raise ValueError(
            f"{PRODUCTION_PROFILE_NAME}/{stage_clean} configuration drift: "
            f"expected {expected}, received {received}"
        )
