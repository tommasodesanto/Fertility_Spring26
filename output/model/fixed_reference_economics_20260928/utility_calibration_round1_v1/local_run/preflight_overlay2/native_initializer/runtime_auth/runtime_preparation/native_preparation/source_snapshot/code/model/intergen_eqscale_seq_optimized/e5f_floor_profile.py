"""Default-off E5F child-room-floor experimental profile.

E5F retains the repaired E5 state space and calibration contract.  It replaces
the two child expenditure-share tilts with one housing-service floor per child
currently at home, leaving nine free parameters against the same twelve target
moments.
"""

from __future__ import annotations

from typing import Any

from .e5_profile import E5_DOMAIN, E5_FIXED, E5_TARGET_SET


E5F_PROFILE_NAME = "eqscale_seq_e5f_child_room_floor_20260805"
E5F_TARGET_SET = E5_TARGET_SET
E5F_DOMAIN: tuple[tuple[str, float, float, str], ...] = tuple(
    row for row in E5_DOMAIN if row[0] not in {"delta_alpha", "delta_alpha_jump"}
) + (("hbar_child_rooms", 0.10, 1.80, "log"),)
E5F_FIXED = {
    **E5_FIXED,
    "delta_alpha": 0.0,
    "delta_alpha_jump": 0.0,
}


def e5f_overrides() -> dict[str, Any]:
    """Activate E5F mechanics; a zero floor remains exactly inert."""
    return {
        "child_room_floor": True,
        "hbar_child_rooms": 0.0,
        "alpha_cons": 0.733,
        "delta_alpha": 0.0,
        "delta_alpha_jump": 0.0,
        "c_bar_0": 0.0,
    }


def e5f_metadata() -> dict[str, Any]:
    return {
        "profile": E5F_PROFILE_NAME,
        "mechanism": "housing_service_floor_per_child_currently_at_home",
        "floor_formula": "hbar(m) = hbar_child_rooms * m",
        "free_parameter_count": 9,
        "target_count": 12,
        "removed_free_parameters": ["delta_alpha", "delta_alpha_jump"],
        "added_free_parameter": "hbar_child_rooms",
        "identifying_room_moments": [
            "housing_increment_0to1",
            "prime30_55_parent_3plus_minus_1to2_mean_rooms",
        ],
        "state_dimension_change": False,
    }
