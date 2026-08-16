"""Small E5F repair: measured permanent income and a first-birth intercept.

The child-room-floor household structure is unchanged.  The income process
uses the already measured E6B permanent component, while one nonnegative
first-birth cost, paid when parenthood begins, separates the extensive margin from continuation
fertility.  The income distribution is externally fixed, so only the
first-birth cost adds a calibrated parameter.
"""

from __future__ import annotations

from typing import Any

from .e5f_floor_profile import E5F_DOMAIN, E5F_FIXED, E5F_TARGET_SET
from .e6b_profile import e6b_metadata, e6b_overrides


E5F_INCOME_ENTRY_PROFILE_NAME = "eqscale_seq_e5f_income_entry_20260816"
E5F_INCOME_ENTRY_TARGET_SET = E5F_TARGET_SET
E5F_INCOME_ENTRY_DOMAIN: tuple[tuple[str, float, float, str], ...] = (
    E5F_DOMAIN
    + (("first_birth_fixed_cost", 0.0, 8.0, "softzero"),)
)
E5F_INCOME_ENTRY_FIXED = dict(E5F_FIXED)


def e5f_income_entry_overrides() -> dict[str, Any]:
    """Return the external income process and a nesting value for the new cost."""
    return {
        **e6b_overrides(),
        "first_birth_fixed_cost": 0.0,
    }


def e5f_income_entry_metadata() -> dict[str, Any]:
    """Describe the two repair objects and their empirical roles."""
    return {
        "profile": E5F_INCOME_ENTRY_PROFILE_NAME,
        "base_profile": "eqscale_seq_e5f_child_room_floor_20260805",
        "free_parameter_count": len(E5F_INCOME_ENTRY_DOMAIN),
        "target_count": 12,
        "added_free_parameter": "first_birth_fixed_cost",
        "first_birth_cost_identification": "completed childlessness, jointly with first-birth timing",
        "permanent_income": e6b_metadata(),
        "state_dimension_change": "income states expand from 5 to 15; no new lifecycle state",
    }
