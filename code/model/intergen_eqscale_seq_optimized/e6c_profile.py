"""E6c default-off readiness-arrival profile.

The binary state reuses the otherwise invalid childless ``child_state=1``
cell: state 0 is unsettled and state 1 is settled. First-child entry remains
the existing logit choice but is available only when settled. Readiness
arrives irreversibly according to a logistic age distribution.
"""

from __future__ import annotations

from typing import Any


E6C_PROFILE_NAME = "eqscale_seq_e6c_readiness_arrival_20260727"
E6C_DOMAIN: tuple[tuple[str, float, float, str], ...] = (
    ("readiness_location_age", 8.0, 32.0, "linear"),
    ("readiness_spread_years", 0.5, 10.0, "log"),
)
E6C_SEED = {
    "readiness_location_age": 14.0,
    "readiness_spread_years": 2.0,
}


def e6c_overrides() -> dict[str, Any]:
    """Enable the E6c state without imposing either free hazard coordinate."""
    return {"readiness_gate_enabled": True}


def e6c_metadata() -> dict[str, Any]:
    return {
        "profile": E6C_PROFILE_NAME,
        "state": "binary_unsettled_settled",
        "transition": "irreversible_logistic_age_arrival",
        "first_child_entry": "settled_only",
        "free_coordinates": [name for name, *_ in E6C_DOMAIN],
        "identifying_moments": [
            "mean_age_first_birth",
            "share_first_births_age30plus",
        ],
        "seed": dict(E6C_SEED),
    }
