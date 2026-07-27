"""Default-off E6a fecundity-tail external for the signed E5 target system."""

from __future__ import annotations

import math
from typing import Any

E6A_PROFILE_NAME = "eqscale_seq_e6a_fecundity_tail_20260727"
E6A_FECUNDITY_OMEGA1 = 0.01331
E6A_FECUNDITY_OMEGA2 = 0.14960
E6A_FECUNDITY_TAIL_START_AGE = 40.0
E6A_FECUNDITY_TERMINAL_AGE = 45.0
E6A_TERMINAL_SUCCESS_PROBABILITY = 0.03


def terminal_decay() -> float:
    """Decay needed to reach the terminal success probability just before 45."""
    base_terminal = 1.0 - E6A_FECUNDITY_OMEGA1 * math.exp(
        E6A_FECUNDITY_OMEGA2 * (E6A_FECUNDITY_TERMINAL_AGE - 18.0)
    )
    return math.log(base_terminal / E6A_TERMINAL_SUCCESS_PROBABILITY) / (
        E6A_FECUNDITY_TERMINAL_AGE - E6A_FECUNDITY_TAIL_START_AGE
    )


def e6a_overrides() -> dict[str, Any]:
    """Return the E6a external schedule without changing the free parameters."""
    return {
        "fecundity_omega1": E6A_FECUNDITY_OMEGA1,
        "fecundity_omega2": E6A_FECUNDITY_OMEGA2,
        "fecundity_tail_start_age": E6A_FECUNDITY_TAIL_START_AGE,
        "fecundity_terminal_decay": terminal_decay(),
        "fecundity_terminal_age": E6A_FECUNDITY_TERMINAL_AGE,
    }


def e6a_metadata() -> dict[str, Any]:
    """Serializable evidence and implementation record for collectors."""
    return {
        "profile": E6A_PROFILE_NAME,
        "omega1": E6A_FECUNDITY_OMEGA1,
        "omega2": E6A_FECUNDITY_OMEGA2,
        "tail_start_age": E6A_FECUNDITY_TAIL_START_AGE,
        "terminal_decay": terminal_decay(),
        "terminal_age": E6A_FECUNDITY_TERMINAL_AGE,
        "terminal_success_probability_before_hard_close": (
            E6A_TERMINAL_SUCCESS_PROBABILITY
        ),
        "identification": (
            "External demographic schedule: four-year conception anchors at "
            "30/35/40 and the late-age sterility tail; no added free parameter."
        ),
    }
