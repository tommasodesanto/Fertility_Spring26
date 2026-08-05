"""Author-directed repair of the E5 sequential child-maturation process.

Parity remains the number of children ever born.  The existing child-state
dimension is reinterpreted as the number currently at home.  Each at-home
child leaves independently, with the common hazard chosen so expected time at
home is 18 years.  No additional state dimension or parameter is introduced.
"""

from __future__ import annotations

from typing import Any


E5_MATURATION_REPAIR_NAME = "e5_independent_child_maturation_20260805"


def e5_maturation_repair_overrides() -> dict[str, Any]:
    return {
        "child_state_mode": "independent_count",
        "A_m": 18.0,
        "use_stochastic_aging": True,
    }


def e5_maturation_repair_metadata() -> dict[str, Any]:
    return {
        "profile": E5_MATURATION_REPAIR_NAME,
        "parity_state": "children_ever_born",
        "child_state": "children_currently_at_home",
        "expected_years_at_home_per_child": 18.0,
        "maturation_process": "independent_binomial_thinning",
        "adds_state_dimension": False,
        "adds_free_parameter": False,
    }
