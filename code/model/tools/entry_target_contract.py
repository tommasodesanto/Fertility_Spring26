"""Empirical entry-target accounting shared by production policy drivers."""

from __future__ import annotations

import math
from typing import Any


def entry_target_from_outside_origin_share(
    solution: Any,
    parameters: Any,
    *,
    outside_origin_share: float,
) -> dict[str, float]:
    """Map an outside-origin entrant share into the model entry probability.

    At normalized baseline scale, total entrants satisfy
    ``E0 = qstar * (M + B0)``. The outside-origin share is
    ``s_out = qstar * M / E0 = 1 - qstar * B0 / E0``. Hence
    ``qstar = (1 - s_out) * E0 / B0``. The calculation is candidate-specific
    because both baseline entrant flow ``E0`` and mature city-born flow ``B0``
    are model outcomes.
    """

    share = float(outside_origin_share)
    if not math.isfinite(share) or not 0.0 < share < 1.0:
        raise ValueError("outside-origin entrant share must lie strictly between zero and one")
    entry_flow = float(solution.entry_rate)
    local_weight = float(getattr(parameters, "local_birth_entry_weight", 1.0))
    mature_cityborn_flow = local_weight * float(solution.entrants_mature_total)
    if not math.isfinite(entry_flow) or entry_flow <= 0.0:
        raise RuntimeError(f"baseline entrant flow must be positive and finite, got {entry_flow}")
    if not math.isfinite(mature_cityborn_flow) or mature_cityborn_flow <= 0.0:
        raise RuntimeError(
            "baseline retained mature city-born flow must be positive and finite, "
            f"got {mature_cityborn_flow}"
        )
    target_qbar = (1.0 - share) * entry_flow / mature_cityborn_flow
    if not 0.0 < target_qbar < 1.0:
        raise RuntimeError(
            "empirical outside-origin share implies an infeasible model entry probability: "
            f"qstar={target_qbar:.12g}, s_out={share:.12g}, "
            f"E0={entry_flow:.12g}, B0={mature_cityborn_flow:.12g}"
        )
    recovered_share = 1.0 - target_qbar * mature_cityborn_flow / entry_flow
    if abs(recovered_share - share) > 1.0e-12:
        raise RuntimeError(
            "entry-target accounting failed: "
            f"requested={share:.12g}, recovered={recovered_share:.12g}"
        )
    return {
        "outside_origin_entrant_share_target": share,
        "target_qbar": float(target_qbar),
        "baseline_entry_flow": entry_flow,
        "baseline_retained_mature_cityborn_flow": mature_cityborn_flow,
        "baseline_mature_cityborn_to_entry_ratio": mature_cityborn_flow / entry_flow,
        "outside_origin_entrant_share_identity_check": float(recovered_share),
    }
