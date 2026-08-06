"""Empirical entry-target accounting shared by production policy drivers."""

from __future__ import annotations

import math
from types import SimpleNamespace
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


def quota_closure_from_outside_origin_share(
    solution: Any,
    parameters: Any,
    *,
    outside_origin_share: float,
) -> dict[str, float | str]:
    """Construct the fixed-retention, fixed-inflow population closure.

    At baseline scale one, ``E0 = Rbar * B0 + Mbar``.  The empirical
    outside-origin entrant share identifies ``Mbar = s_out * E0`` and the
    accounting residual identifies ``Rbar = (1 - s_out) * E0 / B0``.  Neither
    an outside value nor an entry-taste scale enters this closure.
    """

    share = float(outside_origin_share)
    if not math.isfinite(share) or not 0.0 < share < 1.0:
        raise ValueError("outside-origin entrant share must lie strictly between zero and one")
    local_weight = float(getattr(parameters, "local_birth_entry_weight", 1.0))
    if not math.isclose(local_weight, 1.0, rel_tol=0.0, abs_tol=1.0e-14):
        raise RuntimeError(
            "quota closure requires local_birth_entry_weight=1 so Rbar has the "
            f"retention interpretation; got {local_weight:.12g}"
        )
    entry_flow = float(solution.entry_rate)
    mature_cityborn_flow = float(solution.entrants_mature_total)
    if not math.isfinite(entry_flow) or entry_flow <= 0.0:
        raise RuntimeError(f"baseline entrant flow must be positive and finite, got {entry_flow}")
    if not math.isfinite(mature_cityborn_flow) or mature_cityborn_flow <= 0.0:
        raise RuntimeError(
            "baseline mature city-born flow must be positive and finite, "
            f"got {mature_cityborn_flow}"
        )
    retention_rate = (1.0 - share) * entry_flow / mature_cityborn_flow
    if not math.isfinite(retention_rate) or retention_rate <= 0.0:
        raise RuntimeError(f"quota closure implies invalid Rbar={retention_rate:.12g}")
    if retention_rate > 1.0:
        raise RuntimeError(
            "quota closure is infeasible because Rbar exceeds one: "
            f"Rbar={retention_rate:.12g}, s_out={share:.12g}, "
            f"E0={entry_flow:.12g}, B0={mature_cityborn_flow:.12g}"
        )
    outside_inflow = share * entry_flow
    identity_residual = entry_flow - (
        retention_rate * mature_cityborn_flow + outside_inflow
    )
    recovered_share = outside_inflow / entry_flow
    if abs(identity_residual) > 1.0e-12 or abs(recovered_share - share) > 1.0e-12:
        raise RuntimeError(
            "quota baseline accounting failed: "
            f"residual={identity_residual:.12g}, requested_share={share:.12g}, "
            f"recovered_share={recovered_share:.12g}"
        )
    return {
        "closure_mode": "quota",
        "outside_origin_entrant_share_target": share,
        "retention_rate": float(retention_rate),
        "outside_inflow": float(outside_inflow),
        "baseline_entry_flow": entry_flow,
        "baseline_mature_cityborn_flow": mature_cityborn_flow,
        "baseline_scale_factor": 1.0,
        "baseline_identity_residual": float(identity_residual),
        "outside_origin_entrant_share_identity_check": float(recovered_share),
    }


def quota_population_scale(
    solution: Any,
    parameters: Any,
    quota: dict[str, float | str],
) -> tuple[SimpleNamespace, list[dict[str, Any]]]:
    """Solve ``S E0 = Rbar S B(policy) + Mbar`` without reading values.

    The unused ``parameters`` argument keeps the callable compatible with the
    joint price-transfer solver's existing population-scale interface.
    """

    del parameters
    entry_flow = float(quota["baseline_entry_flow"])
    retention_rate = float(quota["retention_rate"])
    outside_inflow = float(quota["outside_inflow"])
    mature_cityborn_flow = float(solution.entrants_mature_total)
    if not math.isfinite(mature_cityborn_flow) or mature_cityborn_flow <= 0.0:
        raise RuntimeError(
            "policy mature city-born flow must be positive and finite, "
            f"got {mature_cityborn_flow}"
        )
    denominator = entry_flow - retention_rate * mature_cityborn_flow
    finite = bool(
        math.isfinite(denominator)
        and denominator > 1.0e-12
        and math.isfinite(outside_inflow)
        and outside_inflow > 0.0
    )
    if finite:
        scale = outside_inflow / denominator
        implied_entry = scale * entry_flow
        retained_cityborn = retention_rate * scale * mature_cityborn_flow
        entry_residual = implied_entry - (retained_cityborn + outside_inflow)
        outside_share = outside_inflow / implied_entry
    else:
        scale = math.inf
        implied_entry = math.inf
        retained_cityborn = math.inf
        entry_residual = math.nan
        outside_share = math.nan
    total_mass = float(getattr(solution, "total_mass", 1.0))
    return (
        SimpleNamespace(
            finite=finite,
            scale_factor=float(scale),
            denominator=float(denominator),
            implied_entry_total=float(implied_entry),
            implied_retained_mature_cityborn_flow=float(retained_cityborn),
            outside_inflow=float(outside_inflow),
            outside_origin_entrant_share=float(outside_share),
            entry_residual=float(entry_residual),
            implied_total_population=total_mass * float(scale),
        ),
        [],
    )
