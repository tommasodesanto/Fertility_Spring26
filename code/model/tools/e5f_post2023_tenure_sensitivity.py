"""Fail-closed post-2023 tenure-taste sensitivity helpers.

The dated E5f history and common 2023 household state are reconstructed under
the selected calibration before this override is applied.  The override is
therefore a conditional post-2023 sensitivity, not a recalibrated E5f model.
"""

from __future__ import annotations

import math
from typing import Any


def apply_post2023_tenure_choice_kappa(
    prepared: Any, requested_value: float | None
) -> dict[str, Any]:
    """Apply and describe a default-off post-2023 tenure-logit scale."""
    baseline = float(getattr(prepared.parameters, "tenure_choice_kappa", 0.0))
    if requested_value is None:
        return {
            "active": False,
            "scope": "none",
            "baseline_tenure_choice_kappa": baseline,
            "post2023_tenure_choice_kappa": baseline,
        }
    value = float(requested_value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(
            "post-2023 tenure_choice_kappa must be finite and strictly positive"
        )
    prepared.base_overrides["tenure_choice_kappa"] = value
    prepared.parameters.tenure_choice_kappa = value
    if float(prepared.base_overrides["tenure_choice_kappa"]) != value:
        raise RuntimeError("Failed to apply tenure_choice_kappa to base overrides")
    if float(prepared.parameters.tenure_choice_kappa) != value:
        raise RuntimeError("Failed to apply tenure_choice_kappa to parameters")
    return {
        "active": True,
        "scope": "post-2023 conditional sensitivity with common reconstructed 2023 state",
        "baseline_tenure_choice_kappa": baseline,
        "post2023_tenure_choice_kappa": value,
    }
