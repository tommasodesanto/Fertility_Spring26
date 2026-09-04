"""Diagnostic E5 target profile with young ownership restored.

This profile is deliberately separate from :mod:`e5_profile`: selecting it is
an explicit sensitivity exercise and cannot alter the production default.  It
appends the clean household-head homeownership rate at ages 25--34 to the live
E5 dated-transition target system.  The row uses the same provisional
five-percent-of-target calibration scale as other E5 rows without a measured
sampling standard error; that scale is not represented as an empirical SE.
"""
from __future__ import annotations

from typing import Any

from .e5_profile import E5_TARGETS, E5_TARGET_SET, E5_WEIGHTS
from .target_system import TargetSystem


E5_YOUNG_OWNERSHIP_PROFILE = "young-ownership-overid-v1"
E5_BASELINE_TARGET_PROFILE = "baseline"
E5_YOUNG_OWNERSHIP_TARGET_SET = (
    f"{E5_TARGET_SET}_plus_young_ownership_v1"
)
E5_YOUNG_OWNERSHIP_TARGET = 0.34116609
E5_YOUNG_OWNERSHIP_CALIBRATION_SCALE = 0.05 * E5_YOUNG_OWNERSHIP_TARGET
E5_YOUNG_OWNERSHIP_WEIGHT = 1.0 / E5_YOUNG_OWNERSHIP_CALIBRATION_SCALE**2
E5_YOUNG_OWNERSHIP_PROVENANCE: dict[str, Any] = {
    "moment": "own_rate_2534",
    "estimate": E5_YOUNG_OWNERSHIP_TARGET,
    "authoritative_builder": (
        "code/data/mms_center_periphery/audit_ownership_targets.R"
    ),
    "authoritative_audit": "docs/model/intergen_target_object_audit.md",
    "sample": (
        "ACS household heads in the DUE-housing sample, ages 25--34; the "
        "household-head object matches the model's independent-household unit"
    ),
    "estimator": "HHWT-weighted homeownership rate",
    "fixed_effects": "none",
    "clustering": "none; the objective uses a declared calibration scale",
    "measurement_definition": (
        "owner mass divided by total household mass for model ages 25--34"
    ),
    "uncertainty": (
        "no promoted empirical sampling SE; diagnostic weight uses five percent "
        "of the point target, matching the live E5 provisional rule"
    ),
    "calibration_scale": E5_YOUNG_OWNERSHIP_CALIBRATION_SCALE,
    "weight": E5_YOUNG_OWNERSHIP_WEIGHT,
    "status": (
        "clean diagnostic overidentifying row; not promoted to the production "
        "target contract"
    ),
}


def e5_young_ownership_target_system() -> TargetSystem:
    """Return the live E5 target block plus the young-ownership row."""

    pending = [name for name, value in E5_TARGETS.items() if value is None]
    if pending:
        raise ValueError(f"E5 target values are pending: {pending}")
    targets = {name: float(value) for name, value in E5_TARGETS.items()}
    weights = {name: float(value) for name, value in E5_WEIGHTS.items()}
    targets["own_rate_2534"] = E5_YOUNG_OWNERSHIP_TARGET
    weights["own_rate_2534"] = E5_YOUNG_OWNERSHIP_WEIGHT
    return TargetSystem.from_mappings(
        E5_YOUNG_OWNERSHIP_TARGET_SET,
        targets,
        weights,
    )


def e5_target_system_for_profile(profile: str) -> TargetSystem:
    """Select an explicit target profile and fail closed on unknown names."""

    from .e5_profile import e5_target_system

    if profile == E5_BASELINE_TARGET_PROFILE:
        return e5_target_system()
    if profile == E5_YOUNG_OWNERSHIP_PROFILE:
        return e5_young_ownership_target_system()
    raise ValueError(f"Unknown target profile: {profile}")
