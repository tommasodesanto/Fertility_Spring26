"""Transparent alternative weights for the E6 calibration diagnostic.

The canonical E5/E6 objective mixes measured, declared, and synthetic
standard errors.  This module leaves every target and model moment unchanged
and replaces only the diagonal weights.  Each gap is divided by the absolute
target, each row receives equal weight within its economic block, and the
three blocks receive equal aggregate weight.

The resulting loss is the sum of the three block mean squared proportional
gaps.  It is a diagnostic calibration objective, not a statistical J test.
"""

from __future__ import annotations

from collections.abc import Mapping


PLAIN_WEIGHT_SCHEME = "target_relative_block_equal"
PLAIN_L1_WEIGHT_SCHEME = "target_relative_block_equal_l1"

MOMENT_BLOCKS: dict[str, tuple[str, ...]] = {
    "fertility": (
        "tfr",
        "childless_rate",
        "mean_age_first_birth",
        "share_first_births_age30plus",
    ),
    "housing_tenure": (
        "housing_increment_0to1",
        "prime30_55_parent_3plus_minus_1to2_mean_rooms",
        "own_family_gap",
        "own_rate",
        "aggregate_mean_occupied_rooms_18_85",
    ),
    "wealth_bequest": (
        "aggregate_wealth_to_annual_gross_labor_earnings",
        "annual_bequest_flow_to_aggregate_wealth",
        "old_total_wealth_to_annual_income_p90_p50_7684",
    ),
}


def target_relative_block_equal_weights(
    targets: Mapping[str, float],
) -> dict[str, float]:
    """Return weights for a sum of block mean squared proportional gaps."""

    expected = {name for names in MOMENT_BLOCKS.values() for name in names}
    actual = set(targets)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(
            "plain E6 weighting requires the exact twelve-row target contract; "
            f"missing={missing}, extra={extra}"
        )

    weights: dict[str, float] = {}
    for names in MOMENT_BLOCKS.values():
        block_size = float(len(names))
        for name in names:
            target = abs(float(targets[name]))
            if target <= 0.0:
                raise ValueError(f"plain E6 weighting requires nonzero target {name!r}")
            weights[name] = 1.0 / (block_size * target**2)
    return {name: weights[name] for name in targets}


def target_relative_block_equal_l1_weights(
    targets: Mapping[str, float],
) -> dict[str, float]:
    """Return coefficients for a sum of block mean absolute proportional gaps."""

    expected = {name for names in MOMENT_BLOCKS.values() for name in names}
    actual = set(targets)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(
            "plain E6 L1 weighting requires the exact twelve-row target contract; "
            f"missing={missing}, extra={extra}"
        )

    weights: dict[str, float] = {}
    for names in MOMENT_BLOCKS.values():
        block_size = float(len(names))
        for name in names:
            target = abs(float(targets[name]))
            if target <= 0.0:
                raise ValueError(f"plain E6 L1 weighting requires nonzero target {name!r}")
            weights[name] = 1.0 / (block_size * target)
    return {name: weights[name] for name in targets}
