#!/usr/bin/env python3
"""Compare the plain-weight E6AB winner with the canonical E6AB rescue."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


BLOCKS: dict[str, tuple[str, ...]] = {
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plain-source", type=Path, required=True)
    parser.add_argument("--canonical-source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def winner(payload: dict[str, Any]) -> dict[str, Any]:
    row = (payload.get("winners") or {}).get("E1")
    if not isinstance(row, dict):
        raise ValueError("source omits winners.E1")
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_relative_gap_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    """Plot signed percentage gaps under the two calibrated objectives."""

    import matplotlib.pyplot as plt
    import numpy as np

    labels = [str(row["moment"]) for row in rows]
    y = np.arange(len(rows), dtype=float)
    canonical = 100.0 * np.array(
        [float(row["canonical_relative_gap"]) for row in rows]
    )
    plain = 100.0 * np.array([float(row["plain_relative_gap"]) for row in rows])
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.barh(y + 0.18, canonical, height=0.34, label="Canonical E6AB rescue")
    ax.barh(y - 0.18, plain, height=0.34, label="Block-equal proportional-gap E6AB")
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Model minus target, percent of target")
    ax.set_title("E6AB target gaps under alternative calibration objectives")
    ax.legend(loc="lower right")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def metrics(
    moments: dict[str, float],
    targets: dict[str, float],
    plain_weights: dict[str, float],
    canonical_weights: dict[str, float],
) -> dict[str, float]:
    result = {
        "plain_loss": sum(
            plain_weights[name] * (float(moments[name]) - target) ** 2
            for name, target in targets.items()
        ),
        "canonical_loss": sum(
            canonical_weights[name] * (float(moments[name]) - target) ** 2
            for name, target in targets.items()
        ),
        "raw_absolute_gap_sum": sum(
            abs(float(moments[name]) - target) for name, target in targets.items()
        ),
        "mean_absolute_percentage_gap": sum(
            abs((float(moments[name]) - target) / target)
            for name, target in targets.items()
        ) / len(targets),
    }
    for block, names in BLOCKS.items():
        result[f"{block}_plain_loss"] = sum(
            plain_weights[name] * (float(moments[name]) - targets[name]) ** 2
            for name in names
        )
        result[f"{block}_canonical_loss"] = sum(
            canonical_weights[name] * (float(moments[name]) - targets[name]) ** 2
            for name in names
        )
        result[f"{block}_mean_absolute_percentage_gap"] = sum(
            abs((float(moments[name]) - targets[name]) / targets[name])
            for name in names
        ) / len(names)
    return result


def main() -> None:
    args = parse_args()
    plain_payload = json.loads(args.plain_source.read_text())
    canonical_payload = json.loads(args.canonical_source.read_text())
    plain_winner, canonical_winner = winner(plain_payload), winner(canonical_payload)
    metadata = plain_payload.get("winner_metadata") or {}
    targets = {name: float(value) for name, value in (metadata.get("targets") or {}).items()}
    plain_weights = {name: float(value) for name, value in (metadata.get("weights") or {}).items()}
    canonical_weights = {
        name: float(value)
        for name, value in (metadata.get("canonical_weights") or {}).items()
    }
    if not targets or set(targets) != set(plain_weights) or set(targets) != set(canonical_weights):
        raise ValueError("plain source lacks a complete dual-weight twelve-row contract")

    canonical_metadata = canonical_payload.get("winner_metadata") or {}
    canonical_targets = canonical_metadata.get("targets") or {}
    if targets != {name: float(value) for name, value in canonical_targets.items()}:
        raise ValueError("plain and canonical sources use different targets")

    plain_moments = {
        name: float(plain_winner["moments"][name]) for name in targets
    }
    canonical_moments = {
        name: float(canonical_winner["moments"][name]) for name in targets
    }
    plain_metrics = metrics(plain_moments, targets, plain_weights, canonical_weights)
    canonical_metrics = metrics(canonical_moments, targets, plain_weights, canonical_weights)

    moment_rows: list[dict[str, Any]] = []
    block_for = {name: block for block, names in BLOCKS.items() for name in names}
    for name, target in targets.items():
        old_model, new_model = canonical_moments[name], plain_moments[name]
        old_gap, new_gap = old_model - target, new_model - target
        moment_rows.append(
            {
                "block": block_for[name],
                "moment": name,
                "target": target,
                "canonical_model": old_model,
                "plain_model": new_model,
                "canonical_gap": old_gap,
                "plain_gap": new_gap,
                "canonical_relative_gap": old_gap / target,
                "plain_relative_gap": new_gap / target,
                "plain_weight": plain_weights[name],
                "canonical_plain_contribution": plain_weights[name] * old_gap**2,
                "plain_plain_contribution": plain_weights[name] * new_gap**2,
                "canonical_weight": canonical_weights[name],
                "canonical_canonical_contribution": canonical_weights[name] * old_gap**2,
                "plain_canonical_contribution": canonical_weights[name] * new_gap**2,
            }
        )

    parameter_rows: list[dict[str, Any]] = []
    old_theta, new_theta = canonical_winner["theta"], plain_winner["theta"]
    for spec in metadata.get("active_domain") or []:
        name = str(spec["name"])
        theta_name = "beta" if name == "beta_annual" else name
        old_value, new_value = float(old_theta[theta_name]), float(new_theta[theta_name])
        if name == "beta_annual":
            old_value, new_value = old_value**0.25, new_value**0.25
        lower, upper = float(spec["lower"]), float(spec["upper"])
        distance = min(new_value - lower, upper - new_value)
        parameter_rows.append(
            {
                "parameter": name,
                "canonical_estimate": old_value,
                "plain_estimate": new_value,
                "plain_minus_canonical": new_value - old_value,
                "lower_bound": lower,
                "upper_bound": upper,
                "plain_near_bound_2pct": distance <= 0.02 * (upper - lower),
            }
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    objective_rows = [
        {"calibration": "canonical_E6AB_rescue", **canonical_metrics},
        {"calibration": "plain_block_equal_E6AB", **plain_metrics},
    ]
    write_csv(args.outdir / "objective_comparison.csv", objective_rows)
    write_csv(args.outdir / "moment_comparison.csv", moment_rows)
    write_csv(args.outdir / "parameter_comparison.csv", parameter_rows)
    write_relative_gap_plot(args.outdir / "relative_gap_comparison.png", moment_rows)
    summary = {
        "plain_source": str(args.plain_source),
        "canonical_source": str(args.canonical_source),
        "plain_metrics": plain_metrics,
        "canonical_metrics": canonical_metrics,
        "plain_loss_improvement_fraction": (
            canonical_metrics["plain_loss"] - plain_metrics["plain_loss"]
        ) / canonical_metrics["plain_loss"],
        "canonical_loss_change_fraction": (
            plain_metrics["canonical_loss"] - canonical_metrics["canonical_loss"]
        ) / canonical_metrics["canonical_loss"],
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    (args.outdir / "README.md").write_text(
        "\n".join(
            [
                "# E6AB plain-weight comparison",
                "",
                "This packet compares the certified block-equal proportional-gap winner",
                "with the canonical E6AB rescue on the same twelve targets.",
                "The plain loss is the sum of the three block mean squared proportional gaps.",
                "It is a transparent calibration criterion, not a statistical J test.",
                "",
                f"- Plain loss: `{plain_metrics['plain_loss']:.12f}` versus `{canonical_metrics['plain_loss']:.12f}` at the canonical winner.",
                f"- Canonical loss: `{plain_metrics['canonical_loss']:.12f}` versus `{canonical_metrics['canonical_loss']:.12f}` at the canonical winner.",
                f"- Mean absolute percentage gap: `{plain_metrics['mean_absolute_percentage_gap']:.6%}` versus `{canonical_metrics['mean_absolute_percentage_gap']:.6%}`.",
                "- `moment_comparison.csv` reports every target, both model moments, both gaps, both weights, and both contribution systems.",
                "- `parameter_comparison.csv` reports all ten free parameters and bounds.",
                "- `relative_gap_comparison.png` plots every signed proportional target gap under both estimates.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
