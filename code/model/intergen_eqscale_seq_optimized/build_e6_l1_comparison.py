#!/usr/bin/env python3
"""Compare certified E6AB winners under canonical, squared, and L1 objectives."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

MODEL_ROOT = Path(__file__).resolve().parents[1]
if str(MODEL_ROOT) not in sys.path:
    sys.path.insert(0, str(MODEL_ROOT))

from intergen_eqscale_seq_optimized.build_e6_plain_comparison import BLOCKS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--l1-source", type=Path, required=True)
    parser.add_argument("--squared-source", type=Path, required=True)
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


def evaluate(
    moments: dict[str, float],
    targets: dict[str, float],
    l1_weights: dict[str, float],
    squared_weights: dict[str, float],
    canonical_weights: dict[str, float],
) -> dict[str, float]:
    result = {
        "l1_loss": sum(
            l1_weights[name] * abs(float(moments[name]) - target)
            for name, target in targets.items()
        ),
        "squared_proportional_loss": sum(
            squared_weights[name] * (float(moments[name]) - target) ** 2
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
        result[f"{block}_mean_absolute_percentage_gap"] = sum(
            abs((float(moments[name]) - targets[name]) / targets[name])
            for name in names
        ) / len(names)
    return result


def parameter_value(theta: dict[str, Any], name: str) -> float:
    value = float(theta["beta" if name == "beta_annual" else name])
    return value**0.25 if name == "beta_annual" else value


def main() -> None:
    args = parse_args()
    payloads = {
        "canonical_E6AB_rescue": json.loads(args.canonical_source.read_text()),
        "squared_proportional_E6AB": json.loads(args.squared_source.read_text()),
        "l1_proportional_E6AB": json.loads(args.l1_source.read_text()),
    }
    winners = {label: winner(payload) for label, payload in payloads.items()}
    l1_metadata = payloads["l1_proportional_E6AB"].get("winner_metadata") or {}
    squared_metadata = payloads["squared_proportional_E6AB"].get("winner_metadata") or {}
    targets = {
        name: float(value) for name, value in (l1_metadata.get("targets") or {}).items()
    }
    l1_weights = {
        name: float(value) for name, value in (l1_metadata.get("weights") or {}).items()
    }
    squared_weights = {
        name: float(value)
        for name, value in (squared_metadata.get("weights") or {}).items()
    }
    canonical_weights = {
        name: float(value)
        for name, value in (l1_metadata.get("canonical_weights") or {}).items()
    }
    if not targets or not all(
        set(targets) == set(weights)
        for weights in (l1_weights, squared_weights, canonical_weights)
    ):
        raise ValueError("sources lack the complete three-objective twelve-row contract")
    for label, payload in payloads.items():
        source_targets = {
            name: float(value)
            for name, value in ((payload.get("winner_metadata") or {}).get("targets") or {}).items()
        }
        if source_targets != targets:
            raise ValueError(f"{label} uses a different target contract")

    moments = {
        label: {name: float(row["moments"][name]) for name in targets}
        for label, row in winners.items()
    }
    objective_rows = [
        {
            "calibration": label,
            **evaluate(row, targets, l1_weights, squared_weights, canonical_weights),
        }
        for label, row in moments.items()
    ]

    block_for = {name: block for block, names in BLOCKS.items() for name in names}
    moment_rows: list[dict[str, Any]] = []
    for name, target in targets.items():
        row: dict[str, Any] = {
            "block": block_for[name],
            "moment": name,
            "target": target,
            "l1_weight": l1_weights[name],
            "squared_weight": squared_weights[name],
            "canonical_weight": canonical_weights[name],
        }
        for label in moments:
            model = moments[label][name]
            gap = model - target
            row[f"{label}_model"] = model
            row[f"{label}_gap"] = gap
            row[f"{label}_relative_gap"] = gap / target
            row[f"{label}_l1_contribution"] = l1_weights[name] * abs(gap)
            row[f"{label}_squared_contribution"] = squared_weights[name] * gap**2
            row[f"{label}_canonical_contribution"] = canonical_weights[name] * gap**2
        moment_rows.append(row)

    parameter_rows: list[dict[str, Any]] = []
    for spec in l1_metadata.get("active_domain") or []:
        name = str(spec["name"])
        lower, upper = float(spec["lower"]), float(spec["upper"])
        row = {
            "parameter": name,
            "lower_bound": lower,
            "upper_bound": upper,
        }
        for label, selected in winners.items():
            estimate = parameter_value(selected["theta"], name)
            row[f"{label}_estimate"] = estimate
            row[f"{label}_near_bound_2pct"] = min(
                estimate - lower, upper - estimate
            ) <= 0.02 * (upper - lower)
        parameter_rows.append(row)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "objective_comparison.csv", objective_rows)
    write_csv(args.outdir / "moment_comparison.csv", moment_rows)
    write_csv(args.outdir / "parameter_comparison.csv", parameter_rows)

    import matplotlib.pyplot as plt
    import numpy as np

    labels = [row["moment"] for row in moment_rows]
    y = np.arange(len(labels), dtype=float)
    width = 0.24
    plot_labels = {
        "canonical_E6AB_rescue": "Canonical rescue",
        "squared_proportional_E6AB": "Squared proportional gaps",
        "l1_proportional_E6AB": "Absolute proportional gaps",
    }
    offsets = (-width, 0.0, width)
    fig, ax = plt.subplots(figsize=(12, 8))
    for offset, label in zip(offsets, moments):
        values = [
            100.0 * float(row[f"{label}_relative_gap"]) for row in moment_rows
        ]
        ax.barh(y + offset, values, height=width, label=plot_labels[label])
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Model minus target, percent of target")
    ax.set_title("E6AB target gaps under three calibration objectives")
    ax.legend(loc="lower right")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.outdir / "relative_gap_comparison.png", dpi=180)
    plt.close(fig)

    summary = {
        "sources": {
            "l1": str(args.l1_source),
            "squared": str(args.squared_source),
            "canonical": str(args.canonical_source),
        },
        "metrics": {row["calibration"]: row for row in objective_rows},
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    (args.outdir / "README.md").write_text(
        "\n".join(
            [
                "# E6AB three-objective comparison",
                "",
                "This packet evaluates three certified E6AB estimates against the same",
                "twelve targets under absolute proportional, squared proportional,",
                "canonical, raw absolute, and mean absolute percentage criteria.",
                "The first three CSV rows are directly comparable across estimates; none",
                "of these scalar calibration criteria is interpreted as a formal J test.",
                "",
                "- `objective_comparison.csv` reports all aggregate and block criteria.",
                "- `moment_comparison.csv` reports every target, model moment, gap, weight, and contribution.",
                "- `parameter_comparison.csv` reports all ten estimates and bounds.",
                "- `relative_gap_comparison.png` plots all signed proportional gaps.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
