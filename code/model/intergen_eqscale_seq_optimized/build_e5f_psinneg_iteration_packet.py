#!/usr/bin/env python3
"""Build the floor-versus-tilt calibration and funded-policy comparison packet."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--floor-report", type=Path, required=True)
    parser.add_argument("--tilt-report", type=Path, required=True)
    parser.add_argument("--floor-policy", type=Path, required=True)
    parser.add_argument("--tilt-policy", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def target_comparison(floor_report: Path, tilt_report: Path) -> list[dict[str, Any]]:
    floor = read_csv(floor_report / "target_fit_full.csv")
    tilt = read_csv(tilt_report / "target_fit_full.csv")
    if len(floor) != 12 or len(tilt) != 12:
        raise ValueError("both reports must contain exactly twelve target rows")
    tilt_by_name = {row["moment"]: row for row in tilt}
    if {row["moment"] for row in floor} != set(tilt_by_name):
        raise ValueError("floor and tilt reports use different target names")
    rows: list[dict[str, Any]] = []
    for floor_row in floor:
        name = floor_row["moment"]
        tilt_row = tilt_by_name[name]
        target = float(floor_row["target"])
        weight = float(floor_row["weight"])
        if not np.isclose(target, float(tilt_row["target"]), rtol=0.0, atol=1e-12):
            raise ValueError(f"target changed for {name}")
        if not np.isclose(weight, float(tilt_row["weight"]), rtol=0.0, atol=1e-12):
            raise ValueError(f"weight changed for {name}")
        rows.append(
            {
                "moment": name,
                "target": target,
                "weight": weight,
                "floor_model": float(floor_row["model"]),
                "floor_gap": float(floor_row["gap"]),
                "floor_relative_gap": float(floor_row["gap"]) / target,
                "floor_loss_contribution": float(floor_row["loss_contribution"]),
                "tilt_model": float(tilt_row["model"]),
                "tilt_gap": float(tilt_row["gap"]),
                "tilt_relative_gap": float(tilt_row["gap"]) / target,
                "tilt_loss_contribution": float(tilt_row["loss_contribution"]),
            }
        )
    return rows


def parameter_comparison(floor_report: Path, tilt_report: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for model, report in (("floor", floor_report), ("tilt", tilt_report)):
        rows.extend(
            {"model": model, **parameter}
            for parameter in read_csv(report / "parameter_table_full.csv")
        )
    return rows


def policy_comparison(floor_policy: Path, tilt_policy: Path) -> list[dict[str, Any]]:
    expected = [
        "rebated_tax1_baseline",
        "rebated_tax2",
        "rebated_tax2_grant0p4_Hge6",
    ]
    rows: list[dict[str, Any]] = []
    for model, path in (("floor", floor_policy), ("tilt", tilt_policy)):
        current = read_csv(path / "summary.csv")
        if [row["case"] for row in current] != expected:
            raise ValueError(f"unexpected policy cases in {path}")
        rows.extend({"model": model, **row} for row in current)
    return rows


def validate_entry_contract(floor_policy: Path, tilt_policy: Path) -> float:
    """Require both policy packets to implement the empirical entry contract."""

    shares: list[float] = []
    for path in (floor_policy, tilt_policy):
        payload = json.loads((path / "benchmark_outside_objects.json").read_text())
        required = {
            "outside_origin_entrant_share_target",
            "outside_origin_entrant_share_recovered",
            "identity_scale_factor",
        }
        missing = required.difference(payload)
        if missing:
            raise ValueError(
                f"policy packet {path} does not carry the empirical entry contract; "
                f"missing {sorted(missing)}"
            )
        target = float(payload["outside_origin_entrant_share_target"])
        recovered = float(payload["outside_origin_entrant_share_recovered"])
        identity_scale = float(payload["identity_scale_factor"])
        if not np.isclose(recovered, target, rtol=0.0, atol=1.0e-10):
            raise ValueError(f"outside-origin entrant share is not reproduced in {path}")
        if not np.isclose(identity_scale, 1.0, rtol=0.0, atol=1.0e-10):
            raise ValueError(f"baseline population identity fails in {path}")
        shares.append(target)
    if not np.isclose(shares[0], shares[1], rtol=0.0, atol=1.0e-12):
        raise ValueError("floor and tilt policy packets use different entry targets")
    return shares[0]


def plot_targets(rows: list[dict[str, Any]], path: Path) -> None:
    labels = [row["moment"] for row in rows]
    y = np.arange(len(rows), dtype=float)
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.barh(
        y - 0.18,
        [100.0 * float(row["floor_relative_gap"]) for row in rows],
        height=0.34,
        label="Child-room floor",
        color="#17365D",
    )
    ax.barh(
        y + 0.18,
        [100.0 * float(row["tilt_relative_gap"]) for row in rows],
        height=0.34,
        label="Expenditure-share tilts",
        color="#C45A3B",
    )
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Model minus target, percent of target")
    ax.set_title("Common 12-target fit")
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_policies(rows: list[dict[str, Any]], path: Path, outside_origin_share: float) -> None:
    cases = ["rebated_tax2", "rebated_tax2_grant0p4_Hge6"]
    labels = ["2% tax", "2% tax + purchase grant"]
    metrics = [
        ("tfr_change_percent", "Completed fertility (%)"),
        ("price_change_percent", "House price (%)"),
        ("population_change_percent", "Population (%)"),
        ("total_births_change_percent", "Total births (%)"),
    ]
    by_model_case = {(row["model"], row["case"]): row for row in rows}
    x = np.arange(len(cases), dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for axis, (metric, title) in zip(axes.reshape(-1), metrics):
        axis.bar(
            x - 0.18,
            [float(by_model_case[("floor", case)][metric]) for case in cases],
            width=0.34,
            color="#17365D",
            label="Child-room floor",
        )
        axis.bar(
            x + 0.18,
            [float(by_model_case[("tilt", case)][metric]) for case in cases],
            width=0.34,
            color="#C45A3B",
            label="Expenditure-share tilts",
        )
        axis.axhline(0.0, color="black", linewidth=0.7)
        axis.set_title(title, loc="left")
        axis.grid(axis="y", alpha=0.25)
        axis.set_xticks(x, labels)
    axes[0, 1].legend(frameon=False, loc="best")
    fig.suptitle(
        "Funded policy effects relative to each model's 1% tax baseline\n"
        f"Baseline outside-origin entrant share: {100.0 * outside_origin_share:.1f}%"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    targets = target_comparison(args.floor_report, args.tilt_report)
    parameters = parameter_comparison(args.floor_report, args.tilt_report)
    outside_origin_share = validate_entry_contract(args.floor_policy, args.tilt_policy)
    policies = policy_comparison(args.floor_policy, args.tilt_policy)
    write_csv(args.outdir / "target_fit_comparison_full.csv", targets)
    write_csv(args.outdir / "parameter_comparison_full.csv", parameters)
    write_csv(args.outdir / "policy_comparison_full.csv", policies)
    plot_targets(targets, args.outdir / "target_relative_gap_comparison.png")
    plot_policies(
        policies,
        args.outdir / "funded_policy_comparison.png",
        outside_origin_share,
    )


if __name__ == "__main__":
    main()
