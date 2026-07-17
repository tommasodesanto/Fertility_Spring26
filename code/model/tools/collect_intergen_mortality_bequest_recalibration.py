#!/usr/bin/env python3
"""Collect the clean mortality-plus-child-blind-bequest calibration."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:  # Torch's batch environment omits plotting packages.
    plt = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--mortality-results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    if not keys:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def eligible_tight(summary: dict[str, Any]) -> dict[str, Any] | None:
    tight = summary.get("best_tight")
    repeat = summary.get("tight_repeat_check") or {}
    if not (
        isinstance(tight, dict)
        and tight.get("strict_converged")
        and repeat.get("both_strict")
        and float(repeat.get("loss_abs_difference", 1.0)) <= 1e-10
        and float(repeat.get("max_abs_moment_difference", 1.0)) <= 1e-10
    ):
        return None
    return tight


def headline_row(arm: str, winner: dict[str, Any]) -> dict[str, Any]:
    life = winner["lifecycle"]
    return {
        "arm": arm,
        "loss": winner["rank_loss"],
        "market_residual": winner["market_residual"],
        "theta0": winner["theta"].get("theta0", 0.0),
        "old_household_mass_share_62plus": life["old_household_mass_share_62plus"],
        "old_owner_rooms_share_62plus": life["old_owner_rooms_share_62plus"],
        "old_occupied_rooms_share_62plus": life["old_occupied_rooms_share_62plus"],
        "ownership_path_acceptance": life["hard_acceptance_pass"],
        "decumulation_ratio_74plus_over_62_74": life["decum_ratio_wealth_74plus_over_62_74"],
    }


def make_plots(outdir: Path, winners: dict[str, dict[str, Any]]) -> None:
    if plt is None:
        return
    colors = {"M1": "#4C78A8", "M2": "#E45756"}
    labels = {"M1": "Mortality, no bequest", "M2": "Mortality + child-blind bequest"}
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    fields = (
        ("ownership_by_age", "Ownership rate"),
        ("age_mass", "Age mass"),
        ("mean_liquid_wealth_by_age", "Mean liquid wealth"),
        ("occupied_rooms_by_age", "Occupied rooms (aggregate)"),
    )
    for arm, winner in winners.items():
        life = winner["lifecycle"]
        ages = life["ages"]
        for ax, (field, title) in zip(axes.flat, fields):
            ax.plot(ages, life[field], marker="o", ms=3, lw=1.7, color=colors[arm], label=labels[arm])
            ax.set_title(title)
            ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("Age")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "lifecycle_comparison.png", dpi=180)
    plt.close(fig)

    moments = [row["moment"] for row in winners["M1"]["target_fit"]]
    contributions = {
        arm: {row["moment"]: float(row["loss_contribution"]) for row in winner["target_fit"]}
        for arm, winner in winners.items()
    }
    order = sorted(moments, key=lambda name: max(contributions[arm][name] for arm in winners), reverse=True)
    y = list(range(len(order)))
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh([v + 0.18 for v in y], [contributions["M1"][name] for name in order], height=0.34, color=colors["M1"], label=labels["M1"])
    ax.barh([v - 0.18 for v in y], [contributions["M2"][name] for name in order], height=0.34, color=colors["M2"], label=labels["M2"])
    ax.set_yticks(y, order)
    ax.invert_yaxis()
    ax.set_xlabel("Weighted loss contribution")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "loss_contributions.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    summaries = [json.loads(path.read_text()) for path in sorted(args.results_dir.glob("task_*/summary.json"))]
    if len(summaries) != 4:
        raise RuntimeError(f"expected 4 completed M2 chain summaries, found {len(summaries)}")

    baseline_payload = json.loads(args.mortality_results.read_text())
    m1 = baseline_payload.get("winners", {}).get("M1")
    if not isinstance(m1, dict):
        raise RuntimeError("mortality result does not contain the strict M1 winner")

    eligible: list[dict[str, Any]] = []
    chain_rows: list[dict[str, Any]] = []
    for summary in summaries:
        metadata = summary["metadata"]
        tight = eligible_tight(summary)
        chain_rows.append(
            {
                "arm": metadata["arm"],
                "seed": metadata["seed"],
                "start_mix": metadata["start_mix"],
                "eligible": tight is not None,
                "search_cases": summary["n_cases_completed"],
                "search_strict": summary["n_strict"],
                "tight_loss": tight["rank_loss"] if tight else None,
                "tight_residual": tight["market_residual"] if tight else None,
            }
        )
        if tight is not None:
            eligible.append(tight)
    if not eligible:
        raise RuntimeError("no strict, exactly repeated M2 winner")
    m2 = min(eligible, key=lambda row: float(row["rank_loss"]))
    if float(m2["rank_loss"]) > float(m1["rank_loss"]) + 1e-8:
        raise RuntimeError("M2 missed its nested M1 zero-bequest seed")

    winners = {"M1": m1, "M2": m2}
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "chain_summary.csv", chain_rows)
    write_csv(
        args.outdir / "target_fit_full.csv",
        [{"arm": arm, **row} for arm, winner in winners.items() for row in winner["target_fit"]],
    )
    write_csv(
        args.outdir / "parameter_table_full.csv",
        [{"arm": arm, **row} for arm, winner in winners.items() for row in winner["parameters"]],
    )
    write_csv(
        args.outdir / "lifecycle_decomposition.csv",
        [
            {
                "arm": arm,
                "age": age,
                "age_mass": mass,
                "ownership_rate": own,
                "mean_liquid_wealth": liquid,
                "owner_rooms": owner_rooms,
                "occupied_rooms": occupied_rooms,
            }
            for arm, winner in winners.items()
            for age, mass, own, liquid, owner_rooms, occupied_rooms in zip(
                winner["lifecycle"]["ages"],
                winner["lifecycle"]["age_mass"],
                winner["lifecycle"]["ownership_by_age"],
                winner["lifecycle"]["mean_liquid_wealth_by_age"],
                winner["lifecycle"]["owner_rooms_by_age"],
                winner["lifecycle"]["occupied_rooms_by_age"],
            )
        ],
    )
    comparison = [headline_row(arm, winner) for arm, winner in winners.items()]
    write_csv(args.outdir / "headline_comparison.csv", comparison)
    (args.outdir / "results.json").write_text(
        json.dumps({"winners": winners, "eligible_M2_chains": len(eligible)}, indent=2, sort_keys=True)
    )
    make_plots(args.outdir, winners)

    delta = float(m2["rank_loss"]) - float(m1["rank_loss"])
    lines = [
        "# Mortality plus child-blind bequest recalibration",
        "",
        "M1 and M2 use the same 15 moments and externally pinned survival schedule. M1 re-estimates all 11 clean-frontier parameters; M2 re-estimates those 11 plus bequest strength, with the luxury shift externally fixed at 0.25.",
        "",
        "| Arm | Loss | Residual | theta0 | Old household mass 62+ | Old share of all rooms | Ownership path passes |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in comparison:
        lines.append(
            f"| {row['arm']} | {float(row['loss']):.6f} | {float(row['market_residual']):.2e} | "
            f"{float(row['theta0']):.6f} | {float(row['old_household_mass_share_62plus']):.4f} | "
            f"{float(row['old_occupied_rooms_share_62plus']):.4f} | {row['ownership_path_acceptance']} |"
        )
    lines.extend(
        [
            "",
            f"M2 minus M1 loss: {delta:+.6f}.",
            "",
            "Complete target, parameter, lifecycle, and visual diagnostic files are adjacent.",
        ]
    )
    (args.outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
