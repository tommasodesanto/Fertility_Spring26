#!/usr/bin/env python3
"""Collect strict winners from the two-arm mortality recalibration."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    summaries = [json.loads(path.read_text()) for path in sorted(args.results_dir.glob("task_*/summary.json"))]
    if len(summaries) != 8:
        raise RuntimeError(f"expected 8 completed chain summaries, found {len(summaries)}")

    eligible: dict[str, list[dict[str, Any]]] = {"M0": [], "M1": []}
    chain_rows: list[dict[str, Any]] = []
    for summary in summaries:
        metadata = summary["metadata"]
        arm = str(metadata["arm"])
        tight = summary.get("best_tight")
        repeat = summary.get("tight_repeat_check") or {}
        ok = bool(
            tight
            and tight.get("strict_converged")
            and repeat.get("both_strict")
            and float(repeat.get("loss_abs_difference", 1.0)) <= 1e-10
            and float(repeat.get("max_abs_moment_difference", 1.0)) <= 1e-10
        )
        chain_rows.append(
            {
                "arm": arm,
                "seed": metadata["seed"],
                "start_mix": metadata["start_mix"],
                "eligible": ok,
                "search_cases": summary["n_cases_completed"],
                "search_strict": summary["n_strict"],
                "tight_loss": tight["rank_loss"] if tight else None,
                "tight_residual": tight["market_residual"] if tight else None,
            }
        )
        if ok:
            eligible[arm].append(tight)

    if any(not eligible[arm] for arm in ("M0", "M1")):
        raise RuntimeError(f"missing eligible tight winner: { {arm: len(rows) for arm, rows in eligible.items()} }")
    winners = {arm: min(rows, key=lambda row: float(row["rank_loss"])) for arm, rows in eligible.items()}

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "chain_summary.csv", chain_rows)
    fit_rows: list[dict[str, Any]] = []
    parameter_rows: list[dict[str, Any]] = []
    lifecycle_rows: list[dict[str, Any]] = []
    for arm, winner in winners.items():
        fit_rows.extend({"arm": arm, **row} for row in winner["target_fit"])
        parameter_rows.extend({"arm": arm, **row} for row in winner["parameters"])
        life = winner["lifecycle"]
        for age, mass, own, owner_rooms, occupied_rooms in zip(
            life["ages"],
            life["age_mass"],
            life["ownership_by_age"],
            life["owner_rooms_by_age"],
            life["occupied_rooms_by_age"],
        ):
            lifecycle_rows.append(
                {
                    "arm": arm,
                    "age": age,
                    "age_mass": mass,
                    "ownership_rate": own,
                    "owner_rooms": owner_rooms,
                    "occupied_rooms": occupied_rooms,
                }
            )
    write_csv(args.outdir / "target_fit_full.csv", fit_rows)
    write_csv(args.outdir / "parameter_table_full.csv", parameter_rows)
    write_csv(args.outdir / "lifecycle_decomposition.csv", lifecycle_rows)

    comparison = []
    for arm, winner in winners.items():
        life = winner["lifecycle"]
        comparison.append(
            {
                "arm": arm,
                "loss": winner["rank_loss"],
                "market_residual": winner["market_residual"],
                "price": winner["price"],
                "old_household_mass_share_62plus": life["old_household_mass_share_62plus"],
                "old_owner_rooms_share_62plus": life["old_owner_rooms_share_62plus"],
                "old_occupied_rooms_share_62plus": life["old_occupied_rooms_share_62plus"],
                "ownership_path_acceptance": life["hard_acceptance_pass"],
            }
        )
    write_csv(args.outdir / "headline_comparison.csv", comparison)
    payload = {"winners": winners, "eligible_chain_counts": {arm: len(rows) for arm, rows in eligible.items()}}
    (args.outdir / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True))

    m0, m1 = comparison
    lines = [
        "# Post-retirement mortality recalibration",
        "",
        "Every arm re-estimates the same 11 parameters against all 15 moments. Only strict, exactly repeated tight winners are eligible.",
        "",
        "| Arm | Loss | Residual | Old household mass 62+ | Old share of owner rooms | Old share of all rooms | Ownership path passes |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in comparison:
        lines.append(
            f"| {row['arm']} | {row['loss']:.6f} | {row['market_residual']:.2e} | "
            f"{row['old_household_mass_share_62plus']:.4f} | {row['old_owner_rooms_share_62plus']:.4f} | "
            f"{row['old_occupied_rooms_share_62plus']:.4f} | {row['ownership_path_acceptance']} |"
        )
    lines.extend(
        [
            "",
            f"Survival minus control loss: {float(m1['loss']) - float(m0['loss']):+.6f}.",
            "",
            "Complete target, parameter, and lifecycle tables are in the adjacent CSV files.",
        ]
    )
    (args.outdir / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
